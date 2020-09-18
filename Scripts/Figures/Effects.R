
# Load packages and data --------------------------------------------------
library(tidyverse)
library(tidytext)
library(glue)
library(ggtext)

morphology = readRDS("Intermediate/Morphology/Effects.rds") %>%
  filter(Trait != "RWC")
photosynthesis = readRDS("Intermediate/Photosynthesis/Effects.rds") %>%
  mutate(Genotype = "Average")
fitness = readRDS("Intermediate/Fitness/Effects.rds") %>% mutate(mu = NA)
medsize = readRDS("Intermediate/Biosorter/Effects.rds")
rates = readRDS("Intermediate/Rates/Effects.rds")

# Combine data for different traits ---------------------------------------

# Merge traits and calculate quantiles
effects = do.call("rbind", list(morphology, photosynthesis, fitness, medsize, rates)) %>%
  select(-mu)

# Add category of trait for coloring in figures
category = c(Angle = "Morphology", DW = "Performance", LeafSize = "Morphology",
             N = "Physiology", Nmol = "Physiology", PetioleRatio = "Morphology", 
             RootLength = "Morphology", ShapeBlade = "Morphology", SLA = "Morphology", 
             A = "Physiology", Chl = "Physiology", gsw = "Physiology", LUE = "Physiology",
             Rd = "Physiology", FT = "Development", Phyllochron = "Development", Yield = "Performance",
             SeedSize = "Morphology", RGR = "Performance", LeafRate = "Development",
             LN = "Development")

category_color = c(Development = "#F8766D", Morphology = "#7CAE00",
                   Performance = "#00BFC4", Physiology = "#C77CFF")
  
effects = mutate(effects, Category = category[Trait], 
                          Color = category_color[Category])

# Given proper names to the genotypes and the traits
pretty_names = c(An1 = "An-1", Bay0 = "Bay-0", Col0 = "Col-0",
                 Average = "Average", Lp26 = "Lp2-6")
effects = mutate(effects, Genotype = pretty_names[Genotype])

# Reorder genotype and effects
effects$Genotype = factor(effects$Genotype, 
                          levels = rev(c("Average", "An-1","Bay-0","Col-0","Lp2-6")))

# Rename the effects and reorder them
effects$Effect = factor(effects$Effect, 
                        levels = c("HT", "HTD", "HT_D", "PS",  "PSD", "PS_D", "S", "D"))
levels(effects$Effect)[3] = "HT:D"
levels(effects$Effect)[6] = "PS:D"

effects = ungroup(effects)

sum_effects = effects  %>%
  group_by(Genotype, Trait, Effect) %>%
  summarise(mu  = median(Value), 
            lwr = quantile(Value, 0.05),
            upr = quantile(Value, 0.95),
            Category = Category[1],
            Color = Color[1])

# Separate proper traits and performance indices
sum_trait_effects = sum_effects

# Include version without photosynthesis for comparison of genotypes to average
sum_trait_effects_nophot = filter(sum_trait_effects, 
                                  !(Trait %in% c("A", "Chl", "gsw", "LUE","Rd")))


abbrev = c(Angle = "Angle", 
                   DW = "DW", 
                   LeafSize = "BA",
                   N = "N", 
                   Nmol = 'N', 
                   PetioleRatio = "PR", 
                   RootLength = "RL",
                   ShapeBlade = "BS", 
                   SLA = "SPA", 
                   A = "A", 
                   Chl = "Chl",
                   gsw = "g", 
                   LUE = "LUE",
                   Rd = "R", 
                   FT = "FT",
                   Phyllochron = "PHY", 
                   Yield = "Y",
                   SeedSize = "SL", 
                   RGR = "RGR",
                   LeafRate = "LAR",
                   LN = "LN")

abbrev_sub = c(Angle = "", 
               DW = "", 
               LeafSize = "",
               N = "", 
               Nmol = 'area', 
               PetioleRatio = "", 
               RootLength = "",
               ShapeBlade = "", 
               SLA = "", 
               A = "", 
               Chl = "",
               gsw = "s", 
               LUE = "",
               Rd = "d", 
               FT = "",
               Phyllochron = "", 
               Yield = "",
               SeedSize = "", 
               RGR = "",
               LeafRate = "",
               LN = "")


sum_trait_effects$Abb = abbrev[sum_trait_effects$Trait]
sum_trait_effects_nophot$Abb = abbrev[sum_trait_effects_nophot$Trait]
sum_trait_effects$Trait_sb = abbrev_sub[sum_trait_effects$Trait]
sum_trait_effects_nophot$Trait_sb = abbrev_sub[sum_trait_effects_nophot$Trait]

sum_trait_effects = mutate(sum_trait_effects, 
      md_name = glue("<i style='color:{Color}'>{Abb}<sub>{Trait_sb}</sub></i>"))
sum_trait_effects_nophot = mutate(sum_trait_effects_nophot, 
      md_name = glue("<i style='color:{Color}'>{Abb}<sub>{Trait_sb}</sub></i>"))

# Version of the scale_x_reordered function to parse the labels in the x axis
# scale_x_reordered_parsed = function (..., sep = "___") {
#   reg <- paste0(sep, ".+$")
#   ggplot2::scale_x_discrete(labels = function(x) parse(text = gsub(reg, "", x)), ...)
# }

# Effects on traits
avg_traits = ggplot(data = filter(sum_trait_effects, Genotype == "Average"), 
                    aes(x = reorder_within(md_name, mu, Effect), y = mu, 
                        col = Category, shape = Category)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(width = 0.5)) +
  geom_segment(aes(xend = reorder_within(md_name, mu, Effect),
                   y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
  scale_shape_manual(values = c(Development = 1, Morphology = 2, 
                                Performance = 6, Physiology = 4)) +
  scale_x_reordered() +
  coord_flip(ylim = c(-1.2,1.0)) + 
  facet_wrap(vars(Effect), scales = "free", ncol = 3) +
  labs(y = "Relative effect on trait", x = "Trait", color = "", shape = "") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "line"))

# Calculate difference of each genotype with respect to average -----------

# Compare Col-0 to average
col0_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Col-0")), 
                     aes(x = reorder_within(md_name, mu, Effect), y = mu, 
                         col = Genotype, shape = Genotype)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(width = 0.5)) +
  geom_segment(aes(xend = reorder_within(md_name, mu, Effect),
                   y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
  scale_x_reordered() +
  scale_shape_manual(values = c(Average = 1, `Col-0` = 4)) + 
  scale_color_manual(values = c(Average = "black", `Col-0` = "gold3")) + 
  coord_flip(ylim = c(-1.2,1.0)) + 
  facet_wrap(vars(Effect), scales = "free", ncol = 3) +
  labs(y = "Relative effect on trait", x = "Trait", color = "", shape = "") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "line"))

# Compare Bay-0 to average
bay0_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Bay-0")), 
                     aes(x = reorder_within(md_name, mu, Effect), y = mu, 
                         col = Genotype, shape = Genotype)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(width = 0.5)) +
  geom_segment(aes(xend = reorder_within(md_name, mu, Effect),
                   y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
  scale_x_reordered() +
  scale_shape_manual(values = c(Average = 1, `Bay-0` = 4)) +  
  scale_color_manual(values = c(Average = "black", `Bay-0` = "gold3")) + 
  coord_flip(ylim = c(-1.2,1.0)) + 
  facet_wrap(vars(Effect), scales = "free", ncol = 3) +
  labs(y = "Relative effect on trait", x = "Trait", color = "", shape = "") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "line"))

# Compare An-1 to average
an1_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "An-1")), 
                    aes(x = reorder_within(md_name, mu, Effect), y = mu, 
                        col = Genotype, shape = Genotype)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(width = 0.5)) +
  geom_segment(aes(xend = reorder_within(md_name, mu, Effect),
                   y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
  scale_x_reordered() +
  scale_shape_manual(values = c(Average = 1, `An-1` = 4)) +  
  scale_color_manual(values = c(Average = "black", `An-1` = "gold3")) + 
  coord_flip(ylim = c(-1.2,1.0)) + 
  facet_wrap(vars(Effect), scales = "free", ncol = 3) +
  labs(y = "Relative effect on trait", x = "Trait", color = "", shape = "") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "line"))


lp26_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Lp2-6")), 
                     aes(x = reorder_within(md_name, mu, Effect), y = mu, 
                         col = Genotype, shape = Genotype)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(width = 0.5)) +
  geom_segment(aes(xend = reorder_within(md_name, mu, Effect),
                   y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
  scale_x_reordered() +
  scale_shape_manual(values = c(Average = 1, `Lp2-6` = 4)) +  
  scale_color_manual(values = c(Average = "black", `Lp2-6` = "gold3")) + 
  coord_flip(ylim = c(-1.2,1.0)) + 
  facet_wrap(vars(Effect), scales = "free", ncol = 3) +
  labs(y = "Relative effect on trait", x = "Trait", color = "", shape = "") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "line"))


# Save the figures to files
ggsave("Results/Effects/Average.png", plot = avg_traits, width = 22, height = 28,
       units = "cm",device = "png", dpi = 300)

ggsave("Results/Effects/Bay0.png", plot = bay0_traits, width = 22, height = 28,
       units = "cm",device = "png", dpi = 300)

ggsave("Results/Effects/Col0.png", plot = col0_traits, width = 22, height = 28,
       units = "cm",device = "png", dpi = 300)

ggsave("Results/Effects/LP26.png", plot = lp26_traits, width = 22, height = 28,
       units = "cm",device = "png", dpi = 300)

ggsave("Results/Effects/An1.png", plot = an1_traits, width = 22, height = 28,
       units = "cm",device = "png", dpi = 300)


# 
# # Load packages and data --------------------------------------------------
# library(tidyverse)
# library(tidytext)
# 
# morphology = readRDS("Intermediate/Morphology/Effects.rds") %>%
#   filter(Trait != "RWC")
# photosynthesis = readRDS("Intermediate/Photosynthesis/Effects.rds") %>%
#   mutate(Genotype = "Average")
# fitness = readRDS("Intermediate/Fitness/Effects.rds") %>% mutate(mu = NA)
# medsize = readRDS("Intermediate/Biosorter/Effects.rds")
# rates = readRDS("Intermediate/Rates/Effects.rds")
# 
# # Combine data for different traits ---------------------------------------
# 
# # Merge traits and calculate quantiles
# effects = do.call("rbind", list(morphology, photosynthesis, fitness, medsize, rates)) %>%
#   select(-mu)
# 
# # Add category of trait for coloring in figures
# category = c(Angle = "Morphology", DW = "Perf", LeafSize = "Morphology",
#              N = "Physiology", Nmol = "Physiology", PetioleRatio = "Morphology", 
#              RootLength = "Morphology", ShapeBlade = "Morphology", SLA = "Morphology", 
#              A = "Physiology", Chl = "Physiology", gsw = "Physiology", LUE = "Physiology",
#              Rd = "Physiology", FT = "Development", Phyllochron = "Development", Yield = "Perf",
#              SeedSize = "Morphology", RGR = "Physiology", LeafRate = "Development",
#              LN = "Development")
# 
# effects = mutate(effects, Category = category[Trait])
# 
# # Given proper names to the genotypes and the traits
# pretty_names = c(An1 = "An-1", Bay0 = "Bay-0", Col0 = "Col-0",
#                  Average = "Average", Lp26 = "Lp2-6")
# effects = mutate(effects, Genotype = pretty_names[Genotype])
# 
# # Reorder genotype and effects
# effects$Genotype = factor(effects$Genotype, 
#                     levels = rev(c("Average", "An-1","Bay-0","Col-0","Lp2-6")))
# 
# # Rename the effects and reorder them
# effects$Effect = factor(effects$Effect, 
#                    levels = c("HT", "HTD", "HT_D", "PS",  "PSD", "PS_D", "S", "D"))
# levels(effects$Effect)[3] = "HT:D"
# levels(effects$Effect)[6] = "PS:D"
# 
# effects = ungroup(effects)
# 
# sum_effects = effects  %>%
#   group_by(Genotype, Trait, Effect) %>%
#   summarise(mu  = median(Value), 
#             lwr = quantile(Value, 0.05),
#             upr = quantile(Value, 0.95),
#             Category = Category[1])
# 
# # Separate proper traits and performance indices
# sum_trait_effects = filter(sum_effects, !(Trait %in% c("Yield", "DW", "RGR")))
# sum_performance_effects = filter(sum_effects, Trait %in% c("Yield", "DW", "RGR"))
# 
# # Include version without photosynthesis for comparison of genotypes to average
# sum_trait_effects_nophot = filter(sum_trait_effects, 
#                             !(Trait %in% c("A", "Chl", "gsw", "LUE","Rd")))
# 
# # Prettify the names of the traits
# pretty_names = c(Angle = "italic(Angle)", 
#                  DW = "italic(DW)", 
#                  LeafSize = "italic(LS)",
#                  N = "italic(N)", 
#                  Nmol = 'italic(N[area])', 
#                  PetioleRatio = "italic(PR)", 
#                  RootLength = "italic(RL)",
#                  ShapeBlade = "italic(SB)", 
#                  SLA = "italic(SPA)", 
#                  A = "italic(A)", 
#                  Chl = "italic(Chl)",
#                  gsw = "italic(g[s])", 
#                  LUE = "italic(LUE)",
#                  Rd = "italic(R[d])", 
#                  FT = "italic(FT)",
#                  Phyllochron = "italic(PHY)", 
#                  Yield = "italic(Y)",
#                  SeedSize = "italic(SL)", 
#                  RGR = "italic(RGR)",
#                  LeafRate = "italic(LAR)",
#                  LN = "italic(LN)")
# 
# sum_trait_effects$Trait = pretty_names[sum_trait_effects$Trait]
# sum_performance_effects$Trait = pretty_names[sum_performance_effects$Trait]
# sum_trait_effects_nophot$Trait = pretty_names[sum_trait_effects_nophot$Trait]
# 
# # Version of the scale_x_reordered function to parse the labels in the x axis
# scale_x_reordered_parsed = function (..., sep = "___") {
#   reg <- paste0(sep, ".+$")
#   ggplot2::scale_x_discrete(labels = function(x) parse(text = gsub(reg, "", x)), ...)
# }
# 
# # Effects on traits
# avg_traits = ggplot(data = filter(sum_trait_effects, Genotype == "Average"), 
#        aes(x = reorder_within(Trait, mu, Effect), y = mu, col = Category)) + 
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.5)) +
#   geom_segment(aes(xend = reorder_within(Trait, mu, Effect),
#                    y = lwr), yend = -1.5, lty = 2, col = "lightgray") +
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = c(-0.3, 0.3), lty = 2, col = "lightgray") +
#   scale_x_reordered_parsed() +
#   coord_flip(ylim = c(-1.2,1.0)) + 
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on trait", x = "Trait", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         legend.text = element_text(size = 16),
#         legend.key.width = unit(2, "line"))
# 
# # Effects on DW and Yield
# avg_perf = ggplot(data = sum_performance_effects, 
#        aes(x = Trait, y = mu, col = Genotype)) + 
#   geom_point(position = position_dodge(width = 0.75)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.75)) +
#   geom_hline(yintercept = 0) + 
#   coord_flip(ylim = c(-1.2,1)) + 
#   scale_x_discrete(labels = parse(text = unique(sum_performance_effects$Trait))) + 
#   scale_color_manual(values = c(Average = "black", `An-1` = "red",
#                                 `Bay-0` = "blue", `Col-0` = "gold",
#                                 `Lp2-6` = "green")) +
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on performance", x = "Performance indicator", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0))
# 
# # Calculate difference of each genotype with respect to average -----------
# 
# # Compare Col-0 to average
# col0_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Col-0")), 
#        aes(x = reorder_within(Trait, mu, Effect), y = mu, col = Genotype)) + 
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 0) + 
#   scale_x_reordered_parsed() +
#   coord_flip(ylim = c(-1.2,1)) + 
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on trait", x = "Trait", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         legend.text = element_text(size = 16),
#         legend.key.width = unit(2, "line"))
# 
# # Compare Bay-0 to average
# bay0_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Bay-0")), 
#        aes(x = reorder_within(Trait, mu, Effect), y = mu, col = Genotype)) + 
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 0) + 
#   scale_x_reordered_parsed() +
#   coord_flip(ylim = c(-1.2,1.5)) + 
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on trait", x = "Trait", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         legend.text = element_text(size = 16),
#         legend.key.width = unit(2, "line"))
# 
# # Compare An-1 to average
# an1_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "An-1")), 
#        aes(x = reorder_within(Trait, mu, Effect), y = mu, col = Genotype)) + 
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 0) + 
#   scale_x_reordered_parsed() +
#   coord_flip(ylim = c(-1.2,2)) + 
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on trait", x = "Trait", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         legend.text = element_text(size = 16),
#         legend.key.width = unit(2, "line"))
# 
# 
# lp26_traits = ggplot(data = filter(sum_trait_effects_nophot, Genotype %in% c("Average", "Lp2-6")), 
#        aes(x = reorder_within(Trait, mu, Effect), y = mu, col = Genotype)) + 
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), 
#                 position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 0) + 
#   scale_x_reordered_parsed() +
#   coord_flip(ylim = c(-1.2,1)) + 
#   facet_wrap(vars(Effect), scales = "free", ncol = 3) +
#   labs(y = "Relative effect on trait", x = "Trait", color = "") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         legend.text = element_text(size = 16),
#         legend.key.width = unit(2, "line"))
# 
# 
# # Save the figures to files
# ggsave("Results/Effects/Average_perf.png", plot = avg_perf, width = 16, height = 16,
#        units = "cm",device = "png", dpi = 300)
# 
# ggsave("Results/Effects/Average_traits.png", plot = avg_traits, width = 22, height = 22,
#        units = "cm",device = "png", dpi = 300)
# 
# ggsave("Results/Effects/Bay0.png", plot = bay0_traits, width = 22, height = 22,
#        units = "cm",device = "png", dpi = 300)
# 
# ggsave("Results/Effects/Col0.png", plot = col0_traits, width = 22, height = 22,
#        units = "cm",device = "png", dpi = 300)
# 
# ggsave("Results/Effects/LP26.png", plot = lp26_traits, width = 22, height = 22,
#        units = "cm",device = "png", dpi = 300)
# 
# ggsave("Results/Effects/An1.png", plot = an1_traits, width = 22, height = 22,
#        units = "cm",device = "png", dpi = 300)
