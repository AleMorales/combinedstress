# Load data and packages --------------------------------------------------
library(tidyverse)
library(readxl)
library(boot)

data = read_excel("Data/ColFrigida.xlsx", sheet = "Raw Data", range = "A1:J121")

# Prepare data for plotting -----------------------------------------------

# Assign relative time to each timepoint
Time = c(LS11 = 0, `Day 4` = 4, `Day 8` = 8)
data = mutate(data, Time = Time[Timepoint])

# Rename treatments for consistency
Treatments = c(C = "C", D = "D", HT = "HT", `HT&D` = "HTD")
data = mutate(data, Treatment = Treatments[Treatment])

# Duplicate initial timepoint for D and HTD
C1  = filter(data, Time == 0, Treatment == "C") %>% mutate(Treatment = "D")
HT1 = filter(data, Time == 0, Treatment == "HT") %>% mutate(Treatment = "HTD")
data = rbind(data, C1, HT1)

# Adjust units to the rest of the dataset
data = mutate(data, DW = DW*1e3, SLA = SLA*10)

# Remove big outlier in angle data
data$Angle[which(data$Angle > 60)] = NA

# Switch to longer format
data = data %>%
        dplyr::select(Genotype, Treatment, Time, PetioleLength, Angle, DW, SLA) %>%
        pivot_longer(cols = PetioleLength:SLA, names_to = "Trait", values_to = "Value")

# Prettier names for the traits
pretty_names = c(Angle = "italic(Angle)~(degree)", 
                 DW = "italic(DW)~(mg)", 
                 PetioleLength = "Petiole~length~(mm)",
                 SLA = "italic(SPA)~(cm^{2}~mg^{-1})")
data = mutate(data, Trait = pretty_names[Trait],
                    Genotype = c(`Col-0` = "Col-0", 
                              `Col-FRIGIDA` = "Col-italic(FRIGIDA)")[Genotype])

# Bootstrap the mean values
boot_data = data %>%
              group_by(Genotype, Treatment, Trait, Time) %>%
              summarise(
                  Value = list(as.numeric(boot(Value, function(x,i) mean(x[i], na.rm = T), 5e3)$t))
              )

# Calculate summary statistics from the distributions
sum_data = boot_data %>%
              mutate(mu = map_dbl(Value, median),
                     upr = map_dbl(Value, quantile, prob = 0.975),
                     lwr = map_dbl(Value, quantile, prob = 0.025)) %>%
          select(-Value)

sum_data_final = filter(sum_data, Time == 8)

# Visualize the data ------------------------------------------------------
frigida_plot = ggplot(data = sum_data_final, aes(x = Treatment, y = mu, 
                      color = Genotype, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.3, width = 0.8)  + 
  geom_errorbar(aes(ymax = upr, ymin = lwr), position = position_dodge(width = 0.8), width = 0.5) + 
  scale_color_discrete(labels = scales::label_parse()) +
  scale_fill_discrete(labels = scales::label_parse()) +
  labs(y = "Average trait value") +
  facet_wrap(~Trait, scales = "free",labeller = "label_parsed") +
  theme_bw() +
  theme(#legend.position = "none",
        panel.grid = element_blank())

ggsave(filename = "Results/Supplementary/FrigidaComparison.png", plot = frigida_plot,
       width = 16, height = 16, units = "cm",device = "png", dpi = 300)

