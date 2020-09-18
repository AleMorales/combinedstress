
# Load packages and data --------------------------------------------------
library(tidyverse)
library(ggcorrplot)
library(glue)
library(ggtext)

corr = readRDS(file = "Intermediate/Correlations.rds")

# Calculate summary statistics --------------------------------------------

calc_pval = function(val) {
  val = val[!is.na(val)]
  if(length(val) == 0) return(NA)
  pval_pos = sum(val > 0)/length(val)
  pval_neg = sum(val < 0)/length(val)
  min(pval_pos, pval_neg)
}

# Calculate summary statistics
sum_corr = corr %>% 
            group_by(Genotype, Mode, Trait1, Trait2) %>%
            summarise(mu = median(Correlation, na.rm = TRUE),
                      upr = quantile(Correlation, 0.95, na.rm = TRUE),
                      lwr = quantile(Correlation, 0.05, na.rm = TRUE),
                      pval = calc_pval(Correlation))

# Create a matrix with the correlations for plotting
get_cor_mat = function(data, genotype = "Average", mode = "All") {
  temp = filter(data, Genotype == genotype, Mode == mode)
  traits = c("DW", "Yield", "Rd", "LUE", "A", "gsw", "Chl", "N", "Nmol",
             "SLA", "LeafSize", "Angle", "PetioleRatio", "ShapeBlade", 
             "Phyllochron", "FT", "RootLength","SeedLength", "RGR", "LeafRate")
  nt = length(traits)
  cor_mat = matrix(NA, nt, nt)
  p_mat   = matrix(NA, nt, nt)
  colnames(p_mat) = colnames(cor_mat) = rownames(p_mat) = rownames(cor_mat) = traits
  for(t1 in traits) {
    for(t2 in traits) {
      if(t1 == t2) {
        cor_mat[t1,t2] = 1
        p_mat[t1,t2]   = 0
      } else {
        val = filter(temp, Trait1 == t1, Trait2 == t2)$mu
        p   = filter(temp, Trait1 == t1, Trait2 == t2)$pval
        if(length(val) == 0) {
          cor_mat[t1,t2] = filter(temp, Trait1 == t2, Trait2 == t1)$mu
          p_mat[t1,t2]   = filter(temp, Trait1 == t2, Trait2 == t1)$pval
        } else {
          cor_mat[t1,t2] = val
          p_mat[t1,t2] = p
        }
      }
    }
  }
  list(cor_mat, p_mat)
}

# The three matrices for the average of genotypes
mat_avg_all = get_cor_mat(sum_corr, genotype = "Average", mode = "All")
cor_mat_avg_all = mat_avg_all[[1]]
p_mat_avg_all   = mat_avg_all[[2]]
mat_avg_htd = get_cor_mat(sum_corr, genotype = "Average", mode = "HTD")
cor_mat_avg_htd = mat_avg_htd[[1]]
p_mat_avg_htd   = mat_avg_htd[[2]]
mat_avg_psd = get_cor_mat(sum_corr, genotype = "Average", mode = "PSD")
cor_mat_avg_psd = mat_avg_psd[[1]]
p_mat_avg_psd   = mat_avg_psd[[2]]

# Calculate the names and colors of each trait

category = c(Angle = "Morphology", DW = "Performance", LeafSize = "Morphology",
             N = "Physiology", Nmol = "Physiology", PetioleRatio = "Morphology", 
             RootLength = "Morphology", ShapeBlade = "Morphology", SLA = "Morphology", 
             A = "Physiology", Chl = "Physiology", gsw = "Physiology", LUE = "Physiology",
             Rd = "Physiology", FT = "Development", Phyllochron = "Development", Yield = "Performance",
             SeedLength = "Morphology", RGR = "Performance", LeafRate = "Development",
             LN = "Development")

category_color = c(Development = "#F8766D", Morphology = "#7CAE00",
                   Performance = "#00BFC4", Physiology = "#C77CFF")

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
           SeedLength = "SL", 
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
               SeedLength = "", 
               RGR = "",
               LeafRate = "",
               LN = "")

trait_names = unique(c(colnames(cor_mat_avg_all), 
                       rownames(cor_mat_avg_all)))
names_traits = tibble(Trait = trait_names) %>%
                mutate(Category = category[Trait],
                       Color    = category_color[Category],
                       abb      = abbrev[Trait],
                       sub      = abbrev_sub[Trait],
                       md_name = glue("<i style='color:{Color}'>{abb}<sub>{sub}</sub></i>"))

pretty_names = names_traits$md_name
names(pretty_names) = names_traits$Trait

colnames(cor_mat_avg_all) = pretty_names[colnames(cor_mat_avg_all)]
rownames(cor_mat_avg_all) = pretty_names[rownames(cor_mat_avg_all)]
colnames(cor_mat_avg_htd) = pretty_names[colnames(cor_mat_avg_htd)]
rownames(cor_mat_avg_htd) = pretty_names[rownames(cor_mat_avg_htd)]
colnames(cor_mat_avg_psd) = pretty_names[colnames(cor_mat_avg_psd)]
rownames(cor_mat_avg_psd) = pretty_names[rownames(cor_mat_avg_psd)]
colnames(p_mat_avg_all) = pretty_names[colnames(p_mat_avg_all)]
rownames(p_mat_avg_all) = pretty_names[rownames(p_mat_avg_all)]
colnames(p_mat_avg_htd) = pretty_names[colnames(p_mat_avg_htd)]
rownames(p_mat_avg_htd) = pretty_names[rownames(p_mat_avg_htd)]
colnames(p_mat_avg_psd) = pretty_names[colnames(p_mat_avg_psd)]
rownames(p_mat_avg_psd) = pretty_names[rownames(p_mat_avg_psd)]


# Generate figures --------------------------------------------------------

# Figures for average genotype and the three calculations
avg_all = ggcorrplot(round(cor_mat_avg_all, 1), method = "circle", type = "lower", 
           outline.col = "white", lab = TRUE, lab_size = 3, p.mat = p_mat_avg_all,
           insig = "blank", show.legend = FALSE, lab_col = "white",
           colors = c("red", "white", "blue"))

avg_all = avg_all + 
  #scale_x_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-1])) +
  #scale_y_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-20])) + 
  theme(axis.text.y = element_markdown(),
        axis.text.x = element_markdown())

ggsave(filename = "Results/Correlations/Average_corr_all.png", plot = avg_all, 
       width = 22, height = 22, units = "cm",device = "png", dpi = 300)


avg_htd = ggcorrplot(round(cor_mat_avg_htd, 1), method = "circle", type = "lower", 
           outline.col = "white", lab = TRUE, lab_size = 3, p.mat = p_mat_avg_htd,
           insig = "blank", show.legend = FALSE, lab_col = "white",
           colors = c("red", "white", "blue"))

avg_htd = avg_htd +  
  #scale_x_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-1])) +
  #scale_y_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-20])) + 
  theme(axis.text.y = element_markdown(),
        axis.text.x = element_markdown())

ggsave(filename = "Results/Correlations/Average_corr_htd.png", plot = avg_htd, 
       width = 22, height = 22, units = "cm",device = "png", dpi = 300)


avg_psd = ggcorrplot(round(cor_mat_avg_psd, 1), method = "circle", type = "lower", 
           outline.col = "white", lab = TRUE, lab_size = 3, p.mat = p_mat_avg_psd,
           insig = "blank", show.legend = FALSE, lab_col = "white",
           colors = c("red", "white", "blue"))

avg_psd = avg_psd +  
  #scale_x_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-1])) +
  #scale_y_discrete(labels = parse(text = unique(colnames(cor_mat_avg_all))[-20])) + 
  theme(axis.text.y = element_markdown(),
        axis.text.x = element_markdown())

ggsave(filename = "Results/Correlations/Average_corr_psd.png", plot = avg_psd, 
       width = 22, height = 22, units = "cm",device = "png", dpi = 300)



# Write the summary of the correlations to csv file -----------------------

write_csv(sum_corr, path = "Results/Correlations/AllCorrelations.csv")


# Correlation networks ----------------------------------------------------
library(igraph)
library(ggraph)
library(patchwork)

category = c(Angle = "Morphology", DW = "Performance", LeafSize = "Morphology",
             N = "Physiology", Nmol = "Physiology", PetioleRatio = "Morphology", 
             RootLength = "Morphology", ShapeBlade = "Morphology", SLA = "Morphology", 
             A = "Physiology", Chl = "Physiology", gsw = "Physiology", LUE = "Physiology",
             Rd = "Physiology", FT = "Development", Phyllochron = "Development", Yield = "Performance",
             SeedSize = "Morphology", RGR = "Performance", LeafRate = "Development",
             LN = "Development", SeedLength = "Morphology")


network_corr = filter(sum_corr, Genotype == "Average", Mode == "All", pval < 0.05)
               
# Prepare vertex metadata
vertices = unique(c(network_corr$Trait1, network_corr$Trait2))

vertex_score = map_dbl(vertices, function(x) {
  data = filter(network_corr, (Trait1 == x | Trait2 == x))
  mean(abs(data$mu))
})

vertices_metadata = data.frame(Vertex   = pretty_names[vertices], 
                               Category = category[vertices], 
                               Score    = vertex_score)

network_corr = network_corr %>%
  ungroup() %>%
  select(Trait1, Trait2, mu) %>%
  mutate(Trait1 = pretty_names[Trait1], Trait2 = pretty_names[Trait2])


network = graph_from_data_frame(network_corr, directed = FALSE, 
                                    vertices = vertices_metadata)

network_plot = ggraph(network, layout = "fr") +
  geom_edge_link(aes(edge_width = abs(mu), color = mu), alpha = 0.3, show.legend = TRUE) +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2"),
                              guide = guide_edge_colorbar(title = "")) +
  scale_edge_width(range = c(0.2,3)) +
  geom_node_point(aes(color = Category, size = exp(Score)), show.legend = (color = TRUE)) +
  scale_size(range = c(2,10)) +
  scale_color_manual(values = c(Development = "green", Morphology = "purple",
                                Performance = "darkgrey", Physiology = "gold")) +
  guides(edge_alpha = "none", edge_width = "none", size = "none", edge_color = "none") + 
  geom_node_text(aes(label = name), repel = TRUE, parse = TRUE) +
  labs(colour = "") +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box="vertical")

#network_plot

ggsave(filename = "Results/Correlations/Network.png", plot = network_plot, 
       width = 20, height = 16, units = "cm",device = "png", dpi = 300)


#




























       
# # Trait correlation network excluding submergence
# cor_htd = filter(network_corr, Mode == "HTD")
# 
# # Prepare vertex metadata
# vertices = unique(c(cor_htd$Trait1, cor_htd$Trait2))
# 
# vertex_score = map_dbl(vertices, function(x) {
#   data = filter(cor_htd, (Trait1 == x | Trait2 == x))
#   mean(abs(data$mu))
# })
# 
# vertices_metadata = data.frame(Vertex   = pretty_names[vertices], 
#                                Category = category[vertices], 
#                                Score    = vertex_score)
# 
# cor_htd = cor_htd %>%
#             ungroup() %>%
#             select(Trait1, Trait2, mu) %>%
#             mutate(Trait1 = pretty_names[Trait1], Trait2 = pretty_names[Trait2])
# 
# network_htd = graph_from_data_frame(cor_htd, directed = FALSE, 
#                                     vertices = vertices_metadata)
# 
# htd_plot = ggraph(network_htd, layout = "auto") +
#   geom_edge_link(aes(edge_width = abs(mu), color = mu), alpha = 0.3, show.legend = TRUE) +
#   guides(edge_alpha = "none", edge_width = "none") + 
#   scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2"),
#                               guide = guide_edge_colorbar(title = "Correlation")) +
#   scale_edge_width(range = c(0.2,3)) +
#   scale_alpha(range = c(0.1,0.3)) + 
#   geom_node_point(aes(color = Category, size = Score), show.legend = TRUE) +
#   scale_size(range = c(1,8)) +
#   geom_node_text(aes(label = name), repel = TRUE, parse = TRUE) +
#   ggtitle("A") +
#   theme_graph(plot_margin = margin(0, 10, 10, 10)) + 
#   theme(legend.position = "none")
# 
# # Trait correlation network excluding high temperature
# cor_psd = filter(network_corr, Mode == "PSD")
# 
# # Prepare vertex metadata
# vertices = unique(c(cor_psd$Trait1, cor_psd$Trait2))
# 
# vertex_score = map_dbl(vertices, function(x) {
#   data = filter(cor_psd, (Trait1 == x | Trait2 == x))
#   mean(abs(data$mu))
# })
# 
# vertices_metadata = data.frame(Vertex   = pretty_names[vertices], 
#                                Category = category[vertices], 
#                                Score    = vertex_score)
# 
# cor_psd = cor_psd %>%
#   ungroup() %>%
#   select(Trait1, Trait2, mu) %>%
#   mutate(Trait1 = pretty_names[Trait1], Trait2 = pretty_names[Trait2])
# 
# 
# network_psd = graph_from_data_frame(cor_psd, directed = FALSE, 
#                                     vertices = vertices_metadata)
# 
# psd_plot = ggraph(network_psd, layout = "auto") +
#   geom_edge_link(aes(edge_width = abs(mu), color = mu), alpha = 0.3, show.legend = TRUE) +
#   scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2"),
#                               guide = guide_edge_colorbar(title = "")) +
#   scale_edge_width(range = c(0.2,2)) +
#   geom_node_point(aes(color = Category, size = Score), show.legend = TRUE) +
#   scale_size(range = c(1,8)) +
#   guides(edge_alpha = "none", edge_width = "none", size = "none") + 
#   geom_node_text(aes(label = name), repel = TRUE, parse = TRUE) +
#   labs(title = "B", colour = "") +
#   theme_graph(plot_margin = margin(15, 10, 0, 10)) + 
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.box="vertical")
# 
# network_plot = htd_plot + psd_plot + plot_layout(ncol = 1, guides = "keep") 
# 
# ggsave(filename = "Results/Correlations/Network.png", plot = network_plot, 
#        width = 20, height = 27, units = "cm",device = "png", dpi = 300)
