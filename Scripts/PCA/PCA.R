
# Load packages and data --------------------------------------------------
library(tidyverse)
library(patchwork)

source("Scripts/PCA/PlotPCA.R")

photosynthesis = readRDS("Intermediate/Photosynthesis/DataForPCA.rds") %>%
                  mutate(Genotype = "Average")
morphology     = readRDS("Intermediate/Morphology/DataForPCA.rds") %>%
                  filter(Trait != "RWC")
fitness        = readRDS("Intermediate/Fitness/DataForPCA.rds")
rates          = readRDS("Intermediate/Rates/DataForPCA.rds")

# Combine traits and create matrix ----------------------------------------
traits = do.call("rbind", list(photosynthesis, morphology, fitness, rates)) #%>%
          #filter(Treatment != "S")

# Prettify the names of the traits
pretty_names = c(Angle = "italic(Angle)", 
                 DW = "italic(DW)", 
                 LeafSize = "italic(LS)",
                 N = "italic(N)", 
                 Nmol = 'italic(N[area])', 
                 PetioleRatio = "italic(PR)", 
                 RootLength = "italic(RL)",
                 ShapeBlade = "italic(SB)", 
                 SLA = "italic(SPA)", 
                 A = "italic(A)", 
                 Chl = "italic(Chl)",
                 gsw = "italic(g[s])", 
                 LUE = "italic(LUE)",
                 Rd = "italic(R[d])", 
                 FT = "italic(FT)",
                 Phyllochron = "italic(PHY)", 
                 Yield = "italic(Y)",
                 SeedSize = "italic(SL)", 
                 RGR = "italic(RGR)",
                 LeafRate = "italic(LAR)")

traits = mutate(traits, Trait = pretty_names[Trait])

# Matrix with trait values for the average of all genotypes, as required by prcomp
avg_traits = filter(traits, Genotype == "Average") %>% 
              select(-Genotype) %>%
              pivot_wider(names_from = Trait, values_from = mu) %>%
              as.data.frame
rownames(avg_traits) = avg_traits$Treatment
avg_traits = as.matrix(avg_traits[,-1])
avg_traits["S",c("italic(FT)", "italic(PHY)", "italic(Y)")] = 
   avg_traits["PS",c("italic(FT)", "italic(PHY)", "italic(Y)")]

bay0_traits = filter(traits, Genotype == "Bay0") %>% 
  select(-Genotype) %>%
  pivot_wider(names_from = Trait, values_from = mu) %>%
  as.data.frame
rownames(bay0_traits) = bay0_traits$Treatment
bay0_traits = as.matrix(bay0_traits[,-1])
bay0_traits["S",c("italic(FT)", "italic(PHY)", "italic(Y)")] = 
  bay0_traits["PS",c("italic(FT)", "italic(PHY)", "italic(Y)")]

col0_traits = filter(traits, Genotype == "Col0") %>% 
  select(-Genotype) %>%
  pivot_wider(names_from = Trait, values_from = mu) %>%
  as.data.frame
rownames(col0_traits) = col0_traits$Treatment
col0_traits = as.matrix(col0_traits[,-1])
col0_traits["S",c("italic(FT)", "italic(PHY)", "italic(Y)")] = 
  col0_traits["PS",c("italic(FT)", "italic(PHY)", "italic(Y)")]

an1_traits = filter(traits, Genotype == "An1") %>% 
  select(-Genotype) %>%
  pivot_wider(names_from = Trait, values_from = mu) %>%
  as.data.frame
rownames(an1_traits) = an1_traits$Treatment
an1_traits = as.matrix(an1_traits[,-1])
an1_traits["S",c("italic(FT)", "italic(PHY)", "italic(Y)")] = 
  an1_traits["PS",c("italic(FT)", "italic(PHY)", "italic(Y)")]


lp26_traits = filter(traits, Genotype == "Lp26") %>% 
  select(-Genotype) %>%
  pivot_wider(names_from = Trait, values_from = mu) %>%
  as.data.frame
rownames(lp26_traits) = lp26_traits$Treatment
lp26_traits = as.matrix(lp26_traits[,-1])
lp26_traits["S",c("italic(FT)", "italic(PHY)", "italic(Y)")] = 
  lp26_traits["PS",c("italic(FT)", "italic(PHY)", "italic(Y)")]


# Perform the PCA ---------------------------------------------------------

biplot = function(pca, treatments = rownames(avg_traits), genotype = "",
                  x = 1, y = 2) {
  my_byplot(object = pca, 
            data = data.frame(Treatment = treatments), 
            shape = 16, 
            label = TRUE, 
            label.label = "Treatment", 
            label.colour = "blue", 
            label.repel = TRUE,
            label.size = 3,
            loadings = TRUE, 
            loadings.colour = 'gray',
            loadings.label = TRUE, 
            loadings.label.size = 3, 
            loadings.label.repel = TRUE,
            loadings.label.colour = "red", 
            x = x, y = y) + 
    theme_bw() + 
    ggtitle(genotype)
}

# Average genotype (including photosynthesis traits)
pca_avg = prcomp(avg_traits, center = TRUE, scale. = TRUE)
avg_biplot_12 = biplot(pca_avg, x = 1, y = 2)
avg_biplot_23 = biplot(pca_avg, x = 3, y = 2)
avg_biplot_13 = biplot(pca_avg, x = 1, y = 3)
# avg_biplot = avg_biplot_12 + avg_biplot_23 + avg_biplot_13 + plot_layout(ncol = 2)
# avg_biplot
avg_biplot2 = avg_biplot_12 + avg_biplot_23

# Biplots 
pca_bay0 = prcomp(bay0_traits, center = TRUE, scale. = TRUE)
bay0_biplot_12 = biplot(pca_bay0, genotype = "Bay-0", x = 1, y = 2)
bay0_biplot_23 = biplot(pca_bay0, genotype = "Bay-0", x = 2, y = 3)

pca_col0 = prcomp(col0_traits, center = TRUE, scale. = TRUE)
col0_biplot_12 = biplot(pca_col0, genotype = "Col-0", x = 1, y = 2)
col0_biplot_23 = biplot(pca_col0, genotype = "Col-0", x = 2, y = 3)

pca_an1 = prcomp(an1_traits, center = TRUE, scale. = TRUE)
an1_biplot_12 = biplot(pca_an1, genotype = "An-1", x = 1, y = 2)
an1_biplot_23 = biplot(pca_an1, genotype = "An-1", x = 2, y = 3)

pca_lp26 = prcomp(lp26_traits, center = TRUE, scale. = TRUE)
lp26_biplot_12 = biplot(pca_lp26, genotype = "Lp2-6", x = 1, y = 2)
lp26_biplot_23 = biplot(pca_lp26, genotype = "Lp2-6", x = 2, y = 3)

all_pca_12 = bay0_biplot_12 + col0_biplot_12 + an1_biplot_12 + lp26_biplot_12 +
          plot_layout(ncol = 2, nrow = 2)

all_pca_23 = bay0_biplot_23 + col0_biplot_23 + an1_biplot_23 + lp26_biplot_23 +
  plot_layout(ncol = 2, nrow = 2)

# Average genotype (excluding photosynthesis traits)
index_Photo = which(colnames(avg_traits) %in% c("italic(R[d])", "italic(A)",
                                             "italic(Chl)", "italic(g[s])",
                                             "italic(LUE)"))
pca_avg_NP = prcomp(avg_traits[,-index_Photo], center = TRUE, scale. = TRUE)
avg_NP_biplot_12 = biplot(pca_avg_NP, x = 1, y = 2)
avg_NP_biplot_23 = biplot(pca_avg_NP, x = 3, y = 2)
avg_NP_biplot = avg_NP_biplot_12 + avg_NP_biplot_23


# Average genotype (excluding performance traits)
index_NG = which(colnames(avg_traits) %in% c("italic(RGR)", "italic(DW)",
                                                "italic(Y)"))
pca_avg_NG = prcomp(avg_traits[,-index_NG], center = TRUE, scale. = TRUE)
avg_NG_biplot_12 = biplot(pca_avg_NG, x = 1, y = 2)
avg_NG_biplot_23 = biplot(pca_avg_NG, x = 3, y = 2)
avg_NG_biplot = avg_NG_biplot_12 + avg_NG_biplot_23
avg_NG_biplot

# Save the figures to files
ggsave(filename = "Results/PCA/avg_pca.png", plot = avg_biplot2,
       device = png(units = 'cm', width = 22, height = 12, res = 300))

ggsave(filename = "Results/PCA/avg_pca_NP.png", plot = avg_NP_biplot,
       device = png(units = 'cm', width = 22, height = 12, res = 300))

ggsave(filename = "Results/PCA/avg_pca_NG.png", plot = avg_NG_biplot,
       device = png(units = 'cm', width = 22, height = 12, res = 300))


ggsave(filename = "Results/PCA/all_pca_12.png", plot = all_pca_12,
       device = png(units = 'cm', width = 22, height = 22, res = 300))

ggsave(filename = "Results/PCA/all_pca_23.png", plot = all_pca_23,
       device = png(units = 'cm', width = 22, height = 22, res = 300))


# 3D PCA plot -------------------------------------------------------------

library(pca3d)

pca3d(pca_avg, show.labels = TRUE)
snapshotPCA3d("Results/PCA/avg_pca_3D.png")



# Using PCA to model treatments -------------------------------------------

treatments = pca_avg2$x[,1:3]

treatments = rbind(treatments[c("C","HT","S","PS"),],
                   `PSD - PS` = treatments["PSD",] - treatments["PS",],
                   `HTD - HT` = treatments["HTD",] - treatments["HT",],
                   `D - C` = treatments["D",] - treatments["C",])

treatments = tibble(Treatment = rep(rownames(treatments), 3),
                    Component = as.factor(rep(1:3, each = 7)),
                    Loading = as.numeric(treatments))

ggplot(data = treatments, aes(x = Treatment, y = Loading)) +
  geom_col(position = "dodge") + facet_wrap(~Component)
  geom_point(position = position_dodge(width = 0.5))

# Associating traits to treatments ----------------------------------------

# Describe each trait/index and treatment as vectors
pca_avg2 = prcomp(avg_traits, center = TRUE, scale. = TRUE)
traits = pca_avg2$rotation
treatments = pca_avg2$x

# Calculate cos(angle) with respect to PC
cos_angle = function(v1, v2) {
  prod = sum(v1*v2)
  mag  = sqrt(sum(v1*v1))*sqrt(sum(v2*v2))
  prod/mag
}

cosangles = apply(traits, 1, function(y) apply(treatments, 1, function(x) cos_angle(y, x)))

library("ggcorrplot")


# Plot correlation of each trait with each PC
performance = which(colnames(cosangles) %in% c("italic(DW)", "italic(Y)", "italic(RGR)"))
#trait_corr_PCA = 
  ggcorrplot(round(cosangles, 1)[c("HT", "PS", "S"), -performance], method = "circle", type = "full", 
           outline.col = "white", lab = TRUE, lab_size = 3, p.mat = NULL,
           insig = "blank", show.legend = FALSE, lab_col = "white", 
           colors = c("red", "white", "blue")) +
  #scale_x_discrete(labels = parse(text = rownames(cosangles))) +
  scale_y_discrete(labels = parse(text = colnames(cosangles)[-performance])) #+
  #geom_vline(xintercept = c(2.5, 4.5, 6.5), linetype = 2) 

trait_corr_PCA

ggsave(filename = "Results/PCA/trait_corr_PC.png", plot = trait_corr_PCA, 
       width = 22, height = 22, units = "cm",device = "png", dpi = 300)
