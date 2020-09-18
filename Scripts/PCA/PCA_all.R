
# Load packages and data --------------------------------------------------
library(tidyverse)
library(patchwork)
library(ggrepel)

source("Scripts/PCA/PlotPCA.R")

photosynthesis = readRDS("Intermediate/Photosynthesis/DataForPCA_all.rds") %>%
  mutate(Genotype = "Average")
morphology     = readRDS("Intermediate/Morphology/DataForPCA_all.rds") %>%
  filter(Trait != "RWC")
fitness        = readRDS("Intermediate/Fitness/DataForPCA_all.rds")
rates          = readRDS("Intermediate/Rates/DataForPCA_all.rds")
seedlength     = readRDS("Intermediate/BioSorter/DataForPCA_all.rds") %>% mutate(Trait = "SeedLength")

# Combine traits and create matrix ----------------------------------------
traits = do.call("rbind", list(photosynthesis, morphology, fitness, rates, seedlength)) 

# Matrix with trait values for the average of all genotypes, as required by prcomp
avg_traits = filter(traits, Genotype == "Average") %>% 
                select(-Genotype) %>%
                pivot_wider(names_from = Trait, values_from = Value) %>%
                as.data.frame

rownames(avg_traits) = paste0(avg_traits$Treatment, avg_traits$ID)
treatments = avg_traits$Treatment
indexS = which(treatments == "S")
indexPS = which(treatments == "PS")
avg_traits = as.matrix(avg_traits[,-(1:2)])
avg_traits[indexS,c("FT", "Phyllochron", "Yield", "SeedLength")] = 
  avg_traits[indexPS,c("FT", "Phyllochron", "Yield", "SeedLength")]

bay0_traits = filter(traits, Genotype == "Bay0") %>% 
                select(-Genotype) %>%
                pivot_wider(names_from = Trait, values_from = Value) %>%
                as.data.frame

rownames(bay0_traits) = paste0(bay0_traits$Treatment, bay0_traits$ID)
treatments = bay0_traits$Treatment
indexS = which(treatments == "S")
indexPS = which(treatments == "PS")
bay0_traits = as.matrix(bay0_traits[,-(1:2)])
bay0_traits[indexS,c("FT", "Phyllochron", "Yield", "SeedLength")] = 
  bay0_traits[indexPS,c("FT", "Phyllochron", "Yield", "SeedLength")]


col0_traits = filter(traits, Genotype == "Col0") %>% 
                select(-Genotype) %>%
                pivot_wider(names_from = Trait, values_from = Value) %>%
                as.data.frame

rownames(col0_traits) = paste0(col0_traits$Treatment, col0_traits$ID)
treatments = col0_traits$Treatment
indexS = which(treatments == "S")
indexPS = which(treatments == "PS")
col0_traits = as.matrix(col0_traits[,-(1:2)])
col0_traits[indexS,c("FT", "Phyllochron", "Yield", "SeedLength")] = 
  col0_traits[indexPS,c("FT", "Phyllochron", "Yield", "SeedLength")]

an1_traits = filter(traits, Genotype == "An1") %>% 
                select(-Genotype) %>%
                pivot_wider(names_from = Trait, values_from = Value) %>%
                as.data.frame

rownames(an1_traits) = paste0(an1_traits$Treatment, an1_traits$ID)
treatments = an1_traits$Treatment
indexS = which(treatments == "S")
indexPS = which(treatments == "PS")
an1_traits = as.matrix(an1_traits[,-(1:2)])
an1_traits[indexS,c("FT", "Phyllochron", "Yield", "SeedLength")] = 
  an1_traits[indexPS,c("FT", "Phyllochron", "Yield", "SeedLength")]


lp26_traits = filter(traits, Genotype == "Lp26") %>% 
              select(-Genotype) %>%
              pivot_wider(names_from = Trait, values_from = Value) %>%
              as.data.frame

rownames(lp26_traits) = paste0(lp26_traits$Treatment, lp26_traits$ID)
treatments = lp26_traits$Treatment
indexS = which(treatments == "S")
indexPS = which(treatments == "PS")
lp26_traits = as.matrix(lp26_traits[,-(1:2)])
lp26_traits[indexS,c("FT", "Phyllochron", "Yield")] = 
  lp26_traits[indexPS,c("FT", "Phyllochron", "Yield")]


index_Photo = which(colnames(avg_traits) %in% c("Rd", "A", "Chl", "gsw", "LUE"))
index_NG = which(colnames(avg_traits) %in% c("RGR", "DW", "Yield"))



# Generate labels with colors ---------------------------------------------

# Prettify the names of the traits
pretty_names = c(Angle = "italic(Angle)", 
                 DW = "italic(DW)", 
                 LeafSize = "italic(BA)",
                 N = "italic(N)", 
                 Nmol = 'italic(N[area])', 
                 PetioleRatio = "italic(PR)", 
                 RootLength = "italic(RL)",
                 ShapeBlade = "italic(BS)", 
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
                 LeafRate = "italic(LAR)",
                 SeedLength = "italic(SL)")

trait_category = c(Angle = "Morphology", DW = "Performance", LeafSize = "Morphology",
             N = "Physiology", Nmol = "Physiology", PetioleRatio = "Morphology", 
             RootLength = "Morphology", ShapeBlade = "Morphology", SLA = "Morphology", 
             A = "Physiology", Chl = "Physiology", gsw = "Physiology", LUE = "Physiology",
             Rd = "Physiology", FT = "Development", Phyllochron = "Development", Yield = "Performance",
             SeedLength = "Morphology", RGR = "Performance", LeafRate = "Development",
             LN = "Development")



# Perform the PCA ---------------------------------------------------------


# Calculating loadings and individual values from PCA
process_pca = function(pca, scale = 0.2) {
  # Scaling factor
  lam = pca$sdev[1:3]
  sqrtn = sqrt(nrow(pca$x))
  lam = (lam*sqrtn)^scale
  
  # Fraction variance
  fvar = pca$sdev^2/sum(pca$sdev^2)
  
  # Extract the individuals
  individuals = tibble(Treatment = treatments, 
                       PC1 = pca$x[,1]/lam[1],
                       PC2 = pca$x[,2]/lam[2],
                       PC3 = pca$x[,3]/lam[3])
  
  # Calculate median centroid of the cloud of points
  sum_individuals = group_by(individuals, Treatment) %>%
    summarise(PC1 = median(PC1),
              PC2 = median(PC2),
              PC3 = median(PC3))
  
  # Extract the loadings of variables
  loadings = as_tibble(pca$rotation)
  loadings = mutate(loadings, 
                    Trait = rownames(pca$rotation),
                    PC1   = PC1*lam[1],
                    PC2   = PC2*lam[2],
                    PC3   = PC3*lam[3],
                    parse_trait = pretty_names[Trait],
                    Category = trait_category[Trait])
  loadings = select(loadings, Trait, PC1, PC2, PC3, parse_trait, Category)
  
  # Add fraction of variance to the PC labels
  PC = paste0("PC", 1:3, " (", round(fvar[1:3]*100,2),'%)')
  
  # Return all info
  list(loadings = loadings, individuals = individuals, 
       sum_individuals = sum_individuals, PC = PC)
}

my_biplot = function(pca, first = TRUE, genotype = "") {

  if(first)
      base = ggplot(data = pca$individuals, aes(x = PC1, y = PC2))
  else
    base = ggplot(data = pca$individuals, aes(x = PC2, y = PC3))
    
    plot = base +
    
    stat_ellipse(aes(group = Treatment), col = "black", alpha = 0.5) +
    
    geom_point(aes(group = Treatment), data = pca$sum_individuals, col = "black") +
    
    geom_segment(data = pca$loadings, aes(xend = 0, yend = 0, color = Category),
                 alpha  =0.5,
                 arrow = arrow(ends = "first",
                               length = unit(0.15, "inches"))) + 
    
    geom_text_repel(data = pca$loadings, aes(label = parse_trait, color = Category), 
                    show.legend = F, parse = TRUE,
                    size = 3)  +
    
    geom_text_repel(data = pca$sum_individuals, 
                    aes(label = Treatment), 
                    show.legend = F, col = "black")  +
      
    scale_color_manual(values = c(Development = "#F8766D",
                                  Morphology = "#7CAE00",
                                  Performance = "#00BFC4",
                                  Physiology = "#C77CFF")) + 
      
    geom_hline(yintercept = 0, linetype = 2, col = "lightgray") +
    
    geom_vline(xintercept = 0, linetype = 2, col = "lightgray") +
      
    ggtitle(genotype) +
      
    theme_bw() +
      
    theme(panel.grid = element_blank(),
          legend.position = "none")
    
    if(first)
      plot + labs(x = pca$PC[1], y = pca$PC[2], color = "")
    else
      plot + labs(x = pca$PC[2], y = pca$PC[3], color = "")
}


# Perform the PCA for the different subsets of the data
pca_avg  = process_pca(prcomp(avg_traits, center = TRUE, scale. = TRUE))
pca_bay0 = process_pca(prcomp(bay0_traits, center = TRUE, scale. = TRUE))
pca_col0 = process_pca(prcomp(col0_traits, center = TRUE, scale. = TRUE))
pca_an1  = process_pca(prcomp(an1_traits, center = TRUE, scale. = TRUE))
pca_lp26 = process_pca(prcomp(lp26_traits, center = TRUE, scale. = TRUE), scale = 0.15)
pca_avg_NP = process_pca(prcomp(avg_traits[,-index_Photo], center = TRUE, scale. = TRUE))
pca_avg_NG = process_pca(prcomp(avg_traits[,-index_NG], center = TRUE, scale. = TRUE))

# Generate all the biplots and combine to generate the final figures
avg_biplot_12    = my_biplot(pca_avg, T)
avg_biplot_23    = my_biplot(pca_avg, F)
avg_NP_biplot_12 = my_biplot(pca_avg_NP, T)
avg_NP_biplot_23 = my_biplot(pca_avg_NP, F)
avg_NG_biplot_12 = my_biplot(pca_avg_NG, T)
avg_NG_biplot_23 = my_biplot(pca_avg_NG, F)
bay0_biplot_12   = my_biplot(pca_bay0, T, genotype = "Bay-0")
bay0_biplot_23   = my_biplot(pca_bay0, F, genotype = "Bay-0")
col0_biplot_12   = my_biplot(pca_col0, T, genotype = "Col-0")
col0_biplot_23   = my_biplot(pca_col0, F, genotype = "Col-0")
an1_biplot_12    = my_biplot(pca_an1, T, genotype = "An-1")
an1_biplot_23    = my_biplot(pca_an1, F, genotype = "An-1")
lp26_biplot_12   = my_biplot(pca_lp26, T, genotype = "Lp2-6")
lp26_biplot_23   = my_biplot(pca_lp26, F, genotype = "Lp2-6")

# Combine the plots with different axes
avg_biplot    = avg_biplot_12    + avg_biplot_23
avg_NP_biplot = avg_NP_biplot_12 + avg_NP_biplot_23
avg_NG_biplot = avg_NG_biplot_12 + avg_NG_biplot_23
all_pca_12    = bay0_biplot_12   + col0_biplot_12 +
                an1_biplot_12    + lp26_biplot_12
all_pca_23    = bay0_biplot_23   + col0_biplot_23 +
                an1_biplot_23    + lp26_biplot_23

# Save the figures to files
ggsave(filename = "Results/PCA/avg_pca.png", plot = avg_biplot,
       device = png(units = 'cm', width = 22, height = 12, res = 300))

ggsave(filename = "Results/PCA/avg_pca_NP.png", plot = avg_NP_biplot,
       device = png(units = 'cm', width = 22, height = 12, res = 300))

ggsave(filename = "Results/PCA/avg_pca_NG.png", plot = avg_NG_biplot,
       device = png(units = 'cm', width = 22, height = 12, res = 300))

ggsave(filename = "Results/PCA/all_pca_12.png", plot = all_pca_12,
       device = png(units = 'cm', width = 22, height = 22, res = 300))

ggsave(filename = "Results/PCA/all_pca_23.png", plot = all_pca_23,
       device = png(units = 'cm', width = 22, height = 22, res = 300))

