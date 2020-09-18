
# Load packages and data --------------------------------------------------
library(tidyverse)
library(multidplyr)

photosynthesis = readRDS(file = "Intermediate/Photosynthesis/DataForCorrelations.rds")
morphology     = readRDS(file = "Intermediate/Morphology/DataForCorrelations.rds")
fitness        = readRDS(file = "Intermediate/Fitness/DataForCorrelations.rds")
medSize        = readRDS(file = "Intermediate/Biosorter/DataForCorrelations.rds")
rates        = readRDS(file = "Intermediate/Rates/DataForCorrelations.rds")

cluster = new_cluster(8)

# Combine all the datasets ------------------------------------------------

corr = full_join(photosynthesis, morphology) %>% 
        full_join(fitness)  %>% full_join(medSize) %>% 
  full_join(rates)

# Calculate correlation coefficients ---------------------------------------

traits = names(corr)[-(1:3)]
all_combs = apply(combn(traits, 2), 2, paste, collapse='_')

trait1 = map_chr(strsplit(all_combs, "_"), ~.x[1])
trait2 = map_chr(strsplit(all_combs, "_"), ~.x[2])
cat(paste0(all_combs, " = cor(", trait1, ", ", trait2,', use = "na.or.complete")'), sep = ",\n")

# Function to calculate correlations among every pair of traits
calc_corrs = function(data) {
  summarise(data, 
            Rd_LUE = cor(Rd, LUE, use = "na.or.complete"),
            Rd_A = cor(Rd, A, use = "na.or.complete"),
            Rd_Chl = cor(Rd, Chl, use = "na.or.complete"),
            Rd_gsw = cor(Rd, gsw, use = "na.or.complete"),
            Rd_DW = cor(Rd, DW, use = "na.or.complete"),
            Rd_RootLength = cor(Rd, RootLength, use = "na.or.complete"),
            Rd_Angle = cor(Rd, Angle, use = "na.or.complete"),
            Rd_PetioleRatio = cor(Rd, PetioleRatio, use = "na.or.complete"),
            Rd_ShapeBlade = cor(Rd, ShapeBlade, use = "na.or.complete"),
            Rd_SLA = cor(Rd, SLA, use = "na.or.complete"),
            Rd_LeafSize = cor(Rd, LeafSize, use = "na.or.complete"),
            Rd_N = cor(Rd, N, use = "na.or.complete"),
            Rd_Nmol = cor(Rd, Nmol, use = "na.or.complete"),
            Rd_FT = cor(Rd, FT, use = "na.or.complete"),
            Rd_Yield = cor(Rd, Yield, use = "na.or.complete"),
            Rd_Phyllochron = cor(Rd, Phyllochron, use = "na.or.complete"),
            Rd_SeedLength = cor(Rd, SeedLength, use = "na.or.complete"),
            Rd_RGR = cor(Rd, RGR, use = "na.or.complete"),
            Rd_LeafRate = cor(Rd, LeafRate, use = "na.or.complete"),
            LUE_A = cor(LUE, A, use = "na.or.complete"),
            LUE_Chl = cor(LUE, Chl, use = "na.or.complete"),
            LUE_gsw = cor(LUE, gsw, use = "na.or.complete"),
            LUE_DW = cor(LUE, DW, use = "na.or.complete"),
            LUE_RootLength = cor(LUE, RootLength, use = "na.or.complete"),
            LUE_Angle = cor(LUE, Angle, use = "na.or.complete"),
            LUE_PetioleRatio = cor(LUE, PetioleRatio, use = "na.or.complete"),
            LUE_ShapeBlade = cor(LUE, ShapeBlade, use = "na.or.complete"),
            LUE_SLA = cor(LUE, SLA, use = "na.or.complete"),
            LUE_LeafSize = cor(LUE, LeafSize, use = "na.or.complete"),
            LUE_N = cor(LUE, N, use = "na.or.complete"),
            LUE_Nmol = cor(LUE, Nmol, use = "na.or.complete"),
            LUE_FT = cor(LUE, FT, use = "na.or.complete"),
            LUE_Yield = cor(LUE, Yield, use = "na.or.complete"),
            LUE_Phyllochron = cor(LUE, Phyllochron, use = "na.or.complete"),
            LUE_SeedLength = cor(LUE, SeedLength, use = "na.or.complete"),
            LUE_RGR = cor(LUE, RGR, use = "na.or.complete"),
            LUE_LeafRate = cor(LUE, LeafRate, use = "na.or.complete"),
            A_Chl = cor(A, Chl, use = "na.or.complete"),
            A_gsw = cor(A, gsw, use = "na.or.complete"),
            A_DW = cor(A, DW, use = "na.or.complete"),
            A_RootLength = cor(A, RootLength, use = "na.or.complete"),
            A_Angle = cor(A, Angle, use = "na.or.complete"),
            A_PetioleRatio = cor(A, PetioleRatio, use = "na.or.complete"),
            A_ShapeBlade = cor(A, ShapeBlade, use = "na.or.complete"),
            A_SLA = cor(A, SLA, use = "na.or.complete"),
            A_LeafSize = cor(A, LeafSize, use = "na.or.complete"),
            A_N = cor(A, N, use = "na.or.complete"),
            A_Nmol = cor(A, Nmol, use = "na.or.complete"),
            A_FT = cor(A, FT, use = "na.or.complete"),
            A_Yield = cor(A, Yield, use = "na.or.complete"),
            A_Phyllochron = cor(A, Phyllochron, use = "na.or.complete"),
            A_SeedLength = cor(A, SeedLength, use = "na.or.complete"),
            A_RGR = cor(A, RGR, use = "na.or.complete"),
            A_LeafRate = cor(A, LeafRate, use = "na.or.complete"),
            Chl_gsw = cor(Chl, gsw, use = "na.or.complete"),
            Chl_DW = cor(Chl, DW, use = "na.or.complete"),
            Chl_RootLength = cor(Chl, RootLength, use = "na.or.complete"),
            Chl_Angle = cor(Chl, Angle, use = "na.or.complete"),
            Chl_PetioleRatio = cor(Chl, PetioleRatio, use = "na.or.complete"),
            Chl_ShapeBlade = cor(Chl, ShapeBlade, use = "na.or.complete"),
            Chl_SLA = cor(Chl, SLA, use = "na.or.complete"),
            Chl_LeafSize = cor(Chl, LeafSize, use = "na.or.complete"),
            Chl_N = cor(Chl, N, use = "na.or.complete"),
            Chl_Nmol = cor(Chl, Nmol, use = "na.or.complete"),
            Chl_FT = cor(Chl, FT, use = "na.or.complete"),
            Chl_Yield = cor(Chl, Yield, use = "na.or.complete"),
            Chl_Phyllochron = cor(Chl, Phyllochron, use = "na.or.complete"),
            Chl_SeedLength = cor(Chl, SeedLength, use = "na.or.complete"),
            Chl_RGR = cor(Chl, RGR, use = "na.or.complete"),
            Chl_LeafRate = cor(Chl, LeafRate, use = "na.or.complete"),
            gsw_DW = cor(gsw, DW, use = "na.or.complete"),
            gsw_RootLength = cor(gsw, RootLength, use = "na.or.complete"),
            gsw_Angle = cor(gsw, Angle, use = "na.or.complete"),
            gsw_PetioleRatio = cor(gsw, PetioleRatio, use = "na.or.complete"),
            gsw_ShapeBlade = cor(gsw, ShapeBlade, use = "na.or.complete"),
            gsw_SLA = cor(gsw, SLA, use = "na.or.complete"),
            gsw_LeafSize = cor(gsw, LeafSize, use = "na.or.complete"),
            gsw_N = cor(gsw, N, use = "na.or.complete"),
            gsw_Nmol = cor(gsw, Nmol, use = "na.or.complete"),
            gsw_FT = cor(gsw, FT, use = "na.or.complete"),
            gsw_Yield = cor(gsw, Yield, use = "na.or.complete"),
            gsw_Phyllochron = cor(gsw, Phyllochron, use = "na.or.complete"),
            gsw_SeedLength = cor(gsw, SeedLength, use = "na.or.complete"),
            gsw_RGR = cor(gsw, RGR, use = "na.or.complete"),
            gsw_LeafRate = cor(gsw, LeafRate, use = "na.or.complete"),
            DW_RootLength = cor(DW, RootLength, use = "na.or.complete"),
            DW_Angle = cor(DW, Angle, use = "na.or.complete"),
            DW_PetioleRatio = cor(DW, PetioleRatio, use = "na.or.complete"),
            DW_ShapeBlade = cor(DW, ShapeBlade, use = "na.or.complete"),
            DW_SLA = cor(DW, SLA, use = "na.or.complete"),
            DW_LeafSize = cor(DW, LeafSize, use = "na.or.complete"),
            DW_N = cor(DW, N, use = "na.or.complete"),
            DW_Nmol = cor(DW, Nmol, use = "na.or.complete"),
            DW_FT = cor(DW, FT, use = "na.or.complete"),
            DW_Yield = cor(DW, Yield, use = "na.or.complete"),
            DW_Phyllochron = cor(DW, Phyllochron, use = "na.or.complete"),
            DW_SeedLength = cor(DW, SeedLength, use = "na.or.complete"),
            DW_RGR = cor(DW, RGR, use = "na.or.complete"),
            DW_LeafRate = cor(DW, LeafRate, use = "na.or.complete"),
            RootLength_Angle = cor(RootLength, Angle, use = "na.or.complete"),
            RootLength_PetioleRatio = cor(RootLength, PetioleRatio, use = "na.or.complete"),
            RootLength_ShapeBlade = cor(RootLength, ShapeBlade, use = "na.or.complete"),
            RootLength_SLA = cor(RootLength, SLA, use = "na.or.complete"),
            RootLength_LeafSize = cor(RootLength, LeafSize, use = "na.or.complete"),
            RootLength_N = cor(RootLength, N, use = "na.or.complete"),
            RootLength_Nmol = cor(RootLength, Nmol, use = "na.or.complete"),
            RootLength_FT = cor(RootLength, FT, use = "na.or.complete"),
            RootLength_Yield = cor(RootLength, Yield, use = "na.or.complete"),
            RootLength_Phyllochron = cor(RootLength, Phyllochron, use = "na.or.complete"),
            RootLength_SeedLength = cor(RootLength, SeedLength, use = "na.or.complete"),
            RootLength_RGR = cor(RootLength, RGR, use = "na.or.complete"),
            RootLength_LeafRate = cor(RootLength, LeafRate, use = "na.or.complete"),
            Angle_PetioleRatio = cor(Angle, PetioleRatio, use = "na.or.complete"),
            Angle_ShapeBlade = cor(Angle, ShapeBlade, use = "na.or.complete"),
            Angle_SLA = cor(Angle, SLA, use = "na.or.complete"),
            Angle_LeafSize = cor(Angle, LeafSize, use = "na.or.complete"),
            Angle_N = cor(Angle, N, use = "na.or.complete"),
            Angle_Nmol = cor(Angle, Nmol, use = "na.or.complete"),
            Angle_FT = cor(Angle, FT, use = "na.or.complete"),
            Angle_Yield = cor(Angle, Yield, use = "na.or.complete"),
            Angle_Phyllochron = cor(Angle, Phyllochron, use = "na.or.complete"),
            Angle_SeedLength = cor(Angle, SeedLength, use = "na.or.complete"),
            Angle_RGR = cor(Angle, RGR, use = "na.or.complete"),
            Angle_LeafRate = cor(Angle, LeafRate, use = "na.or.complete"),
            PetioleRatio_ShapeBlade = cor(PetioleRatio, ShapeBlade, use = "na.or.complete"),
            PetioleRatio_SLA = cor(PetioleRatio, SLA, use = "na.or.complete"),
            PetioleRatio_LeafSize = cor(PetioleRatio, LeafSize, use = "na.or.complete"),
            PetioleRatio_N = cor(PetioleRatio, N, use = "na.or.complete"),
            PetioleRatio_Nmol = cor(PetioleRatio, Nmol, use = "na.or.complete"),
            PetioleRatio_FT = cor(PetioleRatio, FT, use = "na.or.complete"),
            PetioleRatio_Yield = cor(PetioleRatio, Yield, use = "na.or.complete"),
            PetioleRatio_Phyllochron = cor(PetioleRatio, Phyllochron, use = "na.or.complete"),
            PetioleRatio_SeedLength = cor(PetioleRatio, SeedLength, use = "na.or.complete"),
            PetioleRatio_RGR = cor(PetioleRatio, RGR, use = "na.or.complete"),
            PetioleRatio_LeafRate = cor(PetioleRatio, LeafRate, use = "na.or.complete"),
            ShapeBlade_SLA = cor(ShapeBlade, SLA, use = "na.or.complete"),
            ShapeBlade_LeafSize = cor(ShapeBlade, LeafSize, use = "na.or.complete"),
            ShapeBlade_N = cor(ShapeBlade, N, use = "na.or.complete"),
            ShapeBlade_Nmol = cor(ShapeBlade, Nmol, use = "na.or.complete"),
            ShapeBlade_FT = cor(ShapeBlade, FT, use = "na.or.complete"),
            ShapeBlade_Yield = cor(ShapeBlade, Yield, use = "na.or.complete"),
            ShapeBlade_Phyllochron = cor(ShapeBlade, Phyllochron, use = "na.or.complete"),
            ShapeBlade_SeedLength = cor(ShapeBlade, SeedLength, use = "na.or.complete"),
            ShapeBlade_RGR = cor(ShapeBlade, RGR, use = "na.or.complete"),
            ShapeBlade_LeafRate = cor(ShapeBlade, LeafRate, use = "na.or.complete"),
            SLA_LeafSize = cor(SLA, LeafSize, use = "na.or.complete"),
            SLA_N = cor(SLA, N, use = "na.or.complete"),
            SLA_Nmol = cor(SLA, Nmol, use = "na.or.complete"),
            SLA_FT = cor(SLA, FT, use = "na.or.complete"),
            SLA_Yield = cor(SLA, Yield, use = "na.or.complete"),
            SLA_Phyllochron = cor(SLA, Phyllochron, use = "na.or.complete"),
            SLA_SeedLength = cor(SLA, SeedLength, use = "na.or.complete"),
            SLA_RGR = cor(SLA, RGR, use = "na.or.complete"),
            SLA_LeafRate = cor(SLA, LeafRate, use = "na.or.complete"),
            LeafSize_N = cor(LeafSize, N, use = "na.or.complete"),
            LeafSize_Nmol = cor(LeafSize, Nmol, use = "na.or.complete"),
            LeafSize_FT = cor(LeafSize, FT, use = "na.or.complete"),
            LeafSize_Yield = cor(LeafSize, Yield, use = "na.or.complete"),
            LeafSize_Phyllochron = cor(LeafSize, Phyllochron, use = "na.or.complete"),
            LeafSize_SeedLength = cor(LeafSize, SeedLength, use = "na.or.complete"),
            LeafSize_RGR = cor(LeafSize, RGR, use = "na.or.complete"),
            LeafSize_LeafRate = cor(LeafSize, LeafRate, use = "na.or.complete"),
            N_Nmol = cor(N, Nmol, use = "na.or.complete"),
            N_FT = cor(N, FT, use = "na.or.complete"),
            N_Yield = cor(N, Yield, use = "na.or.complete"),
            N_Phyllochron = cor(N, Phyllochron, use = "na.or.complete"),
            N_SeedLength = cor(N, SeedLength, use = "na.or.complete"),
            N_RGR = cor(N, RGR, use = "na.or.complete"),
            N_LeafRate = cor(N, LeafRate, use = "na.or.complete"),
            Nmol_FT = cor(Nmol, FT, use = "na.or.complete"),
            Nmol_Yield = cor(Nmol, Yield, use = "na.or.complete"),
            Nmol_Phyllochron = cor(Nmol, Phyllochron, use = "na.or.complete"),
            Nmol_SeedLength = cor(Nmol, SeedLength, use = "na.or.complete"),
            Nmol_RGR = cor(Nmol, RGR, use = "na.or.complete"),
            Nmol_LeafRate = cor(Nmol, LeafRate, use = "na.or.complete"),
            FT_Yield = cor(FT, Yield, use = "na.or.complete"),
            FT_Phyllochron = cor(FT, Phyllochron, use = "na.or.complete"),
            FT_SeedLength = cor(FT, SeedLength, use = "na.or.complete"),
            FT_RGR = cor(FT, RGR, use = "na.or.complete"),
            FT_LeafRate = cor(FT, LeafRate, use = "na.or.complete"),
            Yield_Phyllochron = cor(Yield, Phyllochron, use = "na.or.complete"),
            Yield_SeedLength = cor(Yield, SeedLength, use = "na.or.complete"),
            Yield_RGR = cor(Yield, RGR, use = "na.or.complete"),
            Yield_LeafRate = cor(Yield, LeafRate, use = "na.or.complete"),
            Phyllochron_SeedLength = cor(Phyllochron, SeedLength, use = "na.or.complete"),
            Phyllochron_RGR = cor(Phyllochron, RGR, use = "na.or.complete"),
            Phyllochron_LeafRate = cor(Phyllochron, LeafRate, use = "na.or.complete"),
            SeedLength_RGR = cor(SeedLength, RGR, use = "na.or.complete"),
            SeedLength_LeafRate = cor(SeedLength, LeafRate, use = "na.or.complete"),
            RGR_LeafRate = cor(RGR, LeafRate, use = "na.or.complete"))
}

cluster_copy(cluster, "calc_corrs")

# Correlation across all treatments
corr_all = corr %>% 
  group_by(Genotype, ID) %>%
  partition(cluster) %>%
  calc_corrs %>%
  collect %>%
  pivot_longer(c(-Genotype, -ID), 
               names_to = c("Trait1", "Trait2"), 
               values_to = "Correlation",names_sep = "_") %>%
  mutate(Mode = "All")

# Correlation across treatments of the HTD experiments
corr_HTD = filter(corr, Effect %in% c("C", "D", "HT", "HTD")) %>% 
  group_by(Genotype, ID) %>%
  partition(cluster) %>%
  calc_corrs %>%
  collect %>%
  pivot_longer(c(-Genotype, -ID), 
               names_to = c("Trait1", "Trait2"), 
               values_to = "Correlation",names_sep = "_") %>%
  mutate(Mode = "HTD")

# Correlation across treatments of the PSD experiments
corr_PSD = filter(corr, Effect %in% c("C", "D", "PS", "PSD")) %>% 
  group_by(Genotype, ID) %>%
  partition(cluster) %>%
  calc_corrs %>%
  collect %>%
  pivot_longer(c(-Genotype, -ID), 
               names_to = c("Trait1", "Trait2"), 
               values_to = "Correlation",names_sep = "_") %>%
  mutate(Mode = "PSD")


corr = rbind(corr_all, corr_HTD, corr_PSD)

# Save results ------------------------------------------------------------

saveRDS(corr, file = "Intermediate/Correlations.rds")

