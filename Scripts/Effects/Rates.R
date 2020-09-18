
# Load packages and data --------------------------------------------------
library(tidyverse)

coefs = readRDS(file = "Intermediate/Rates/LinearFits.rds")

# Combine into single table
effects = do.call("rbind", coefs)

# Compute the different combinations to get the effects on RGR and LeafRate
effects = effects %>% group_by(Trait) %>%
  mutate(# An-1
    An1_HT   = `Time:TreatmentHT`/Time,
    An1_S    = `Time:TreatmentS`/Time,
    An1_PS   = `Time:TreatmentPS`/Time,
    An1_D    = `Time:TreatmentD`/Time,
    An1_HTD  = `Time:TreatmentHTD`/Time,
    An1_PSD  = `Time:TreatmentPSD`/Time,
    An1_HT_D = (`Time:TreatmentHTD` - `Time:TreatmentHT` - `Time:TreatmentD`)/Time,
    An1_PS_D = (`Time:TreatmentPSD` - `Time:TreatmentPS` - `Time:TreatmentD`)/Time,
    # Bay-0
    Bay0_HT  = (`Time:TreatmentHT` + `Time:GenotypeBay-0:TreatmentHT`)/(Time + `Time:GenotypeBay-0`),
    Bay0_S    = (`Time:TreatmentS` + `Time:GenotypeBay-0:TreatmentS`)/(Time + `Time:GenotypeBay-0`),
    Bay0_PS  = (`Time:TreatmentPS` + `Time:GenotypeBay-0:TreatmentPS`)/(Time + `Time:GenotypeBay-0`),
    Bay0_D   = (`Time:TreatmentD`  + `Time:GenotypeBay-0:TreatmentD`)/(Time + `Time:GenotypeBay-0`),
    Bay0_HTD   = (`Time:TreatmentHTD`  + `Time:GenotypeBay-0:TreatmentHTD`)/(Time + `Time:GenotypeBay-0`),
    Bay0_PSD   = (`Time:TreatmentPSD`  + `Time:GenotypeBay-0:TreatmentPSD`)/(Time + `Time:GenotypeBay-0`),
    Bay0_HT_D = (`Time:TreatmentHTD` + `Time:GenotypeBay-0:TreatmentHTD` -
                  `Time:TreatmentHT`  - `Time:GenotypeBay-0:TreatmentHT`  -
                  `Time:TreatmentD`   - `Time:GenotypeBay-0:TreatmentD`)/(Time + `Time:GenotypeBay-0`),
    Bay0_PS_D = (`Time:TreatmentPSD` + `Time:GenotypeBay-0:TreatmentPSD` -
                  `Time:TreatmentPS`  - `Time:GenotypeBay-0:TreatmentPS`  -
                  `Time:TreatmentD`   - `Time:GenotypeBay-0:TreatmentD`)/(Time + `Time:GenotypeBay-0`),
    # Col-0
    Col0_HT  = (`Time:TreatmentHT` + `Time:GenotypeCol-0:TreatmentHT`)/(Time + `Time:GenotypeCol-0`),
    Col0_S    = (`Time:TreatmentS` + `Time:GenotypeCol-0:TreatmentS`)/(Time + `Time:GenotypeCol-0`),
    Col0_PS  = (`Time:TreatmentPS` + `Time:GenotypeCol-0:TreatmentPS`)/(Time + `Time:GenotypeCol-0`),
    Col0_D   = (`Time:TreatmentD`  + `Time:GenotypeCol-0:TreatmentD`)/(Time + `Time:GenotypeCol-0`),
    Col0_HTD   = (`Time:TreatmentHTD`  + `Time:GenotypeCol-0:TreatmentHTD`)/(Time + `Time:GenotypeCol-0`),
    Col0_PSD   = (`Time:TreatmentPSD`  + `Time:GenotypeCol-0:TreatmentPSD`)/(Time + `Time:GenotypeCol-0`),
    Col0_HT_D = (`Time:TreatmentHTD` + `Time:GenotypeCol-0:TreatmentHTD` -
                  `Time:TreatmentHT`  - `Time:GenotypeCol-0:TreatmentHT`  -
                  `Time:TreatmentD`   - `Time:GenotypeCol-0:TreatmentD`)/(Time + `Time:GenotypeCol-0`),
    Col0_PS_D = (`Time:TreatmentPSD` + `Time:GenotypeCol-0:TreatmentPSD` -
                  `Time:TreatmentPS`  - `Time:GenotypeCol-0:TreatmentPS`  -
                  `Time:TreatmentD`   - `Time:GenotypeCol-0:TreatmentD`)/(Time + `Time:GenotypeCol-0`),
    # Lp2-6
    Lp26_HT  = (`Time:TreatmentHT` + `Time:GenotypeLp2-6:TreatmentHT`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_S    = (`Time:TreatmentS` + `Time:GenotypeLp2-6:TreatmentS`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_PS  = (`Time:TreatmentPS` + `Time:GenotypeLp2-6:TreatmentPS`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_D   = (`Time:TreatmentD`  + `Time:GenotypeLp2-6:TreatmentD`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_HTD   = (`Time:TreatmentHTD`  + `Time:GenotypeLp2-6:TreatmentHTD`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_PSD   = (`Time:TreatmentPSD`  + `Time:GenotypeLp2-6:TreatmentPSD`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_HT_D = (`Time:TreatmentHTD` + `Time:GenotypeLp2-6:TreatmentHTD` -
                  `Time:TreatmentHT`  - `Time:GenotypeLp2-6:TreatmentHT`  -
                  `Time:TreatmentD`   - `Time:GenotypeLp2-6:TreatmentD`)/(Time + `Time:GenotypeLp2-6`),
    Lp26_PS_D = (`Time:TreatmentPSD` + `Time:GenotypeLp2-6:TreatmentPSD` -
                  `Time:TreatmentPS`  - `Time:GenotypeLp2-6:TreatmentPS`  -
                  `Time:TreatmentD`   - `Time:GenotypeLp2-6:TreatmentD`)/(Time + `Time:GenotypeLp2-6`)) %>%
  dplyr::select(Trait, An1_HT, An1_S, An1_PS, An1_D, An1_HTD, An1_PSD, An1_HT_D, An1_PS_D,
                Bay0_HT, Bay0_S, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD, Bay0_HT_D, Bay0_PS_D,
                Col0_HT, Col0_S, Col0_PS, Col0_D, Col0_HTD, Col0_PSD, Col0_HT_D, Col0_PS_D,
                Lp26_HT, Lp26_S, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD, Lp26_HT_D, Lp26_PS_D)

# Create average of genotypes
effects = mutate(effects, 
                 Average_HT   = (An1_HT  + Bay0_HT  + Col0_HT  + Lp26_HT )/4,
                 Average_S    = (An1_S   + Bay0_S   + Col0_S   + Lp26_S  )/4,
                 Average_PS   = (An1_PS  + Bay0_PS  + Col0_PS  + Lp26_PS )/4,
                 Average_D    = (An1_D   + Bay0_D   + Col0_D   + Lp26_D  )/4,
                 Average_HTD  = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
                 Average_PSD  = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4,
                 Average_HT_D = (An1_HT_D + Bay0_HT_D + Col0_HT_D + Lp26_HT_D)/4,
                 Average_PS_D = (An1_PS_D + Bay0_PS_D + Col0_PS_D + Lp26_PS_D)/4)

# Split dataset across genotypes, add genotype, rename and rbind it
effects = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(effects, Trait, contains(x))
  names(out) = c("Trait", "HT", "S", "PS", "D", "HTD", "PSD", "HT_D", "PS_D")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to long format
effects = pivot_longer(effects, c(-Trait, -Genotype), 
                       names_to = "Effect", values_to = "Value")

effects = mutate(effects, mu = NA)


# Save results -------------------------------------------------------------

saveRDS(object = effects, file   = "Intermediate/Rates/Effects.rds")
