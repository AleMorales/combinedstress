
# Load data and packages --------------------------------------------------
library(tidyverse)
library(boot)

# Load all the data
data = readRDS("Intermediate/Rates/AllTraits.rds")

# Basic modifications and selection for this analysis
HTDtime = c(0,5,9)
PSDtime = c(0,2,5,9,12)
data = data %>% 
          mutate(XPBlock = paste0(XP, Block),
                 Time = ifelse(XP == "HTD", HTDtime[Timepoint], PSDtime[Timepoint])) %>%
          select(Genotype, XP, Block, XPBlock, Treatment, 
                 Time, DW, LeafNumber, Angle, PetioleRatio, RootLength)

# Count number of distinct timepoints with each batch and remove the obvious choices
to_remove = data %>% 
        group_by(Genotype, XPBlock) %>%
        summarise(n = length(unique(Time))) %>%
        filter(n == 1) %>%
        mutate(GenXPBlock = paste0(Genotype, XPBlock))
        
data = filter(data, !(paste0(Genotype, XPBlock) %in% to_remove$GenXPBlock))

saveRDS(object = data, file = "Intermediate/Rates/DataForLinearModel.rds") 

# Calculate the average trait value for each Gen, XP, Block, TP and trait
means = data %>%
          group_by(Genotype, XP, Block, Time, Treatment) %>%
          summarise(DW = list(boot(DW, function(x,i) mean(x[i], na.rm = T), 1e3)$t[1,]),
                    LeafNumber = list(boot(LeafNumber, function(x,i) mean(x[i], na.rm = T), 1e3)$t),
                    Angle = list(boot(Angle, function(x,i) mean(x[i], na.rm = T), 1e3)$t),
                    PetioleRatio = list(boot(PetioleRatio, function(x,i) mean(x[i], na.rm = T), 1e3)$t),
                    RootLength = list(boot(RootLength, function(x,i) mean(x[i], na.rm = T), 1e3)$t))

# Associate each measurement to each initial Time
TP1 = means %>% filter(Time == 0, Treatment == "C") %>% ungroup() %>% select(-Time)
  
means = full_join(means, TP1, by = c("Genotype", "XP", "Block"),
                 suffix = c("","1")) 

# Normalize the means by the initial values (just a ratio)
means = means %>%
        mutate(DW = map2(DW, DW1, function(x,y) x/y),
               LeafNumber = map2(LeafNumber, LeafNumber1, function(x,y) x/y),
               Angle = map2(Angle, Angle1, function(x,y) x/y),
               PetioleRatio = map2(PetioleRatio, PetioleRatio1, function(x,y) x/y),
               RootLength = map2(RootLength, RootLength1, function(x,y) x/y)) %>%
        select(-DW1, -LeafNumber1, -Angle1, -PetioleRatio1, -RootLength1, -Treatment1)

# Elongate
means = means %>% ungroup() %>% pivot_longer(DW:RootLength, names_to = "Trait", values_to = "Value")

# Average across blocks
means = means %>%
        group_by(Genotype, XP, Time, Treatment, Trait) %>%
          summarise(Value = list(rowMeans(do.call("cbind", Value), na.rm = T)))

# Average across genotypes
avg_means = means %>%
              group_by(XP, Time, Treatment, Trait) %>%
              summarise(Value = list(rowMeans(do.call("cbind", Value), na.rm = T)))

# Compute summary stats for both datasets
sum_means = means %>%
              group_by(Genotype, XP, Time, Treatment, Trait) %>%
              summarise(mu = mean(Value[[1]], na.rm = T),
                        upr = quantile(Value[[1]], probs = 0.975, na.rm = T),
                        lwr = quantile(Value[[1]], probs = 0.025, na.rm = T))

sum_avg_means = avg_means %>%
                    mutate(mu  = map_dbl(Value, median, na.rm = T),
                           upr = map_dbl(Value, quantile, probs = 0.975, na.rm = T),
                           lwr = map_dbl(Value, quantile, probs = 0.025, na.rm = T))



# Save results ------------------------------------------------------------

saveRDS(object = list(sum_means = sum_means, sum_avg_means = sum_avg_means), 
        file   = "Intermediate/Rates/AvgTraits.rds") 

