
# Load packages and data --------------------------------------------------
library(tidyverse)
library(robustlmm)
library(furrr)

plan(multiprocess)

data = readRDS(file = "Intermediate/Rates/DataForLinearModel.rds") 

# Calculate relative growth rate of DW and rate of leaf appearance --------
rates_data = data %>%  
  mutate(Batch = paste0(Genotype, XP, Block),
         lDW = log(DW)) %>%
  arrange(Batch)

# Robust linear mixed models with random slope on lDW and LeafNumber vs time
fit_lDW = rlmerRcpp(lDW~Time*Genotype*Treatment + (0 + Time | Batch), 
                    data = rates_data, method = "DASvar",
                    rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                    rho.sigma.b = psi2propII(smoothPsi, k = 2.28))

fit_LN = rlmerRcpp(LeafNumber~Time*Genotype*Treatment + (0 + Time | Batch), 
                   data = rates_data, method = "DASvar",
                   rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                   rho.sigma.b = psi2propII(smoothPsi, k = 2.28))

# Create the function for refitting the model to a bootstrap sample
create_refit_model = function(fit) {
  data = fit@frame
  res = residuals(fit)
  ranef = ranef(fit)$Batch[[1]]
  blocks = as.numeric(table(data$Batch))
  refit_model  = function(trait) {
    if(trait == 'lDW')
      new_data = mutate(data, lDW = predict(fit, re.form = NA) + 
                          sample(res, replace = TRUE) + 
                          rep(sample(ranef, replace = TRUE), blocks)*Time)
    else
      new_data = mutate(data, LeafNumber = predict(fit, re.form = NA) + 
                          sample(res, replace = TRUE) + 
                          rep(sample(ranef, replace = TRUE), blocks)*Time)      
    new_fit = update(fit, data = new_data)
    fixef(new_fit)
  }
}

N = 5e3

# DW
DW_fun = create_refit_model(fit_lDW)
DW_coef = future_map_dfr(1:N, ~DW_fun("lDW"), .progress = T, 
                         .options = future_options(packages = "robustlmm"))

# LeafNumber
LeafNumber_fun = create_refit_model(fit_LN)
LN_coef = future_map_dfr(1:N, ~LeafNumber_fun("LeafNumber"), .progress = T, 
                         .options = future_options(packages = "robustlmm"))


# Calculate the slope of the linear model for each Genotype x Treatment
DW_coef = mutate(DW_coef, Trait = "RGR")
LN_coef = mutate(LN_coef, Trait = "LeafRate")
coefs = rbind(DW_coef, LN_coef)


# Save the results of fitting the linear models ---------------------------

saveRDS(list(DW_coef,LN_coef), file = "Intermediate/Rates/LinearFits.rds")
