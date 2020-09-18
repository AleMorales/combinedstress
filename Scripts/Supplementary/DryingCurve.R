
# Load packages and data --------------------------------------------------
library(tidyverse)
library(readxl)
library(manipulate)

DryingCurve = read_excel("Data/DryingCurve.xlsx", col_names = FALSE, skip = 1) %>% as.matrix()
Timepoints = c(0,1,2,3,4,5,6,7,8,9,12)


# Prepare data for plotting -----------------------------------------------

DryingCurve = tibble(RelativeWeight = as.numeric(t(DryingCurve)),
                     Time = rep(Timepoints, times = 370),
                     ID = rep(1:370, each = 11))

FinalRW = filter(DryingCurve, Time == 12)
Control = filter(FinalRW, RelativeWeight > 0.8)$ID
Drought = filter(FinalRW, RelativeWeight < 0.3)$ID

DryingCurve = mutate(DryingCurve, Treatment = ifelse(ID %in% Control, "Control",
                                              ifelse(ID %in% Drought, "Drought", "Other")))

# Check the individual traces
ggplot(data = filter(DryingCurve, Treatment %in% c("Control", "Drought")), 
       aes(x = Time, y = RelativeWeight, group = ID, col = Treatment)) + 
  geom_line()

# Calculate median of the individual traces
sum_DryingCurve = group_by(DryingCurve, Treatment, Time) %>%
                    summarise(mu = median(RelativeWeight, na.rm = T),
                              upr = quantile(RelativeWeight, probs = 0.90, na.rm = T),
                              lwr = quantile(RelativeWeight, probs = 0.10, na.rm = T))

# Check the individual traces
pellet_weights = ggplot(data = filter(sum_DryingCurve, Treatment %in% c("Control", "Drought")), 
       aes(x = Time, y = mu, col = Treatment, linetype = Treatment, shape = Treatment)) + 
  geom_line() + 
  geom_point() +
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = Treatment), col = NA, alpha = 0.2) +
  ylim(0,1.05) +
  labs(y = "Relative pellet weight", x = "Time (days)", fill = "", col = "", shape = "", linetype = "") + 
  scale_color_manual(values = c(Control = "Blue", Drought = "Red")) + 
  scale_fill_manual(values = c(Control = "Blue", Drought = "Red")) + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.2,0.2))

ggsave(filename = "Results/Supplementary/DryingCurve.png", plot = pellet_weights,
       width = 8, height = 8, units = "cm",device = "png", dpi = 300)
