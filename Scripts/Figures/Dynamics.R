# Load data and packages --------------------------------------------------
library(tidyverse)

time_series = readRDS("Intermediate/Rates/AvgTraits.rds") 

sum_means = time_series$sum_means
sum_avg_means = time_series$sum_avg_means

# Prepare for plotting ----------------------------------------------------

# Rename traits
pretty_names = c(Angle = "italic(Angle)", 
                 DW = "italic(DW)", 
                 LeafNumber = "italic(LN)",  
                 PetioleRatio = "italic(PR)", 
                 RootLength = "italic(RL)")

sum_means = mutate(sum_means, Trait = pretty_names[Trait])
sum_avg_means = mutate(sum_avg_means, Trait = pretty_names[Trait])


XP_names = c(PSD = "Submergence~experiments",
             HTD = "High~temperature~experiments")

sum_means = mutate(sum_means, XP = XP_names[XP])
sum_avg_means = mutate(sum_avg_means, XP = XP_names[XP])

# Plot time series average across genotypes -------------------------------

avg_plot = ggplot(data = sum_avg_means, 
       aes(x = Time, y = mu, col = Treatment, 
           linetype = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,lty = 2, col = "lightgray") +
  geom_vline(data = data.frame(x0 = 5, XP = "Submergence~experiments"), 
             aes(xintercept = x0), lty = 2, col = "lightgray") +
  labs(x = "Time (days from 10LS)", y = "Normalized trait value") +
  scale_color_manual(values = c(C = "black", D = "gold", HT = "red1", HTD = "sienna",
                                S = "blue", PS = "blue", PSD = "darkgreen")) +
  scale_shape_manual(values = c(C = 1, D = 2, HT = 3, HTD = 4,S = 5, PS = 6, PSD = 7)) +
  facet_grid(Trait~XP, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(plot = avg_plot, filename = "Results/Dynamics/avg_dynamics.png",
       width = 20, height = 16,
       units = "cm",device = "png", dpi = 300)

lp26_plot = ggplot(data = filter(sum_means, Genotype == "Lp2-6"), 
                   aes(x = Time, y = mu, col = Treatment, 
                       linetype = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,lty = 2, col = "lightgray") +
  geom_vline(data = data.frame(x0 = 5, XP = "Submergence~experiments"), 
             aes(xintercept = x0), lty = 2, col = "lightgray") +
  labs(x = "Time (days from 10LS)", y = "Normalized trait value") +
  scale_color_manual(values = c(C = "black", D = "gold", HT = "red1", HTD = "sienna",
                                S = "blue", PS = "blue", PSD = "darkgreen")) +
  scale_shape_manual(values = c(C = 1, D = 2, HT = 3, HTD = 4,S = 5, PS = 6, PSD = 7)) +
  facet_grid(Trait~XP, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(panel.grid = element_blank())


ggsave(plot = lp26_plot, filename = "Results/Dynamics/lp26_dynamics.png",
       width = 20, height = 16,
       units = "cm",device = "png", dpi = 300)


col0_plot = ggplot(data = filter(sum_means, Genotype == "Col-0"), 
                   aes(x = Time, y = mu, col = Treatment, 
                       linetype = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,lty = 2, col = "lightgray") +
  geom_vline(data = data.frame(x0 = 5, XP = "Submergence~experiments"), 
             aes(xintercept = x0), lty = 2, col = "lightgray") +
  labs(x = "Time (days from 10LS)", y = "Normalized trait value") +
  scale_color_manual(values = c(C = "black", D = "gold", HT = "red1", HTD = "sienna",
                                S = "blue", PS = "blue", PSD = "darkgreen")) +
  scale_shape_manual(values = c(C = 1, D = 2, HT = 3, HTD = 4,S = 5, PS = 6, PSD = 7)) +
  facet_grid(Trait~XP, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(panel.grid = element_blank())


ggsave(plot = col0_plot, filename = "Results/Dynamics/col0_dynamics.png",
       width = 20, height = 16,
       units = "cm",device = "png", dpi = 300)


bay0_plot = ggplot(data = filter(sum_means, Genotype == "Bay-0"), 
                   aes(x = Time, y = mu, col = Treatment, 
                       linetype = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,lty = 2, col = "lightgray") +
  geom_vline(data = data.frame(x0 = 5, XP = "Submergence~experiments"), 
             aes(xintercept = x0), lty = 2, col = "lightgray") +
  labs(x = "Time (days from 10LS)", y = "Normalized trait value") +
  scale_color_manual(values = c(C = "black", D = "gold", HT = "red1", HTD = "sienna",
                                S = "blue", PS = "blue", PSD = "darkgreen")) +
  scale_shape_manual(values = c(C = 1, D = 2, HT = 3, HTD = 4,S = 5, PS = 6, PSD = 7)) +
  facet_grid(Trait~XP, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(plot = bay0_plot, filename = "Results/Dynamics/bay0_dynamics.png",
       width = 20, height = 16,
       units = "cm",device = "png", dpi = 300)

an1_plot = ggplot(data = filter(sum_means, Genotype == "An-1"), 
                  aes(x = Time, y = mu, col = Treatment, 
                      linetype = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,lty = 2, col = "lightgray") +
  geom_vline(data = data.frame(x0 = 5, XP = "Submergence~experiments"), 
             aes(xintercept = x0), lty = 2, col = "lightgray") +
  labs(x = "Time (days from 10LS)", y = "Normalized trait value") +
  scale_color_manual(values = c(C = "black", D = "gold", HT = "red1", HTD = "sienna",
                                S = "blue", PS = "blue", PSD = "darkgreen")) +
  scale_shape_manual(values = c(C = 1, D = 2, HT = 3, HTD = 4,S = 5, PS = 6, PSD = 7)) +
  facet_grid(Trait~XP, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = an1_plot, filename = "Results/Dynamics/an1_dynamics.png",
       width = 20, height = 16,
       units = "cm",device = "png", dpi = 300)
