# simulate_transitivity_violations_t_invariance.R
#
# Simulates how a violation of t-invariance affects transitivity (Figure 3,
# panel "T-invariance violations"). For the chain l -> s -> b, the indirect
# match x_lsb is compared with the direct match x_lb from the same loudness
# standard, without noise.
#
# t-invariance for the middle modality vibration requires
# rho_s<-l = rho_s->b. Both references are varied over a grid, all other
# conditions hold: s-invariance (rho_l->b = rho_l->s), v-invariance
# (rho_b<-s = rho_b<-l) and W(1) = 1 (omega_1 = 1).
# Points on the diagonal mean transitivity holds; the colour shows
# rho_s<-l - rho_s->b.
#
# The _r_i suffix of the object names stems from role-independence; the
# simulated condition is t-invariance.
#
# Used by the combined script (p_r_i, params_base_r_i, loud_std_r_i,
# rho_grid_r_i).
# Output: transitivity_violations_t_invariance.png

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(latex2exp)

# Provides simulate_gpm()
if(!require(globalpsychophysics)){
  devtools::install_github("https://github.com/Kaanwoj/globalpsychophysics.git")
}
library(globalpsychophysics)

# Constants -------------------------------------------------------------------
loud_std_r_i   <- 30   # loudness standard of the chain
sigma      <- 0          # no response noise
ntrials    <- 1          # one deterministic production per step

# Base parameters (fixed across sweep) ----------------------------------------
params_base_r_i <- list(
  alpha_l = 50,
  alpha_s = 5,
  alpha_b = 1,
  beta_s  = 0.65,
  beta_b  = 0.65,
  beta_l  = 0.65,
  omega_1 = 1, # W(1) = 1
  # fixed rhos (not swept)
  rho_bfroms = 50, # v-invariance holds (rho_b<-s = rho_b<-l)
  rho_bfroml = 50,  # v-invariance
  rho_ltob = 15,  # s-invariance holds (rho_l->b = rho_l->s)
  rho_ltos = 15   # s-invariance
)

# Subject reference points ----------------------------------------------------
# Example parameter sets for an optional overlay (commented out in the plot)
subject_params <- data.frame(
  subj     = 1:4,
  rho_sfroml = c(5, 5, 5, 5),
  rho_stob = c(5, 10, 30, 50)
)

# Sweep grid ------------------------------------------------------------------
# Target and standard reference of vibration from 0 to 50 dB
rho_grid_r_i <- expand.grid(
  rho_sfroml = seq(0, 50, by = 1),
  rho_stob = seq(0, 50, by = 1)
)

# Simulate --------------------------------------------------------------------
# step1: l -> s from the loudness standard
# step2: s -> b with the result of step1 as standard (indirect match x_lsb)
# step3: l -> b from the same loudness standard (direct match x_lb)
# Parameter combinations that give invalid (NaN) predictions are dropped.
results <- pmap_dfr(rho_grid_r_i, function(rho_sfroml, rho_stob) {
  p <- params_base_r_i
  p$rho_sfroml <- rho_sfroml
  p$rho_stob <- rho_stob
  
  step1 <- simulate_gpm(ntrials, data.frame(task = "loud_strong",
                                            std = loud_std_r_i, p = 1, sigma = sigma), p)
  step2 <- simulate_gpm(ntrials, data.frame(task = "strong_bright",
                                            std = step1$tgt, p = 1, sigma = sigma), p)
  step3 <- simulate_gpm(ntrials, data.frame(task = "loud_bright",
                                            std = loud_std_r_i, p = 1, sigma = sigma), p)
  
  if (any(is.nan(c(step1$mu, step2$mu, step3$mu)))) return(NULL)
  
  data.frame(rho_sfroml, rho_stob, lbs = step2$tgt, ls = step3$tgt)
})

# Compute subject outputs for overlay -----------------------------------------
# Same three steps for the example parameter sets
subject_outputs <- pmap_dfr(subject_params, function(subj, rho_sfroml, rho_stob) {
  p <- params_base_r_i
  p$rho_sfroml <- rho_sfroml
  p$rho_stob <- rho_stob
  
  step1 <- simulate_gpm(ntrials, data.frame(task = "loud_strong",
                                            std = loud_std_r_i, p = 1, sigma = sigma), p)
  step2 <- simulate_gpm(ntrials, data.frame(task = "strong_bright",
                                            std = step1$tgt, p = 1, sigma = sigma), p)
  step3 <- simulate_gpm(ntrials, data.frame(task = "loud_bright",
                                            std = loud_std_r_i, p = 1, sigma = sigma), p)
  
  data.frame(subj, rho_sfroml, rho_stob, lbs = step2$tgt, ls = step3$tgt)
})

# Plot ------------------------------------------------------------------------
set_colors <- c("#A51E37", "#164467", "#76B3AD", "#B6BD86")

results <- results %>%
  mutate(rho_diff = (rho_sfroml - rho_stob))

subject_outputs <- subject_outputs %>%
  mutate(rho_diff = (rho_sfroml - rho_stob))

rho_range <- range(results$rho_diff, na.rm = TRUE)

# Square-root colour scale: stretches the colours around a difference of 0
narrow_mid_rescaler <- function(x, to = c(0, 1), from = rho_range) {
  x_transformed <- sign(x) * abs(x)^0.5
  from_transformed <- sign(from) * abs(from)^0.5
  scales::rescale(x_transformed, to = to, from = from_transformed)
}

p_r_i <- ggplot(results |> arrange(rho_diff), 
                aes(x = lbs, y = ls, color = rho_diff)) +
  geom_point(size = 2, alpha = 1) +
  geom_abline(slope = 1, intercept = 0, lty = 2, linewidth = 0.6) +
  # geom_point(data = subject_outputs,
  #            aes(x = lbs, y = ls, fill = factor(subj)),
  #            inherit.aes = FALSE,
  #            shape = 21, size = 4, color = "black", stroke = .5) +
  coord_cartesian(xlim = c(32, 64), ylim = c(32,64)) +
  # geom_text(data = subject_outputs,
  #           aes(x = lbs, y = ls, label = subj),
  #           inherit.aes = FALSE,
  #           size = 3, vjust = -1) +
  scale_color_gradientn(
    colors = c("#164467", "#9B7EA6", "#A51E37"),
    limits = rho_range,
    rescaler = narrow_mid_rescaler,
    name = TeX(r"($\rho_{s \leftarrow l} - \rho_{s \to b}$)")) +
  scale_fill_manual(values = set_colors, name = "Subject") +
  labs(
    x = TeX(r"($x_{lsb;11}$)"),
    y = TeX(r"($x_{lb;1}$)"))+
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5)
  )+
  ggtitle("T-invariance violations")
p_r_i

ggsave(plot = p_r_i, filename = "transitivity_violations_t_invariance.png", 
       width = 8, height = 6, 
       dpi = 300)
