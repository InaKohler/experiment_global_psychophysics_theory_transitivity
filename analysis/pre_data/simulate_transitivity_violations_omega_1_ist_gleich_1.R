# simulate_transitivity_violations_omega_1_ist_gleich_1.R
#
# Simulates how a violation of W(1) = 1 affects transitivity (Figure 3,
# panel "W(1)=1 violations"). For the chain l -> s -> b, the indirect match
# x_lsb is compared with the direct match x_lb from the same loudness
# standard, without noise.
#
# With p = 1 the weighting function is W(1) = omega_1, so W(1) = 1 requires
# omega_1 = 1. omega_1 is varied from 0 to 10, all invariances hold:
# s-invariance (rho_l->b = rho_l->s), v-invariance (rho_b<-s = rho_b<-l)
# and t-invariance (rho_s<-l = rho_s->b).
# Points on the diagonal mean transitivity holds; the colour shows omega_1.
#
# Used by the combined script (p_w, params_base_w, loud_std_w, omega_grid).
# Output: transitivity_violations_omega_1_euqals_1.png

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
loud_std_w   <- 30   # loudness standard of the chain
sigma      <- 0        # no response noise
ntrials    <- 1        # one deterministic production per step

# Base parameters (fixed across sweep) ----------------------------------------
params_base_w <- list(
  alpha_l = 50,
  alpha_s = 5,
  alpha_b = 1,
  beta_s  = 0.65,
  beta_b  = 0.65,
  beta_l  = 0.65,
  # omega_1 is swept, see below
  # fixed rhos (not swept)
  rho_sfroml = 5, # rho_g_from_f = rho_g_to_h: t-invariance holds
  rho_stob   = 5, # rho_g_to_h
  rho_bfroms = 50, # v-invariance holds (rho_b<-s = rho_b<-l)
  rho_bfroml = 50,  # v-invariance
  rho_ltob = 15,  # s-invariance holds (rho_l->b = rho_l->s)
  rho_ltos = 15   # s-invariance
)

# Subject reference points ----------------------------------------------------
# Example parameter sets for an optional overlay (commented out in the plot)
subject_params <- data.frame(
  subj     = 1:4,
  omega_1 = c(0, 0.5, 1, 2)
)

# Sweep grid ------------------------------------------------------------------
# omega_1 from 0 to 10 in steps of 0.1
omega_grid <- expand.grid(
  omega_1 = seq(0, 10, by = .1)
)

# Simulate --------------------------------------------------------------------
# step1: l -> s from the loudness standard
# step2: s -> b with the result of step1 as standard (indirect match x_lsb)
# step3: l -> b from the same loudness standard (direct match x_lb)
# Parameter combinations that give invalid (NaN) predictions are dropped.
results <- pmap_dfr(omega_grid, function(omega_1) {
  p <- params_base_w
  p$omega_1 <- omega_1
  
  step1 <- simulate_gpm(ntrials, data.frame(task = "loud_strong",
                                            std = loud_std_w, p = 1, sigma = sigma), p)
  step2 <- simulate_gpm(ntrials, data.frame(task = "strong_bright",
                                            std = step1$tgt, p = 1, sigma = sigma), p)
  step3 <- simulate_gpm(ntrials, data.frame(task = "loud_bright",
                                            std = loud_std_w, p = 1, sigma = sigma), p)
  
  if (any(is.nan(c(step1$mu, step2$mu, step3$mu)))) return(NULL)
  
  data.frame(omega_1, lbs = step2$tgt, ls = step3$tgt)
})

# Compute subject outputs for overlay -----------------------------------------
# Same three steps for the example parameter sets
subject_outputs <- pmap_dfr(subject_params, function(subj, omega_1) {
  p <- params_base_w
  p$omega_1 <- omega_1
  
  step1 <- simulate_gpm(ntrials, data.frame(task = "loud_strong",
                                            std = loud_std_w, p = 1, sigma = sigma), p)
  step2 <- simulate_gpm(ntrials, data.frame(task = "strong_bright",
                                            std = step1$tgt, p = 1, sigma = sigma), p)
  step3 <- simulate_gpm(ntrials, data.frame(task = "loud_bright",
                                            std = loud_std_w, p = 1, sigma = sigma), p)
  
  data.frame(subj, omega_1, lbs = step2$tgt, ls = step3$tgt)
})

# Plot ------------------------------------------------------------------------
set_colors <- c("#A51E37", "#164467", "#76B3AD", "#B6BD86")

p_w <- ggplot(results |> arrange(omega_1), 
              aes(x = lbs, y = ls, color = omega_1)) +
  geom_point(size = 2, alpha = 1) +
  geom_abline(slope = 1, intercept = 0, lty = 2, linewidth = 0.6) +
  # geom_point(data = subject_outputs,
  #            aes(x = lbs, y = ls, fill = factor(subj)),
  #            inherit.aes = FALSE,
  #            shape = 21, size = 4, color = "black", stroke = .5) +
  coord_cartesian(xlim = c(49, 85), ylim = c(49,85)) +
  # geom_text(data = subject_outputs,
  #           aes(x = lbs, y = ls, label = subj),
  #           inherit.aes = FALSE,
  #           size = 3, vjust = -1) +
  # Colour scale centred on omega_1 = 1 (purple): blue below 0.5,
  # red above 1.5
  scale_color_gradientn(
    colors = c("#164467", "#164467", "#9B7EA6", "#A51E37", "#A51E37"),
    values  = scales::rescale(c(0, 0.5, 1, 1.5, 10), to = c(0, 1)),
    name    = TeX(r"($\omega_1$)")
  )+
  scale_fill_manual(values = set_colors, name = "Subject") +
  labs(
    x = TeX(r"($x_{lsb;111}$)"),
    y = TeX(r"($x_{lb;11}$)"))+
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5)
  )+
  ggtitle("W(1)=1 violations")
p_w

ggsave(plot = p_w, filename = "transitivity_violations_omega_1_euqals_1.png", 
       width = 8, height = 6, 
       dpi = 300)
