#' ---
#' title: Producer for graph overviewing a single prior set 
#' input: model_building/prior_transitivity.csv (column set_1)
#' output: output/prior_set_density_plot.pdf
#' ---
#'
#' Plots the prior density of each model parameter in one panel per
#' parameter (Figure A7 in the paper). References with
#' no trials in the design (l -> l, v -> v) are grayed out.
#' The standard references rho_<f>to<g> are truncated in the model: they
#' must lie below the lowest standard of their modality. Their priors are
#' drawn as truncated normal densities, with the bound marked by a dashed
#' line.

# Setup---------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggh4x)

output_path <- "output/"

# Prior parameters-----------------------------------------------------------
# Named list: <parameter>_logmu/_logsigma (lognormal) or _mu/_sigma (normal)
prior_set <- read.table("model_building/prior_transitivity.csv", header = TRUE, sep = ";")
priors <- setNames(as.list(prior_set$set_1), prior_set$parname)

# Upper bounds of the standard references in dB (lowest standard of each
# modality), the same as in the prior predictive check, SBC and sensitivity
# analysis. In the fit to the real data the bounds are the lowest standards
# in each participant's data.
truncation_bounds <- c(
  rho_ltol = 34, rho_ltob = 34, rho_ltos = 34,   # loudness, dB SPL
  rho_btol = 54, rho_btob = 54, rho_btos = 54,   # brightness, dB Lambert
  rho_stol = 10, rho_stob = 10, rho_stos = 10    # vibration, dB displacement
)

# Function------------------------------------------------------------------
plot_single_prior_set <- function(priors, upper = truncation_bounds) {
  
  n_points <- 2000
  plot_data <- data.frame()
  
  # Parameters to gray out (self-loop rhos that are redundant/excluded:
  # no l -> l or v -> v trials in the design)
  grayed_params <- c("rho_stos", "rho_sfroms", "rho_ltol", "rho_lfroml")
  
  # Parameter definitions: name -> distribution type and parameter keys
  # (alpha_s is not included)
  param_info <- list(
    alpha_l    = list(dist = "lognormal", mu_name = "alpha_l_logmu",  sd_name = "alpha_l_logsigma"),
    beta_b     = list(dist = "lognormal", mu_name = "beta_b_logmu",   sd_name = "beta_b_logsigma"),
    beta_l     = list(dist = "lognormal", mu_name = "beta_l_logmu",   sd_name = "beta_l_logsigma"),
    beta_s     = list(dist = "lognormal", mu_name = "beta_s_logmu",   sd_name = "beta_s_logsigma"),
    omega      = list(dist = "lognormal", mu_name = "omega_logmu",    sd_name = "omega_logsigma"),
    omega1     = list(dist = "lognormal", mu_name = "omega_1_logmu",  sd_name = "omega_1_logsigma"),
    # l rhos
    rho_ltol   = list(dist = "normal", mu_name = "rho_ltol_mu",   sd_name = "rho_ltol_sigma"),
    rho_ltob   = list(dist = "normal", mu_name = "rho_ltob_mu",   sd_name = "rho_ltob_sigma"),
    rho_ltos   = list(dist = "normal", mu_name = "rho_ltos_mu",   sd_name = "rho_ltos_sigma"),
    rho_lfroml = list(dist = "normal", mu_name = "rho_lfroml_mu", sd_name = "rho_lfroml_sigma"),
    rho_lfromb = list(dist = "normal", mu_name = "rho_lfromb_mu", sd_name = "rho_lfromb_sigma"),
    rho_lfroms = list(dist = "normal", mu_name = "rho_lfroms_mu", sd_name = "rho_lfroms_sigma"),
    # b rhos
    rho_btob   = list(dist = "normal", mu_name = "rho_btob_mu",   sd_name = "rho_btob_sigma"),
    rho_btol   = list(dist = "normal", mu_name = "rho_btol_mu",   sd_name = "rho_btol_sigma"),
    rho_btos   = list(dist = "normal", mu_name = "rho_btos_mu",   sd_name = "rho_btos_sigma"),
    rho_bfromb = list(dist = "normal", mu_name = "rho_bfromb_mu", sd_name = "rho_bfromb_sigma"),
    rho_bfroml = list(dist = "normal", mu_name = "rho_bfroml_mu", sd_name = "rho_bfroml_sigma"),
    rho_bfroms = list(dist = "normal", mu_name = "rho_bfroms_mu", sd_name = "rho_bfroms_sigma"),
    # s rhos
    rho_stos   = list(dist = "normal", mu_name = "rho_stos_mu",   sd_name = "rho_stos_sigma"),
    rho_stol   = list(dist = "normal", mu_name = "rho_stol_mu",   sd_name = "rho_stol_sigma"),
    rho_stob   = list(dist = "normal", mu_name = "rho_stob_mu",   sd_name = "rho_stob_sigma"),
    rho_sfroms = list(dist = "normal", mu_name = "rho_sfroms_mu", sd_name = "rho_sfroms_sigma"),
    rho_sfroml = list(dist = "normal", mu_name = "rho_sfroml_mu", sd_name = "rho_sfroml_sigma"),
    rho_sfromb = list(dist = "normal", mu_name = "rho_sfromb_mu", sd_name = "rho_sfromb_sigma")
  )
  
  # x ranges per parameter (references in dB)
  x_range <- function(param_name) {
    if (param_name %in% c("beta_b", "beta_l", "beta_s", "omega")) {
      seq(0, 1.1, length.out = n_points)
    } else if (param_name == "omega1") {
      seq(0, 3, length.out = n_points)
    } else if (param_name == "sig") {
      seq(0, 20, length.out = n_points)
    } else {
      seq(0, 100, length.out = n_points)
    }
  }
  
  # Density on the x grid, one row per parameter and grid point
  for (param_name in names(param_info)) {
    info   <- param_info[[param_name]]
    x_vals <- x_range(param_name)
    mu     <- priors[[info$mu_name]]
    sigma  <- priors[[info$sd_name]]
    
    density_vals <- if (info$dist == "normal") {
      dnorm(x_vals, mu, sigma)
    } else if (info$dist == "lognormal") {
      dlnorm(x_vals, mu, sigma)
    }
    
    # Truncated normal for the standard references: zero above the bound,
    # rescaled below it so the density integrates to 1
    if (param_name %in% names(upper)) {
      u <- upper[[param_name]]
      density_vals <- ifelse(x_vals <= u,
                             density_vals / pnorm(u, mu, sigma), 0)
    }
    
    plot_data <- rbind(plot_data, data.frame(
      parameter = param_name,
      value     = x_vals,
      density   = density_vals,
      grayed    = param_name %in% grayed_params
    ))
  }
  
  # Labels (plotmath; v = vibration, s in the parameter names)
  greek_labels <- c(
    "alpha_l"    = "alpha[l]",
    "beta_b"     = "beta[b]",
    "beta_l"     = "beta[l]",
    "beta_s"     = "beta[v]",
    "omega"      = "omega",
    "omega1"     = "omega[1]",
    "rho_ltol"   = "rho[l %->% l]",
    "rho_lfroml" = "rho[l %<-% l]",
    "rho_btob"   = "rho[b %->% b]",
    "rho_bfromb" = "rho[b %<-% b]",
    "rho_stos"   = "rho[v %->% v]",
    "rho_sfroms" = "rho[v %<-% v]",
    "rho_ltob"   = "rho[l %->% b]",
    "rho_lfromb" = "rho[l %<-% b]",
    "rho_btol"   = "rho[b %->% l]",
    "rho_bfroml" = "rho[b %<-% l]",
    "rho_stol"   = "rho[v %->% l]",
    "rho_sfroml" = "rho[v %<-% l]",
    "rho_ltos"   = "rho[l %->% v]",
    "rho_lfroms" = "rho[l %<-% v]",
    "rho_btos"   = "rho[b %->% v]",
    "rho_bfroms" = "rho[b %<-% v]",
    "rho_stob"   = "rho[v %->% b]",
    "rho_sfromb" = "rho[v %<-% b]"
  )
  
  plot_data$parameter_label <- factor(greek_labels[plot_data$parameter],
                                      levels = greek_labels)
  
  # Truncation bounds for the dashed lines
  bound_data <- data.frame(
    parameter_label = factor(greek_labels[names(upper)], levels = greek_labels),
    bound           = unname(upper)
  )
  
  single_color <- "#164467"
  gray_color   <- "#aaaaaa"
  
  # Build per-facet strip colors: gray for grayed params, default otherwise
  strip_colors <- ifelse(
    names(greek_labels) %in% grayed_params,
    "#cccccc",   # gray strip background
    "white"      # default strip background
  )
  names(strip_colors) <- greek_labels
  
  # Split into active and grayed data for layered geoms
  plot_data_active <- plot_data[!plot_data$grayed, ]
  plot_data_grayed <- plot_data[plot_data$grayed, ]
  
  # Build strip theme: one element_rect per label, in factor level order
  strip_fills <- strip_colors[levels(plot_data$parameter_label)]
  strip_theme <- strip_themed(
    background_x = elem_list_rect(fill = unname(strip_fills))
  )
  
  # Plot
  comparison_plot <- ggplot(plot_data, aes(x = value, y = density)) +
    # Gray overlay for excluded panels (drawn first, fills full panel area)
    geom_rect(
      data = plot_data_grayed %>% group_by(parameter_label) %>% slice(1),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "#e8e8e8", color = NA, inherit.aes = FALSE
    ) +
    # Active density lines and fills
    geom_line(
      data = plot_data_active,
      linewidth = 0.5, color = single_color
    ) +
    geom_area(
      data = plot_data_active,
      alpha = 0.3, fill = single_color
    ) +
    # Grayed density lines and fills (muted)
    geom_line(
      data = plot_data_grayed,
      linewidth = 0.5, color = gray_color
    ) +
    geom_area(
      data = plot_data_grayed,
      alpha = 0.3, fill = gray_color
    ) +
    # Truncation bounds of the standard references
    geom_vline(
      data = bound_data, aes(xintercept = bound),
      linetype = "dashed", linewidth = 0.4, color = "grey40"
    ) +
    facet_wrap2(~ parameter_label,
                scales = "free",
                ncol = 6,
                labeller = labeller(parameter_label = label_parsed),
                strip = strip_theme) +
    labs(x = "Parameter value", y = "Density") +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(angle = 0, hjust = 1, size = 14, face = "bold"),
      strip.text       = element_text(size = 16, face = "bold"),
      axis.title       = element_text(size = 16),
      axis.text        = element_text(size = 14),
      legend.position  = "none",
      panel.background = element_rect(colour = "black", linewidth = 0.5),
      panel.spacing    = unit(0.5, "lines")
    )
  
  return(comparison_plot)
}

# Usage----------------------------------------------------------------------
prior_density_plot <- plot_single_prior_set(priors)
prior_density_plot

ggsave(filename = paste0(output_path, "prior_set_density_plot.pdf"),
       plot = prior_density_plot, width = 14, height = 7)
