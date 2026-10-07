# Simulated transitivity violations (Introduction)
#
# Combines the four simulations of how violating one condition affects
# transitivity in the case l -> s -> b vs l -> b (Figure 3 in the paper),
# and creates the table of the simulation parameters (Table A1).
#
# Each sourced script varies one set of parameters over a grid while all
# others are held constant, and provides the objects used here:
#   ..._s_invariance.R         p_s,   params_base_s,   loud_std_s,   rho_grid_s
#   ..._v_invariance.R         p_v,   params_base_v,   loud_std_v,   rho_grid_v
#   ..._t_invariance.R         p_r_i, params_base_r_i, loud_std_r_i, rho_grid_r_i
#   ..._omega_1_ist_gleich_1.R p_w,   params_base_w,   loud_std_w,   omega_grid
#
# Output: transitivity_violations_combined.pdf      (Figure 3)
#         appendix_simulation_parameters_table.tex  (Table A1)

library(patchwork)
library(xtable)
library(dplyr)

# Figure ----------------------------------------------------------------------
# 2 x 2 panels: s-invariance, v-invariance / t-invariance, W(1) = 1
source("pre_data/simulate_transitivity_violations_s_invariance.R")
source("pre_data/simulate_transitivity_violations_v_invariance.R")
source("pre_data/simulate_transitivity_violations_t_invariance.R")
source("pre_data/simulate_transitivity_violations_omega_1_ist_gleich_1.R")

combined <- (p_s + p_v) / (p_r_i + p_w) +
  plot_layout(guides = "keep")

ggsave("transitivity_violations_combined.pdf", combined,
       width = 14, height = 12, dpi = 300)


# Table -----------------------------------------------------------------------
# One column per simulation, one row per parameter. Swept parameters are
# shown as their range, constant ones as their value.

# helpers

# Format a number with up to 4 significant digits
fmt <- function(x) {
  if (is.numeric(x)) formatC(x, format = "g", digits = 4) else as.character(x)
}

# LaTeX symbol for each parameter
param_label <- function(name) {
  map <- c(
    alpha_l      = "$\\alpha_l$",
    alpha_s      = "$\\alpha_s$",
    alpha_b      = "$\\alpha_b$",
    beta_l       = "$\\beta_l$",
    beta_s       = "$\\beta_s$",
    beta_b       = "$\\beta_b$",
    omega_1      = "$\\omega_1$",
    rho_sfroml   = "$\\rho_{s \\leftarrow l}$",
    rho_stob     = "$\\rho_{s \\to b}$",
    rho_bfroms   = "$\\rho_{b \\leftarrow s}$",
    rho_bfroml   = "$\\rho_{b \\leftarrow l}$",
    rho_ltob     = "$\\rho_{l \\to b}$",
    rho_ltos     = "$\\rho_{l \\to s}$"
  )
  unname(map[name])
}

# Category shown in the first row of each parameter group
category_label <- function(name) {
  map <- c(
    alpha_l    = "Scale ($\\psi$)",
    alpha_s    = "",
    alpha_b    = "",
    beta_l     = "Exponent ($\\psi$)",
    beta_s     = "",
    beta_b     = "",
    omega_1    = "Weighing function (W)",
    rho_sfroml = "References ($\\rho$)",
    rho_stob   = "",
    rho_bfroms = "",
    rho_bfroml = "",
    rho_ltob   = "",
    rho_ltos   = ""
  )
  unname(map[name])
}

# "[min, max]" if the parameter is varied in the simulation grid, else NA
swept_label <- function(name, grid) {
  if (name %in% names(grid)) {
    r <- range(grid[[name]])
    paste0("[", fmt(r[1]), ", ", fmt(r[2]), "]")
  } else {
    NA_character_
  }
}

# build table

# Row order; loud_std is the loudness standard used in the simulation
param_order <- c("alpha_l", "alpha_s", "alpha_b",
                 "beta_l",  "beta_s",  "beta_b",
                 "omega_1", 
                 "rho_sfroml", "rho_stob", "rho_bfroms", "rho_bfroml",
                 "rho_ltob",  "rho_ltos",
                 "loud_std")

all_params_1 <- c(params_base_s,   list(loud_std = loud_std_s))
all_params_2 <- c(params_base_v,   list(loud_std = loud_std_v))
all_params_3 <- c(params_base_r_i, list(loud_std = loud_std_r_i))
all_params_4 <- c(params_base_w,   list(loud_std = loud_std_w))

# Range if swept, value if constant, "---" if not used in the simulation
cell_value <- function(params, grid, name) {
  sv <- swept_label(name, grid)
  if (!is.na(sv)) return(sv)
  v <- params[[name]]
  if (is.null(v) || is.na(v)) return("---")
  fmt(v)
}

tab <- data.frame(
  Category              = sapply(param_order, category_label),
  Parameter             = sapply(param_order, param_label),
  "S-invariance"        = sapply(param_order, cell_value, params = all_params_1, grid = rho_grid_s),
  "V-invariance"        = sapply(param_order, cell_value, params = all_params_2, grid = rho_grid_v),
  "Role-independence"   = sapply(param_order, cell_value, params = all_params_3, grid = rho_grid_r_i),
  "W(1)=1"              = sapply(param_order, cell_value, params = all_params_4, grid = omega_grid),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# render with xtable
print(
  xtable(tab,
         caption = 'Parameter values for the s-invariance, v-invariance, t-invariance, and W(1)\\,=\\,1 simulations illustrated in Figure~\\ref{fig:simulation_transitivity_invariances}.\\newline Parameters listed as ranges were varied systematically across the sweep grid; all remaining parameters were held constant at the values shown. Sweeping ranges are chosen to give a wide range of outputs while\\ avoiding invalid outputs due to unrealistically large reference points.',
         label   = "app:simulation_table",
         align   = c("l", "l", "l", "c", "c", "c", "c")),
  include.rownames       = FALSE,
  sanitize.text.function = identity,
  booktabs               = TRUE,
  caption.placement      = "top",
  comment                = FALSE,
  file                   = "appendix_simulation_parameters_table.tex"
)
