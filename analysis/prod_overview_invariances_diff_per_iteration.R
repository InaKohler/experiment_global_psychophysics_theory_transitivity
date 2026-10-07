# prod_overview_invariances_diff_per_iteration.R
#
# Overview of all participants and chains f -> g -> h (Figure 5 in the
# paper). Per panel:
#   x_fh - x_fgh : mean difference between direct and indirect match with
#                  its 95% t-interval
#   s-, v-, t-invariance : 95% posterior interval of the difference between
#                  the two references compared
#   omega[1] - 1 : 95% posterior interval of omega_1 - 1 (W(1) = 1)
# A point is blue when its interval includes 0 (transitivity / condition
# holds) and red otherwise. The panel border is red when the difference and
# the conditions disagree (see "Panel border colors").
#
# Input:  data via get_analysis_data()
#         output/rho_samples_<subj>.csv (rho_* and omega_1 draws)
# Output: output/invariance_diff_overview_with_w_1.pdf

# Setup ---------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)

source("get_analysis_data.R")

uni_red <- "#951A36"
modal_blue <- "#164467"
modal_tuerkis <- "#76B3AD"
modal_gruen <- "#B6BD86"

jnd_lookup <- c(visual = 0.7, auditory = 1.1, tactile = 2.2)

arrow_sep <- "->"
task_levels <- c(
  paste("auditory", arrow_sep, "visual"),  
  paste("auditory", arrow_sep, "tactile"),
  paste("visual", arrow_sep, "auditory"),  
  paste("visual", arrow_sep, "tactile"),
  paste("tactile", arrow_sep, "auditory"), 
  paste("tactile", arrow_sep, "visual")
)

# Strip labels: chain f -> g -> h for each task f -> h,
# l = loudness, b = brightness, v = vibration
task_strip_labels <- c(
  "auditory -> visual"  = "l \u2192 v \u2192 b",
  "auditory -> tactile" = "l \u2192 b \u2192 v",
  "visual -> auditory"  = "b \u2192 v \u2192 l",
  "visual -> tactile"   = "b \u2192 l \u2192 v",
  "tactile -> auditory" = "v \u2192 b \u2192 l",
  "tactile -> visual"   = "v \u2192 l \u2192 b"
)

modalities <- c("auditory", "visual", "tactile")
# Letters used in the rho_ parameter names (s = tactile/strong)
modality_letter <- c(auditory = "l", visual = "b", tactile = "s")

if (!exists("dat")) {
  dat <- get_analysis_data()
}

# Direct (basic) and indirect (successive) match of each chain f -> g -> h,
# see get_transitivity_pairs() in get_analysis_data.R
dat_wide <- get_transitivity_pairs(dat) %>%
  mutate(JND = jnd_lookup[target_modality])

# Task f -> h (start and end of the chain) --------------------------------
dat_wide <- dat_wide %>%
  mutate(task = factor(paste(standard_modality, arrow_sep, target_modality),
                       levels = task_levels))

# Mean difference x_fh - x_fgh with 95% t-interval per subject and task ----
diff_summary <- dat_wide %>%
  group_by(subj_id, task) %>%
  summarise(
    n         = n(),
    mean_diff = mean(diff),
    sd_diff   = sd(diff),
    se_diff   = sd_diff / sqrt(n),
    ci_lower  = mean_diff - qt(0.975, df = n - 1) * se_diff,
    ci_upper  = mean_diff + qt(0.975, df = n - 1) * se_diff,
    .groups = "drop"
  ) %>%
  mutate(ci_includes_zero = ci_lower <= 0 & ci_upper >= 0) %>%
  arrange(subj_id, task)

# Load rho posterior samples -------------------------------------------------
# (full MCMC draws saved by prod_rds_file_parameter_estimation.R; needed to
# compute posterior distributions of pairwise differences for invariance plots)
sample_files <- list.files("output", pattern = "^rho_samples_.*\\.csv$",
                           full.names = TRUE)

rho_samples <- lapply(sample_files, read.csv, stringsAsFactors = FALSE)
names(rho_samples) <- gsub(".*rho_samples_|\\.csv$", "", sample_files)

# Returns mean + 95% CI of (param1 - param2) from joint posterior samples
get_diff_ci <- function(code, param1, param2) {
  samp <- rho_samples[[code]]
  if (is.null(samp) || !all(c(param1, param2) %in% names(samp))) {
    return(tibble(mean = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  d <- samp[[param1]] - samp[[param2]]
  tibble(mean     = mean(d),
         ci_lower = quantile(d, 0.025),
         ci_upper = quantile(d, 0.975))
}

# Invariance holds when the CI of the difference includes zero
diff_includes_zero <- function(code, param1, param2) {
  ci <- get_diff_ci(code, param1, param2)
  if (anyNA(ci)) return(NA)
  ci$ci_lower <= 0 & ci$ci_upper >= 0
}

# Invariance judgments per subj x task ----------------------------------------
# For standard F, target H, mediator G:
#   s-invariance(F): rho_FtoG   vs rho_FtoH
#   v-invariance(H): rho_HfromF vs rho_HfromG
#   t-invariance(G): rho_GfromF vs rho_GtoH
#   W(1) = 1:        omega_1 - 1 vs 0 (same for all tasks of a subject)
invariance_by_task <- dat_wide %>%
  distinct(subj_id, subj_code = subj, standard_modality, target_modality) %>%
  rowwise() %>%
  mutate(
    mediator = setdiff(modalities, c(standard_modality, target_modality)),
    F = modality_letter[standard_modality],
    H = modality_letter[target_modality],
    G = modality_letter[mediator],
    s_invariance = diff_includes_zero(subj_code, paste0("rho_", F, "to", G),
                                      paste0("rho_", F, "to", H)),
    v_invariance = diff_includes_zero(subj_code, paste0("rho_", H, "from", F),
                                      paste0("rho_", H, "from", G)),
    t_invariance = diff_includes_zero(subj_code, paste0("rho_", G, "from", F),
                                      paste0("rho_", G, "to", H)),
    w_invariance = {
      samp <- rho_samples[[subj_code]]
      if (is.null(samp) || !"omega_1" %in% names(samp)) {
        NA
      } else {
        q <- quantile(samp$omega_1 - 1, c(0.025, 0.975))
        q[1] <= 0 & q[2] >= 0
      }
    },
    task = factor(paste(standard_modality, arrow_sep, target_modality),
                  levels = task_levels)
  ) %>%
  ungroup() %>%
  select(subj_id, task, s_invariance, v_invariance, t_invariance, w_invariance)

diff_summary <- diff_summary %>%
  left_join(invariance_by_task, by = c("subj_id", "task"))

# Build per-panel plotting data ------------------------------------------------
build_invariance_rows <- function(subj_id, subj_code, standard_modality,
                                  target_modality, diff_row) {
  mediator <- setdiff(modalities, c(standard_modality, target_modality))
  F <- modality_letter[standard_modality]
  H <- modality_letter[target_modality]
  G <- modality_letter[mediator]
  task_lab <- paste(standard_modality, arrow_sep, target_modality)
  
  diff_tbl <- tibble(
    subj_id = subj_id, 
    task = task_lab, 
    group = "diff",
    y_label = "x[list(fh)] - x[list(fgh)]", 
    estimate = diff_row$mean_diff,
    ci_lower = diff_row$ci_lower, 
    ci_upper = diff_row$ci_upper,
    holds = NA
  )
  
  # Row for W(1) = 1: CI of omega_1 - 1
  make_w_row <- function() {
    samp <- rho_samples[[subj_code]]
    if (is.null(samp) || !"omega_1" %in% names(samp)) {
      return(tibble(subj_id = subj_id, task = task_lab,
                    group = "w-invariance", y_label = "omega[1] - 1",
                    estimate = NA_real_, ci_lower = NA_real_,
                    ci_upper = NA_real_, holds = NA))
    }
    v <- samp$omega_1 - 1
    q <- quantile(v, c(0.025, 0.975))
    tibble(subj_id = subj_id, task = task_lab, group = "w-invariance",
           y_label = "omega[1] - 1", estimate = mean(v),
           ci_lower = unname(q[1]), ci_upper = unname(q[2]),
           holds = q[1] <= 0 & q[2] >= 0)
  }
  
  # One row per invariance: CI of the posterior difference (param1 - param2)
  make_pair <- function(group_name, param1, param2) {
    r     <- get_diff_ci(subj_code, param1, param2)
    holds <- if (anyNA(r)) NA else r$ci_lower <= 0 & r$ci_upper >= 0
    tibble(subj_id = subj_id, task = task_lab, group = group_name,
           y_label = group_name, estimate = r$mean,
           ci_lower = r$ci_lower, ci_upper = r$ci_upper, holds = holds)
  }
  
  bind_rows(
    diff_tbl,
    make_pair("s-invariance", paste0("rho_", F, "to", G), paste0("rho_", F, "to", H)),
    make_pair("v-invariance", paste0("rho_", H, "from", F), paste0("rho_", H, "from", G)),
    make_pair("t-invariance", paste0("rho_", G, "from", F), paste0("rho_", G, "to", H)),
    make_w_row()
  )
}

fitted_subj <- names(rho_samples)

task_combos <- dat_wide %>%
  distinct(subj_id, subj_code = subj, standard_modality, target_modality) %>%
  filter(subj_code %in% fitted_subj) %>%
  mutate(task = factor(paste(standard_modality, arrow_sep, target_modality),
                       levels = task_levels)) %>%
  left_join(diff_summary, by = c("subj_id", "task"))

# Sanity check: mean_diff/ci_lower/ci_upper should be populated, not NA
# task_combos %>% select(subj_id, task, mean_diff, ci_lower, ci_upper) %>% head()

plot_data <- pmap_dfr(
  task_combos,
  function(subj_id, subj_code, standard_modality, target_modality,
           mean_diff, ci_lower, ci_upper, ...) {
    diff_row <- tibble(mean_diff = mean_diff,
                       ci_lower = ci_lower, ci_upper = ci_upper)
    build_invariance_rows(subj_id, subj_code, standard_modality,
                          target_modality, diff_row)
  }
)

plot_data <- plot_data %>%
  mutate(
    y_label = factor(y_label, levels = rev(c(
      "x[list(fh)] - x[list(fgh)]",
      "s-invariance",
      "v-invariance",
      "t-invariance",
      "omega[1] - 1"
    ))),
    # Color by whether the CI includes zero (for the W(1) row: omega_1 - 1):
    # modal_blue = yes (holds / consistent with transitivity)
    # uni_red    = no  (violated)
    # grey60     = NA  (not yet fitted)
    includes_zero = case_when(
      group == "diff" ~ ci_lower <= 0 & ci_upper >= 0,
      TRUE            ~ holds
    ),
    point_color = case_when(
      is.na(includes_zero) ~ "grey60",
      includes_zero        ~ modal_blue,
      TRUE                 ~ uni_red
    )
  )

# Panel border colors --------------------------------------------------------
# Red if diff and the conditions (s-, v-, t-invariance, W(1) = 1) are
# inconsistent with each other:
#   - diff CI includes zero (no violation) but any invariance fails
#   - diff CI excludes zero (violation) but all invariances hold
# Grey otherwise (consistent picture).
border_df <- diff_summary %>%
  filter(subj_id %in% unique(plot_data$subj_id)) %>%
  mutate(
    all_inv_hold = s_invariance %in% TRUE & v_invariance %in% TRUE &
      t_invariance %in% TRUE & w_invariance %in% TRUE,
    any_inv_fail = s_invariance %in% FALSE | v_invariance %in% FALSE |
      t_invariance %in% FALSE | w_invariance %in% FALSE,
    border_color = case_when(
      ci_includes_zero  & any_inv_fail  ~ uni_red,
      !ci_includes_zero & all_inv_hold  ~ uni_red,
      TRUE                              ~ "grey80"
    )
  ) %>%
  select(subj_id, task, border_color)

# Plot --------------------------------------------------------------------
p_invariance <- ggplot(plot_data, aes(y = y_label, x = estimate, color = point_color)) +
  # Colored panel borders (drawn first so points sit on top)
  geom_rect(data = border_df,
            aes(color = border_color),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = NA, linewidth = 1.2, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.7) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  orientation = "y", linewidth = 0.5, size = 0.3) +
  facet_grid2(subj_id ~ task, scales = "free_x", independent = "x",
              labeller = labeller(task = task_strip_labels)) +
  scale_color_identity() +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(x = "Estimate (dB)", y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border     = element_blank(),   # replaced by geom_rect
    axis.text.y      = element_text(size = 16),
    axis.text.x      = element_text(size = 14),
    strip.text       = element_text(size = 18, face = "bold"),
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 16),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

p_invariance

ggsave("output/invariance_diff_overview_with.pdf", p_invariance,
       width = 13, height = 15, device = cairo_pdf)

