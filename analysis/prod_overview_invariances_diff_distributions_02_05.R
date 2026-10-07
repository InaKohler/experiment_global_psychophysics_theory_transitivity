# prod_overview_invariances_diff_distributions_02_04.R
#
# For selected participants (target_subjs) and every chain f -> g -> h
# (Figure 6 in the paper). Per panel:
#   x_fh - x_fgh : distribution of the trial-wise differences between direct
#                  and indirect match; blue when the 95% t-interval of the
#                  mean difference includes 0, red otherwise
#   s-, v-, t-invariance : the posterior distributions of the two references
#                  compared, as overlapping ridges; blue when the 95%
#                  posterior interval of their difference includes 0, red
#                  otherwise
#
# Input:  data via get_analysis_data()
#         output/rho_samples_<subj>.csv
# Output: output/invariance_posterior_dist_<subj_id>.pdf

# Setup ---------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggridges)

source("get_analysis_data.R")

jnd_lookup      <- c(visual = 0.7, auditory = 1.1, tactile = 2.2)
modalities      <- c("auditory", "visual", "tactile")
# Letters used in the rho_ parameter names (s = tactile/strong)
modality_letter <- c(auditory = "l", visual = "b", tactile = "s")

target_subjs <- c("02", "05")   # subj_id values to plot

arrow_sep   <- "->"

# Fixed panel order (f -> g -> h, i.e. standard -> mediator -> target),
# l = loudness, b = brightness, v = vibration
combo_levels <- c(
  "l \u2192 v \u2192 b",   # auditory -> visual   (mediator: tactile)
  "l \u2192 b \u2192 v",   # auditory -> tactile  (mediator: visual)
  "b \u2192 v \u2192 l",   # visual   -> auditory (mediator: tactile)
  "b \u2192 l \u2192 v",   # visual   -> tactile  (mediator: auditory)
  "v \u2192 b \u2192 l",   # tactile  -> auditory (mediator: visual)
  "v \u2192 l \u2192 b"    # tactile  -> visual   (mediator: auditory)
)

# Display letters (vibration = v). Keep modality_letter above untouched —
# it uses "s" because the rho_ parameter names were built with "s" for tactile.
modality_display_letter <- c(auditory = "l", visual = "b", tactile = "v")

if (!exists("dat")) dat <- get_analysis_data()

# Direct (basic) and indirect (successive) match of each chain f -> g -> h,
# see get_transitivity_pairs() in get_analysis_data.R. task = f -> h
dat_wide <- get_transitivity_pairs(dat) %>%
  mutate(
    task  = paste(standard_modality, arrow_sep, target_modality),
    combo = factor(
      paste(modality_display_letter[standard_modality], "\u2192",
            modality_display_letter[mediator],          "\u2192",
            modality_display_letter[target_modality]),
      levels = combo_levels
    )
  )

# Per-subject x task t-distribution CI for the diff row color ---------------
diff_summary <- dat_wide %>%
  group_by(subj_id, task) %>%
  summarise(
    n         = n(),
    mean_diff = mean(diff),
    se_diff   = sd(diff) / sqrt(n),
    ci_lower  = mean_diff - qt(0.975, df = n - 1) * se_diff,
    ci_upper  = mean_diff + qt(0.975, df = n - 1) * se_diff,
    .groups   = "drop"
  ) %>%
  mutate(ci_includes_zero = ci_lower <= 0 & ci_upper >= 0)

# Load rho posterior samples -------------------------------------------------
sample_files <- list.files("output", pattern = "^rho_samples_.*\\.csv$",
                           full.names = TRUE)
rho_samples <- lapply(sample_files, read.csv, stringsAsFactors = FALSE)
names(rho_samples) <- gsub(".*rho_samples_|\\.csv$", "", sample_files)

# Build long distribution data for one subject x task combination ------------
# diff row: empirical trial differences, colored by t-CI holds/fails
# invariance rows: two raw rho posteriors shown as overlapping ridges
# For standard F, target H, mediator G:
#   s-invariance(F): rho_FtoG   vs rho_FtoH
#   v-invariance(H): rho_HfromF vs rho_HfromG
#   t-invariance(G): rho_GfromF vs rho_GtoH
build_panel_data <- function(sid, std_mod, tgt_mod, task_lab, combo_lab) {
  subj_code <- sid
  mediator  <- setdiff(modalities, c(std_mod, tgt_mod))
  F <- modality_letter[std_mod]
  H <- modality_letter[tgt_mod]
  G <- modality_letter[mediator]
  
  samp <- rho_samples[[subj_code]]
  
  get_param_vec <- function(p) {
    if (is.null(samp) || !(p %in% names(samp))) return(NULL)
    samp[[p]]
  }
  
  diff_vals      <- dat_wide %>%
    filter(subj_id == sid, task == task_lab) %>% pull(diff)
  diff_incl_zero <- diff_summary %>%
    filter(subj_id == sid, task == task_lab) %>% pull(ci_includes_zero)
  
  make_inv_rows <- function(group_name, p1, p2) {
    v1 <- get_param_vec(p1)
    v2 <- get_param_vec(p2)
    # Check whether 95% CI of (p1 - p2) includes zero — same criterion as CI plot
    holds <- if (!is.null(v1) && !is.null(v2)) {
      d <- v1 - v2
      q <- quantile(d, c(0.025, 0.975))
      q[1] <= 0 & q[2] >= 0
    } else NA
    bind_rows(
      if (!is.null(v1)) tibble(group = group_name, param = "param1",
                               value = v1, inv_holds = holds) else NULL,
      if (!is.null(v2)) tibble(group = group_name, param = "param2",
                               value = v2, inv_holds = holds) else NULL
    )
  }
  
  bind_rows(
    tibble(group = "x[list(f,h)] - x[list(f,g,h)]", param = "diff",
           value = diff_vals, inv_holds = NA,
           includes_zero = if (length(diff_incl_zero)) diff_incl_zero else NA),
    make_inv_rows("s-invariance",
                  paste0("rho_", F, "to", G), paste0("rho_", F, "to", H)),
    make_inv_rows("v-invariance",
                  paste0("rho_", H, "from", F), paste0("rho_", H, "from", G)),
    make_inv_rows("t-invariance",
                  paste0("rho_", G, "from", F), paste0("rho_", G, "to", H))
  ) %>%
    mutate(subj_id = sid, task = task_lab, combo = combo_lab)
}

# Build full data for all target subjects ------------------------------------
task_modality_map <- dat_wide %>%
  filter(subj_id %in% target_subjs) %>%
  distinct(subj_id, task, combo, standard_modality, target_modality)

dist_data <- pmap_dfr(task_modality_map, function(subj_id, task, combo,
                                                  standard_modality,
                                                  target_modality) {
  build_panel_data(subj_id, standard_modality, target_modality,
                   task, as.character(combo))
})

# Factor levels and colors ---------------------------------------------------
group_levels <- rev(c(
  "x[list(f,h)] - x[list(f,g,h)]",
  "s-invariance",
  "v-invariance",
  "t-invariance"
))

dist_data <- dist_data %>%
  mutate(
    group = factor(group, levels = group_levels),
    combo = factor(combo, levels = combo_levels),
    fill_color = case_when(
      param == "diff" & is.na(includes_zero)             ~ "grey60",
      param == "diff" & includes_zero                    ~ modal_blue,
      param == "diff" & !includes_zero                   ~ uni_red,
      param %in% c("param1","param2") & is.na(inv_holds) ~ "grey60",
      param %in% c("param1","param2") & inv_holds        ~ modal_blue,
      param %in% c("param1","param2") & !inv_holds       ~ uni_red,
      TRUE                                               ~ "grey60"
    )
  )

# Plot function for a single subject -----------------------------------------
make_subject_dist <- function(sid, dist_data) {
  d <- dist_data %>% filter(subj_id == sid)
  if (nrow(d) == 0) return(NULL)
  
  # With ncol=3, dir="v", the left column holds combo positions 1 and 2
  left_col_combos <- levels(d$combo)[c(1, 2)]
  
  # y-label renderer: parses plotmath so arrow and subscripts render
  parse_labels <- function(x) parse(text = x)
  blank_labels  <- function(x) rep("", length(x))
  
  ggplot(d, aes(x = value, y = group,
                fill = fill_color, color = fill_color,
                group = interaction(group, param))) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.7) +
    geom_density_ridges(
      alpha = 0.4, linewidth = 0.5,
      scale = 0.9,
      rel_min_height = 0.01
    ) +
    facet_wrap(~ combo, ncol = 3, dir = "v", scales = "free") +
    ggh4x::facetted_pos_scales(
      y = lapply(levels(d$combo), function(cmb) {
        if (cmb %in% left_col_combos) {
          scale_y_discrete(labels = parse_labels)
        } else {
          scale_y_discrete(labels = blank_labels)
        }
      })
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(xlim = c(-25, 75)) +
    labs(x = "Parameter value / Difference (dB)", y = NULL) +
    theme_minimal() +
    theme(
      strip.text       = element_text(size = 20, face = "bold"),
      axis.text.y      = element_text(size = 16),
      axis.text.x      = element_text(size = 16),
      axis.title       = element_text(size = 14),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# Save one PDF per subject ---------------------------------------------------
for (sid in target_subjs) {
  p <- make_subject_dist(sid, dist_data)
  # print(p)
  if (!is.null(p)) {
    ggsave(paste0("output/invariance_posterior_dist_", sid, ".pdf"),
           p, width = 10, height = 7, device = cairo_pdf)
  }
}
