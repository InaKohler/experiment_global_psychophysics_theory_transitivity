# prod_transitivity_plots_descriptive.R
#
# Per participant: direct match x_fh against indirect match x_fgh for each of
# the six chains f -> g -> h (Figures 4 and A8-A14 in the paper). Points on
# the dashed diagonal mean x_fh = x_fgh; the band marks +/- 1 JND of the
# target modality around the diagonal.
#
# make_subject_diff() plots the difference x_fh - x_fgh instead; its output
# is currently switched off in the loop at the end.
#
# Input:  data via get_analysis_data()
# Output: output/transitivity_scatter_<subj_id>.pdf

# Setup -------------
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggh4x)

source("get_analysis_data.R")

jnd_lookup <- c(visual = 0.7, auditory = 1.1, tactile = 2.2)
modalities <- c("auditory", "visual", "tactile")

if (!exists("dat")) {
  dat <- get_analysis_data()
}

# -------------------
# Mapping to the simulation script:
#   step_2 (indirect estimate)  <-> matchDB where trial_type == "successive"
#   step_3 (direct estimate)   <-> matchDB where trial_type == "basic"
# Exactly 6 panels per page, one per chain f -> g -> h:
#   column 1 = chain starts with loudness (f = auditory),
#   column 2 = chain starts with brightness (f = visual),
#   column 3 = chain starts with vibration (f = tactile)

# Fixed panel order (f -> g -> h, i.e. standard -> mediator -> target),
# same as invariance plot. Labels: l = loudness, b = brightness, v = vibration
combo_levels <- c(
  "l \u2192 v \u2192 b",   # auditory -> visual   (mediator: tactile)
  "l \u2192 b \u2192 v",   # auditory -> tactile  (mediator: visual)
  "b \u2192 v \u2192 l",   # visual   -> auditory (mediator: tactile)
  "b \u2192 l \u2192 v",   # visual   -> tactile  (mediator: auditory)
  "v \u2192 b \u2192 l",   # tactile  -> auditory (mediator: visual)
  "v \u2192 l \u2192 b"    # tactile  -> visual   (mediator: auditory)
)

modality_letter <- c(auditory = "l", visual = "b", tactile = "v")

# Direct (basic) and indirect (successive) match of each chain f -> g -> h,
# see get_transitivity_pairs() in get_analysis_data.R
dat_wide <- get_transitivity_pairs(dat) %>%
  mutate(
    combo = factor(
      paste(modality_letter[standard_modality], "\u2192",
            modality_letter[mediator],          "\u2192",
            modality_letter[target_modality]),
      levels = combo_levels
    ),
    JND = jnd_lookup[target_modality]
  )

# Fixed axis limits per facet (same for x and y so the diagonal stays a diagonal)
combo_limits <- list(
  "l \u2192 v \u2192 b" = c(50, 90),   # auditory -> visual
  "l \u2192 b \u2192 v" = c(20, 85),   # auditory -> tactile
  "b \u2192 v \u2192 l" = c(20, 90),   # visual   -> auditory
  "b \u2192 l \u2192 v" = c(20, 85),   # visual   -> tactile
  "v \u2192 b \u2192 l" = c(20, 90),   # tactile  -> auditory
  "v \u2192 l \u2192 b" = c(50, 90)    # tactile  -> visual
)

# Scatter plot (direct vs indirect match) for a single subject ------------
make_subject_scatter <- function(subj_code, dat_wide) {
  d <- dat_wide %>% filter(subj == subj_code)
  if (nrow(d) == 0) return(NULL)
  subj_label <- unique(d$subj_id)
  
  # +/- 1 JND band around the diagonal. It spans slightly beyond the fixed
  # limits; the panel clips it
  ribbon_df <- d %>%
    distinct(combo, JND) %>%
    mutate(
      lo = sapply(as.character(combo), function(cmb) combo_limits[[cmb]][1]) - 10,
      hi = sapply(as.character(combo), function(cmb) combo_limits[[cmb]][2]) + 10
    ) %>%
    tidyr::pivot_longer(cols = c(lo, hi), values_to = "successive") %>%
    select(combo, successive, JND)
  
  ggplot(d, aes(x = successive, y = basic)) +
    geom_ribbon(
      data = ribbon_df,
      aes(x = successive, ymin = successive - JND, ymax = successive + JND),
      inherit.aes = FALSE, fill = "#164467", alpha = 0.3
    ) +
    geom_abline(slope = 1, intercept = 0, lty = 2, color = "grey40", linewidth = .8) +
    geom_point(size = 2, color = "#A51E37", alpha = 0.8) +
    facet_wrap(~ combo, ncol = 3, dir = "v", scales = "free") +
    ggh4x::facetted_pos_scales(
      x = lapply(levels(d$combo), function(cmb) {
        scale_x_continuous(limits = combo_limits[[cmb]], oob = scales::oob_keep)
      }),
      y = lapply(levels(d$combo), function(cmb) {
        scale_y_continuous(limits = combo_limits[[cmb]], oob = scales::oob_keep)
      })
    ) +
    labs(
      x = "Indirect match (dB)",
      y = "Direct match (dB)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 20, face = "bold"),
      axis.text  = element_text(size = 16),
      axis.title = element_text(size = 18)
    )
}

# Difference plot (basic - successive) for a single subject -----------------
# Differences against the indirect match, with the +/- 1 JND band, the mean
# difference +/- 1 SE per chain and the number of differences within 1 JND.
make_subject_diff <- function(subj_code, dat_wide) {
  d <- dat_wide %>% filter(subj == subj_code)
  if (nrow(d) == 0) return(NULL)
  subj_label <- unique(d$subj_id)
  
  band_df <- d %>%
    distinct(combo, JND) %>%
    mutate(lower = -JND, upper = JND)
  
  agg_dif <- d %>%
    group_by(combo) %>%
    summarise(
      mean_diff = mean(diff),
      sd_diff   = sd(diff),
      n_diff    = n(),
      x_pos     = mean(range(successive, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(se_diff = sd_diff / sqrt(n_diff))
  
  count_df <- d %>%
    group_by(combo) %>%
    summarise(
      n_in  = sum(abs(diff) < JND),
      n_tot = n(),
      label = paste0(n_in, "/", n_tot),
      .groups = "drop"
    )
  
  ggplot(d, aes(x = successive, y = diff)) +
    geom_rect(
      data = band_df,
      aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = "#164467", alpha = 0.3
    ) +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    geom_point(size = 1.8, color = "#A51E37", alpha = 0.8) +
    geom_point(data = agg_dif, aes(x = x_pos, y = mean_diff),
               shape = 18, color = "grey30", size = 3, inherit.aes = FALSE) +
    geom_errorbar(
      data = agg_dif,
      aes(x = x_pos, ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
      color = "grey30", width = 1, inherit.aes = FALSE
    ) +
    geom_text(
      data = count_df, aes(label = label), x = Inf, y = Inf,
      hjust = 1.2, vjust = 1.5, size = 3, inherit.aes = FALSE
    ) +
    facet_wrap(~ combo, ncol = 3, dir = "v", scales = "free") +
    labs(
      # title = paste("Subject:", subj_label),
      x = "Indirect match (dB)",
      y = "Difference (Basic \u2212 Successive)"
    ) +
    theme_minimal()
}

# Generate one PDF per subject, named with subj number ---------------------
# (uncomment the p_diff lines to also save the difference plots)
subj_order <- dat_wide %>% distinct(subj, subj_id) %>% arrange(subj_id)

for (i in seq_len(nrow(subj_order))) {
  subj_code  <- subj_order$subj[i]
  subj_label <- subj_order$subj_id[i]
  
  p_scatter <- make_subject_scatter(subj_code, dat_wide)
  # p_diff    <- make_subject_diff(subj_code, dat_wide)
  
  # print(p_scatter)
  if (!is.null(p_scatter)) {
    ggsave(paste0("output/transitivity_scatter_", subj_label, ".pdf"),
           p_scatter, width = 10, height = 7, device = cairo_pdf)
  }
  
  # if (!is.null(p_diff)) {
  #   ggsave(paste0("output/transitivity_diff_", subj_label, ".pdf"),
  #          p_diff, width = 10, height = 7, device = cairo_pdf)
  # }
}
