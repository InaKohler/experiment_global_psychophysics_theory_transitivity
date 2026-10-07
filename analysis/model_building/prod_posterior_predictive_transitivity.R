#' ---
#' title: Posterior predictive checks for global psychophysics model (transitivity)
#' author: ""
#' input: output/gpm_fits.rds (fitted stanfit objects per subject, from model
#'        fitting script), get_analysis_data.R
#' output: PDF-documents in output/
#' ---
#'
#' Perform posterior predictive checks by drawing predicted target intensities
#' from each subject's fitted posterior and evaluating whether predictions fall
#' within psychophysically plausible ranges and exhibit expected monotonic
#' relationships with standard intensity, mirroring the prior predictive checks
#' (same plausible ranges and slope criterion).
#'
#' Output: output/posterior_predictive_quantiles.pdf (Figure A5)
#'         output/posterior_predictive_reg_std.pdf   (Figure A6)
#'
#' The rows of `dat` must match the trials the model was fitted to, in the
#' same order: column j of tgt_pred is the prediction for row j of a
#' subject's data.

# Setup ----
library(dplyr)
library(ggplot2)
library(rstan)

uni_red        <- "#951A36"
modal_blue     <- "#164467"
modal_tuerkis  <- "#76B3AD"

source("get_analysis_data.R")

output_path <- "output"

# All fitted participants, including 06
subj_codes <- c("01", "02", "03", "04", "05", "06", "07", "08")

modality_to_name <- c(visual = "bright", auditory = "loud", tactile = "strong")

nsample <- 5000   # posterior predictive draws per subject

# Rebuild the same observed data frame used for fitting ----
# (the fitting script doesn't save `dat` separately, only the stanfit
# objects, so we reconstruct it here to line up rows with each subject's
# tgt_pred draws)
dat_raw <- get_analysis_data(subj_ids = subj_codes)

dat <- dat_raw %>%
  mutate(
    std_modality_name = modality_to_name[standard_modality],
    tgt_modality_name  = modality_to_name[target_modality],
    task = paste0(std_modality_name, "_", tgt_modality_name),
    std  = standardDB,
    tgt  = matchDB
  ) %>%
  filter(task != "bright_bright") %>%   
  select(subj, task, std, p, tgt)

# Load fitted stanfit objects (single RDS holding all subjects) ----
results <- readRDS(file.path(output_path, "gpm_fits.rds"))

# Extract posterior predictive draws per subject, freeing memory as we go ----
# One row per trial and draw: observed trial data plus tgt_pred and the
# draw index (sampleidx)
pred_smpl_dat <- do.call(rbind, lapply(subj_codes, function(subj_code) {
  res <- results[[subj_code]]
  if (is.null(res) || is.null(res$fit)) {
    message(subj_code, ": no fit available, skipping.")
    return(NULL)
  }
  
  dat_subj <- dat %>% filter(subj == subj_code) %>% select(-subj)
  
  tgt_pred_mat <- rstan::extract(res$fit)$tgt_pred
  n_draws <- nrow(tgt_pred_mat)
  prdidx  <- sample(seq_len(n_draws), min(nsample, n_draws))
  
  tmp <- cbind(
    dat_subj,
    tgt_pred = as.vector(t(tgt_pred_mat[prdidx, ]))
  )
  tmp$sampleidx <- rep(seq_along(prdidx), each = nrow(dat_subj))
  tmp$subj <- subj_code
  tmp
}))
rm(results)
gc()

# Plausible range bounds per target modality (same as prior predictive check) ----
# A draw counts as "inside" for a subject and task when the 2.5% and 97.5%
# quantiles of its predicted matches lie within the bounds of the target
# modality
range_bounds <- list(
  loud   = c(10, 100),   # dB SPL
  bright = c(10, 100),   # dB Lambert
  strong = c(10, 90)     # dB displacement
)

pred_smpl_dat <- pred_smpl_dat %>%
  mutate(tgt_modality = sub(".*_", "", as.character(task)))

# inrange is computed per subject x task x sampleidx (crossing subj into the
# split keeps subjects' draws from being pooled together when checking coverage)
pred_smpl_dat$inrange <- unsplit(
  lapply(split(pred_smpl_dat, list(pred_smpl_dat$subj, pred_smpl_dat$sampleidx,
                                   pred_smpl_dat$task)), function(grp) {
                                     if (nrow(grp) == 0) return(logical(0))
                                     tgt_mod <- grp$tgt_modality[1]
                                     bounds  <- range_bounds[[tgt_mod]]
                                     qu <- quantile(grp$tgt_pred, prob = c(0.025, 0.975), na.rm = TRUE)
                                     rep((qu[1] >= bounds[1]) & (qu[2] <= bounds[2]), nrow(grp))
                                   }),
  list(pred_smpl_dat$subj, pred_smpl_dat$sampleidx, pred_smpl_dat$task)
)

pred_smpl_dat$inrange <- factor(pred_smpl_dat$inrange,
                                levels = c(TRUE, FALSE),
                                labels = c("inside", "outside"))

# Task factor — fixed order + arrow labels (v = vibration) ----
task_levels <- c(
  "bright_loud", "bright_strong",
  "loud_bright", "loud_strong",
  "strong_bright", "strong_loud"
)
task_labels <- c(
  "b \u2192 l",
  "b \u2192 v",
  "l \u2192 b",
  "l \u2192 v",
  "v \u2192 b",
  "v \u2192 l"
)

# Slope of a linear regression of predicted match on standard, per subject,
# task and draw; positive slopes mean matches rise with the standard.
# Build reg_std_df BEFORE converting task to factor (task still has original strings)
reg_std_df <- pred_smpl_dat %>%
  group_by(subj, task, sampleidx) %>%
  summarise(betacoef = coef(lm(tgt_pred ~ std))[2], .groups = "drop") %>%
  mutate(task = factor(task, levels = task_levels, labels = task_labels))

# NOW convert pred_smpl_dat$task
pred_smpl_dat$task <- factor(pred_smpl_dat$task,
                             levels = task_levels,
                             labels = task_labels)

# Proportion inside, aggregated across subjects ----
prop_inside <- pred_smpl_dat %>%
  group_by(task) %>%
  summarise(
    prop  = mean(inrange == "inside"),
    label = sprintf("%.1f%%", prop * 100),
    .groups = "drop"
  )

# Plot: posterior predictive quantile intervals ----
# 90% interval (5%-95%) of the predicted matches of each draw and subject
pp_quantile_plot <- ggplot(pred_smpl_dat,
                           aes(x = sampleidx, y = tgt_pred,
                               group = interaction(subj, sampleidx),
                               color = inrange)) +
  stat_summary(
    fun.data = function(x) {
      data.frame(ymin = quantile(x, 0.05), ymax = quantile(x, 0.95))
    },
    geom = "linerange", linewidth = 0.3
  ) +
  geom_text(data = prop_inside,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5,
            inherit.aes = FALSE, size = 3) +
  scale_color_manual(values = c("inside" = modal_tuerkis, "outside" = uni_red),
                     name = "Plausible range") +
  facet_wrap(~task, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16),
        panel.grid.minor = element_blank()) +
  labs(x = "MCMC sample",
       y = "Predicted target intensity (dB)")

# pp_quantile_plot
ggsave(filename = file.path(output_path, 
                            "posterior_predictive_quantiles.pdf"),
       plot = pp_quantile_plot, width = 12, height = 8,
       device = cairo_pdf)

# Plot: regression slope tgt_pred ~ std, per subject then summarized across subjects ----
pp_reg_std_plot <- ggplot(reg_std_df, aes(y = betacoef, x = task)) +
  geom_boxplot(color = modal_tuerkis, outliers = FALSE) +
  geom_hline(aes(yintercept = 0), linetype = 2, color = "grey70", linewidth = 0.3) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = "task", y = "slope (adjustment ~ standard)")

ggsave(file.path(output_path, "posterior_predictive_reg_std.pdf"),
       pp_reg_std_plot, width = 10, height = 5,
       device = cairo_pdf)

# Regression coefficient summaries ----
reg_std_df$betacoef |> summary()
reg_std_df$betacoef |> mean()
reg_std_df$betacoef |> sd()
