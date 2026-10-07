# prod_transitivity_ci_diff_table.R
#
# Table of the transitivity test per participant and chain f -> g -> h
# (Table A2 in the paper, counts in Table 2):
#   - mean difference between direct and indirect match (x_fh - x_fgh) with
#     its 95% t-interval; transitivity holds when the interval includes 0
#   - 95% posterior intervals of the reference differences for s-, v- and
#     t-invariance; an invariance holds when the interval includes 0
#   - 95% posterior interval of omega_1 = W(1); W(1) = 1 holds when the
#     interval includes 1
#
# Input:  data via get_analysis_data()
#         output/rho_samples_<subj>.csv (rho_* and omega_1 draws)
# Output: output/transitivity_diff_ci_summary.csv
#         output/transitivity_diff_ci_summary.tex 

# Setup -------------
library(dplyr)
library(tidyr)
library(xtable)

source("get_analysis_data.R")

jnd_lookup <- c(visual = 0.7, auditory = 1.1, tactile = 2.2)
modalities <- c("auditory", "visual", "tactile")

if (!exists("dat")) {
  dat <- get_analysis_data()
}

# Chain labels f -> g -> h (standard -> mediator -> target),
# l = loudness, b = brightness, v = vibration
task_levels <- c(
  "l \u2192 v \u2192 b",   # auditory -> visual   (mediator: tactile)
  "l \u2192 b \u2192 v",   # auditory -> tactile  (mediator: visual)
  "b \u2192 v \u2192 l",   # visual   -> auditory (mediator: tactile)
  "b \u2192 l \u2192 v",   # visual   -> tactile  (mediator: auditory)
  "v \u2192 b \u2192 l",   # tactile  -> auditory (mediator: visual)
  "v \u2192 l \u2192 b"    # tactile  -> visual   (mediator: auditory)
)

modality_display <- c(auditory = "loudness", 
                      visual = "brightness", 
                      tactile = "vibration")

# Direct (basic) and indirect (successive) match of each chain f -> g -> h,
# see get_transitivity_pairs() in get_analysis_data.R
dat_wide <- get_transitivity_pairs(dat) %>%
  mutate(
    task = factor(
      paste(substr(modality_display[standard_modality], 1, 1), "\u2192",
            substr(modality_display[mediator],          1, 1), "\u2192",
            substr(modality_display[target_modality],   1, 1)),
      levels = task_levels),
    JND = jnd_lookup[target_modality]
  )

# Per-subject x task summary: mean diff, 95% CI (t-distribution), and
# whether that CI includes zero -----------------------------------------
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

# Invariance testing via posterior difference CIs ---------------------------
# modality_letter matches the l/b/s scheme used in the rho_ parameter names
# (l = auditory/loud, b = visual/bright, s = tactile/strong)
modality_letter <- c(auditory = "l", visual = "b", tactile = "s")

# Load rho posterior samples (full MCMC draws, one CSV per subject)
sample_files <- list.files("output", pattern = "^rho_samples_.*\\.csv$",
                           full.names = TRUE)
rho_samples <- lapply(sample_files, read.csv, stringsAsFactors = FALSE)
names(rho_samples) <- gsub(".*rho_samples_|\\.csv$", "", sample_files)

# Returns named vector (mean, ci_lower, ci_upper) of (param1 - param2).
# unname() on quantile() is essential: quantile() returns "2.5%" / "97.5%"
# as names which overwrite the names in c(), making indexing by name fail.
get_diff_ci <- function(code, param1, param2) {
  samp <- rho_samples[[code]]
  if (is.null(samp) || !all(c(param1, param2) %in% names(samp))) {
    return(c(mean = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  d <- samp[[param1]] - samp[[param2]]
  c(mean     = mean(d),
    ci_lower = unname(quantile(d, 0.025)),
    ci_upper = unname(quantile(d, 0.975)))
}

# For standard F, target H, mediator G:
#   s-invariance(F): rho_FtoG  vs rho_FtoH    (standard's own ref, two targets)
#   v-invariance(H): rho_HfromF vs rho_HfromG (target's own ref, two standards)
#   t-invariance(G): rho_GfromF vs rho_GtoH   (mediator's target-ref vs standard-ref)
#   W(1) = 1:        omega_1 vs 1             (weighting at p = 1, same for all
#                                              chains of a subject)
invariance_by_task <- dat_wide %>%
  distinct(subj_id, subj_code = subj, standard_modality, target_modality) %>%
  rowwise() %>%
  mutate(
    mediator = setdiff(modalities, c(standard_modality, target_modality)),
    F = modality_letter[standard_modality],
    H = modality_letter[target_modality],
    G = modality_letter[mediator],
    # compute each invariance once and index into the named vector
    s_ci     = list(get_diff_ci(subj_code, paste0("rho_", F, "to", G),
                                paste0("rho_", F, "to", H))),
    v_ci     = list(get_diff_ci(subj_code, paste0("rho_", H, "from", F),
                                paste0("rho_", H, "from", G))),
    t_ci     = list(get_diff_ci(subj_code, paste0("rho_", G, "from", F),
                                paste0("rho_", G, "to", H))),
    # W(1) = omega_1 * 1^omega = omega_1; w_diff holds the posterior mean
    # of omega_1 itself, not a difference
    w_ci     = list({
      samp <- rho_samples[[subj_code]]
      if (is.null(samp) || !"omega_1" %in% names(samp)) {
        c(mean = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
      } else {
        v <- samp$omega_1
        c(mean     = mean(v),
          ci_lower = unname(quantile(v, 0.025)),
          ci_upper = unname(quantile(v, 0.975)))
      }
    }),
    w_diff   = w_ci[["mean"]],
    w_lower  = w_ci[["ci_lower"]],
    w_upper  = w_ci[["ci_upper"]],
    s_diff   = s_ci[["mean"]],
    s_lower  = s_ci[["ci_lower"]],
    s_upper  = s_ci[["ci_upper"]],
    v_diff   = v_ci[["mean"]],
    v_lower  = v_ci[["ci_lower"]],
    v_upper  = v_ci[["ci_upper"]],
    t_diff   = t_ci[["mean"]],
    t_lower  = t_ci[["ci_lower"]],
    t_upper  = t_ci[["ci_upper"]],
    task = factor(
      paste(substr(modality_display[standard_modality], 1, 1), "\u2192",
            substr(modality_display[mediator],          1, 1), "\u2192",
            substr(modality_display[target_modality],   1, 1)),
      levels = task_levels)
  ) %>%
  ungroup() %>%
  mutate(
    # Invariance holds when the posterior difference CI includes zero
    s_holds = s_lower <= 0 & s_upper >= 0,
    v_holds = v_lower <= 0 & v_upper >= 0,
    t_holds = t_lower <= 0 & t_upper >= 0,
    # W(1) = 1 holds when the posterior CI of omega_1 includes 1
    w_holds = w_lower <= 1 & w_upper >= 1
  ) %>%
  select(subj_id, task,
         s_diff, s_lower, s_upper, s_holds,
         v_diff, v_lower, v_upper, v_holds,
         t_diff, t_lower, t_upper, t_holds,
         w_diff, w_lower, w_upper, w_holds)

diff_summary <- diff_summary %>%
  left_join(invariance_by_task, by = c("subj_id", "task"))

print(tibble::as_tibble(diff_summary), n = Inf)

write.csv(diff_summary, "output/transitivity_diff_ci_summary.csv", row.names = FALSE)

# Format invariance CI column: [lower, upper] or -- when not fitted ---------
fmt_inv <- function(lower, upper, holds) {
  case_when(
    is.na(holds) ~ "--",
    TRUE         ~ sprintf("[%.1f, %.1f]", lower, upper)
  )
}

# Build LaTeX-ready data frame -----------------------------------------------
library(kableExtra)

tex_table <- diff_summary %>%
  mutate(
    `Mean diff.`   = sprintf("%.2f", mean_diff),
    `95\\% CI`     = sprintf("[%.2f, %.2f]", ci_lower, ci_upper),
    `S-inv.`       = fmt_inv(s_lower, s_upper, s_holds),
    `V-inv.`       = fmt_inv(v_lower, v_upper, v_holds),
    `T-inv.`       = fmt_inv(t_lower, t_upper, t_holds),
    `$W(1)$`       = fmt_inv(w_lower, w_upper, w_holds)
  ) %>%
  arrange(subj_id, task) %>%
  select(subj_id, Task = task, `$n$` = n,
         `Mean diff.`, `95\\% CI`,
         `S-inv.`, `V-inv.`, `T-inv.`, `$W(1)$`)

# Caption of Table A2
shared_caption <- "Basic vs.\\ successive matching differences by task, with 95\\% posterior credibility intervals of the invariance parameter differences. Matching directions are shown as standard $\\rightarrow$ mediator $\\rightarrow$ target with l = loudness, b = brightness, v = vibration."

# One longtable, rows grouped by subject with pack_rows()
# pack_rows wants a named integer vector: name = section header, value = rows
subj_counts <- table(tex_table$subj_id)
pack_idx    <- setNames(as.integer(subj_counts),
                        paste("Subject", names(subj_counts)))

k <- tex_table %>%
  select(-subj_id) %>%
  kbl(format    = "latex",
      booktabs  = TRUE,
      longtable = TRUE,
      escape    = FALSE,
      caption   = shared_caption,
      label     = "transitivity_diff_ci",
      align     = c("l", "c", "c", "c", "c", "c", "c", "c")) %>%
  kable_styling(latex_options = c("repeat_header"),
                font_size     = 9) %>%
  pack_rows(index = pack_idx)

writeLines(as.character(k),
           "output/transitivity_diff_ci_summary.tex")
