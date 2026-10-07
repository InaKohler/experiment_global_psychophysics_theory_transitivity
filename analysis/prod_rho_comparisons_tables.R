# prod_rho_comparisons_tables.R
#
# Directional comparisons of the posterior means of the reference
# parameters, counted across participants (reported in the paper):
#   1. rho_{g->f} < rho_{g<-f}: standard references smaller than target
#      references for the same modality pair (Heller, 2022b)
#   2. v-invariance direction: rho_{h<-f1} > rho_{h<-f2} for the two
#      standards of each target modality h
#   3. rho_{b<-s} and rho_{b<-l} > rho_{b<-b}: cross-modal brightness target
#      references larger than the intra-modal one (Kohler & Heller, 2026)
#
# Input:  output/rho_samples_<subj>.csv
# Output: summary_table and v_summary_table are printed to the console;
#         intra_counts holds the counts for comparison 3

library(dplyr)
library(tidyr)
library(xtable)

# Load rho posterior samples -------------------------------------------------
sample_files <- list.files("output", pattern = "^rho_samples_.*\\.csv$",
                           full.names = TRUE)
rho_samples <- lapply(sample_files, read.csv, stringsAsFactors = FALSE)
names(rho_samples) <- gsub(".*rho_samples_|\\.csv$", "", sample_files)

# Exclude participant 06
rho_samples  <- rho_samples[names(rho_samples) != "06"]

modality_letter  <- c(auditory = "l", visual = "b", tactile = "s")
modality_display <- c(l = "loudness", b = "brightness", s = "vibration")
modalities_raw   <- names(modality_letter)

# All ordered pairs (g, f) with g != f --------------------------------------
pairs <- expand.grid(g = names(modality_letter),
                     f = names(modality_letter),
                     stringsAsFactors = FALSE) %>%
  filter(g != f) %>%
  mutate(G = modality_letter[g],
         F = modality_letter[f],
         param_to   = paste0("rho_", G, "to",   F),  # rho_{g->f}
         param_from = paste0("rho_", G, "from",  F),  # rho_{g<-f}
         pair_label = paste0(modality_display[G], " (g=", G, "), f=", F))

# For each subject x pair: get posterior means and check direction -----------
results <- lapply(names(rho_samples), function(subj_code) {
  samp <- rho_samples[[subj_code]]
  
  lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i, ]
    if (!all(c(p$param_to, p$param_from) %in% names(samp))) return(NULL)
    tibble(
      subj       = subj_code,
      g_mod      = p$g,
      f_mod      = p$f,
      G          = p$G,
      F          = p$F,
      mean_to    = mean(samp[[p$param_to]]),
      mean_from  = mean(samp[[p$param_from]]),
      to_lt_from = mean(samp[[p$param_to]]) < mean(samp[[p$param_from]])
    )
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(pair_label = ifelse(
    to_lt_from,
    paste0("\u03c1(", G, "\u2192", F, ") < \u03c1(", G, "\u2190", F, ")"),
    paste0("\u03c1(", G, "\u2192", F, ") > \u03c1(", G, "\u2190", F, ")")
  ))

# Per modality pair: count subjects where rho_{g->f} < rho_{g<-f} -----------
pair_counts <- results %>%
  group_by(G, F) %>%
  summarise(
    n_subj     = n(),
    n_asymm    = sum(to_lt_from, na.rm = TRUE),
    g_modality = modality_display[first(G)],
    .groups = "drop"
  ) %>%
  mutate(pair_label = paste0("\u03c1(", G, "\u2192", F,
                             ") < \u03c1(", G, "\u2190", F, ")")) %>%
  arrange(G, F)

# Per modality g: aggregate across both f partners --------------------------
mod_counts <- pair_counts %>%
  group_by(G, g_modality) %>%
  summarise(
    pair_label = paste0(modality_display[first(G)], " (total)"),
    n_subj     = sum(n_subj),
    n_asymm    = sum(n_asymm),
    .groups = "drop"
  )

# Overall -------------------------------------------------------------------
overall <- tibble(
  G          = "all",
  g_modality = "Overall",
  pair_label = "Overall",
  n_subj     = sum(pair_counts$n_subj),
  n_asymm    = sum(pair_counts$n_asymm)
)

# Combined table ------------------------------------------------------------
summary_table <- bind_rows(
  pair_counts %>% select(g_modality, pair_label, n_subj, n_asymm),
  mod_counts  %>% select(g_modality, pair_label, n_subj, n_asymm),
  overall     %>% select(g_modality, pair_label, n_subj, n_asymm)
) %>%
  mutate(prop = sprintf("%d / %d", n_asymm, n_subj))

print(summary_table, n = Inf)

# V-invariance directional check ---------------------------------------------
# For each target modality h and its two possible standards f1, f2:
# is rho_{h<-f1} > rho_{h<-f2} consistently across subjects?
# One unordered pair {f1, f2} per h; canonical order: alphabetical by letter.

v_pairs <- expand.grid(h = names(modality_letter),
                       f1 = names(modality_letter),
                       f2 = names(modality_letter),
                       stringsAsFactors = FALSE) %>%
  filter(h != f1, h != f2, f1 != f2) %>%
  mutate(H  = modality_letter[h],
         F1 = modality_letter[f1],
         F2 = modality_letter[f2]) %>%
  filter(F1 < F2) %>%   # keep each unordered pair once
  mutate(
    param_f1   = paste0("rho_", H, "from", F1),  # rho_{h<-f1}
    param_f2   = paste0("rho_", H, "from", F2)   # rho_{h<-f2}
  )

v_results <- lapply(names(rho_samples), function(subj_code) {
  samp <- rho_samples[[subj_code]]
  
  lapply(seq_len(nrow(v_pairs)), function(i) {
    p <- v_pairs[i, ]
    if (!all(c(p$param_f1, p$param_f2) %in% names(samp))) return(NULL)
    tibble(
      subj        = subj_code,
      h_mod       = p$h,
      H           = p$H,
      F1          = p$F1,
      F2          = p$F2,
      mean_f1     = mean(samp[[p$param_f1]]),
      mean_f2     = mean(samp[[p$param_f2]]),
      f1_gt_f2    = mean(samp[[p$param_f1]]) > mean(samp[[p$param_f2]])
    )
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(pair_label = ifelse(
    f1_gt_f2,
    paste0("\u03c1(", H, "\u2190", F1, ") > \u03c1(", H, "\u2190", F2, ")"),
    paste0("\u03c1(", H, "\u2190", F1, ") < \u03c1(", H, "\u2190", F2, ")")
  ))

v_pair_counts <- v_results %>%
  group_by(H, F1, F2) %>%
  summarise(
    n_subj     = n(),
    n_asymm    = sum(f1_gt_f2, na.rm = TRUE),
    h_modality = modality_display[first(H)],
    .groups = "drop"
  ) %>%
  mutate(pair_label = paste0("\u03c1(", H, "\u2190", F1,
                             ") > \u03c1(", H, "\u2190", F2, ")")) %>%
  arrange(H)

v_mod_counts <- v_pair_counts %>%
  group_by(H, h_modality) %>%
  summarise(
    pair_label = paste0(modality_display[first(H)], " (total)"),
    n_subj     = sum(n_subj),
    n_asymm    = sum(n_asymm),
    .groups = "drop"
  )

v_overall <- tibble(
  H          = "all",
  h_modality = "Overall",
  pair_label = "Overall",
  n_subj     = sum(v_pair_counts$n_subj),
  n_asymm    = sum(v_pair_counts$n_asymm)
)

v_summary_table <- bind_rows(
  v_pair_counts %>% select(h_modality, pair_label, n_subj, n_asymm),
  v_mod_counts  %>% select(h_modality, pair_label, n_subj, n_asymm),
  v_overall     %>% select(h_modality, pair_label, n_subj, n_asymm)
) %>%
  mutate(prop = sprintf("%d / %d", n_asymm, n_subj))

print(v_summary_table, n = Inf)

# Cross-modal vs intra-modal brightness target reference --------------------
# rho_{b<-s} vs rho_{b<-b}  and  rho_{b<-l} vs rho_{b<-b}
# Following Kohler & Heller (2026): cross-modal target references are
# expected to exceed the intra-modal reference (rho_{b<-x} > rho_{b<-b})

intra_comparisons <- list(
  list(cross = "rho_bfroms", label_cross = "s", intra = "rho_bfromb"),
  list(cross = "rho_bfroml", label_cross = "l", intra = "rho_bfromb")
)

intra_results <- lapply(names(rho_samples), function(subj_code) {
  samp <- rho_samples[[subj_code]]
  lapply(intra_comparisons, function(comp) {
    if (!all(c(comp$cross, comp$intra) %in% names(samp))) return(NULL)
    cross_mean <- mean(samp[[comp$cross]])
    intra_mean <- mean(samp[[comp$intra]])
    tibble(
      subj         = subj_code,
      cross_gt_intra = cross_mean > intra_mean,
      pair_label   = ifelse(
        cross_mean > intra_mean,
        paste0("\u03c1(b\u2190", comp$label_cross, ") > \u03c1(b\u2190b)"),
        paste0("\u03c1(b\u2190", comp$label_cross, ") < \u03c1(b\u2190b)")
      ),
      comparison   = paste0("b\u2190", comp$label_cross, " vs b\u2190b")
    )
  }) %>% bind_rows()
}) %>% bind_rows()

intra_counts <- intra_results %>%
  group_by(comparison) %>%
  summarise(
    n_subj       = n(),
    n_asymm      = sum(cross_gt_intra, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pair_label = paste0("\u03c1(", gsub(" vs .*", "", comparison),
                        ") > \u03c1(", gsub(".* vs ", "", comparison), ")"),
    prop = sprintf("%d / %d", n_asymm, n_subj)
  )
