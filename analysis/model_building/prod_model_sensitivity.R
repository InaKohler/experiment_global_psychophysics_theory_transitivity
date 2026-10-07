#' ---
#' title: Model sensitivity analysis for global psychophysics model (transitivity)
#' author: ""
#' output: html_document
#' ---
#'
#' Checks how well the model recovers its parameters (Figure A4 in the
#' paper): parameters and data sets are sampled from the priors, the model is
#' fitted to N_sim of these data sets, and for each parameter
#'   z-score     = (posterior mean - true value) / posterior SD
#'   contraction = 1 - posterior variance / prior variance
#' are computed. z-scores near 0 mean unbiased recovery, contraction near 1
#' means the data are informative about the parameter.
#'
#' The fits take long, so every result is saved right away in a cumulative
#' RDS file; rerunning the script skips simulations that are already saved.
#'
#' Run from analysis/:
#'   source("model_building/prod_model_sensitivity.R")
#'
#' Input:  model_building/gpm_generic.stan, model_building/prior_transitivity.csv,
#'         ../experiment/lookuptable/lookuptab_mavo.txt,
#'         data_transitivity.csv (trial structure of participant 01,
#'         session 01 only), gain_to_decibel.R,
#'         model_building/displacement_to_decibel.R
#' Output: output/sensitivity_transitivity.rds  (prior draws, data, metrics)
#'         output/sensitivity_transitivity.csv  (metrics per simulation and
#'                                               parameter)
#'         output/sensitivity_transitivity.pdf  (z-score vs contraction,
#'                                               Figure A4)
#'         output/sensitivity_contraction_transitivity.pdf
#'         output/sensitivity_zscores_transitivity.pdf

# Setup ----
library(dplyr)
library(ggplot2)
library(ggh4x)
library(tictoc)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# globalpsychophysics provides simulate_gpm() and make_datlist_generic()
if (!require(globalpsychophysics)) {
  devtools::install_github("https://github.com/Kaanwoj/globalpsychophysics.git")
}
library(globalpsychophysics)
source("gain_to_decibel.R")
source("model_building/displacement_to_decibel.R")

output_path <- "output"
model       <- "model_building/gpm_generic.stan"
N_sim       <- 500   # number of simulated data sets

# Simulated matches are rounded to values that can occur in the experiment:
# brightness to the nearest dB value of the lookup table, vibration to the
# nearest dB value of the 255 tactor gains, loudness to whole dB
luRGBdB <- read.csv("../experiment/lookuptable/lookuptab_mavo.txt", sep = ";")
get_db_category <- function(tgt_vec) {
  sapply(tgt_vec, \(dbvisual)
         luRGBdB[which.min(abs(luRGBdB$db - dbvisual)), "db"])
}
get_displacement_db_category <- function(db_vec) {
  valid_db <- displacement_to_decibel(gain_to_displacement(1:255))
  sapply(db_vec, \(db_val) valid_db[which.min(abs(valid_db - db_val))])
}

prior_set    <- read.table("model_building/prior_transitivity.csv", header = TRUE, sep = ";")
prior_params <- setNames(as.list(prior_set$set_1), prior_set$parname)

# ── Build data structure ──────────────────────────────────────────────────────
# Same simulation set-up as in the SBC script: trial structure from one
# session file, a dummy data set only provides the structure of datlist.

params_dummy <- list(
  alpha_l    = 50, alpha_s    = 40,
  beta_s     = 0.65, beta_b   = 0.65, beta_l     = 0.65,
  omega_1    = 1,  omega      = 0.8,
  rho_stob   = 5,  rho_sfroml = 15, rho_ltob   = 10,
  rho_ltos   = 20, rho_bfroml = 30, rho_bfroms = 50,
  rho_stol   = 5,  rho_sfromb = 20, rho_lfroms = 10,
  rho_lfromb = 10, rho_btos   = 15, rho_btol   = 40,
  rho_btob   = 40, rho_bfromb = 40
)

df <- read.csv("data_transitivity.csv",
               colClasses = c(subj = "character", session = "character",
                              id = "character", block = "character")) %>%
  filter(subj == "01", session == "01") %>%
  filter(run == min(run)) %>%
  filter(block != "training") %>%
  mutate(tgt_modality = case_when(target_modality == "visual"   ~ "bright",
                                  target_modality == "auditory" ~ "loud",
                                  target_modality == "tactile"  ~ "strong"),
         std_modality = case_when(standard_modality == "visual"   ~ "bright",
                                  standard_modality == "auditory" ~ "loud",
                                  standard_modality == "tactile"  ~ "strong")) %>%
  mutate(task = glue::glue("{std_modality}_{tgt_modality}"))

# Basic conditions; standards converted from device units to dB
# (visual 43 -> 76 dB Lambert, auditory 67 -> 34 dB SPL)
cond <- df %>%
  filter(trial_type == "basic") %>%
  mutate(standard = case_when(standard == 43 ~ 76,
                              standard == 67 ~ 34,
                              TRUE ~ standard)) %>%
  select(task, p, std = standard) %>%
  unique()

ntrials    <- 75   # trials per basic condition
cond$sigma <- 3    # response SD (dB) of the simulated data

# Simulates one data set: ntrials basic matches per condition, then one
# successive match per cross-modal p = 1 basic match, with that match as the
# standard and the third modality as target (f -> g becomes g -> h).
# sigidx groups the trials into the same low/mid/high standard bins as the
# model fit.
simulate_transitivity_dat <- function(params) {
  dat_basic <- simulate_gpm(ntrials, cond, params) %>%
    group_by(task, p, std) %>%
    mutate(idx = cur_group_id()) %>%
    ungroup() %>%
    as.data.frame() %>%
    mutate(tgt = case_when(
      sub("_.*", "", as.character(task)) == "bright" ~ get_db_category(tgt),
      sub("_.*", "", as.character(task)) == "loud"   ~ round(tgt, 0),
      sub("_.*", "", as.character(task)) == "strong" ~ get_displacement_db_category(tgt)
    ))
  
  cond_successive <- dat_basic %>%
    filter(p != 2, task != "bright_bright") %>%
    rename(std_old = std, std = tgt) %>%
    mutate(task = case_when(
      task == "loud_bright"   ~ "bright_strong",
      task == "loud_strong"   ~ "strong_bright",
      task == "bright_loud"   ~ "loud_strong",
      task == "bright_strong" ~ "strong_loud",
      task == "strong_loud"   ~ "loud_bright",
      task == "strong_bright" ~ "bright_loud"
    )) %>%
    select(task, std, p, sigma)
  
  dat_successive <- simulate_gpm(1, cond_successive, params) %>%
    select(task, std, p, sigma, tgt)
  
  rbind(dat_basic %>% select(task, std, p, sigma, tgt), dat_successive) %>%
    mutate(
      task_prefix = sub("_.*", "", as.character(task)),
      sigidx = case_when(
        task_prefix == "bright" & std < 63 ~ 1,
        task_prefix == "bright" & std > 76 ~ 2,
        task_prefix == "bright"            ~ 3,
        task_prefix == "loud"   & std < 34 ~ 4,
        task_prefix == "loud"   & std > 61 ~ 6,
        task_prefix == "loud"              ~ 5,
        task_prefix == "strong" & std < 30 ~ 7,
        task_prefix == "strong" & std > 60 ~ 9,
        task_prefix == "strong"            ~ 8
      )
    ) %>%
    select(-task_prefix)
}

dat_dummy <- suppressWarnings(simulate_transitivity_dat(params = params_dummy))
datlist   <- make_datlist_generic(dat          = dat_dummy,
                                  prior_params = prior_params,
                                  onlyprior    = 1)
datlist[["tgt"]]               <- numeric(datlist$ntotal)
# Fixed lowest standards in dB (upper bounds of the standard references)
datlist[["lowestsoundstd"]]     <- 34
datlist[["lowestlightstd"]]     <- 54
datlist[["lowestvibrationstd"]] <- 10

# ── 1. Sample from prior ──────────────────────────────────────────────────────
# With onlyprior = 1 the model samples parameters and predicted data
# (tgt_pred) from the priors; N_sim draws are used as true values and
# simulated data sets.

m0 <- stan(file = model,
           data = datlist,
           iter = 2000,
           control = list(adapt_delta = 0.99))

idx_sample <- sample(1:length(rstan::extract(m0)[[1]]), N_sim)
tgt_list   <- lapply(idx_sample, \(i) rstan::extract(m0)$tgt_pred[i, ])

parsnames <- c(
  "alpha_l", "alpha_s",
  "beta_b", "beta_l", "beta_s",
  "omega_1", "omega",
  "rho_ltob", "rho_ltos", "rho_ltol",
  "rho_lfromb", "rho_lfroms", "rho_lfroml",
  "rho_btob", "rho_btol", "rho_btos",
  "rho_bfromb", "rho_bfroml", "rho_bfroms",
  "rho_stob", "rho_stol", "rho_stos",
  "rho_sfromb", "rho_sfroml", "rho_sfroms"
)

param0_samples <- rstan::extract(m0)[parsnames] %>%
  as.data.frame() %>%
  .[idx_sample, ]

# ── 2. Prior SDs ──────────────────────────────────────────────────────────────
# Needed for the contraction. Lognormal priors: SD from logmu and logsigma;
# normal priors of the references: sigma (truncation of the standard
# references is not taken into account)

lognormal_sd <- function(logmu, logsigma) {
  sqrt(exp(2 * logmu + logsigma^2) * (exp(logsigma^2) - 1))
}

prior_sds <- c(
  setNames(
    sapply(c("alpha_l", "alpha_s", "beta_b", "beta_l", "beta_s", "omega_1", "omega"),
           \(p) lognormal_sd(prior_params[[paste0(p, "_logmu")]],
                             prior_params[[paste0(p, "_logsigma")]])),
    c("alpha_l", "alpha_s", "beta_b", "beta_l", "beta_s", "omega_1", "omega")
  ),
  setNames(
    sapply(c("rho_ltob", "rho_ltos", "rho_ltol",
             "rho_lfromb", "rho_lfroms", "rho_lfroml",
             "rho_btob", "rho_btol", "rho_btos",
             "rho_bfromb", "rho_bfroml", "rho_bfroms",
             "rho_stob", "rho_stol", "rho_stos",
             "rho_sfromb", "rho_sfroml", "rho_sfroms"),
           \(p) prior_params[[paste0(p, "_sigma")]]),
    c("rho_ltob", "rho_ltos", "rho_ltol",
      "rho_lfromb", "rho_lfroms", "rho_lfroml",
      "rho_btob", "rho_btol", "rho_btos",
      "rho_bfromb", "rho_bfroml", "rho_bfroms",
      "rho_stob", "rho_stol", "rho_stos",
      "rho_sfromb", "rho_sfroml", "rho_sfroms")
  )
)

# ── 3. Cumulative RDS (resume-safe) ───────────────────────────────────────────
# An existing file is reused, so the same draws are fitted after a restart

cumulative_file <- file.path(output_path, "sensitivity_transitivity.rds")
temp_file       <- file.path(output_path, "sensitivity_transitivity_TEMP.rds")

if (file.exists(cumulative_file)) {
  cumulative <- readRDS(cumulative_file)
  if (N_sim > length(cumulative$results)) {
    n_old <- length(cumulative$results)
    cumulative$results <- c(cumulative$results, vector("list", N_sim - n_old))
    message("Resumed: extended from ", n_old, " to ", N_sim, " slots")
  }
  idx_sample     <- cumulative$idx_sample
  tgt_list       <- cumulative$tgt_list
  param0_samples <- cumulative$param0_samples
} else {
  cumulative <- list(
    idx_sample     = idx_sample,
    param0_samples = param0_samples,
    tgt_list       = tgt_list,
    results        = vector("list", N_sim)
  )
  saveRDS(cumulative, cumulative_file)
}

# ── 4. Fit posterior and compute metrics ──────────────────────────────────────
# z-score and contraction per parameter, see header

compute_metrics <- function(param_name, true_value, post_values, prior_sd) {
  post_mean   <- mean(post_values)
  post_sd     <- sd(post_values)
  z_score     <- (post_mean - true_value) / post_sd
  contraction <- 1 - (post_sd^2 / prior_sd^2)
  data.frame(parameter   = param_name,
             true_value  = true_value,
             post_mean   = post_mean,
             post_sd     = post_sd,
             prior_sd    = prior_sd,
             z_score     = z_score,
             contraction = contraction)
}

datlist[["onlyprior"]] <- 0

tic()
for (i in seq_len(N_sim)) {
  if (!is.null(cumulative$results[[i]])) {
    message("Skipping simulation ", i, " (already saved)")
    next
  }
  
  cat(sprintf("\n=== Simulation %d of %d ===\n", i, N_sim))
  
  true_params      <- param0_samples[i, ]
  datlist[["tgt"]] <- tgt_list[[i]]
  
  m_post <- stan(
    file    = model,
    data    = datlist,
    iter    = 2000,
    chains  = 2,
    refresh = 0,
    control = list(max_treedepth = 14)
  )
  
  post_samples <- rstan::extract(m_post)
  
  metrics <- do.call(rbind, lapply(parsnames, \(p)
                                   compute_metrics(p, true_params[[p]],
                                                   post_samples[[p]], prior_sds[[p]])
  ))
  metrics$simulation  <- i
  metrics$n_divergent <- get_num_divergent(m_post)
  
  saveRDS(metrics, temp_file)
  cumulative$results[[i]] <- metrics
  saveRDS(cumulative, cumulative_file)
  file.remove(temp_file)
  
  rm(m_post, post_samples, metrics)
  gc()
}
toc()

# ── 5. Combine and visualize ──────────────────────────────────────────────────
# Dotted lines at z = +/-2 mark the range of unbiased recovery

sensitivity_df <- do.call(rbind, cumulative$results)
write.csv(sensitivity_df,
          file.path(output_path, "sensitivity_transitivity.csv"),
          row.names = FALSE)

# Plotmath labels — same set as SBC script
param_labels <- c(
  alpha_l   = "alpha[l]",
  alpha_s   = "alpha[v]",
  beta_b    = "beta[b]",
  beta_l    = "beta[l]",
  beta_s    = "beta[v]",
  omega_1   = "omega[1]",
  omega     = "omega",
  rho_ltob  = "rho[l %->% b]",
  rho_ltos  = "rho[l %->% v]",
  rho_ltol  = "rho[l %->% l]",
  rho_lfromb = "rho[l %<-% b]",
  rho_lfroms = "rho[l %<-% v]",
  rho_lfroml = "rho[l %<-% l]",
  rho_btob  = "rho[b %->% b]",
  rho_btol  = "rho[b %->% l]",
  rho_btos  = "rho[b %->% v]",
  rho_bfromb = "rho[b %<-% b]",
  rho_bfroml = "rho[b %<-% l]",
  rho_bfroms = "rho[b %<-% v]",
  rho_stob  = "rho[v %->% b]",
  rho_stol  = "rho[v %->% l]",
  rho_stos  = "rho[v %->% v]",
  rho_sfromb = "rho[v %<-% b]",
  rho_sfroml = "rho[v %<-% l]",
  rho_sfroms = "rho[v %<-% v]"
)

sensitivity_df <- sensitivity_df %>%
  mutate(parameter = factor(parameter,
                            levels = parsnames,
                            labels = param_labels[parsnames]))

# z-score vs contraction
sensitivity_plot <- ggplot(sensitivity_df, aes(x = contraction, y = z_score)) +
  geom_point(alpha = 0.6, color = "#76B3AD", size = 2) +
  geom_hline(yintercept = 0,        linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "#951A36") +
  facet_wrap(~parameter, scales = "free",
             labeller = label_parsed) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-10, 10)) +
  labs(x        = "Posterior contraction",
       y        = "Posterior z-score") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 16),
        axis.text        = element_text(size = 12),
        axis.title       = element_text(size = 14),
        panel.grid.minor = element_blank())

ggsave(file.path(output_path, "sensitivity_transitivity.pdf"),
       sensitivity_plot, width = 10, height = 10)

# Contraction boxplots
contraction_plot <- ggplot(sensitivity_df, aes(x = parameter, y = contraction)) +
  geom_boxplot(fill = "#76B3AD", alpha = 0.7, outlier.color = "#B4A069") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_x_discrete(labels = label_parsed) +
  labs(title    = "Posterior Contraction by Parameter",
       subtitle = "Higher values indicate more information gained from data",
       x = "Parameter", y = "Posterior contraction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_path, "sensitivity_contraction_transitivity.pdf"),
       contraction_plot, width = 10, height = 5)

# Z-score boxplots
zscore_plot <- ggplot(sensitivity_df, aes(x = parameter, y = z_score)) +
  geom_boxplot(fill = "#164467", alpha = 0.7, outlier.color = "#951A36") +
  geom_hline(yintercept = 0,        linetype = "solid",  color = "gray50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#951A36") +
  scale_x_discrete(labels = label_parsed) +
  scale_y_continuous(limits = c(-20, 20)) +
  labs(title    = "Posterior Z-scores by Parameter",
       subtitle = "Values near 0 indicate unbiased recovery; |z| > 2 suggests problems",
       x = "Parameter", y = "Posterior z-score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_path, "sensitivity_zscores_transitivity.pdf"),
       zscore_plot, width = 10, height = 5)

