# Simulation-based calibration (SBC) of gpm_generic.stan
#
# Checks whether the model recovers its parameters without bias (Figure A3
# in the paper):
#   1. builds the trial structure of one session from a real data file and
#      samples parameters and data sets from the priors (onlyprior = 1)
#   2. fits the model to N_sim = 500 of these simulated data sets and
#      records the rank of the true (prior) value among the posterior draws
#   3. plots the rank histograms; uniform ranks mean unbiased recovery
#   4. summarises divergences, Rhat and n_eff of the fits
#
# The fits take long, so every result is saved right away in a cumulative
# RDS file; rerunning the script skips simulations that are already saved.
#
# Run from analysis/:
#   source("model_building/prod_simulation_based_calibrations.R")
#
# Input:  model_building/gpm_generic.stan, model_building/prior_transitivity.csv,
#         ../experiment/lookuptable/lookuptab_mavo.txt,
#         data_transitivity.csv (trial structure of participant 01,
#         session 01 only), gain_to_decibel.R,
#         model_building/displacement_to_decibel.R
# Output: output/sbc_transitivity.rds  (prior draws, simulated data, ranks)
#         output/sbc_transitivity.pdf  (rank histograms)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(dplyr)
library(tictoc)
library(bayesplot)
# globalpsychophysics provides simulate_gpm() and make_datlist_generic()
if (!require(globalpsychophysics)) {
  devtools::install_github("https://github.com/Kaanwoj/globalpsychophysics.git")
}
library(globalpsychophysics)
source("gain_to_decibel.R")
source("model_building/displacement_to_decibel.R")

output_path <- "output"
model       <- "model_building/gpm_generic.stan"
offset      <- 94

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

# ── Data simulation helper ────────────────────────────────────────────────────
# Trial structure (tasks, standards, production factors) from one session
# file; modality names converted to the bright/loud/strong naming of the
# model

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

# ── 1. Build datlist and sample from prior ────────────────────────────────────
# A dummy data set with arbitrary parameter values only provides the
# structure of datlist (conditions, indices); with onlyprior = 1 the model
# ignores tgt and samples parameters and predicted data (tgt_pred) from the
# priors.

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

dat_dummy <- suppressWarnings(simulate_transitivity_dat(params = params_dummy))
datlist   <- make_datlist_generic(dat          = dat_dummy,
                                  prior_params = prior_params,
                                  onlyprior    = 1)
datlist[["tgt"]] <- numeric(datlist$ntotal)
# Fixed lowest standards in dB — from the experimental design
datlist[["lowestsoundstd"]]     <- 34   # lowest auditory standard in dB SPL re offset
datlist[["lowestlightstd"]]     <- 54    # lowest visual standard in dB Lambert  
datlist[["lowestvibrationstd"]] <- 10    # lowest tactile standard in dB displacement

tic()
m0 <- stan(file    = model,
           data    = datlist,
           iter    = 2000,
           control = list(adapt_delta = 0.99))
toc()
cat("Divergent transitions (prior):", get_num_divergent(m0), "\n")

# Parameters checked in the SBC (all 25 free parameters except sig)
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

# ── 2. Set up cumulative RDS (resume-safe) ────────────────────────────────────
# N_sim prior draws are selected once and stored with their simulated data
# sets (tgt_list). An existing file is reused, so the same draws are fitted
# after a restart; an invalid file is rebuilt.

N_sim           <- 500
cumulative_file <- file.path(output_path, "sbc_transitivity.rds")
temp_file       <- file.path(output_path, "sbc_transitivity_TEMP.rds")

if (file.exists(cumulative_file)) {
  cumulative <- readRDS(cumulative_file)
  
  if (is.null(cumulative$tgt_list) ||
      length(cumulative$tgt_list) == 0 ||
      any(sapply(cumulative$tgt_list, is.null))) {
    message("Cumulative file invalid — rebuilding from m0")
    file.remove(cumulative_file)
    cumulative <- NULL
  }
}

if (!file.exists(cumulative_file)) {
  idx_sim     <- sample(1:length(rstan::extract(m0)[[1]]), N_sim)
  tgt_list    <- lapply(idx_sim, \(i) rstan::extract(m0)$tgt_pred[i, ])
  param_prior <- rstan::extract(m0)[parsnames] %>%
    as.data.frame() %>%
    .[idx_sim, ]
  
  cumulative <- list(
    idx_sim     = idx_sim,
    param_prior = param_prior,
    tgt_list    = tgt_list,
    results     = vector("list", N_sim)
  )
  saveRDS(cumulative, cumulative_file)
} else {
  if (N_sim > length(cumulative$results)) {
    n_old <- length(cumulative$results)
    cumulative$results <- c(cumulative$results, vector("list", N_sim - n_old))
    message("Resumed: extended from ", n_old, " to ", N_sim, " slots")
  }
  idx_sim     <- cumulative$idx_sim
  param_prior <- cumulative$param_prior
  tgt_list    <- cumulative$tgt_list
}

# ── 3. Fit posterior for each simulated dataset ───────────────────────────────
# For each simulated data set: fit the model, thin the posterior draws to
# every 8th, and count how many draws lie below the true value (rank).
# Each result is written to the cumulative file immediately; failed fits are
# reported and left empty.

N_s               <- 5000
datlist[["onlyprior"]] <- 0

tic()
for (i in seq_len(N_sim)) {
  if (!is.null(cumulative$results[[i]])) {
    message("Skipping simulation ", i, " (already saved)")
    next
  }
  
  cat(sprintf("\n=== SBC simulation %d / %d ===\n", i, N_sim))
  
  datlist[["tgt"]] <- tgt_list[[i]]
  
  tryCatch({
    m_post <- stan(
      file    = model,
      data    = datlist,
      iter    = N_s / 2,
      chains  = 4,
      refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    
    thinner    <- seq(from = 1, to = N_s / 2, by = 8)
    post_draws <- rstan::extract(m_post)[parsnames] %>%
      as.data.frame() %>%
      .[thinner, ]
    
    ranks_i <- setNames(
      sapply(parsnames, \(p) sum(post_draws[, p] < param_prior[i, p])),
      parsnames
    )
    
    result_i <- list(
      ranks       = ranks_i,
      n_divergent = get_num_divergent(m_post),
      rhat        = summary(m_post)$summary[, "Rhat"] %>% mean(na.rm = TRUE),
      n_eff       = summary(m_post)$summary[, "n_eff"] %>% mean(na.rm = TRUE)
    )
    
    saveRDS(result_i, temp_file)
    cumulative$results[[i]] <- result_i
    saveRDS(cumulative, cumulative_file)
    file.remove(temp_file)
    
    rm(m_post, post_draws)
    gc()
    
  }, error = function(e) {
    message(sprintf("Simulation %d failed: %s", i, conditionMessage(e)))
  })
}
toc()

# ── 4. SBC histograms ─────────────────────────────────────────────────────────
# Rank histograms with B bins; the dashed lines mark the 99% interval of the
# bin counts expected under uniform ranks.
library(ggplot2)

thinner <- seq(from = 1, to = N_s / 2, by = 8)

rank_param <- do.call(rbind, lapply(cumulative$results, \(r) r$ranks)) %>%
  as.data.frame()

B    <- 16
ci_l <- qbinom(0.005, size = N_sim, prob = 1 / B)
ci_u <- qbinom(0.995, size = N_sim, prob = 1 / B)

# Plotmath labels for each parameter (v = vibration, s in the parameter names)
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

# Reshape to long format
rank_long <- rank_param %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to  = "parameter",
                      values_to = "rank") %>%
  mutate(parameter = factor(parameter,
                            levels = parsnames,
                            labels = param_labels[parsnames]))

sbc_plot <- ggplot(rank_long, aes(x = rank)) +
  geom_histogram(breaks = seq(0, length(thinner), length.out = B + 1),
                 fill = "#76B3AD", color = "white", linewidth = 0.3) +
  geom_hline(yintercept = c(ci_l, ci_u), linetype = 2, color = "gray40") +
  facet_wrap(~parameter, ncol = 5,
             labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Rank", y = "Count") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 20),
    axis.text        = element_text(size = 16),
    axis.title       = element_text(size = 18),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_path, "sbc_transitivity.pdf"),
       sbc_plot, width = 14, height = 10)

# ── 5. Convergence summary ────────────────────────────────────────────────────
# Divergent transitions, mean Rhat and mean n_eff across all SBC fits

div_trans <- sapply(cumulative$results, \(r) r$n_divergent)
rhat      <- sapply(cumulative$results, \(r) r$rhat)
n_eff     <- sapply(cumulative$results, \(r) r$n_eff)

cat("Divergences table:\n"); print(table(div_trans))
cat("Mean Rhat:", mean(rhat, na.rm = TRUE), "\n")
cat("Mean n_eff:", mean(n_eff, na.rm = TRUE), "\n")

