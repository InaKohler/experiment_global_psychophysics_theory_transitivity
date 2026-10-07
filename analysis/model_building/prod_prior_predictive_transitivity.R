#' ---
#' title: Prior predictive checks for global psychophysics model
#' author: ""
#' input: model_building/gpm_generic.stan, model_building/prior_transitivity.csv,
#'        ../experiment/lookuptable/lookuptab_mavo.txt, data_transitivity.csv
#'        (trial structure of participant 01, session 01 only)
#' output: PDF-documents in output/
#' ---
#' 
#' Perform prior predictive checks for the prior set in prior_transitivity.csv
#' (column set_1): sample parameters and data sets from the priors and
#' evaluate whether predicted adjustments fall within psychophysically
#' plausible ranges (Figure A1 in the paper) and rise with the standard
#' intensity, using linear regression slopes (Figure A2).
#'
#' Output: output/prior_predictive_quantiles_transitiviy.pdf  (Figure A1)
#'         output/prior_predictive_reg_std_transitiviy.pdf    (Figure A2)

# Setup ----
library(dplyr)
library(ggplot2)
library(ggh4x)
library(rstan)

uni_red        <- "#951A36"
modal_blue     <- "#164467"
modal_tuerkis  <- "#76B3AD"

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
model <- "model_building/gpm_generic.stan"
offset <- 94

# Simulated matches are rounded to values that can occur in the experiment:
# brightness to the nearest dB value of the lookup table, vibration to the
# nearest dB value of the 255 tactor gains (loudness to whole dB, see below)
luRGBdB <- read.csv("../experiment/lookuptable/lookuptab_mavo.txt", sep = ";") 
get_db_category <- function(tgt_vec) {
  sapply(tgt_vec, \(dbvisual)
         luRGBdB[which.min(abs(luRGBdB$db - dbvisual)), "db"]
  )}

get_displacement_db_category <- function(db_vec) {
  # Precompute all valid dB levels for gain 1:255
  valid_db <- displacement_to_decibel(gain_to_displacement(1:255))
  
  # Snap each input dB value to the nearest valid level
  sapply(db_vec, \(db_val) {
    valid_db[which.min(abs(valid_db - db_val))]
  })
}

prior_set <- read.table("model_building/prior_transitivity.csv", header = TRUE, sep = ";")
prior_params <- setNames(as.list(prior_set$set_1), prior_set$parname)

# Parameter values for simulating the example data set. It only provides the
# structure of datlist and the plot of the data structure; the prior
# predictive samples come from the priors.
params <- list(
  alpha_l    = 50,
  alpha_s    = 40,
  beta_s     = 0.65,
  beta_b     = 0.65,
  beta_l     = 0.65,
  omega_1    = 1,
  omega      = 0.8,
  rho_stob   = 5,
  rho_sfroml = 5,
  rho_ltob   = 10,
  rho_ltos   = 10,
  rho_bfroml = 30,
  rho_bfroms = 30,
  rho_stol   = 5,
  rho_sfromb = 5,
  rho_lfroms = 10,
  rho_lfromb = 10,
  rho_btos   = 40,
  rho_btol   = 40,
  rho_btob   = 40,
  rho_bfromb = 40
)

# ----
#' ## Build data frame containing the experimental setup
#' Trial structure (tasks, standards, production factors) from one session
#' file; modality names converted to the bright/loud/strong naming of the
#' model.
df <- read.csv("data_transitivity.csv",
               colClasses = c(subj = "character", session = "character",
                              id = "character", block = "character")) %>%
  filter(subj == "01", session == "01") %>%
  filter(run == min(run)) %>%
  filter(block != "training") %>%
  mutate(tgt_modality = case_when(target_modality == "visual" ~ "bright",
                                  target_modality == "auditory" ~ "loud",
                                  target_modality == "tactile" ~ "strong"),
         std_modality = case_when(standard_modality == "visual" ~ "bright",
                                  standard_modality == "auditory" ~ "loud",
                                  standard_modality == "tactile" ~ "strong")) %>%
  mutate(task = glue::glue("{std_modality}_{tgt_modality}")) %>%
  mutate(match_db = case_when(target_modality == "tactile" ~
                                gain_to_decibel(match),
                              target_modality == "visual" ~
                                get_db_category(match),
                              target_modality == "auditory" ~
                                match - offset),
         std_db   = case_when(standard_modality == "tactile" ~
                                gain_to_decibel(standard),
                              standard_modality == "visual" ~
                                get_db_category(standard),
                              standard_modality == "auditory" ~
                                standard - offset))

# Basic conditions; standards converted from device units to dB
# (visual 43 -> 76 dB Lambert, auditory 67 -> 34 dB SPL)
cond <- df %>%
  filter(trial_type == "basic") %>%
  mutate(standard = case_when(standard == 43 ~ 76,
                              standard == 67 ~ 34,
                              TRUE ~ standard)) %>%
  select(task, p, std = standard) %>%
  unique()

# simulate basic trials
ntrials <- 75      # trials per basic condition
cond$sigma <- 3    # response SD (dB) of the simulated data

# simulate adjustements (tgt) for basic trials
dat_basic <- simulate_gpm(ntrials, cond, params)

# index for the conditions
dat_basic <- dat_basic |>
  group_by(task, p, std) |>
  mutate(idx = cur_group_id()) |>
  ungroup() |> as.data.frame()

dat_basic <- dat_basic %>%
  mutate(tgt = case_when(
    sub("_.*", "", as.character(task)) == "bright" ~ get_db_category(tgt),
    sub("_.*", "", as.character(task)) == "loud" ~ round(tgt, 0),
    sub("_.*", "", as.character(task)) == "strong" ~ get_displacement_db_category(tgt)
  ))

# use tgt from basic trials as standards for successive trials
# (cross-modal p = 1 trials only; f -> g becomes g -> h)
cond_successive <- dat_basic %>%
  filter(p != 2,
         task != "bright_bright") %>%
  rename(std_old = std, std = tgt) %>%
  mutate(task = case_when(
    task == "loud_bright"   ~ "bright_strong",
    task == "loud_strong"   ~ "strong_bright",
    task == "bright_loud"   ~ "loud_strong",
    task == "bright_strong" ~ "strong_loud",
    task == "strong_loud"   ~ "loud_bright",
    task == "strong_bright" ~ "bright_loud"
  )) |>
  select(task, std, p, sigma)
dat_successive <- simulate_gpm(1, cond_successive, params) |>
  select(task, std, p, sigma, tgt)

dat <- rbind(dat_basic |> select(task, std, p, sigma, tgt), dat_successive) 

# sigidx: same low/mid/high standard bins as the model fit
dat <- dat %>%
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
  select(-task_prefix)  # drop helper if you don't need it

# dat <- dat_basic
# With onlyprior = 1 the model ignores tgt and samples parameters and
# predicted data (tgt_pred) from the priors
datlist <- make_datlist_generic(dat = dat, 
                                prior_params = prior_params, 
                                onlyprior = 1)

# Fixed lowest standards in dB
datlist[["lowestsoundstd"]]     <- 34   # lowest auditory standard in dB SPL re offset
datlist[["lowestlightstd"]]     <- 54    # lowest visual standard in dB Lambert  
datlist[["lowestvibrationstd"]] <- 10    # lowest tactile standard in dB displacement


# STAN call ----
m_prior <- stan(file = model,
                data = datlist,
                iter = 3000,
                control = list(adapt_delta = 0.99))

# Extract prior predictive samples ----
tgt_pred_mat <- rstan::extract(m_prior)$tgt_pred
ndivergent   <- get_num_divergent(m_prior)
cat("Divergent transitions:", ndivergent, "\n")


# Take a sample ----
# nsample prior predictive data sets, one row per trial and sample
nsample <- 5000
prdidx  <- sample(1:nrow(tgt_pred_mat), nsample)

pred_smpl_dat <- cbind(
  dat,
  tgt_pred = as.vector(t(tgt_pred_mat[prdidx, ]))
)
pred_smpl_dat$sampleidx <- rep(seq_len(nsample), each = nrow(dat))

# Plausible range bounds per target modality ----
# A data set counts as "inside" for a task when the 2.5% and 97.5% quantiles
# of its predicted matches lie within the bounds of the target modality
range_bounds <- list(
  loud   = c(10, 100),   # dB SPL
  bright = c(10, 100),   # dB Lambert
  strong = c(10, 90)     # dB displacement 
)

pred_smpl_dat <- pred_smpl_dat %>%
  mutate(tgt_modality = sub(".*_", "", as.character(task)))

pred_smpl_dat$inrange <- unsplit(
  lapply(split(pred_smpl_dat, list(pred_smpl_dat$sampleidx, 
                                   dat$task)), function(grp) {
                                     if (nrow(grp) == 0) return(logical(0))
                                     tgt_mod  <- grp$tgt_modality[1]
                                     bounds   <- range_bounds[[tgt_mod]]
                                     qu <- quantile(grp$tgt_pred, prob = c(0.025, 0.975), na.rm = TRUE)
                                     rep((qu[1] >= bounds[1]) & (qu[2] <= bounds[2]), nrow(grp))
                                   }),
  list(pred_smpl_dat$sampleidx, pred_smpl_dat$task)
)

pred_smpl_dat$inrange <- factor(pred_smpl_dat$inrange,
                                levels = c(TRUE, FALSE),
                                labels = c("inside", "outside"))

# Task factor — fixed order + arrow labels (v = vibration) ----
task_levels <- c(
  "bright_loud", "bright_strong",
  "loud_bright", "loud_strong",
  "strong_bright", "strong_loud",
  "bright_bright"
)
task_labels <- c(
  "b \u2192 l",
  "b \u2192 v",
  "l \u2192 b",
  "l \u2192 v",
  "v \u2192 b",
  "v \u2192 l",
  "b \u2192 b"
)

pred_smpl_dat$task <- factor(pred_smpl_dat$task,
                             levels = task_levels,
                             labels = task_labels)

# Proportion inside ----
# Percentage of data sets inside the plausible range per task
prop_inside <- pred_smpl_dat %>%
  group_by(task) %>%
  summarise(
    prop  = mean(inrange == "inside"),
    label = sprintf("%.1f%%", prop * 100),
    .groups = "drop"
  )

# Plot: quantile ranges ----
# 90% interval (5%-95%) of the predicted matches of each data set
pp_quantile_plot <- ggplot(pred_smpl_dat,
                           aes(x = sampleidx, y = tgt_pred,
                               group = sampleidx, color = inrange)) +
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

pp_quantile_plot


ggsave(filename = file.path(output_path, 
                            "prior_predictive_quantiles_transitiviy.pdf"),
       plot = pp_quantile_plot, width = 12, height = 8)

# Plot: slope ~ standard intensity ----
# Slope of a linear regression of predicted match on standard, per data set
# and task; positive slopes mean matches rise with the standard
reg_std_df <- tapply(pred_smpl_dat,
                     list(pred_smpl_dat$task, pred_smpl_dat$sampleidx),
                     function(x) coef(lm(tgt_pred ~ std, x))[2]) |>
  array2DF() |>
  setNames(c("task", "sampleidx", "betacoef")) |>
  mutate(task = factor(task, levels = task_labels))

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

pp_reg_std_plot

ggsave(filename = file.path(output_path, "prior_predictive_reg_std_transitiviy.pdf"),
       plot = pp_reg_std_plot, width = 10, height = 5)


