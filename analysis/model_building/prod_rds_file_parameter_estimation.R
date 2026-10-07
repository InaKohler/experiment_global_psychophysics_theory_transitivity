# prod_rds_file_parameter_estimation.R
#
# Fits the global psychophysics model (gpm_generic.stan) separately to the
# data of each participant, using all experimental trials (basic and
# successive) with standard and match in dB.
#
# Run from analysis/:
#   source("model_building/prod_rds_file_parameter_estimation.R")
#
# Input:  data via get_analysis_data()
#         model_building/gpm_generic.stan
#         model_building/prior_transitivity.csv
# Output: output/rho_estimates_<subj>.csv   (posterior summary of all rho)
#         output/rho_samples_<subj>.csv     (posterior draws of all rho and
#                                            omega_1, read by the invariance
#                                            scripts)
#         output/rho_estimates_<subj>.pdf   (caterpillar plot)
#         output/rho_estimates_all_subjects.csv
#         output/gpm_fits.rds               (stanfit objects)

# Setup -------------
library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("get_analysis_data.R")
# Provides make_datlist_generic()
if(!require(globalpsychophysics)){
  devtools::install_github("https://github.com/Kaanwoj/globalpsychophysics.git")
}
library(globalpsychophysics)

output_path <- "output"
model       <- "model_building/gpm_generic.stan"

prior_set    <- read.table("model_building/prior_transitivity.csv", header = TRUE, sep = ";")
prior_params <- setNames(as.list(prior_set$set_1), prior_set$parname)

# get_analysis_data() modality names -> the loud/bright/strong naming the
# Stan model and make_datlist_generic() expect
modality_to_name <- c(visual = "bright", auditory = "loud", tactile = "strong")

# Build the model data frame from real trial data (basic + successive) ------
# std = standardDB, tgt = matchDB: both basic and successive trials already
# carry their own actual standard/match values in dB, so no chain simulation
# is needed here (unlike the prior predictive script, which had to simulate
# successive standards because it was testing priors, not fitting real data)
# All participants are fitted, including 06.
dat_raw <- get_analysis_data(subj_ids = c("01", "02", "03", "04", "05", "06",
                                          "07", "08"))

dat <- dat_raw %>%
  mutate(
    std_modality_name = modality_to_name[standard_modality],
    tgt_modality_name  = modality_to_name[target_modality],
    task = paste0(std_modality_name, "_", tgt_modality_name),
    std  = standardDB,
    tgt  = matchDB
  ) %>% mutate(
    # same low/mid/high standard-intensity bins as the prior predictive script
    sigidx = case_when(
      std_modality_name == "bright" & std < 63 ~ 1,
      std_modality_name == "bright" & std > 76 ~ 2,
      std_modality_name == "bright"            ~ 3,
      std_modality_name == "loud"   & std < 34 ~ 4,
      std_modality_name == "loud"   & std > 61 ~ 6,
      std_modality_name == "loud"              ~ 5,
      std_modality_name == "strong" & std < 30 ~ 7,
      std_modality_name == "strong" & std > 60 ~ 9,
      std_modality_name == "strong"            ~ 8
    )
  ) %>%
  select(subj, task, std, p, tgt, sigidx)

# Fit the model for one subject and extract the rho estimates ---------------
fit_subject_rho <- function(subj_code, dat, model, prior_params, output_path) {
  dat_subj <- dat %>% filter(subj == subj_code) %>% select(-subj)
  
  if (nrow(dat_subj) == 0) {
    warning(subj_code, ": no trials found after filtering, skipping.")
    return(NULL)
  }
  
  datlist <- make_datlist_generic(
    dat = dat_subj, prior_params = prior_params, onlyprior = 0
  )
  
  datlist[["lowestsoundstd"]]     <- min(dat_subj$std[grepl("^loud_", dat_subj$task)],   na.rm = TRUE)
  datlist[["lowestlightstd"]]     <- min(dat_subj$std[grepl("^bright_", dat_subj$task)], na.rm = TRUE)
  datlist[["lowestvibrationstd"]] <- min(dat_subj$std[grepl("^strong_", dat_subj$task)], na.rm = TRUE)
  
  # Safety net: the Stan data block requires `tgt` regardless of the
  # onlyprior flag (onlyprior only changes whether the model block *uses*
  # it). If make_datlist_generic() didn't include it, set it explicitly
  # here instead of letting Stan fail with "variable does not exist".
  if (is.null(datlist$tgt)) {
    message(subj_code, ": 'tgt' missing from datlist, setting it explicitly.")
    datlist$tgt <- dat_subj$tgt
  }
  if (is.null(datlist$ntotal)) {
    datlist$ntotal <- length(datlist$tgt)
  }
  
  cat(subj_code, "- datlist fields:\n")
  str(datlist)
  
  m_fit <- stan(
    file    = model,
    data    = datlist,
    iter    = 4000,
    warmup  = 1000,
    chains  = 4,
    seed    = 1234,
    control = list(adapt_delta = 0.99, max_treedepth = 12)
  )
  
  if (length(m_fit@sim$samples) == 0) {
    warning(subj_code, ": sampling failed, no samples returned.")
    return(list(fit = m_fit, rho_table = NULL))
  }
  
  ndivergent <- get_num_divergent(m_fit)
  cat(subj_code, "- divergent transitions:", ndivergent, "\n")
  
  post_summary <- rstan::summary(m_fit)$summary
  
  rho_table <- as.data.frame(post_summary) %>%
    tibble::rownames_to_column("parameter") %>%
    filter(grepl("^rho_", parameter)) %>%
    transmute(
      subj             = subj_code,
      parameter,
      mean             = mean,
      ci_lower         = `2.5%`,
      ci_upper         = `97.5%`,
      n_eff            = n_eff,
      Rhat             = Rhat,
      ci_excludes_zero = ci_lower > 0 | ci_upper < 0
    )
  
  write.csv(
    rho_table,
    file.path(output_path, paste0("rho_estimates_", subj_code, ".csv")),
    row.names = FALSE
  )
  
  # Extract and save full posterior samples for rho_* and omega_1 ----------
  # (needed to compute posterior distributions of pairwise differences for
  # invariance plots; marginal CIs alone are insufficient for that purpose.
  # omega_1 is kept so the invariance plot can test whether its 95% CI
  # includes 1 (W(1) = 1 condition).)
  all_draws <- rstan::extract(m_fit, permuted = TRUE)
  keep      <- grep("^rho_|^omega_1$", names(all_draws), value = TRUE)
  rho_samples <- as.data.frame(all_draws[keep])
  
  write.csv(
    rho_samples,
    file.path(output_path, paste0("rho_samples_", subj_code, ".csv")),
    row.names = FALSE
  )
  
  p_caterpillar <- ggplot(rho_table, aes(x = reorder(parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                    color = "#164467", linewidth = 0.6) +
    coord_flip() +
    labs(
      title = paste("Rho estimates -", subj_code),
      x = NULL, y = "Posterior estimate (95% CI)"
    ) +
    theme_minimal()
  
  ggsave(
    file.path(output_path, paste0("rho_estimates_", subj_code, ".pdf")),
    p_caterpillar, width = 8, height = 6
  )
  
  list(fit = m_fit, rho_table = rho_table)
}

# Run for all subjects --------------------------------------------------
subj_codes <- unique(dat$subj)

results <- lapply(subj_codes, function(subj_code) {
  tryCatch(
    fit_subject_rho(subj_code, dat = dat, model = model,
                    prior_params = prior_params, output_path = output_path),
    error = function(e) {
      message(subj_code, " failed: ", conditionMessage(e))
      NULL
    }
  )
})
names(results) <- subj_codes

rho_table_all <- bind_rows(lapply(results, function(x) {
  if (is.null(x) || is.null(x$rho_table)) return(NULL)
  x$rho_table
}))
print(rho_table_all, max = nrow(rho_table_all) * ncol(rho_table_all))

write.csv(
  rho_table_all,
  file.path(output_path, "rho_estimates_all_subjects.csv"),
  row.names = FALSE
)

# Save the fitted stanfit objects for later use (posterior predictive
# checks, other parameters, etc.)
saveRDS(results, file.path(output_path, "gpm_fits.rds"))
