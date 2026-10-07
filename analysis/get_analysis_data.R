# get_analysis_data.R
#
# Loads the trial data of the transitivity experiment, converts all
# intensities from device units to dB and pairs direct and indirect matches
# for the transitivity test. Sourced by all analysis scripts.
#
#   get_analysis_data()      : trial data of the experimental blocks
#   get_transitivity_pairs() : direct and indirect matches side by side
#
# Units after conversion (columns standardDB, startDB, matchDB):
#   visual   : dB Lambert, via the monitor lookup table
#   auditory : device units are used unchanged
#   tactile  : dB displacement, via gain_to_decibel()
#
# Input:  data_transitivity.csv
#         ../experiment/lookuptable/lookuptab_mavo.txt
#         matching_gain_to_displacement_250_Hz.txt (read in gain_to_decibel.R)
#
# All analysis scripts are run with analysis/ as working directory.

library(dplyr)

# Setup -------------
source("gain_to_decibel.R")

# Monitor calibration: visualscale (device units) -> dB Lambert
luRGBdB <- read.csv("../experiment/lookuptable/lookuptab_mavo.txt", sep = ";")

# Lab colour palette
uni_red <- "#951A36"
modal_blue <- "#164467"
modal_tuerkis <- "#76B3AD"
modal_gruen <- "#B6BD86"

output_path <- "output"
dir.create(output_path, showWarnings = FALSE)

# Read in data ----------------------------
# Returns all trials of the experimental blocks (training excluded) with the
# columns of the PsychoPy data files plus
#   subj, session, run : participant, session and run within participant
#   subj_id            : identical to subj, used for labels in plots/tables
#   standardDB, startDB, matchDB : intensities in dB (see header)
#
# Participant 06 is excluded by default (only one adjustment in more than
# 50% of trials, see paper).
get_analysis_data <- function(subj_ids = c("01", "02", "03", "04", "05",
                                           "07", "08")){
  dat <- read.csv("data_transitivity.csv",
                  colClasses = c(subj = "character", session = "character",
                                 id = "character", block = "character"),
                  stringsAsFactors = FALSE) %>%
    filter(subj %in% subj_ids) %>%
    mutate(subj_id = subj)

  # Helper: apply a function only to the subset where `condition` is TRUE,
  # leaving everything else NA. Avoids case_when() evaluating gain_to_decibel()
  # on rows it was never meant to receive (visual/auditory values), which is
  # what was producing the NaN warnings.
  safe_apply <- function(x, condition, fun) {
    out <- rep(NA_real_, length(x))
    out[condition] <- fun(x[condition])
    out
  }

  dat <- dat %>%
    filter(block != "training") %>%
    mutate(
      standardDB = case_when(
        standard_modality == "visual"  ~ luRGBdB[match(standard, luRGBdB$visualscale), "db"],
        standard_modality == "tactile" ~ safe_apply(standard, standard_modality == "tactile", gain_to_decibel),
        TRUE ~ standard
      ),
      matchDB = case_when(
        target_modality == "visual"  ~ luRGBdB[match(match, luRGBdB$visualscale), "db"],
        target_modality == "tactile" ~ safe_apply(match, target_modality == "tactile", gain_to_decibel),
        TRUE ~ match
      ),
      startDB = case_when(
        target_modality == "visual"  ~ luRGBdB[match(start, luRGBdB$visualscale), "db"],
        target_modality == "tactile" ~ safe_apply(start, target_modality == "tactile", gain_to_decibel),
        TRUE ~ start
      )
    )
  }


# Pair direct and indirect matches for the transitivity test ----------------
# Returns one row per indirect match x_fgh with the corresponding direct match
# x_fh from the same subject, session, run and block:
#   standard_modality : f, start of the chain
#   mediator          : g, middle modality
#   target_modality   : h, end of the chain
#   successive        : indirect match x_fgh (matchDB)
#   basic             : direct match x_fh (matchDB)
#   diff              : basic - successive
#
# A successive trial <block>.2.<n> continues the basic trial <block>.1.<n>
# (f -> g), so its own standard modality is g. The start modality f is taken
# from that basic trial, and the direct match is the basic trial with p = 1
# from f to h in the same block.
#
# Input: dat from get_analysis_data()
get_transitivity_pairs <- function(dat) {
  basic_trials <- dat %>%
    filter(trial_type == "basic", p == 1)

  dat %>%
    filter(trial_type == "successive") %>%
    mutate(basic_id = sub("\\.2\\.", ".1.", id)) %>%
    # start modality f of the chain
    left_join(basic_trials %>%
                select(subj, session, run, basic_id = id,
                       origin = standard_modality),
              by = c("subj", "session", "run", "basic_id")) %>%
    # direct match f -> h in the same block
    left_join(basic_trials %>%
                select(subj, session, run, block, origin = standard_modality,
                       target_modality, basic = matchDB),
              by = c("subj", "session", "run", "block", "origin",
                     "target_modality")) %>%
    rename(successive = matchDB, mediator = standard_modality) %>%
    rename(standard_modality = origin) %>%
    filter(!is.na(basic) & !is.na(successive)) %>%
    mutate(diff = basic - successive) %>%
    select(subj, subj_id, session, run, block, standard_modality, mediator,
           target_modality, basic, successive, diff)
}
