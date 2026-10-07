# createSession.R
#
# Generates session files for the cross-modal transitivity experiment on
# loudness, brightness and vibration strength, as reported in:
#   Kohler, D., & Heller, J. (in prep.). Do internal references determine
#   transitivity? A cross-modal study on loudness, brightness, and vibration
#   strength.
#
# Each file lists all trials of one session and is read by the PsychoPy
# experiment code to control trial order, standard intensities, starting
# values and production factors.
#
# Each block contains 19 trials:
#   7 basic matches (p = 1): all six cross-modal pairs plus visual -> visual
#   6 successive matches (p = 1): the target of each cross-modal basic match
#     becomes the standard for a match to the third modality. Together with
#     the basic matches these give the indirect (x_fgh) and direct (x_fh)
#     productions needed to test transitivity.
#   6 cross-modal matches with p = 2, needed to identify the internal
#     references in the Bayesian model
#
# This script only defines parameters and create_session_file(). Session
# files are written by session/create_next_session_file_call_from_python.R,
# which PsychoPy calls after each session to create the file for the NEXT
# session. Files for the first session are created once by hand (see
# "Generate first-session files" at the bottom).
#
# Helper functions are defined in session/create_session_utils.R.
#
# All paths are relative to the experiment folder (the parent of session/),
# which must be the working directory.
#
# Input:  session/subj_<id>/standard_match.txt  (sessions > 1 only; matches
#         from earlier sessions, used to set target starting values)
# Output: session/subj_<id>/subj_<id>session_<ses>.txt
#
# Created:  19/11/2024, Dorina Kohler
#           22/12/2024, AO
# Last mod: 30/09/2026, Ina Kohler

library(dplyr, warn.conflicts = FALSE)

source("session/create_session_utils.R")

# ── Experiment parameters ─────────────────────────────────────────────────────
nsubj   <- 8  # number of participants
nblocks <- 7  # blocks per session (= repetitions per standard and factor)

# Placeholder value marking the second step of a successive trial, whose
# standard is only known once the participant has completed the first step.
successive_placeholder <- 999

# ── Standard stimuli ──────────────────────────────────────────────────────────
# One standard per modality in the main experiment. All intensities in this
# script (standards, target ranges, limits) are in device units, which the
# PsychoPy code converts to physical intensities. In physical units the
# standards correspond to 76 dB Lambert (visual), 34 dB SPL (auditory) and
# 56 dB displacement (vibrotactile).
#
# NOTE: The *_train vectors are grouped by TARGET modality, not by standard
# modality. vis_train holds the standards for training trials with a visual
# target, which are auditory, tactile and visual standards (see
# standard_modality in the training block of create_session_file()). The
# same holds for aud_train and tac_train.

# Visual (luminance)
vis_train     <- c(36, 40, 48)
vis_standards <- 34

# Auditory (sound level)
aud_train     <- c(11, 22, 33)
aud_standards <- 67

# Vibrotactile (vibration intensity)
tac_train     <- c(11, 22, 33)
tac_standards <- 90

n_standards_per_modality <- length(vis_standards)

# Standards for which successive trials are created
successive_standards <- c(aud_standards, vis_standards, tac_standards)

# ── Target starting values ────────────────────────────────────────────────────
# Range from which the randomised starting intensity of the target is drawn
# in the first session (and for successive trials in every session).
vis_target_range <- 9:44
aud_target_range <- 40:70
tac_target_range <- 30:200

# Adjustment limits per modality. Starting values derived from earlier
# sessions are clamped to these ranges.
vis_limits <- c(9, 49)
aud_limits <- c(25, 85)
tac_limits <- c(10, 255)

# ── Production factors ────────────────────────────────────────────────────────
# production_factors: factors used for basic trials in every block
# successive_factors: factors of the two steps in a successive trial
# Trials with p = 2 are added separately in each block (see below).
production_factors <- 1
n_prod_factors     <- length(production_factors)
successive_factors <- c(1, 1)

# ── Modality combinations ─────────────────────────────────────────────────────
# All six cross-modal pairs plus visual -> visual as an intramodal condition.
modality_mapping_base <- data.frame(
  standards = c("visual",   "visual",
                "auditory", "auditory",
                "tactile",  "tactile",
                "visual"),
  target    = c("auditory", "tactile",
                "visual",   "tactile",
                "visual",   "auditory",
                "visual")
)

# Crossed with every production factor
modality_mapping <- data.frame(
  standards         = rep(modality_mapping_base$standards, n_prod_factors),
  target            = rep(modality_mapping_base$target, n_prod_factors),
  production_factor = rep(production_factors,
                          each = nrow(modality_mapping_base))
)

# Number of basic trials per block
n_basic_trials <- length(vis_standards) * nrow(modality_mapping)


# ── Main function ─────────────────────────────────────────────────────────────
# Creates and writes the session file for one participant and session.
#
#   subj     : participant ID as string, e.g. "01"
#   session  : session number
#   nblocks  : number of blocks in the session (always pass explicitly, e.g.
#              nblocks = nblocks)
create_session_file <- function(subj, session, nblocks = nblocks) {

  # ── Training trials ───────────────────────────────────────────────────────
  # Only used to familiarise participants with the task. They are labelled
  # "training" so they can be excluded easily from analysis.
  df_train <- data.frame(matrix(NA,
                                nrow = length(vis_train) + length(aud_train) +
                                       length(tac_train),
                                ncol = 0)) %>%
    mutate(
      standard          = c(vis_train, aud_train, tac_train),
      target_modality   = c(rep("visual",   length(vis_train)),
                            rep("auditory", length(aud_train)),
                            rep("tactile",  length(tac_train))),
      standard_modality = c("auditory", "tactile", "visual",    # visual target
                            "visual",   "tactile", "visual",    # auditory target
                            "auditory", "visual",  "visual"),   # tactile target
      trial_type        = "basic",
      start_target      = case_when(
        grepl("visual",   target_modality) ~ sample(vis_target_range, n(),
                                                    replace = TRUE),
        grepl("auditory", target_modality) ~ sample(aud_target_range, n(),
                                                    replace = TRUE),
        grepl("tactile",  target_modality) ~ sample(tac_target_range, n(),
                                                    replace = TRUE),
        TRUE ~ NA_integer_
      ),
      presentation_order = "training",
      prod_factor        = sample(production_factors, size = n(),
                                  replace = TRUE),
      id                 = "training",
      block              = "training"
    )

  # ── Experimental blocks ───────────────────────────────────────────────────
  df_session <- data.frame()

  for (b in 1:nblocks) {

    # Basic trials: one per standard and modality combination, created in
    # fixed order first and shuffled together with the successive trials.
    df_basic <- data.frame(matrix(NA, nrow = n_basic_trials, ncol = 0)) %>%
      mutate(
        standard_modality = rep(modality_mapping$standards,
                                each = n_standards_per_modality),
        target_modality   = rep(modality_mapping$target,
                                each = n_standards_per_modality),
        standard          = case_when(
          standard_modality == "visual"   ~ rep(vis_standards,
                                                length.out = n()),
          standard_modality == "auditory" ~ rep(aud_standards,
                                                length.out = n()),
          standard_modality == "tactile"  ~ rep(tac_standards,
                                                length.out = n())
        ),
        trial_type = "basic",
        # Unique trial ID, links basic and successive trials
        id         = sprintf("%02d", 1:n())
      ) %>%
      # prod_factor must exist before add_start_target() is called
      group_by(standard, target_modality) %>%
      mutate(prod_factor = rep(production_factors, length.out = n())) %>%
      ungroup() %>%
      add_start_target(., session = session, subj = subj)

    # Trials with p = 2 for parameter estimation: one per cross-modal pair,
    # using the first standard of each modality.
    df_p_2 <- modality_mapping %>%
      filter(standards != target) %>%
      mutate(
        standard    = case_when(
          standards == "visual"   ~ vis_standards[1],
          standards == "auditory" ~ aud_standards[1],
          standards == "tactile"  ~ tac_standards[1]
        ),
        trial_type  = "basic",
        prod_factor = 2,
        block       = as.character(b),
        id          = paste0(b, ".1.",
                             sprintf("%02d", n_basic_trials + row_number()))
      ) %>%
      rename(standard_modality = standards, target_modality = target) %>%
      add_start_target(., session = session, subj = subj)

    # Successive trials are added to the basic trials, the order is
    # randomised, and each successive pair is kept in the correct order.
    # IDs get the form <block>.<step>.<trial>.
    df_block <- add_successive_trials(df_basic, subj = subj,
                                      session = session) %>%
      sample_frac(., 1L) %>%
      reorder_pairs() %>%
      mutate(id    = case_when(!grepl("2\\.", id) ~ paste0(b, ".1.", id),
                               TRUE               ~ paste0(b, ".", id)),
             block = as.character(b))

    df_session <- bind_rows(df_session, df_block, df_p_2)
  }

  # ── Combine and write ─────────────────────────────────────────────────────
  df_session <- bind_rows(df_train, df_session) %>%
    mutate(presentation_order = 2) %>%
    select(id, block, trial_type, standard_modality, standard,
           target_modality, start_target, prod_factor, presentation_order)

  # Create participant subfolder if it doesn't exist
  dirname <- paste0("session/subj_", subj)
  if (!dir.exists(dirname)) {
    dir.create(dirname)
  }

  write.table(df_session,
              paste0(dirname, "/subj_", subj,
                     "session_", sprintf("%02d", session), ".txt"),
              quote = FALSE, row.names = FALSE)
}


# ── Generate first-session files ──────────────────────────────────────────────
# Run once for all participants before the experiment starts, then comment
# out again. Later sessions are created by
# create_next_session_file_call_from_python.R.
# for (s in sprintf("%02d", 1:nsubj)) {
#   create_session_file(s, session = 1, nblocks = nblocks)
# }


# ── Sanity check ──────────────────────────────────────────────────────────────
# Counts trials per condition in a generated file. Every combination of
# production factor, trial type and modalities should appear equally often.
if (FALSE) {
  read.table("session/subj_01/subj_01session_01.txt", header = TRUE) %>%
    filter(block != "training") %>%
    count(prod_factor, trial_type, standard_modality, target_modality)
}
