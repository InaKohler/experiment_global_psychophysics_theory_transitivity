# create_next_session_file_call_from_python.R
#
# Creates the session file for one participant's next session. Called by
# PsychoPy at the end of each session, but can also be run by hand.
#
# Usage (from the experiment folder, the parent of session/):
#   Rscript session/create_next_session_file_call_from_python.R <subj> <ses>
#   e.g. Rscript session/create_next_session_file_call_from_python.R 01 3
#
# Without arguments, subject ID and session number are asked for
# interactively.
#
# Input:  session/subj_<id>/standard_match.txt  (sessions > 1 only)
# Output: session/subj_<id>/subj_<id>session_<ses>.txt
#
# Created:  2024, Dorina Kohler
# Last mod: 30/09/2026, Ina Kohler

# Helper functions, parameters and create_session_file()
source("session/create_session_utils.R")
source("session/createSession.R")

# Number of blocks per session
nblocks <- 7

# ── Subject and session ───────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  # Called from Python with arguments
  subj <- args[1]
  ses <- as.integer(args[2])
} else {
  # Fallback for running manually
  cat("Subject ID (e.g., 01, 02, ...): ")
  subj <- readLines(file("stdin"), 1)
  cat("Number of NEXT session (2, 3, ...): ")
  ses <- as.integer(readLines(file("stdin"), 1))
}

# ── Create file ───────────────────────────────────────────────────────────────
create_session_file(subj, session = ses, nblocks = nblocks)
