# create_session_utils.R
#
# Helper functions for createSession.R. They rely on the global parameters
# defined there (successive_standards, successive_factors,
# successive_placeholder, *_target_range, *_limits), so this file must be
# sourced from createSession.R and not run on its own.
#
#   add_successive_trials() : creates the second step of each indirect match
#   add_start_target()      : sets the starting intensity of the target
#   find_closest_match()    : looks up a starting value from earlier sessions
#   reorder_pairs()         : keeps each successive trial after its basic trial
#
# Created:  19/11/2024, Dorina Kohler
# Last mod: 30/09/2026, Ina Kohler

modalities <- c("visual", "auditory", "tactile")


# ── Successive trials ─────────────────────────────────────────────────────────
# For every cross-modal basic trial f -> g, adds a successive trial g -> h,
# where h is the third modality. The standard of the successive trial is the
# participant's final adjustment in the basic trial, which is only known
# during the experiment, so it is stored as successive_placeholder and
# replaced by PsychoPy at runtime.
#
# Successive trials are linked to their basic trial by the id prefix "2.".
add_successive_trials <- function(df, session, subj) {

  new_rows <- df %>%
    filter(prod_factor %in% successive_factors,
           standard %in% successive_standards,
           standard_modality != target_modality) %>%
    mutate(
      # Third modality, neither standard nor target of the basic trial
      third_modality = mapply(function(s, t) setdiff(modalities, c(s, t)),
                              standard_modality, target_modality),
      # Swap to the other production factor (both are 1 in this design)
      prod_factor    = if_else(prod_factor == successive_factors[1],
                               successive_factors[2],
                               successive_factors[1]),
      trial_type     = "successive",
      standard       = successive_placeholder,
      id             = paste0("2.", id)
    ) %>%
    # Target of the basic trial becomes the standard of the successive trial
    select(-standard_modality) %>%
    rename(standard_modality = target_modality,
           target_modality   = third_modality) %>%
    select(standard_modality, target_modality, standard,
           trial_type, id, prod_factor, start_target) %>%
    add_start_target(., subj = subj, session = session)

  rbind(df, new_rows)
}


# ── Target starting values ────────────────────────────────────────────────────
# Session 1: the starting intensity of the target is drawn at random from
# the modality's target range.
# Later sessions: it is drawn from +/- 5 around the mean match the
# participant produced for the closest standard in earlier sessions, read
# from session/subj_<id>/standard_match.txt.
add_start_target <- function(df, session, subj) {

  if (session == 1) {
    df <- df %>%
      mutate(start_target = case_when(
        target_modality == "auditory" ~ sample(aud_target_range, n(),
                                               replace = TRUE),
        target_modality == "visual"   ~ sample(vis_target_range, n(),
                                               replace = TRUE),
        target_modality == "tactile"  ~ sample(tac_target_range, n(),
                                               replace = TRUE),
        TRUE ~ NA_integer_
      ))

  } else {
    standard_match <- read.table(
      paste0("session/subj_", subj, "/standard_match.txt"),
      header = TRUE, as.is = TRUE
    ) %>%
      type.convert(as.is = TRUE)

    # Mean match per standard, modality pair and production factor
    standard_match <- aggregate(
      match ~ standard + target_modality + standard_modality + prod_factor,
      standard_match, mean
    )

    df <- df %>%
      mutate(start_target = round(mapply(
        find_closest_match,
        standard, target_modality, standard_modality, prod_factor,
        MoreArgs = list(standards           = standard_match$standard,
                        target_modalities   = standard_match$target_modality,
                        standard_modalities = standard_match$standard_modality,
                        prod_factors        = standard_match$prod_factor,
                        matches             = standard_match$match)
      )))
  }

  df
}


# Returns the starting value for one trial in sessions > 1.
#
# Successive trials have no fixed standard (successive_placeholder), so their
# starting value is drawn at random from the target range as in session 1.
# For all other trials, the match for the closest standard with the same
# modality pair and production factor is looked up (falling back to the same
# target modality if the pair is missing), a value within +/- 5 of it is
# drawn, and the result is clamped to the modality's adjustment limits.
find_closest_match <- function(x_standard, x_target_modality,
                               x_standard_modality, x_prod_factor,
                               standards, target_modalities,
                               standard_modalities, prod_factors, matches) {

  target_range <- switch(x_target_modality,
                         visual   = vis_target_range,
                         auditory = aud_target_range,
                         tactile  = tac_target_range)
  limits       <- switch(x_target_modality,
                         visual   = vis_limits,
                         auditory = aud_limits,
                         tactile  = tac_limits)

  if (x_standard == successive_placeholder) {
    return(sample(target_range, size = 1))
  }

  same_pair_idx <- which(target_modalities   == x_target_modality &
                         standard_modalities == x_standard_modality &
                         prod_factors        == x_prod_factor)

  if (length(same_pair_idx) == 0) {
    same_pair_idx <- which(target_modalities == x_target_modality &
                           prod_factors      == x_prod_factor)
  }

  sub_standards <- standards[same_pair_idx]
  sub_matches   <- matches[same_pair_idx]
  closest_match <- sub_matches[which.min(abs(sub_standards - x_standard))]

  out <- sample((closest_match - 5):(closest_match + 5), size = 1)
  max(min(out, limits[2]), limits[1])
}


# ── Trial order ───────────────────────────────────────────────────────────────
# After shuffling, a successive trial may come before the basic trial it
# depends on. Such pairs are swapped so the basic trial is always first.
reorder_pairs <- function(df) {

  for (i in 1:(nrow(df) - 1)) {
    current_id <- df$id[i]

    if (grepl("^2\\.", current_id)) {
      pair_idx <- which(df$id == sub("^2\\.", "", current_id))

      if (length(pair_idx) > 0 && pair_idx > i) {
        temp          <- df[i, ]
        df[i, ]       <- df[pair_idx, ]
        df[pair_idx, ] <- temp
      }
    }
  }

  df
}
