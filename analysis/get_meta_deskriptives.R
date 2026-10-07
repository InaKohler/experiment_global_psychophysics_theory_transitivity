# get_meta_deskriptives.R
#
# Trial-level descriptives per participant (Table A3) and key press
# frequencies per participant and target modality:
#   percentage of final adjustments at ceiling/floor, percentage of trials
#   with only one adjustment, mean number of adjustments, mean session
#   duration (minutes), frequency of each key
#
# The key-percentages table is split across two pages (loudness and
# brightness on the first, vibration on the second, separated by \clearpage)
# so it doesn't overflow the page.
#
# Input:  data via get_analysis_data() (or `dat` if it already exists)
# Output: output/trial_descriptives.tex
#         output/key_percentages.tex
#
# last mod: 2026/10/01

library(dplyr)
library(knitr)
library(kableExtra)

source("get_analysis_data.R")
output_table <- "output"

if (!exists("dat")) {
  dat <- get_analysis_data()
}


# ---------------------------------------------------------------------------
upper_limit_spl         <- 85   # auditory, raw `match` scale
lower_limit_spl         <- 15
upper_limit_visualscale <- 49   # visual, raw `match` scale
lower_limit_visualscale <- 0
upper_limit_tactile     <- 81.6 # tactile, dB scale -> checked against matchDB
lower_limit_tactile     <- 0

## Subject numbers, filter out subject 06
dat <- dat %>%
  mutate(subj_number = subj_id) %>%
  filter(!is.na(subj_number), subj_number != "06")

## Trial-specific descriptives
dat <- dat %>%
  mutate(
    ceiling = case_when(
      target_modality == "visual"   ~ match == upper_limit_visualscale,
      target_modality == "auditory" ~ match == upper_limit_spl,
      target_modality == "tactile"  ~ matchDB == upper_limit_tactile,
      TRUE ~ NA
    ),
    floor = case_when(
      target_modality == "visual"   ~ match == lower_limit_visualscale,
      target_modality == "auditory" ~ match == lower_limit_spl,
      target_modality == "tactile"  ~ matchDB == lower_limit_tactile,
      TRUE ~ NA
    )
  )
dat$minadjust <- nchar(dat$keys) == 6   # 6 characters = "space" + 1
dat$nadjust <- nchar(dat$keys) - 5      # remove 5 characters "space"

dftrial <- sapply(split(dat, dat$subj_number), \(subj) {
  c(
    Subject = subj$subj_number |> unique(),
    percceiling = 100 * mean(subj$ceiling) |> round(2),
    percfloor = 100 * mean(subj$floor) |> round(4),
    perconeadj = 100 * mean(subj$minadjust) |> round(2),
    meanadj = paste0(mean(subj$nadjust) |> round(2), " (",
                     sd(subj$nadjust) |> round(2), ")"),
    meansestime = paste0(mean(
      sapply(split(subj, subj$session), \(ses)
             diff(range(ses$exptime)) / 60)) |> round(2), " (",
      sd(sapply(split(subj, subj$session), \(ses)
                diff(range(ses$exptime)) / 60)) |> round(2), ")")
  )
}) |> t() |> as.data.frame()

make_header <- function(top, bottom) {
  paste0("\\makecell[c]{", top, "\\\\", bottom, "}")
}

col_names <- c(
  "Subject",
  make_header("Ceiling", "(in \\%)"),
  make_header("Floor", "(in \\%)"),
  make_header("One Adjustment", "(in \\%)"),
  make_header("Number of", "Adjustment (M, SD)"),
  make_header("Session Duration", "(M, SD)")
)

writeLines(
  kbl(dftrial,
      format = "latex",
      row.names = FALSE,
      col.names = col_names,
      align = c("l", "r", "l", "c", "c", "c"),
      booktabs = TRUE,
      escape = FALSE,
      caption = "Trial-Level Descriptive Statistics by Subject",
      label = "trial_descriptives") |>
    footnote(
      general = "Ceiling and floor indicate the percentage of trials in which the final adjustment reached the upper or lower limit of the respective modality's observable range. One Adjustment indicates the percentage of trials in which the participant confirmed their initial value after one mandatory adjustement. Number of Adjustments and Session Duration (minutes) are given as mean (SD).",
      general_title = "",
      threeparttable = TRUE,
      escape = FALSE
    ),
  file.path(output_table, "trial_descriptives.tex")
)

## Key-specific descriptives
# Key mapping to their step size
# down: small (a), medium (s), large (d)
# up:   small (w), medium (e), large (r)
keymap <- c(a = "down:small (a)", s = "down:medium (s)", d = "down:large (d)",
            w = "up:small (w)",   e = "up:medium (e)",   r = "up:large (r)")

modalities <- c("auditory", "visual", "tactile")

dfkey <- sapply(split(dat, dat$subj_number), \(subj) {
  keyprop <- sapply(split(subj, subj$target_modality), \(mod) {
    100 * (strsplit(mod$keys, "space") |> unlist() |> strsplit("") |>
             unlist() |> table() |> prop.table() |>
             _[c("a", "s", "d", "r", "e", "w")] |> round(2))
  })
  unlist(lapply(modalities, \(m)
                setNames(keyprop[, m], paste0(m, ":", keymap[rownames(keyprop)]))))
}) |> as.data.frame()

dfkey$"Step Size (key)" <- sub("^([^:]*:){2}", "", rownames(dfkey))
dfkey <- dfkey[, c(ncol(dfkey), 1:(ncol(dfkey) - 1))]

# Split by modality rows so the table can be stretched across two pages:
# Auditory + Visual (rows 1-12) on page 1, Tactile (rows 13-18) on page 2.
n_subj <- ncol(dfkey) - 1

make_key_table <- function(df_part, mod_groups, caption_suffix, label_suffix,
                           add_footnote) {
  tbl <- kbl(df_part, format = "latex", row.names = FALSE, booktabs = TRUE,
             escape = FALSE,
             caption = paste0("Key Press Frequencies by Step Size, Modality, and Subject",
                              caption_suffix),
             label = paste0("key_percentages", label_suffix)) |>
    add_header_above(c(" " = 1, "Subject" = n_subj))
  
  for (grp in mod_groups) {
    tbl <- tbl |>
      pack_rows(grp$name, grp$start, grp$end) |>
      pack_rows("down", grp$start,     grp$start + 2) |>
      pack_rows("up",   grp$start + 3, grp$end)
  }
  
  if (add_footnote) {
    tbl <- tbl |>
      footnote(
        general = "Values are percentages of key presses per subject and target modality, out of all adjustments in that modality. Step sizes are ordered from largest (top) to smallest (bottom) within each direction (down, up).",
        general_title = "",
        threeparttable = TRUE,
        escape = FALSE
      )
  }
  tbl
}

table_part1 <- make_key_table(
  dfkey[1:12, ],
  list(list(name = "Loudness",   start = 1, end = 6),
       list(name = "Brightness", start = 7, end = 12)),
  "", "", FALSE
)
table_part2 <- make_key_table(
  dfkey[13:18, ],
  list(list(name = "Vibration", start = 1, end = 6)),
  " (continued)", "_cont", TRUE
)

writeLines(
  c(table_part1, "\\clearpage", table_part2),
  file.path(output_table, "key_percentages.tex")
)
