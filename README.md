# Cross-Modal Transitivity: Experiment and Analysis Code

Experiment and analysis code for a psychophysics experiment testing the transitivity of cross-modal magnitude production for loudness, brightness and vibration strength, as reported in:

> Kohler, D., & Heller, J. (in preparation). *Do internal references determine transitivity? A cross-modal study on loudness, brightness, and vibration strength.*

Supplementary material is archived on OSF: [osf.io/8f7ed](https://osf.io/8f7ed/)

---

## Background

If a tone is matched to a light, and a vibration is then matched to that tone, is the vibration equally strong as one matched to the light directly? Transitivity requires that the indirect match x<sub>fgh</sub> (from modality f via g to h) equals the direct match x<sub>fh</sub>. With three modalities there are six such chains:

| Indirect match | Direct match |
| -------------- | ------------ |
| l → b → s      | l → s        |
| l → s → b      | l → b        |
| b → l → s      | b → s        |
| b → s → l      | b → l        |
| s → l → b      | s → b        |
| s → b → l      | s → l        |

l = loudness, b = brightness, s = vibration strength

In the global psychophysical theory of Heller (2021), each match depends on internal references ρ that may differ with the role of the stimulus (standard or target) and the modality it is paired with. Transitivity follows when three invariances hold for these references:

- **s-invariance**: the reference of the standard does not depend on the modality of the target (ρ<sub>f→g</sub> = ρ<sub>f→h</sub>)
- **v-invariance**: the reference of the target does not depend on the modality of the standard (ρ<sub>h←f</sub> = ρ<sub>h←g</sub>)
- **t-invariance**: the reference of the middle modality is the same as target and as standard (ρ<sub>g←f</sub> = ρ<sub>g→h</sub>)

The experiment measures direct and indirect matches for all six chains, and a Bayesian model estimates the internal references of each participant to test the three invariances.

---

## Repository Structure

```
├── experiment/                        # PsychoPy experiment code, see experiment/README.md
│   ├── run_experiment.py              # Runs one session
│   ├── global_setup.py                # Window, stimuli, parameters, tutorial
│   ├── trialFunction.py               # A single matching trial
│   ├── conversions.py                 # Visual units: visualscale ↔ RGB ↔ dB Lambert
│   ├── vib_exp_v2_successive.py       # Connects to the tactor controller
│   ├── lookuptable/
│   │   └── lookuptab_mavo.txt         # Monitor calibration
│   └── session/                       # Creation of the session files (R)
│
└── analysis/
    ├── data_transitivity.csv          # Trial-level data of all participants
    ├── get_analysis_data.R            # Loads data, converts to dB, pairs direct/indirect matches
    ├── gain_to_decibel.R              # Tactor gain → dB displacement
    ├── matching_gain_to_displacement_250_Hz.txt
    ├── get_meta_deskriptives.R        # Table A3, key press frequencies
    ├── prod_transitivity_plots_descriptive.R              # Figures 4, A8–A14
    ├── prod_transitivity_ci_diff_table.R                  # Table A2 (counts in Table 2)
    ├── prod_overview_invariances_diff_per_iteration.R     # Figure 5
    ├── prod_overview_invariances_diff_distributions_02_05.R  # Figure 6
    ├── prod_rho_comparisons_tables.R  # Directional comparisons of the references
    │
    ├── pre_data/                      # Before data collection
    │   ├── simulate_transitivity_violations_s_invariance.R      # Figure 3 panels
    │   ├── simulate_transitivity_violations_v_invariance.R
    │   ├── simulate_transitivity_violations_t_invariance.R
    │   ├── simulate_transitivity_violations_omega_1_ist_gleich_1.R
    │   ├── prod_invariance_simulations.R                       # Figure 3, Table A1
    │   └── samplesize_calculations_confidence_interval.R       # Sample size (CI width)
    │
    └── model_building/                # Bayesian model, validation and fit
        ├── gpm_generic.stan           # Bayesian global psychophysics model
        ├── prior_transitivity.csv     # Priors
        ├── displacement_to_decibel.R  # Gain ↔ displacement fit for the tactor
        ├── prod_overview_priors.R                       # Figure A7
        ├── prod_prior_predictive_transitivity.R         # Figures A1, A2
        ├── prod_simulation_based_calibrations.R         # Figure A3
        ├── prod_model_sensitivity.R                     # Figure A4
        ├── prod_rds_file_parameter_estimation.R         # Model fit per participant
        └── prod_posterior_predictive_transitivity.R     # Figures A5, A6
```

---

## Design

Each block contains 19 trials in random order:

| Trials | p | Modalities |
| ------ | - | ---------- |
| 7 basic | 1 | all six cross-modal pairs plus brightness → brightness |
| 6 successive | 1 | target of each cross-modal basic trial → third modality |
| 6 basic | 2 | all six cross-modal pairs |

In a successive trial, the standard is the participant's final adjustment in the basic trial f → g of the same block, so together the two trials give the indirect match x<sub>fgh</sub>. The direct match x<sub>fh</sub> is the basic trial f → h of the same block. The trials with p = 2 make the internal references identifiable, and the brightness → brightness trials allow a comparison of cross-modal and intra-modal references.

---

## Analysis

All analysis scripts are written in R and run with `analysis/` as working directory, each in a fresh R session, e.g. `source("model_building/prod_model_sensitivity.R")`. Scripts that use data read them through `get_analysis_data.R`. Results are written to `analysis/output/`; the simulations in `pre_data/` save their figures to `analysis/`.

1. **Simulations (`pre_data/`):** the four `simulate_transitivity_violations_*` scripts show how a violation of s-, v-, t-invariance or W(1) = 1 affects transitivity in the chain l → s → b; `prod_invariance_simulations.R` combines them into Figure 3 and writes the parameter table (Table A1). `samplesize_calculations_confidence_interval.R` simulates the width of the confidence interval of the matches for different numbers of trials.
2. **Model validation (`model_building/`):** prior overview (Figure A7), prior predictive checks (Figures A1, A2), simulation-based calibration (Figure A3) and model sensitivity (Figure A4). These scripts simulate data with the trial structure of one session (participant 01, session 01 in `data_transitivity.csv`). SBC and sensitivity save every fit to a cumulative RDS file and resume where they stopped.
3. **Model fit:** `prod_rds_file_parameter_estimation.R` fits `gpm_generic.stan` separately to the data of each participant and saves the posterior draws of all references and ω<sub>1</sub> (`output/rho_samples_<subj>.csv`) and the fits (`output/gpm_fits.rds`). `prod_posterior_predictive_transitivity.R` uses the fits for the posterior predictive checks (Figures A5, A6).
4. **Transitivity and invariances:** the `prod_*` scripts in `analysis/` compare direct and indirect matches and the posterior differences of the references for s-, v- and t-invariance, and check W(1) = 1. Transitivity holds when the 95% interval of the mean difference between direct and indirect match includes zero; an invariance holds when the 95% posterior interval of the difference between the two references includes zero; W(1) = 1 holds when the 95% posterior interval of ω<sub>1</sub> includes 1.
5. **Descriptives:** `get_meta_deskriptives.R` creates the trial-level descriptives per participant (Table A3).

### Pairing of direct and indirect matches

A successive trial g → h does not carry the start modality f of its chain. `get_transitivity_pairs()` in `get_analysis_data.R` therefore links each successive trial `<block>.2.<n>` to its basic trial `<block>.1.<n>` to find f, and pairs it with the direct match f → h of the same block. All transitivity scripts use this function.

### Key R packages

```
dplyr, tidyr, purrr, tidyverse, glue  # Data wrangling
ggplot2, ggh4x, ggridges, patchwork, latex2exp, bayesplot  # Visualization
knitr, kableExtra, xtable            # LaTeX tables
rstan, tictoc           # Bayesian model
globalpsychophysics     # simulate_gpm(), make_datlist_generic(), github.com/Kaanwoj/globalpsychophysics
```

---

## Data

`analysis/data_transitivity.csv` contains all trials of all participants (training trials included) with pseudonymised participant IDs 01–08. Participant 06 is excluded from all analyses because they made only one adjustment in more than 50% of trials. `get_analysis_data()` excludes 06 by default; the model is fitted to all eight participants.

| Column | Content |
| ------ | ------- |
| `subj` | participant ID |
| `session`, `run` | session number, run number per participant |
| `id` | trial id `<block>.<step>.<trial>`, step 1 = basic, 2 = successive |
| `block` | block number or `training` |
| `trial_type` | `basic` or `successive` |
| `standard_modality`, `target_modality` | `visual`, `auditory` or `tactile` |
| `standard` | intensity of the standard |
| `start` | starting intensity of the target |
| `match` | final adjustment |
| `p` | production factor |
| `keys` | all keys pressed in the trial |
| `rt` | time until confirmation (s) |
| `exptime` | time since start of the session (s) |
| `lightRGB` | RGB value of the circle at the end of the trial |

All intensities are in device units of the experiment (see `experiment/README.md`). `get_analysis_data()` converts them to dB: visual via the monitor lookup table to dB Lambert, tactile via `gain_to_decibel()` to dB displacement; auditory values are used unchanged.

---

## Experiment Code

The experiment was run in PsychoPy 2024 with an OLED monitor, headphones and a C-2 tactor (Engineering Acoustics, Inc.) at the University of Tübingen. Setup, session flow and output files are described in [`experiment/README.md`](experiment/README.md).

---

## Citation

If you use this code, please cite:

```
Kohler, D., & Heller, J. (in preparation). Do internal references determine
transitivity? A cross-modal study on loudness, brightness, and vibration
strength.
```

---

## License

MIT, see [LICENSE](LICENSE)
