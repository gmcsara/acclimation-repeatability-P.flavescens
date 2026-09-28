# README — RawData21.xlsx and BP_Analysis_Complete.R

Data and code accompanying the manuscript:

**Effect of acclimation duration on the repeatability of boldness in *Pomatoschistus flavescens***
Martins-Cardoso, S. & Faria, A.M.
Submitted to *Behavioural Processes*.

## General information

| | |
|---|---|
| Species | *Pomatoschistus flavescens* (two-spotted goby) |
| Life stage | Wild-caught adults |
| Collection site | Arrábida Marine Park, Portugal (38°28'N; 8°59'W) |
| Collection date | April 2021 |
| N individuals | 31 (17 males, 14 females); 2 min: n = 10; 5 min: n = 11; 15 min: n = 10 |
| N observations | 93 (31 individuals × 3 trials) |

## Files

- **`RawData21.xlsx`** — raw behavioural data (single sheet, 1 header row + 93 observations, 10 columns).
- **`BP_Analysis_Complete.R`** — complete analysis pipeline. Place it in the same folder as `RawData21.xlsx` and run it in R. Tables are written as `.csv` files and figures to the `Figures_Behavioural_Processes` subfolder (created automatically).

Required R packages: `tidyverse`, `readxl`, `lme4`, `lmerTest`, `rptR`, `survival`, `coxme`, `car`, `patchwork`.

Note: the bootstrap and power simulations are computationally intensive; a full run can take more than 30 minutes.

## Variable descriptions (`RawData21.xlsx`)

| Variable | Description |
|---|---|
| `ind_ID` | Individual identifier (integer). Each fish was housed in a uniquely numbered compartment, enabling identification without marking. |
| `acclimation_time` | Acclimation treatment group (2, 5, or 15 minutes). Each individual was assigned to one treatment and tested under that condition in all trials. |
| `trial_number` | Trial number (1, 2, or 3). Trials 1 and 2 were conducted one week apart; Trial 3 was conducted approximately six weeks after Trial 1. |
| `Sex` | Sex of the individual (1 = male, 2 = female), determined by the number of diagnostic lateral spots (males: two spots; females: one spot; Miller, 1986). |
| `SL` | Standard length (cm), measured before the start of the experiment. |
| `Latency_shelter` | Latency to leave the shelter (seconds): time from door opening until the fish first entered the test arena, as recorded in BORIS. Values close to 600 s correspond to fish that did not emerge (see *Censoring* below). |
| `Shelter` | Time spent in the shelter compartment (seconds) after first emergence. |
| `TB1` | Time spent in zone 1 (seconds): the 6 cm strip adjacent to the walls of the test arena. |
| `TB2` | Time spent in zone 2 (seconds): the intermediate 6 cm strip. |
| `TB3` | Time spent in zone 3 (seconds): the central, most exposed area of the test arena. |

`Shelter`, `TB1`, `TB2` and `TB3` are only scored after the first emergence; fish that never emerged have 0 in all four.

## Variables used in the analysis (built in the script)

**Latency to exit shelter.** Analysed with a Cox proportional hazards mixed model (`coxme`), with right-censoring at 600 s.

*Censoring:* a trial is treated as an emergence only if the fish left the shelter before 600 s **and** spent time in at least one zone afterwards. In 15 trials the video ended a fraction of a second before 600 s (recorded latency 599.3–599.9 s) with zero time in all zones; these are non-emergences and are censored at 600 s. Genuine emergences reach at most 524.7 s. In total, 52 of 93 trials (56%) are censored. Repeatability is reported as the ICC on the latent scale (McCune et al., 2025).

**Open-zone use.** (`TB2` + `TB3`) / 600, logit-transformed (floor/ceiling adj = 0.001).

**Near-wall use.** (`Shelter` + `TB1`) / 600, logit-transformed (adj = 0.001).

Space-use repeatability was estimated with `rptR` (Gaussian LMM), with trial number and standardised SL as fixed effects and individual ID as a random effect. Sensitivity to the logit adjustment (adj = 0.001, 0.01, 0.05) is reported in Figure S1.

Fixed effects (trial number and SL) are taken from the same models used for the variance components (Table 4) and are written to `Table_fixed_effects.csv`.
