# Egg Neglect Behaviour in Little Auk (*Alle alle*) — Analysis Repository

This repository contains the R scripts and data files used to analyse
egg neglect behaviour during incubation in the Little Auk, a small
Arctic seabird breeding in Hornsund, Svalbard. The analyses are
organised into two stages: data preparation scripts (prefixed `DP`) and
analysis scripts (prefixed `EDA`). The main results are compiled into a
single R Markdown report (`egg_neglect_behaviour.Rmd`).

---

## How to Run

Run scripts in this order:

```
DP1 → DP2 → DP3 → EDA scripts → Rmd report
```

All scripts assume the working directory is set to the folder containing
the raw data files. Raw input files (`.xlsx`, `.csv`) are read from
absolute paths — **update these paths** to match your local setup before
running.

---

## Stage 1 — Data Preparation Scripts

### `DP1` — Nest-level egg neglect index

**Input:** `DPA_d05_status_activity data.rds`,
`SHO_LIAK_phenology_productivity.xlsx`, `SI_PS.csv`

**What it does:**
- Filters activity data to incubation sessions only
  (`incubation1`, `incubation2`, `incubation3`)
- Excludes known problematic nests (hardcoded list)
- Detects and calculates **egg neglect gaps** — periods when no parent
  was present in the nest — by merging overlapping attendance intervals
  and finding spaces within the session window
- Calculates **negotiation overlaps** — periods when both partners were
  simultaneously present in the nest
- Joins phenology data (hatch date, lay day) to compute days prior to
  hatch (`day_prior_hatch`)
- Joins pair status and similarity data from `SI_PS.csv`
- Creates hourly grid of exits per nest for diurnal analysis
- Generates **Gantt verification plots** per nest saved to `VERIFICATION/`

**Outputs:**
- `EDA_egg_neglect_index_gaps_table.RDS` — raw gap list per nest/session
- `EDA_egg_neglect_index_overlaps_table.RDS` — raw overlap list
- `EDA_egg_neglect_index_neglect_table.RDS` — **main nest-level analysis
  table** (one row per season × session × nest)
- `EDA_exit_hourly.RDS` — hourly exit counts per nest for diurnal models
- `VERIFICATION/*.jpg` — Gantt plots for visual verification of gaps

---

### `DP2` — Sex-separated egg neglect index

**Input:** `DPA_d05_status_activity data.rds`,
`SHO_LIAK_phenology_productivity.xlsx`, `SI_PS.csv`

**What it does:**
- Applies identical filtering and exclusions as DP1
- Fixes missing sex labels using nest-level context
- Detects gaps using the same merging logic as DP1
- **Attributes each gap to a responsible sex** using a tiered decision
  rule:
  - *Beginning rule*: gaps at session start attributed to the partner
    of whoever arrives first
  - *Short gaps* (< 1 hour): attributed based on who left and who
    arrived around the gap
  - *Long gaps* (> 1 hour): attributed using cumulative workload and
    dynamic shift-length thresholds; hard-coded overrides for known
    problem nests
- Summarises neglect per sex per session (frequency, total duration,
  mean duration)
- Creates sex-separated wide table (one row per season × session × nest
  × sex) with ring numbers, phenology, pair status and similarity
- Calculates nest-level neglect bias (male vs female contribution)
- Generates **sex-separated Gantt verification plots** saved to
  `VERIFICATION_SEX/`

**Outputs:**
- `EDA_sex_neglect_wide.RDS` — **main sex-level analysis table**
  (one row per season × session × nest × sex)
- `EDA_egg_neglect_index_sex.RDS` — nest-level summary with sex-specific
  neglect columns
- `EDA_sex_gaps_hourly.RDS` — hourly gap records with responsible sex
  attributed, used for diurnal models
- `VERIFICATION_SEX/*.jpg` — sex-separated Gantt plots

---

### `DP3` — Pair similarity and pair status

**Input:** `SHO_LIAK_capturing.xlsx`, `CapturingDt2025.xlsx`

**What it does:**
- Loads and combines capturing data across all seasons including 2025
- Filters out first-year birds (`Age != "1"`)
- Fills missing wing/head/tarsus measurements within individuals across
  years using `fill()`
- Calculates **three morphological similarity indices** for each pair
  using a bootstrap-normalised absolute difference approach:
  - `Similarity_Index` — based on wing length
  - `Similarity_Index_Head` — based on total head length (THL)
  - `Similarity_Index_Tarsus` — based on tarsus length
  - All indices scaled 0–1, where 1 = identical morphology
- Determines **pair status** (New vs Old) based on whether partners
  appeared together in previous seasons at the same nest

**Output:**
- `SI_PS.csv` — pair similarity indices and pair status per nest per
  season, used as input by DP1 and DP2

---

## Stage 2 — Analysis Scripts (EDA) and R Markdown Report

All EDA analyses are compiled in `egg_neglect_behaviour.Rmd`. The
sections are:

| Section | Research question | Data used | Method |
|---------|------------------|-----------|--------|
| Diurnal pattern of exits | Is there a sex-specific daily rhythm in nest exits? | `EDA_exit_hourly.RDS`, `EDA_sex_neglect_wide.RDS`, `EDA_sex_gaps_hourly.RDS` | Bayesian Poisson GLM with B-splines (`quap`) |
| Neglect across incubation | How do gap frequency, probability, and duration change toward hatching, by sex? | `EDA_sex_neglect_wide.RDS` | Poisson GLM + Binomial GLM + Normal GLM on log scale (`quap`) |
| Pair status & similarity | Does pair experience or morphological similarity predict neglect? | `EDA_egg_neglect_index_neglect_table.RDS` | Poisson GLM + Normal GLM on log scale (`quap`) |
| Seasonal differences | Does neglect vary across breeding seasons 2019–2025? | `EDA_egg_neglect_index_neglect_table.RDS` | Poisson + Binomial + Normal GLMs with season intercepts (`quap`) |
| Repeatability | Are individuals consistent in neglect tendency across sessions? | `EDA_sex_neglect_wide.RDS` | Binomial repeatability (`rptR`) |
| Effect on incubation duration | Does total neglect predict incubation length? | `EDA_egg_neglect_index_neglect_table.RDS` | Normal GLM with log-transformed predictor (`quap`) |

To reproduce the report, knit `egg_neglect_behaviour.Rmd` in RStudio
after ensuring all `.RDS` files are present in the working directory.

---

## Data Files — Variable Definitions

### `EDA_egg_neglect_index_neglect_table.RDS`
One row per season × session × nest.

| Variable | Description |
|----------|-------------|
| `season` | Year of study (character) |
| `session` | Incubation session (`incubation1`, `incubation2`, `incubation3`) |
| `nest` | Nest identity code |
| `n_gaps` | Number of neglect gaps (both sexes combined) |
| `mean_gaps` | Mean gap duration in seconds (gaps > 0 only) |
| `sum_gaps` | Total neglect duration in seconds |
| `n_overlaps` | Number of negotiation overlaps (both partners present) |
| `mean_overlaps` | Mean overlap duration in seconds |
| `sum_overlaps` | Total overlap duration in seconds |
| `session_start` | Start datetime of the 48h recording window |
| `session_end_trimmed` | End datetime of the recording window |
| `mean_date` | Midpoint datetime of the recording window |
| `Julian_Date_point` | Julian day of the recording midpoint |
| `hatch_date` | Hatching date |
| `hatch_succ` | Hatching success (1 = hatched) |
| `Lay_day` | Julian day of laying |
| `J_hatchday` | Julian day of hatching |
| `Incubation_Duration` | Incubation duration in days (`J_hatchday - Lay_day`) |
| `day_prior_hatch` | Days before hatching at time of recording |
| `ID1`, `ID2` | Ring numbers of the two pair members |
| `Pair_status` | Pair experience (`New` or `Old`) |
| `Similarity_Index` | Wing-based morphological similarity (0–1) |
| `Similarity_Index_Head` | Head-based morphological similarity (0–1) |
| `Similarity_Index_Tarsus` | Tarsus-based morphological similarity (0–1) |

---

### `EDA_sex_neglect_wide.RDS`
One row per season × session × nest × sex.

| Variable | Description |
|----------|-------------|
| `season` | Year of study |
| `session` | Incubation session |
| `nest` | Nest identity |
| `sx` | Sex (`f` = female, `m` = male) |
| `ringno` | Individual ring number |
| `n_gaps` | Number of gaps attributed to this individual |
| `sum_neglect_sec` | Total neglect duration attributed to this individual (seconds) |
| `mean_neglect_sec` | Mean gap duration attributed to this individual (seconds) |
| `n_incubations` | Number of incubation bouts recorded for this individual |
| `mean_incubation_duration` | Mean incubation bout duration (seconds) |
| `day_prior_hatch` | Days before hatching at time of recording |
| `hatch_date` | Hatching date |
| `Lay_day` | Julian day of laying |
| `J_hatchday` | Julian day of hatching |
| `Incubation_Duration` | Incubation duration in days |
| `Julian_Date_point` | Julian day of recording midpoint |
| `ID1`, `ID2` | Ring numbers of both pair members |
| `Pair_status` | Pair experience (`New` or `Old`) |
| `Similarity_Index` | Wing-based morphological similarity (0–1) |
| `Similarity_Index_Head` | Head-based morphological similarity (0–1) |
| `Similarity_Index_Tarsus` | Tarsus-based morphological similarity (0–1) |

---

### `EDA_exit_hourly.RDS`
One row per season × session × nest × date × hour.

| Variable | Description |
|----------|-------------|
| `season` | Year of study |
| `session` | Incubation session |
| `nest` | Nest identity |
| `date` | Calendar date of the hour |
| `hour` | Hour of day (0–23) |
| `n_gaps_per_hour` | Number of nest exits in that hour (0 if none) |
| `day_prior_hatch` | Days before hatching |
| `Pair_status` | Pair experience |
| `Similarity_Index` | Wing-based morphological similarity |
| `ID1`, `ID2` | Ring numbers of pair members |

---

### `EDA_sex_gaps_hourly.RDS`
One row per detected gap, with hour and responsible sex attributed.

| Variable | Description |
|----------|-------------|
| `season` | Year of study |
| `session` | Incubation session |
| `nest` | Nest identity |
| `gap_start` | Start datetime of the gap |
| `gap_end` | End datetime of the gap |
| `dur` | Gap duration in seconds |
| `responsible_sex` | Sex attributed as responsible for the gap (`f` or `m`) |
| `hour` | Hour of day when gap started (0–23) |

---

## Verification Outputs

Two folders of Gantt plots are generated during data preparation for
manual quality control.

| Folder | Contents |
|--------|----------|
| `VERIFICATION/` | One plot per nest — green = nest occupied, red shading = egg neglect gap, with duration labels |
| `VERIFICATION_SEX/` | One plot per nest — green = nest attended, blue = gap attributed to male, red = gap attributed to female |

---

## Software Requirements

```r
library(rethinking)   # Bayesian model fitting — requires cmdstanr or rstan
library(rptR)         # Repeatability estimation
library(tidyverse)    # Data manipulation and plotting
library(lubridate)    # Datetime handling
library(splines)      # B-spline basis functions
library(patchwork)    # Plot composition
library(readxl)       # Reading Excel files
library(plyr)         # ldply for list-to-dataframe conversion
library(zoo)          # fill() for morphometric data
library(knitr)        # Table rendering in Rmd
```

> Installing `rethinking`: see https://github.com/rmcelreath/rethinking

---


- **Start with the Rmd report** (`egg_neglect_behaviour.Rmd`) for a
  full narrative of the analyses with inline results and plots
- **Check the VERIFICATION plots** for a sample of nests to confirm
  gap detection looks correct, and **VERIFICATION_SEX plots** to verify
  sex attribution of gaps
- The sex attribution algorithm in DP2 uses a tiered logic — the most
  complex cases (long gaps, ambiguous sex) have hard-coded overrides
  for specific nests that were manually verified.
- The `Similarity_Index` used in pair similarity analyses is based on
  **wing length** only; head and tarsus indices are also available in
  the data but were not the primary predictor
- All Bayesian models use **quadratic approximation** (`quap`) from the
  `rethinking` package 
