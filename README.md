# IEEE S&P 2026 Submission: “In Grid We Trust: Electric Network Frequency Signatures for Chip Geolocation”

## 1. Introduction

This repository contains the artifacts for the paper titled **“In Grid We Trust: Electric Network Frequency Signatures for Chip Geolocation,”** submitted to the **47th IEEE Symposium on Security and Privacy (IEEE S&P 2026)**. The artifacts provided here allow for independent verification of our main experimental results on extracting and using Electric Network Frequency (ENF) signatures from the local environment of **DC-powered PCBs / FPGA boards**.

The repository includes:

- The raw ENF datasets used in the paper:
  - Baseline ambient EM feasibility measurements.
  - FPGA-based ambient EM traces.
  - Temporal reliability measurements at a single grid location (Experiment 3).
  - Server-room robustness measurements (Experiment 4).
  - Multi-location cross-grid validation measurements (Experiment 5).
- MATLAB scripts that implement the ENF extraction, correlation analysis, and statistical summaries.
- A compiled MATLAB analysis core (`proc_enf_analysis.p`) used by the top-level scripts.

*Note:* Our core analysis function, originally developed under authors' university policy that restricts source release, is provided here as a compiled MATLAB function (`proc_enf_analysis.p`) rather than in source form. The top-level scripts in this repository call that compiled function so that all analyses in the paper can be reproduced without exposing proprietary internals. 


## 2. Folder and File Hierarchy

At a high level, the repository is organized as follows:

- **`exp_inputs/`** – All ENF input traces (ambient EM + mains references) used in the experiments.

  - **`BASE/`** – Baseline ambient EM feasibility experiment (no FPGA workload / simple setup).
    - **`NO_FB/`**  
      - `fpga_em_trace_dc.wav` – Ambient EM trace captured with only picoscope.  
      - `mains_pow_trace_ac.wav` – Ground-truth mains reference trace.
    - **`W_FB/`**  
      - Same file naming as `NO_FB`, but with the setup placed inside a shielding (e.g., Faraday bag) to test isolation effects.

  - **`FPGA/`** – FPGA-board ambient EM experiment (pair-wise ENF feasibility).
    - **`SAKU/`**  
      - `fpga_em_trace_dc.wav` – Ambient EM near the Sakura-G (or equivalent) FPGA board.  
      - `mains_pow_trace_ac.wav` – Co-recorded mains reference.

  - **`MULTI/`** – **Experiment 5 – ENF Multi-Location Validation.**  
    ENF measurements comparing a single local EM trace against mains references from multiple grid regions.

    - **`US_60/`** – 60 Hz (US) multi-location dataset.
      - **`AUG/`**, **`OCT/`** – Months of collection.
        - `MON/`, `WED/`, `TUE/`, etc. – Calendar day.
          - `T01/`, `T02/`, `T03/` – Measurement folders (each one 200-s trace pair set).
            - `fpga_em_trace_dc_egrid_citya_lab.wav` – Local ambient EM trace near the FPGA.
            - `mains_pow_trace_ac_egrid_citya_lab.wav` – Local mains reference (same grid as EM sensor).
            - `mains_pow_trace_ac_egrid_citya_home.wav` – Additional local mains (home) reference.
            - `mains_pow_trace_ac_egrid_worcester.wav` – US Eastern grid reference.
            - `mains_pow_trace_ac_tgrid_richardson.wav` – US Texas grid reference.
            - `mains_pow_trace_ac_wgrid_tucson.wav` – US Western grid reference.

    - **`DE_50/`** – 50 Hz (European) multi-location dataset.
      - **`AUG/`, `OCT/`** – Months of collection.
        - `WED/`, `THU/`, `TUE/`, … – Calendar day.
          - `T01/`, `T02/`, `T03/` – Measurement folders.
            - `fpga_em_trace_dc_citya_lab.wav` – Local 50 Hz EM trace near FPGA.
            - `mains_pow_trace_ac_citya_lab.wav` – Local 50 Hz mains reference.
            - `mains_pow_trace_ac_dresden.wav` – Remote European mains reference.

    Each `Txx` folder thus contains one EM trace and multiple mains references from different grid regions, enabling the cross-grid correlation matrices used in Experiment 5.

  - **`TREND/`** – **Experiment 3 – ENF Temporal Reliability.**  
    Designed to measure how stable ENF correlation is over **Week**, **Day-of-Week**, and **Time-of-Day** at a single grid location.

    - **`WK01/`, `WK02/`** – Two consecutive weeks.
      - **`WDAY/`**, **`WEND/`** – Weekday vs weekend partition.
        - **`WED/`, `THUR/`, `SAT/`, `SUN/`** – Specific days.
          - **`EMRN/`, `MORN/`, `AFTN/`, `EVEN/`** – Time-of-day slots:
            - EMRN ≈ early morning (4 AM).  
            - MORN ≈ morning (9 AM).  
            - AFTN ≈ afternoon (2 PM).  
            - EVEN ≈ evening (7 PM).
          - **`T01/`–`T05/`** – Five repeated trials per (Week, Day, Time-of-Day) condition.
            - `fpga_em_trace_dc.wav` – Ambient EM trace near the FPGA board.
            - `mains_pow_trace_ac.wav` – Mains reference trace.

    This design yields 2 weeks × 4 days × 4 times × 5 trials = 160 measurements.

  - **`SRV_L/`** – **Experiment 4 – ENF Server-Room Robustness.**  
    ENF measurements in a high-density computing environment, studying robustness to sensor placement on the server PSU and day-of-week.

    - **`SBOX/`**, **`SPSU/`** – Two sensor locations:
      - For each site:
        - **`MON/`, `WED/`** – Two days.
          - **`T01/`–`T05/`** – Five repeated trials per (Site, Day) condition.
            - `fpga_em_trace_dc.wav` – Ambient EM trace near server PSU/BOX.
            - `mains_pow_trace_ac.wav` – Co-recorded mains reference.

- **`exp_scripts/`** – MATLAB scripts required to run the ENF extraction, correlation, and statistical analysis.

  - `enf_analysis_top_pair_wise.m`  
    – Top-level script for basic pair-wise ENF signature analysis (for singular analysis).

  - `enf_analysis_top_reliability_stat.m`  
    – Top-level script for **Experiment 3 – Temporal Reliability** (TREND).

  - `enf_analysis_top_reliability_stat_server.m`  
    – Top-level script for **Experiment 4 – Server-Room Robustness** (SRV_L).

  - `enf_analysis_top_multi_loc_stat.m`  
    – Top-level script for **Experiment 5 – Multi-Location Validation** (MULTI).

  - `proc_enf_analysis.p`  
    – Compiled MATLAB function that implements the core ENF extraction and correlation pipeline shared across all scripts.

- **`README.md`** – This file.


## 3. Experiment Scripts Overview

This section summarizes what each top-level MATLAB script does and how it relates to the paper’s experiments.

### 3.1 Pair-Wise ENF Signature Analysis (Feasibility Pipeline)

**Script:** `enf_analysis_top_pair_wise.m`  
**Data roots:** `exp_inputs/*`

**Goal:** Validate the feasibility of extracting ENF from **DC-powered hardware** by comparing a sensed ambient EM trace against a mains reference.

**Overview:**

1. **Inputs**
   - `fpga_em_trace_dc*.wav`: Ambient EM trace captured near the FPGA/PCB.
   - `mains_pow_trace_ac*.wav`: Co-recorded mains reference trace.

2. **Spectrogram Generation**
   - Applies Short-Time Fourier Transform (STFT) to both traces to obtain time–frequency spectrograms around the nominal grid frequency and its harmonics.

3. **ENF Estimation**
   - Uses a weighted-harmonic ENF estimator to extract the instantaneous ENF trajectory from each spectrogram.

4. **Correlation & Outputs**
   - Computes the **Pearson correlation coefficient** between the sensed ENF and the mains ENF.
   - Optionally produces aligned ENF plots and summary correlation statistics for the baseline and FPGA experiments.

This script is the simplest entry point if you want to understand the core ENF extraction and correlation pipeline on a single pair of traces.

---

### 3.2 Experiment 3 – ENF Temporal Reliability

**Script:** `enf_analysis_top_reliability_stat.m`  
**Data root:** `exp_inputs/TREND/`

**Goal:** Quantify how stable the ENF correlation is over **time-of-day, day-of-week, and week-to-week** at a fixed grid location.

**Overview:**

1. **Inputs**
   - Walks the `TREND/WKxx/…` hierarchy described above.
   - For each (Week, Day, Time-of-Day, Trial) condition, loads:
     - `mains_pow_trace_ac.wav` – Ground-truth mains reference.
     - `fpga_em_trace_dc.wav` – Ambient EM trace near the board.

2. **ENF Extraction & Correlation**
   - Uses the same STFT+weighted-harmonic estimator as the pair-wise script.
   - Computes correlation **r** per file pair and applies the **Fisher transform** `z = atanh(r)` for parametric inference.

3. **Statistics**
   - Computes descriptive statistics for **r** and **z** (N, mean, SD, SEM, 95% CI).
   - Performs one-factor-at-a-time analyses on Fisher-z:
     - Time-of-Day (EMRN, MORN, AFTN, EVEN).
     - Day-of-Week (WED, THUR, SAT, SUN).
     - Week-to-Week (WK01 vs WK02) via paired t-tests on matched conditions.

4. **Visualization**
   - Dot-and-whisker (mean ± 95% CI) plots for the factors above.
   - Optional histograms and z-distribution diagnostics.

---

### 3.3 Experiment 4 – ENF Server-Room Robustness

**Script:** `enf_analysis_top_reliability_stat_server.m`  
**Data root:** `exp_inputs/SRV_L/`

**Goal:** Evaluate how robust ENF extraction is inside a **high-density server-room environment**, focusing on:

- Sensor position on the server PSU (SBOX vs SPSU).
- Day-of-week effects (MON vs WED).

**Overview:**

1. **Inputs**
   - Traverses `SRV_L/SBOX/` and `SRV_L/SPSU/`, each with `MON/` and `WED/` subfolders.
   - Within each (Site, Day) folder, it loads:
     - `fpga_em_trace_dc.wav` – Ambient EM at the server PSU.
     - `mains_pow_trace_ac.wav` – Mains reference.
   - Five trials (`T01`–`T05`) per condition.

2. **ENF Extraction & Correlation**
   - STFT around 60 Hz harmonics.
   - Weighted-harmonic ENF estimator.
   - Pearson correlation **r** per trial, then Fisher-z for inference.

3. **Statistics**
   - Descriptive statistics on **r** and **z**.
   - One-way repeated-measures ANOVA on **z** for:
     - **Day**: MON vs WED (subjects = repeated trials, averaged over Site).
     - **Site**: SBOX vs SPSU (subjects = repeated trials, averaged over Day).
   - Reports effect sizes and 95% CIs, then back-transforms model estimates to **r**.

4. **Visualization**
   - Dot-and-whisker plots for Site and Day.
   - Boxplots and histograms for **r** and **z**.
   - Diagnostic z-distribution plots.

---

### 3.4 Experiment 5 – ENF Multi-Location Validation

**Script:** `enf_analysis_top_multi_loc_stat.m`  
**Data root:** `exp_inputs/MULTI/`

**Goal:** Quantify **cross-location separability** of ENF signatures by correlating a single local ambient-EM ENF against mains references drawn from multiple grid regions (e.g., US East / Texas / West 60 Hz, German 50 Hz).

**Overview:**

1. **Inputs**
   - Root directory contains trial folders (`T01/`, `T02/`, `T03/`, …) for each date.
   - Each `Txx` folder holds:
     - One local EM trace: `fpga_em_trace_dc_*lab.wav`.
     - Multiple mains references whose filenames encode **grid/region**:
       - `*_egrid_citya_lab.wav`, `*_egrid_citya_home.wav`, `*_egrid_worcester.wav`, `*_tgrid_richardson.wav`, `*_wgrid_tucson.wav`, `*_dresden.wav`, etc.
   - Subtrees `US_60/` and `DE_50/` correspond to 60 Hz and 50 Hz grids.

2. **ENF Extraction & Correlation**
   - STFT centered on the appropriate nominal frequency (50 or 60 Hz) and harmonics.
   - ENF estimation via weighted-harmonic approach for each trace.
   - For every EM trace, computes correlation **r** with every available mains reference in that trial folder.
   - Aggregates pairwise correlations into:
     - Per-folder lower-triangle correlation matrices.
     - Per-grid average correlation matrices and summaries.

3. **Statistics**
   - Descriptive statistics on **r** and **z = atanh(r)** for:
     - Intra-grid (same grid) vs inter-grid (cross-grid) comparisons.
   - Two-sample t-tests on **z** (intra-grid > inter-grid, Welch unequal-variance).
   - Optional parametric **Equal-Error Rate (EER)** estimates under a Gaussian model for the intra-grid vs inter-grid distributions, including “conservative” variants using reliability-derived variance.

4. **Visualization**
   - Lower-triangle heatmaps with annotated correlation coefficients.
   - Per-grid average heatmaps (e.g., `US_60` and `DE_50` summaries).
   - Histograms and boxplots for **r** and **z**.
   - Optional plots of z-distributions and ROC/EER operating points.

---

### 3.5 Core ENF Extraction Function

**File:** `proc_enf_analysis.p`

This compiled MATLAB function encapsulates the shared pipeline used by all top-level scripts:

1. Pre-processing & optional resampling.
2. STFT/spectrogram computation around the nominal ENF harmonics.
3. Weighted-harmonic ENF estimation per trace.
4. Temporal alignment and correlation computation.
5. Optional plot generation.

The top-level scripts configure and call this function; you do not need to invoke it directly.


## 4. Reproducing the Analysis

### 4.1 Requirements

To run the analysis end-to-end, you will need:

- **MATLAB R2018b or newer**.
- The **Signal Processing Toolbox**.
- Sufficient disk space to load the WAV files and write figures / logs.

No additional toolboxes are required beyond what the top-level scripts are already using.

### 4.2 Quick Start

1. **Clone or download this repository** to your local machine.

2. **Open MATLAB** and set the current folder to the repository root:
   ```matlab
   cd('PATH/TO/SP2026_OPENSCI_GRIDTRUST');
