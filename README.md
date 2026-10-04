# DBN_BrainNet_Schizophrenia_Connections

Dynamic Bayesian network (DBN) analysis of fMRI connectivity in schizophrenia,
comparing patients (S) and healthy controls (H) across encoding (E) and
retrieval (R) task phases.

## Repository structure

| Folder | Contents |
|---|---|
| `data/` | Raw/source data: `240621_Final_Data_for_Sparsity_Manuscript.csv` (main dataset), `Demographics.csv` |
| `clinical_data/` | Clinical variables, BPRS subscales, CPZ dose data |
| `analysis/` | Prepped `.rds` files ready for modeling (`prepped_data.rds`, `prepped_data_nonsmoker.rds`) |
| `scripts_full_resolution/` | Full analysis pipeline (R scripts, numbered in run order) — see its own README files for details |
| `boot_results_full/`, `boot_results_nonsmoker/`, `tier1_boot/` | Saved bootstrap fit outputs |
| `results/` | Final outputs: summary stats, posteriors, stability/permutation results |
| `run_all.R` | Master script that sources the full pipeline end to end |

## Getting started

1. Clone the repo:
   ```bash
   git clone https://github.com/MatthieuVignes/DBN_BrainNet_Schizophrenia_Connections.git
   ```
2. Open `run_all.R` in RStudio and run it — it sources each pipeline step in
   `scripts_full_resolution/` in order (data prep → DBN fitting → bootstrapping →
   statistical tests → clinical correlations → simulation validation).
3. See `scripts_full_resolution/README.md` for a full walkthrough of each
   analysis step and how it maps to results in `results/`. Additional
   topic-specific notes are in `README_CLINICAL_CORRELATION.md`,
   `README_FULL_RESOLUTION.md`, `README_SIMULATION_VALIDATION.md`, and
   `README_SMOKING_SENSITIVITY.md` in the same folder.

## Requirements

- R (with `bnlearn`, `dbnR`, `rstudioapi`, and standard tidyverse packages)
- Python (only needed to rerun `18_clinical_prep.py` from scratch — not
  required if using the provided `clinical_data/` files)

## icence

- MIT (see `LICENSE`).
