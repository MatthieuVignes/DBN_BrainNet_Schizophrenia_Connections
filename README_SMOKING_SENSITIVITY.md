# Smoking Sensitivity Analysis (R2 Major Comment 4)

Addresses the smoking half of Reviewer 2's Major Comment 4 (unmodelled
demographic/clinical confounds). Restricts the full-resolution sample to
non-smokers only in both groups (34 HC, 17 PDS vs. 45/45 in the full
sample), and re-runs the density comparison to check whether the group difference 
seen in the full sample survives, once smoking is removed.

The medication half of this comment is NOT addressed here: no
per-participant antipsychotic dose/medication data has been provided. If
such a file becomes available, an equivalent script can be written
following the same pattern. Until then, this remains a stated limitation
in the response letter, not an analysis.

## Prerequisites

- `00_prep.R` already run (produces `analysis/prepped_data.rds`)
- `data/Demographics.csv` available (edit the path near the top of
  `15_prep_smoking_subsample.R`)
- Same R packages as the main pipeline (`bnlearn`, `dbnR`, `data.table`,
  `parallel`)
  
## Steps

```bash
Rscript 15_prep_smoking_subsample.R 2>&1 | tee smoking_prep.log
```
Prints smoking status counts by group and saves `analysis/prepped_data_nonsmoker.rds`.
Confirm the group sizes look right (34 HC / 17 PDS non-smokers) before
continuing.

```bash
Rscript 16_bootstrap_parallel_nonsmoker.R 2>&1 | tee smoking_bootstrap.log
```
**Update**: this script now sets single-threaded BLAS/data.table internally
(via `RhpcBLASctl::blas_set_num_threads(1)` and `setDTthreads(1)`) (the 
manual `export OPENBLAS_NUM_THREADS=1` / `export OMP_NUM_THREADS=1` shell
steps should not be necessary; the script handles it itself now). Harmless to 
still set them if you're in the habit, just redundant. This script runs 50 bootstraps ×
4 groups = 200 fits (half the primary analysis's bootstrap count, as a secondary 
sensitivity check). Safe to interrupt and re-run; completed fits are skipped.

```bash
Rscript 17_aggregate_bayesian_nonsmoker.R 2>&1 | tee smoking_bayesian.log
```
Prints the non-smoker-only posterior density estimates and the two
group-difference contrasts (encoding, retrieval). Compare these directly
against the full-sample numbers: if the direction and credibility of the effect 
hold up in the non-smoker subsample, that is good evidence the group difference 
is not just a proxy for the smoking confound.

## Valid inference step (use this, not script 17's intervals)

```bash
Rscript 24_extract_per_bootstrap_counts.R boot_results_nonsmoker 50 per_bootstrap_edge_counts_nonsmoker.csv
Rscript 25_bootstrap_difference_test.R per_bootstrap_edge_counts_nonsmoker.csv bootstrap_difference_nonsmoker.csv
```
(`ADDR_SPACE` for the non-smoker run can be passed as a third argument to 25)

## Results so far, and caveat

The Beta-Binomial contrasts from script 17 were: encoding PDS - HC +0.0013 density, 95% CrI [0.0000, 0.0025]; retrieval +0.0140 [0.0127, 0.0153]. **These intervals are not valid evidence**: the same test has a 30% false-positive rate under the null (Tier 2), and in the full sample the corrected test finds no significant contrast. Treat the numbers as descriptive only until the two commands above have been run on the non-smoker fits.

The medication half of R2 Major Comment 4 is handled in `README_CLINICAL_CORRELATION.md` (CPZ-equivalent dose as a correlate; no medicated/unmedicated group-level sensitivity analysis yet).
