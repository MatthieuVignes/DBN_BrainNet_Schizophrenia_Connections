# Smoking Sensitivity Analysis (R2 Major Comment 4) — Run Instructions

Addresses the smoking half of Reviewer 2's Major Comment 4 (unmodelled
demographic/clinical confounds). Restricts the full-resolution sample to
non-smokers only in both groups (34 HC, 17 PDS — vs. 45/45 in the full
sample), and re-runs the same Bayesian density comparison, to check whether
the group difference found in the full sample survives once the
smoking/diagnosis confound is removed.

**The medication half of this comment is NOT addressed here** — no
per-participant antipsychotic dose/medication data has been provided. If
such a file becomes available, an equivalent script can be written
following the same pattern. Until then, this remains a stated limitation
in the response letter, not an analysis.

## Prerequisites

- `00_prep.R` already run (produces `prepped_data.rds`)
- `Demographics.csv` available — edit the path near the top of
  `15_prep_smoking_subsample.R` if it isn't at
  `/mnt/user-data/uploads/Demographics.csv` on your machine
- Same R packages as the main pipeline (`bnlearn`, `dbnR`, `data.table`,
  `parallel`) — already installed if you've run the full-resolution pipeline

## Steps

```bash
Rscript 15_prep_smoking_subsample.R 2>&1 | tee smoking_prep.log
```
Prints smoking status counts by group and saves `prepped_data_nonsmoker.rds`.
Confirm the group sizes look right (34 HC / 17 PDS non-smokers) before
continuing.

```bash
Rscript 16_bootstrap_parallel_nonsmoker.R 2>&1 | tee smoking_bootstrap.log
```
**Update**: this script now sets single-threaded BLAS/data.table internally
(via `RhpcBLASctl::blas_set_num_threads(1)` and `setDTthreads(1)`), so the
manual shell `export OPENBLAS_NUM_THREADS=1` / `export OMP_NUM_THREADS=1`
steps that were previously required before this command are no longer
necessary — the script handles it itself now, more reliably than the shell
env vars did. Harmless to still set them if you're in the habit, just
redundant. This script runs 50 bootstraps ×
4 groups = 200 fits (half the primary analysis's bootstrap count, since this
is a secondary sensitivity check — documented as a deliberate choice, not
a shortcut). Safe to interrupt and re-run; completed fits are skipped.

```bash
Rscript 17_aggregate_bayesian_nonsmoker.R 2>&1 | tee smoking_bayesian.log
```
Prints the non-smoker-only posterior density estimates and the two
group-difference contrasts (encoding, retrieval). Compare these directly
against the full-sample numbers already in the response letter: if the
direction and credibility of the effect hold up in the non-smoker
subsample, that's good evidence the group difference isn't just a proxy
for the smoking confound.

## What to do with the result

Send me the three log files (or just the printed summary lines) once this
finishes, and I'll write up the comparison against the full-sample result
and drop it into the response letter under R2 Major Comment 4.
