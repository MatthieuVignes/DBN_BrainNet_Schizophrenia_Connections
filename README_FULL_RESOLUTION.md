# Full-Resolution (116-Region) DBN Pipeline

This is the main methodology run at the manuscript's actual 116-region resolution. 
(run on an 8-core/32GB Ubuntu machine)

## Prerequisites

```r
install.packages(c("data.table","bnlearn"))
# dbnR isn't on CRAN -> install from GitHub:
install.packages("remotes")
remotes::install_github("dkesada/dbnR")
```

Confirm OpenBLAS is active (`sessionInfo()`, check the `BLAS:` line) before
running anything computational.

## What changed vs. the macro-ROI pilot

- **All 116 AAL regions used directly** without aggregation. `Scan_no` is 
  included as a model node; `DVARS` too.
- **`maxp` cap added** to `bnlearn::hc()`. It bounds parents-per-node, which is
  what makes fitting 118 nodes (236 once unrolled into two time slices)
  computationally feasible.
- **Parallelized bootstrapping** via `parallel::mclapply` (Linux/Mac only)

## Step-by-step

**1. Copy files onto your machine**, keeping them in one directory, and copy
your main CSV data to `data/240621_Final_Data_for_Sparsity_Manuscript.csv` 
(or edit the path at the top of `00_prep.R`).

**2. Prep the data (once, fast):**
```bash
Rscript 00_prep.R
```
Produces `analysis/prepped_data.rds`.

**3. Benchmark** (before committing to anything computationally expensive):
```bash
Rscript 13_benchmark_full_res.R
```
This fits ONE group ONE time at a few different `maxp` values and prints:
- actual fit time on your hardware
- an extrapolated estimate for the full 400-fit bootstrap run at that `maxp`

Use this to choose `maxp`. Lower `maxp` = faster but more constrained search;
if timings look reasonable, try a higher `maxp` (e.g. 12+) for a search closer
to the manuscript's unbounded hill-climbing. If a run takes many minutes for
one fit, that maxp is too high to bootstrap 400 times in reasonable time.

**4. Set `MAXP` in `14_bootstrap_parallel_runner.R`** (top of the file) to
whatever you decided in step 3, then run:
```bash
Rscript 14_bootstrap_parallel_runner.R
```
This fits 100 bootstraps × 4 groups = 400 total, distributed across
`N_CORES` (defaults to 7; adjust if you like).
**It's safe to interrupt and re-run**: completed fits are skipped
automatically, so you can stop and resume freely (e.g. overnight in
chunks, or if something crashes partway through).

**5. Aggregate results and get the optimal threshold:**
```bash
Rscript 06_aggregate_strength_full.R
```
Note the "Addressable arc space" number it prints at the end. **You need it
for step 6.**

**6. Group-difference inference:**
```bash
Rscript 24_extract_per_bootstrap_counts.R # edge count of each of the 400 stored bootstrap fits
Rscript 25_bootstrap_difference_test.R    # pairwise-difference interval, z and p for the 4 contrasts
```
`07_bayesian_density_full.R` is kept only as a descriptive posterior-density summary; 
do not use its credible intervals as a test (see below).

**7. The remaining analyses** (set `MAXP` at the top of each to match what
you used in step 4):
```bash
Rscript 08_directed_edges_full.R          # directed/temporal vs intra-slice edges
Rscript 09_correlation_comparison_full.R  # vs. conventional Pearson correlation network
Rscript 10_stability_check_full.R         # split-half overfitting/stability check (fits 8 new networks)
Rscript 11_permutation_runner_full.R      # formal frequentist test (parallelized, 104 new fits)
Rscript 12_permutation_pvalues_full.R     # p-values from the above
```

## Run order

| Order | Script | Purpose | Time (8 cores) |
|---|---|---|---|
| 1 | 00_prep.R, 01_helpers.R | data prep, helpers | seconds |
| 2 | 13_benchmark_full_res.R | pick `maxp` (4/6/8/12: 14/25/38/77 s per fit) | ~5 min |
| 3 | 14_bootstrap_parallel_runner.R | 400 bootstrap fits (B=100, maxp=8) | ~hours (theoretical); 1-2 days (reality) |
| 4 | 06 - 10 | threshold, density, directed edges, Pearson comparison, split-half | minutes |
| 5 | 11, 12 | permutation test | long (104 fits) / seconds |
| 6 | 24, 25 | bootstrap-based group inference | seconds |

## Results (maxp = 8, B = 100)

- Optimal threshold 0.50 in all groups; arcs retained: H.E 466, H.R 472, S.E 513, S.R 496 (addressable arc space 13,688).
- Directed (transition) share of retained edges: 21.0 / 22.9 / 21.2 / 21.0 %.
- Split-half Jaccard 0.110 - 0.116: individual edges are unstable; only network-level summaries should be interpreted.
- Pearson network (FDR) retains ~6,700-6,900 edges vs ~180-200 for the DBN static layer; overlap Jaccard 0.03.
- Permutation tests: encoding |dEdges| = 2, p = 1.00; retrieval |dEdges| = 154, p = 0.24.
- Per-bootstrap edge counts (mean, SD): H.E 2280 (181), H.R 2272 (178), S.E 2225 (189), S.R 2553 (117).

| Contrast | Mean diff (edges) | 95% interval | p |
|---|---|---|---|
| Encoding PDS - HC | -55 | [-527, 489] | 0.83 |
| Retrieval PDS - HC | +281 | [-153, 659] | 0.19 |
| HC retrieval - encoding | -7 | [-492, 483] | 0.98 |
| PDS retrieval - encoding | +328 | [-150, 681] | 0.14 |

**No group or phase contrast is statistically significant.** The retrieval PDS > HC difference is in the direction reported in the manuscript but is not distinguishable from bootstrap variability.

## Why two earlier tests were replaced

1. The Beta-Binomial density test pooled bootstrap*arc pairs as independent observations; in the Tier 2 null simulation it had a 30% false-positive rate (9/30).
2. A first "fix" (Wilcoxon / mean of the 100 replicates) is also anti-conservative: 84-86% false positives in `26_toy_bootstrap_calibration.py`. Bootstrap replicates are not independent samples, so their count is not a sample size.
3. The pairwise-difference interval (`25_`) uses the bootstrap distribution as the sampling distribution (~7% false positives in the same toy check).
