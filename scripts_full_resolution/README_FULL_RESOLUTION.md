# Full-Resolution (116-Region) DBN Pipeline — Run on Your Machine

This is the same methodology as the macro-ROI pilot, adapted to run at the
manuscript's actual 116-region resolution, on your 8-core/32GB Ubuntu machine.

## Prerequisites

```r
install.packages(c("data.table","bnlearn"))
# dbnR isn't on CRAN -- install from GitHub:
install.packages("remotes")
remotes::install_github("dkesada/dbnR")
```

Confirm OpenBLAS is active (`sessionInfo()`, check the `BLAS:` line) before
running anything computational — you already did this.

## What changed vs. the macro-ROI pilot

- **All 116 AAL regions used directly** — no aggregation. `Scan_no` is now
  included as a model node (it was dropped in the pilot); `DVARS` was already
  included.
- **`maxp` cap added** to `bnlearn::hc()` — bounds parents-per-node, which is
  what makes fitting 118 nodes (236 once unrolled into two time slices)
  computationally feasible at all. This is a documented methodological choice,
  not present in the pilot (where the small node count made it unnecessary).
- **Parallelized bootstrapping** via `parallel::mclapply` (fork-based —
  Linux/Mac only, which is fine for your Ubuntu box) instead of my sandbox's
  serial checkpoint-loop, which existed only to work around a 1-core, 250-second-per-call
  environment you don't have.

## Step-by-step

**1. Copy files onto your machine**, keeping them in one directory, and copy
your CSV to `/mnt/user-data/uploads/...` — or just edit the path at the top of
`00_prep.R` to point at wherever you put the CSV.

**2. Prep the data (once, fast):**
```bash
Rscript 00_prep.R
```
Produces `prepped_data.rds`.

**3. Benchmark before committing to anything expensive:**
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
`N_CORES` (defaults to 7, leaving one core free — adjust if you like).
**It's safe to interrupt and re-run** — completed fits are skipped
automatically, so you can stop and resume freely (e.g. overnight in
chunks, or if something crashes partway through).

**5. Aggregate results and get the optimal threshold:**
```bash
Rscript 06_aggregate_strength_full.R
```
Note the "Addressable arc space" number it prints at the end — you need it
for step 6.

**6. Formal Bayesian inference:** open `07_bayesian_density_full.R`, paste the
addressable-arc-space number from step 5 into the `ADDR_SPACE <- NULL` line,
then run it.

**7. The remaining analyses** (set `MAXP` at the top of each to match what
you used in step 4):
```bash
Rscript 08_directed_edges_full.R          # directed/temporal vs intra-slice edges
Rscript 09_correlation_comparison_full.R  # vs. conventional Pearson correlation network
Rscript 10_stability_check_full.R         # split-half overfitting/stability check (fits 8 new networks)
Rscript 11_permutation_runner_full.R      # formal frequentist test (parallelized, 104 new fits)
Rscript 12_permutation_pvalues_full.R     # p-values from the above
```

## What I'd genuinely expect

With OpenBLAS + 8 cores, I'd guess (not tested on your actual hardware —
that's exactly what the benchmark step is for) something in the range of
minutes-to tens-of-minutes for the full bootstrap step at a moderate `maxp`
(6-10), rather than the multi-hour-plus territory unbounded search hit in my
1-core sandbox. If the benchmark comes back much slower than that, the most
likely lever left is dropping `maxp` further, or reducing `B` from 100
bootstraps to something like 50 as a documented compromise.

## One thing worth deciding before you run the full 400-fit job

Whether to report the full-resolution results as your primary numbers for the
response letter, replacing the macro-ROI pilot, or to report both (pilot as a
robustness/sensitivity check, full-resolution as primary). I'd lean toward the
latter — showing the pattern holds (or doesn't) across resolutions is itself a
useful robustness argument against R2's overfitting/instability comment.
