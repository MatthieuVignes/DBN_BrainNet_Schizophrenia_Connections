# Simulation-Based Method Validation (Tier 1 + Tier 2)

Answers a question no amount of analysis on the real data can answer on its
own: **Does this method actually recover true network structure, and are its
formal statistical tests well-calibrated?** On real data you never know the
true network, so this can only be checked with simulated data, where the
ground truth is known by construction.

## What is simulated, and why it is meant to be "realistic"

- **Modular, sparse topology** (a stochastic block model), as real brain networks 
  are modular.
- **AR(1)-autocorrelated innovations**, mimicking hemodynamic-response
  smoothing in real BOLD signal. If `analysis/prepped_data.rds` (from the real
  pipeline's `00_prep.R`) is present, the AR(1) coefficient is **estimated
  from the real data's autocorrelation** rather than picked arbitrarily: check the 
  printed "AR(1) params" line to see whether this happened, or whether it fell back 
  to a literature-typical default.
- **Between-participant heterogeneity**: every true coefficient gets
  participant-specific random-effect jitter, so no two simulated
  participants share identical dynamics.
- **Both edge types deliberately embedded**: contemporaneous (intra-slice)
  and lagged (transition, past→present) edges with known, non-zero
  coefficients, so recovery can be evaluated separately for each.
- **A "known-null" node**: `Scan_no` is included exactly as it is in the
  real pipeline (a genuine model node, per the manuscript's stated use of it
  as a control variable), but with **zero true edges**. Any edge the
  recovery pipeline finds involving `Scan_no` is, by construction, a false
  positive. Built-in specificity check.

## Tier 1 - does the method recover a known network at all?

`21_tier1_recovery.R`: simulates one group (default: 45 participants × 171
scans, matching H.E's dimensions) from a single known ground truth, runs the
exact same bootstrap + optimal-threshold pipeline as the real analysis
(reusing `01_helpers.R` unmodified), and reports:
- Precision/recall/F1 for intra-slice edges, on the **undirected skeleton**
- Precision/recall/F1 for transition edges, on the **directed (ordered)
  pair**. Notice that direction accuracy is folded directly into this, since a true
  A→B edge recovered as B→A counts as both a false negative for A→B and a
  false positive for B→A
- A supplementary "direction accuracy given detection" number
- A full precision-recall curve across every candidate threshold
  (`results/tier1_pr_curve.csv`), so the optimal-threshold procedure's chosen point
  can be shown against the best achievable trade-off
- The Scan_no specificity check described above

**Why intra-slice edges are scored on the undirected skeleton, not
direction**: purely contemporaneous linear-Gaussian relationships are not
always identifiable from observational data alone (Markov equivalence; see
Spirtes, Glymour & Scheines, 2000). A correctly-working method can recover
the right undirected relationship while "getting the direction wrong" simply
because both directions are statistically indistinguishable given the data.
Scoring intra-slice edges on direction would penalise the method for a
limitation of the data-generating process, not a flaw in the method. Hence
direction is only scored for transition edges, where time itself breaks the
symmetry.

```bash
Rscript 21_tier1_recovery.R 2>&1 | tee tier1.log
```
Expected runtime: minutes to 1h (50 bootstrap fits).

## Tier 2 - are the formal statistical tests well-calibrated, and do they have power?

Two scenarios, each repeated `N_REPLICATES` times (default 30) with fresh
simulated data:
- **null**: group A and group B are simulated from the *identical* ground
  truth. Any replicate where the Bayesian density-comparison test (same
  method as `07_bayesian_density_full.R`) finds a credible difference is a
  **false positive**. Across replicates, this rate should sit near the
  nominal ~5% for a well-calibrated test.
- **effect**: group B's ground truth has `N_EXTRA_EFFECT_EDGES` (default 15)
  additional true edges layered on top of group A's. This reflects a known, designed
  density difference. The fraction of replicates that correctly detect a
  credible difference **in the correct direction** is the method's empirical
  **power** at this effect size and this sample size (45/group, matching the
  real data).

`B_BOOT` (bootstraps per group per replicate) and `N_REPLICATES` are both
reduced from the primary analyses, as total cost multiplies across
replicates × scenarios × groups × bootstraps (documented trade-off). 
Increase `N_REPLICATES` for tighter calibration/power estimates if you have time to 
let this run longer; each replicate is independent, so cost scales linearly.

```bash
Rscript 22_tier2_runner.R 2>&1 | tee tier2_run.log
```
This prints an estimated wall-time before starting, based on your earlier
pilot's per-fit timing at this node count. Safe to interrupt and
resume. Completed fits are skipped on rerun, same pattern as the main
pipeline's bootstrap runners. Sets BLAS/data.table to single-threaded
internally (no shell `export` commands needed).

```bash
Rscript 23_tier2_aggregate.R 2>&1 | tee tier2_aggregate.log
```
Prints the false-positive rate (null scenario) and power (effect scenario),
each with a binomial confidence interval reflecting `N_REPLICATES`'
precision, plus how often a credible difference was found in the *wrong*
direction under the effect scenario (should be rare/zero for a trustworthy
test).

## Results

**Tier 1** (one simulated group, 45 x 171, AR(1) phi = 0.518 estimated from the real data, optimal threshold 0.50): intra-slice skeleton TP 28, FP 46, FN 0 (precision 0.38, recall 1.00, F1 0.55); transition edges TP 47, FP 78, FN 0 (precision 0.38, recall 1.00, F1 0.55); direction accuracy given detection 47/47; 3 spurious edges involving the null node `Scan_no`. 
The method finds every true edge but over-predicts: roughly 60% of retained edges are false, so individual edges should not be interpreted or very cautiously/validated.

**Tier 2** (30 replicates per scenario):
- Null: 9/30 replicates flagged a "credible" difference, i.e. a **30% false-positive rate** for the Beta-Binomial density test (nominal 5%); 7 of those were in the "correct" direction by chance.
- Effect (+15 true edges): 30/30 detected in the right direction (power 100%), which says little because an anti-conservative test detects everything.

## Consequence for the main analysis

The Beta-Binomial test is not valid and was replaced by the pairwise bootstrap-difference interval (`25_bootstrap_difference_test.R`). 
Under it no group contrast is significant in the real data. The Tier 2 null/effect fits (`tier2_boot/`) can be re-scored with the new test (no new fitting) to check its calibration and power at this effect size.

## To be done

- Re-scoring Tier 2 with the replacement test.
- Simulations at higher `maxp` or with larger effects; power curves.
