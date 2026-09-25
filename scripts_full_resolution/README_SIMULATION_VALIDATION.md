# Simulation-Based Method Validation (Tier 1 + Tier 2)

Answers the question no amount of analysis on the real data can answer on its
own: does this method actually recover true network structure, and are its
formal statistical tests well-calibrated? On real data you never know the
true network, so this can only be checked with simulated data where the
ground truth is known by construction.

## What's simulated, and why it's meant to be "realistic," not just convenient

- **Modular, sparse topology** (a stochastic block model), not a uniformly
  random graph — real brain networks are modular, and a random-graph target
  would be an unrealistically easy test that invites the objection that it
  doesn't resemble what the method is actually used on.
- **AR(1)-autocorrelated innovations**, mimicking hemodynamic-response
  smoothing in real BOLD signal. If `prepped_data.rds` (from the real
  pipeline's `00_prep.R`) is present, the AR(1) coefficient is **estimated
  from your own real data's autocorrelation** rather than picked arbitrarily
  — check the printed "AR(1) params" line to see whether this happened, or
  whether it fell back to a literature-typical default.
- **Between-participant heterogeneity**: every true coefficient gets
  participant-specific random-effect jitter, so no two simulated
  participants share identical dynamics.
- **Both edge types deliberately embedded**: contemporaneous (intra-slice)
  and lagged (transition, past→present) edges with known, non-zero
  coefficients, so recovery can be evaluated separately for each.
- **A "known-null" node**: `Scan_no` is included exactly as it is in the
  real pipeline (a genuine model node, per the manuscript's stated use of it
  as a control variable), but with **zero true edges** — any edge the
  recovery pipeline finds involving `Scan_no` is, by construction, a false
  positive. Built-in specificity check.

## Tier 1 — does the method recover a known network at all?

`21_tier1_recovery.R`: simulates one group (default: 45 participants × 171
scans, matching H.E's dimensions) from a single known ground truth, runs the
exact same bootstrap + optimal-threshold pipeline as the real analysis
(reusing `01_helpers.R` unmodified), and reports:
- Precision/recall/F1 for intra-slice edges, on the **undirected skeleton**
- Precision/recall/F1 for transition edges, on the **directed (ordered)
  pair** — direction accuracy is folded directly into this, since a true
  A→B edge recovered as B→A counts as both a false negative for A→B and a
  false positive for B→A
- A supplementary "direction accuracy given detection" number
- A full precision-recall curve across every candidate threshold
  (`tier1_pr_curve.csv`), so the optimal-threshold procedure's chosen point
  can be shown against the best achievable trade-off
- The Scan_no specificity check described above

**Why intra-slice edges are scored on the undirected skeleton, not
direction**: purely contemporaneous linear-Gaussian relationships are not
always identifiable from observational data alone (Markov equivalence — see
Spirtes, Glymour & Scheines, 2000). A correctly-working method can recover
the right undirected relationship while "getting the direction wrong" simply
because both directions are statistically indistinguishable given the data.
Scoring intra-slice edges on direction would penalize the method for a
limitation of the data-generating process, not a flaw in the method — hence
direction is only scored for transition edges, where time itself breaks the
symmetry.

```bash
Rscript 21_tier1_recovery.R 2>&1 | tee tier1.log
```
Expected runtime: a few minutes (50 bootstrap fits at macro-ROI-scale node
count, similar to your earlier pilot's ~6s/fit).

## Tier 2 — are the formal statistical tests well-calibrated, and do they have power?

Two scenarios, each repeated `N_REPLICATES` times (default 30) with fresh
simulated data:
- **null**: group A and group B are simulated from the *identical* ground
  truth. Any replicate where the Bayesian density-comparison test (same
  method as `07_bayesian_density_full.R`) finds a credible difference is a
  **false positive**. Across replicates, this rate should sit near the
  nominal ~5% for a well-calibrated test.
- **effect**: group B's ground truth has `N_EXTRA_EFFECT_EDGES` (default 15)
  additional true edges layered on top of group A's — a known, designed
  density difference. The fraction of replicates that correctly detect a
  credible difference **in the correct direction** is the method's empirical
  **power** at this effect size and this sample size (45/group, matching the
  real data).

`B_BOOT` (bootstraps per group per replicate) and `N_REPLICATES` are both
reduced from the primary analyses — total cost multiplies across
replicates × scenarios × groups × bootstraps, so this is a deliberate,
documented trade-off. Increase `N_REPLICATES` for tighter calibration/power
estimates if you have time to let this run longer; each replicate is
independent, so cost scales linearly, not combinatorially.

```bash
Rscript 22_tier2_runner.R 2>&1 | tee tier2_run.log
```
This prints an estimated wall-time before starting, based on your earlier
macro-ROI pilot's per-fit timing at this node count. Safe to interrupt and
resume — completed fits are skipped on rerun, same pattern as the main
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

## What to do with the results

If the false-positive rate sits near 5% and power is reasonably high at the
embedded effect size, that's strong, concrete evidence for a rebuttal
section arguing the method works and the real-data findings are worth
taking seriously. If either number looks off (inflated false-positive rate,
very low power), that's honestly still valuable to know before publication
rather than after — send me the results either way and I'll help interpret
them and, if needed, figure out what in the pipeline needs adjusting.
