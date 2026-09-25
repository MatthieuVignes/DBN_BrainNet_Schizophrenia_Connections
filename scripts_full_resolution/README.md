# Corrected DBN Analysis — Results Memo

This package contains a from-scratch R re-implementation of the manuscript's dynamic
Bayesian network (DBN) pipeline, run on the uploaded participant-level data, extended
with the analyses requested by Reviewers 1 and 2. All code is in `scripts/`, all
numeric outputs in `results/`.

## 1. What was replicated, and how

- **Trial filtering** (`00_prep.R`): encoding groups (H.E, S.E) keep `TASK` trials
  only; retrieval groups (H.R, S.R) keep `CORRECTLY` + `INCORRECTLY` trials;
  `CONTROL`/`Firstob` dropped. This exactly reproduces the manuscript's reported
  171 encoding / 214 retrieval scans per participant.
- **Transform**: BOLD region values are signed (min ≈ −717), so a literal `sqrt()`
  is undefined for much of the data. We used the signed square-root transform,
  `sign(x)·√|x|` — variance-stabilising in the same spirit as the manuscript's
  stated `sqrt()` step, but the manuscript text doesn't specify how negative
  values were handled, so **this is a documented assumption**, not a confirmed match.
- **DBN structure learning** (`01_helpers.R`): two-slice network — a static
  (intra-slice) network learned once via hill-climbing, applied to both time
  slices, joined with a transition network restricted (via blacklist) to
  past→present arcs only — reusing `dbnR`'s own blacklist/merge internals, with
  `bnlearn::hc()` doing the greedy score-based search. This mirrors the
  manuscript's stated "hill-climbing algorithm" directly. Score used: `bic-g`
  (matches the manuscript) for the main results; `bge` (Bayesian Gaussian
  equivalent score) is wired in as a drop-in option (`score="bge"` in
  `fit_group_dbn_hc()`) for a fully Bayesian structure-learning criterion.
- **Bootstrapping** (`05_bootstrap_runner.R`): 100 bootstrap resamples per group
  (participants resampled with replacement), refitting the full DBN each time —
  400 fits total, all completed.
- **Optimal threshold** (`06_aggregate_strength.R`): applied `bnlearn`'s own
  implementation of the Scutari & Nagarajan (2013) L1-optimal threshold — the
  exact method the manuscript cites — to the bootstrap edge-confidence scores.

### One necessary deviation: region resolution

The manuscript's 116-region network (236 nodes once unrolled into two time
slices) is **computationally intractable in an interactive session** — a single
static-network fit did not complete in over 4 minutes, and scaling tests
(`03_scaling_test.R`-style benchmarking, not included as it was exploratory)
showed fit time growing roughly with the 4th–5th power of node count. To make
bootstrapped structure learning tractable, the 116 AAL regions were aggregated
into **25 anatomically-grounded bilateral macro-ROIs** (`00b_macro_roi.R`) —
standard lobar/system groupings (Frontal, Frontal-Orbital, Motor, Parietal,
Occipital, Temporal, Cingulate, Limbic, Insula, Basal Ganglia, Thalamus,
Cerebellum × L/R, plus midline Vermis), by simple averaging. Every one of the
116 original regions is assigned to exactly one macro-ROI (verified
programmatically — nothing dropped or double-counted). `Scan_no` was **not**
included as a model node in this pilot (unlike the manuscript, which used it as
a control node) — a minor further simplification, noted for transparency.

**This is a scoped pilot, not a substitute for the full analysis.** All
functions in `01_helpers.R` accept an arbitrary region set, so the identical
pipeline can be rerun at full 116-region resolution given adequate (HPC/cluster)
compute. Absolute network density at 25-node resolution is almost certainly
inflated relative to the 116-region networks (averaging regions reduces noise
and inflates apparent connectivity — see density values below, which are much
higher than typically reported at full resolution), so **only the relative,
between-group comparisons here should be treated as informative; the absolute
density numbers should not be quoted as replicating the original 116-region
results.**

## 2. Results, mapped to specific reviewer comments

### R2 Major Comment #1 — no formal statistical inference

Two independent formal tests were added, deliberately using different
statistics so the results can be read against each other:

**(a) Bayesian Beta-Binomial network-density comparison**
(`07_bayesian_density.R`, `results/threshold_and_bayesian_density_summary.csv`) —
treats each of the 624 addressable arcs, across all 100 bootstrap fits per
group, as a Bernoulli(θ_group) trial; flat Beta(1,1) prior; posterior credible
intervals and pairwise contrasts computed by Monte Carlo.

| Contrast | Mean density difference | 95% CrI | P(difference > 0) |
|---|---|---|---|
| Encoding: PDS − HC | −0.0002 | [−0.0037, 0.0033] | 0.45 (no credible difference) |
| Retrieval: PDS − HC | +0.0126 | [0.0095, 0.0157] | 1.00 (credible PDS > HC) |

**(b) Frequentist permutation test** (`11_permutation_runner.R` /
`12_permutation_pvalues.R`, `results/permutation_test_results.csv`) — group
labels shuffled 50× per contrast (encoding, retrieval), a single (non-bootstrap)
DBN fit per pseudo-group per permutation, compared against the true-label
observed fit on edge-count difference and Jaccard overlap.

| Contrast | Observed \|edge diff\| | Permutation p | Observed Jaccard | Permutation p |
|---|---|---|---|---|
| Encoding (H.E vs S.E) | 22 | 0.46 | 0.334 | 0.12 |
| Retrieval (H.R vs S.R) | 9 | 0.78 | 0.394 | 0.38 |

**Important tension, worth reporting honestly in the response letter:** the
bootstrap-aggregated Bayesian test finds a credible retrieval-phase group
difference; the permutation test on single (non-aggregated) fits does not reach
significance in either phase. This is not a contradiction so much as a
reminder that single DBN fits are noisy (see stability results below) and that
the manuscript's own instinct to bootstrap and threshold, rather than trust one
fit, is the right one — but it also means claims of a group difference should
rest on the bootstrap-aggregated test, stated as such, not on a single edge
count comparison.

### R2 Major Comment #2 — DBN's directional/temporal structure unused

`08_directed_edges.R` / `results/directed_edge_summary.csv` — breaks the averaged
(thresholded) network into intra-slice vs. genuine transition (past→present)
edges:

| Group | Total edges | Transition edges | % transition |
|---|---|---|---|
| H.E | 506 | 22 | 4.3% |
| H.R | 473 | 25 | 5.3% |
| S.E | 477 | 23 | 4.8% |
| S.R | 480 | 26 | 5.4% |

Only a small, fairly consistent minority of edges are genuinely temporal in
this dataset/resolution — worth stating plainly rather than glossing over, since
it tempers how much the "DBN's directionality" analysis can carry the paper's
main story, but it directly answers the reviewer's request for this breakdown.

### R2 Major Comment #3 — overfitting / instability

`10_stability_check.R` / `results/split_half_stability.csv` — each group's
participants split into two independent halves (n≈22 each), DBN fit
independently on each half, Jaccard overlap computed:

| Group | Jaccard (half1 vs half2) |
|---|---|
| H.E | 0.241 |
| H.R | 0.326 |
| S.E | 0.408 |
| S.R | 0.347 |

Substantial instability at n≈22 — a legitimate finding that supports R2's
concern, and simultaneously the strongest argument *for* keeping the
manuscript's bootstrap-and-threshold approach rather than reporting any single
fit's edge list.

### R2 Major Comment #6 — comparison to a conventional method

`09_correlation_comparison.R` / `results/correlation_vs_dbn_comparison.csv` —
FDR-corrected (BH, q<.05) Pearson correlation network vs. the DBN's static
(intra-slice) undirected skeleton:

| Group | Correlation edges | DBN static edges | Overlap | Jaccard |
|---|---|---|---|---|
| H.E | 300 | 242 | 239 | 0.789 |
| H.R | 325 | 224 | 224 | 0.689 |
| S.E | 325 | 227 | 227 | 0.698 |
| S.R | 300 | 227 | 225 | 0.745 |

The DBN skeleton is almost entirely a subset of the FDR-significant correlation
edges (near-total overlap, DBN always sparser) — supports framing the DBN as a
more parsimonious representation of largely the same pairwise structure, which
is a defensible, non-overclaiming way to answer this comment.

### R1 — optimal threshold detail

`06_aggregate_strength.R` — using `bnlearn`'s own implementation of the cited
Scutari & Nagarajan (2013) method, the optimal threshold **independently
converged to 0.50 for all four groups**, matching the manuscript's reported
50/100 value. This is useful corroborating evidence for the response letter's
methodological-transparency section.

## 2b. Simulation-based method validation

Everything above evaluates the method on real data, where the true network
is unknown. `README_SIMULATION_VALIDATION.md` covers a complementary
validation using simulated data with a KNOWN ground truth — precision/
recall/F1 for network recovery (Tier 1), plus false-positive-rate and power
estimates for the formal statistical tests (Tier 2), under a realistic
(modular, autocorrelated, participant-heterogeneous) data-generating
process. This is the strongest available evidence that the method works in
principle, independent of anything about the real dataset's specific
findings.

## 3. Suggested next steps

1. Decide, with the reviewers'/editor's expectations in mind, whether the
   macro-ROI pilot results are sufficient to cite directly in the response
   letter (with the resolution caveat stated plainly), or whether the full
   116-region version should be run on proper compute before resubmission.
2. If proceeding to full resolution: the exact same `01_helpers.R` functions
   run unmodified on the 116-region data (just skip `00b_macro_roi.R`); budget
   for substantially longer per-fit runtime and consider a compute cluster or
   at minimum an overnight/unattended batch job rather than an interactive
   session.
3. Consider whether to also run the `bge` (Bayesian score) variant of the
   structure learning at full resolution for the main manuscript results,
   given reviewers asked for a Bayesian framework where possible — the code
   supports this as a one-argument change.
