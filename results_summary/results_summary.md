# Results summary: DBN analysis of PDS vs HC (116 regions)

Full-resolution two-slice dynamic Bayesian network, 100 participant-level bootstraps per group and phase (`maxp = 8`). Code and run instructions: `README_FULL_RESOLUTION.md` in the repository.

## Headline

**The data do not provide statistically significant evidence that network density differs between PDS and HC, in either task phase.** 
The retrieval difference (PDS > HC) has the direction reported in the manuscript, but its interval includes zero. 
Earlier drafts of this summary reported "credible" opposite-direction differences; 
that came from a test that was miscalibrated (30% false positives under the null) and is withdrawn.

![Fig: Forest plot of contrasts](./fig3_forest_plot.png)

| Contrast | Mean difference (edges) | 95% interval | p | Non-smokers only: difference [interval], p |
|---|---|---|---|---|
| Encoding, PDS − HC | −55 | [−527, 489] | 0.83 | +17 [−499, 501], 0.95 |
| Retrieval, PDS − HC | +281 | [−153, 659] | 0.19 | +192 [−220, 571], 0.36 |
| HC, retrieval − encoding | −7 | [−492, 483] | 0.98 | +8 [−441, 435], 0.97 |
| PDS, retrieval − encoding | +328 | [−150, 681] | 0.14 | +183 [−281, 649], 0.46 |

Interval = percentile interval of all pairwise bootstrap differences (`25_bootstrap_difference_test.R`). 
Divide edge differences by 13,688 (addressable arc space) for density: 
retrieval PDS - HC is about +0.020. 
The label-permutation test on single fits agrees (encoding p = 1.00, retrieval p = 0.24).

![Edge counts per bootstrap fit](./fig2_edge_count_distribution.png)

*The boxes show spread across bootstrap fits. Spread across refits is the sampling variability being tested, so separation of the boxes alone is not evidence of a group effect.*

## Threshold

The L1-optimal threshold (Scutari and Nagarajan 2013) is computed per group and equals 0.50 in all four. 
Arcs retained: H.E 466, H.R 472, S.E 513, S.R 496.

![Edge strength distribution](./fig1_edge_strength.png)

*(log y-axis; most observed arcs have confidence ≤ 0.05, with a small cluster near 1.0.)*

## Network structure and robustness

| Check | Result |
|---|---|
| Directed (transition) edges | 21.0-22.9% of retained edges |
| Split-half Jaccard | 0.110-0.116: individual edges are unstable |
| DBN vs FDR-corrected Pearson network | 184-202 vs 6,692-6,902 edges; DBN edges nested within the correlation set; overall Jaccard 0.03 |
| Simulation (known truth) | all true edges recovered (recall 1.00) but precision ≈ 0.38; direction of transition edges 47/47 correct |

Interpret network-level summaries only, not individual edges (about 60% are not correct).

## Clinical correlation (R1 Comment 4)

Each PDS participant's own data were projected on the PDS group's template edges (network-expression score) and correlated (Spearman, n = 45) with BPRS subscales and CPZ-equivalent dose, separately for encoding and retrieval. 
Ten tests; none passed BH-FDR (smallest FDR p = 0.56; raw |ρ| ≤ 0.26). With n = 45 only |ρ| ≳ 0.4 is reliably detectable, so this is a low-power null. 
CPZ-equivalent dose: mean 398.9 mg (SD 445.3) vs 370.4 (SD 299.3) in the manuscript.

## Smoking sensitivity (R2 Major Comment 4)

Non-smokers only (34 HC, 17 PDS; 50 bootstraps per group). 
No contrast is significant (table above). 
The retrieval PDS - HC difference keeps its sign (+192 edges, p = 0.36) with an interval similar in width to the full sample, so this analysis neither supports nor refutes a smoking explanation.

## Simulation validation

Tier 1 (known network, AR(1) φ = 0.518 from real data): recall 1.00, precision ≈ 0.38 for both intra-slice (skeleton) and transition edges, 3 spurious edges on the null node `Scan_no`. 
Tier 2: the original Beta-Binomial density test flagged a difference in 9/30 null replicates (30% false positives), which is why it was replaced. 
The replacement has not yet been re-scored on the Tier 2 fits.

## What this supports

A sparse, bootstrap-stabilised network pipeline that behaves sensibly in simulation and shows far greater sparsity than correlation networks, applied to data in which a group difference in density is suggested at retrieval but not demonstrated. 
Claims should be framed as exploratory.
