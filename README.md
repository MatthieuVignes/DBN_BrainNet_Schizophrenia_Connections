# Dynamic Bayesian network analysis of PDS vs HC: supplementary code

Code behind the revision of "Brain networks for schizophrenia show differences in 
functional connections compared to healthy case controls for memory encoding and \
retrieval tasks" (*Journal of the Royal Society of New Zealand*). 
It rebuilds the two-slice DBN analysis from the participant-level data and adds the 
relevant analyses the reviewers requested.

## Headline result

Under formal inference there is **no statistically significant group difference in network density**, 
in either encoding or retrieval (bootstrap-difference intervals and permutation tests). 
Retrieval PDS > HC (+281 edges, 95% interval [-153, 659], p = 0.19) matches the direction 
of the original claim but is not distinguishable from resampling variability. 
Individual edges are unstable (split-half Jaccard ≈ 0.11) and the pipeline over-calls 
edges in simulation (precision ≈ 0.38), so results should be read at network level and as 
exploratory.

An earlier Beta-Binomial density test reported "credible" differences. A simulation 
showed it is miscalibrated (30% false positives under the null), so it has been replaced. 
See `README_FULL_RESOLUTION.md`.

## Where to look

| Question | Read |
|---|---|
| How to run the main 116-region pipeline, run order, timings, results | `README_FULL_RESOLUTION.md` (scripts `00`-`14`, `24`, `25`, `27`) |
| Smoking (non-smoker) sensitivity | `README_SMOKING_SENSITIVITY.md` (scripts `15`-`17`, then `24`, `25`) |
| BPRS and antipsychotic dose correlations | `README_CLINICAL_CORRELATION.md` (scripts `18`, `19`) |
| Simulation with known ground truth | `README_SIMULATION_VALIDATION.md` (scripts `20`-`23`, `26`) |

## Mapping to reviewer comments

| Analysis | Scripts | Comment |
|---|---|---|
| Data prep, two-slice DBN, `maxp` cap | `00`, `01`, `13`, `14` | R1-1 |
| Bootstrap and L1-optimal threshold (0.50) | `14`, `06` | R1-2 |
| Group inference (bootstrap-difference interval; permutation) | `24`, `25`; `11`, `12` | R2-1 |
| Directed vs intra-slice edges | `08` | R2-2 |
| Split-half stability; simulation | `10`; `20`–`23` | R2-3 |
| Smoking sensitivity; CPZ dose | `15`–`17`; `18`, `19` | R2-4 |
| Pearson-network comparison | `09` | R2-6 |
| BPRS correlations | `18`, `19` | R1-4 |

## Folders

- `scripts_full_resolution/`: primary pipeline (use this).
- `clinical_data/`: derived BPRS and CPZ-equivalent files.
- `results/`: main outputs.

## Assumptions and deviations

- Signed square-root transform `sign(x)·√|x|` of BOLD values.
- Parents-per-node cap (`maxp = 8`) in hill-climbing for tractability at 118 nodes.
- `Scan_no` is included as a model node.
- Single-threaded BLAS and data.table inside the parallel runners (required to avoid fork deadlocks).

## Limitations not yet addressed

Medicated vs unmedicated group-level comparison; education as a covariate; second-order DBN; re-scoring the simulation fits with the replacement test.
