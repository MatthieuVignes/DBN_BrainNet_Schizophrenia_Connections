# Clinical Correlation Analysis (R1 Comment 4 + medication confound for R2 Major Comment 4)

Two analyses here:
1. **CPZ-equivalent antipsychotic dose** is computed from `medication.tsv`
   (completing the medication half of R2 Major Comment 4, alongside the
   smoking sensitivity analysis already built).
2. **BPRS subscales are correlated against a per-participant network
   score**, addressing R1 Comment 4.

## Step 1 - build the clinical variables file

```bash
python3 18_clinical_prep.py
```

This reads `medication.tsv` and `bprs.tsv`, and produces
`clinical_data/clinical_variables_pds.csv` (45 PDS participants × BPRS subscales + total
CPZ-equivalent dose). Read the printed output carefully:

- **Flagged dose adjustments** are printed explicitly.
  Two categories occur in this dataset:
  - A Risperdal **Consta** (long-acting injectable) dose was converted from
    its biweekly injection dose to an oral daily-equivalent before applying
    the standard oral CPZ factor. Applying the oral factor directly to an
    injectable dose would have overstated exposure roughly 6-fold.
  - Two Olanzapine doses (60mg, 200mg/day) exceeded 1.5× the recognised
    maximum oral daily dose and were winsorised to that maximum (30mg/day). These 
    are almost certainly data-entry artifacts, but capped rather than 
    dropped, and reported in `cpz_dose_flags.csv` for transparency.
- **Validation against the manuscript**: the resulting mean CPZ-equivalent
  dose (398.9mg, SD 445.3) is reasonably close to the manuscript's reported
  370.4mg (SD 299.3). Notice that full agreement is not expected, since published
  chlorpromazine-equivalence tables commonly differ by 10-30% depending on
  source, but this closeness is good corroborating evidence the
  classification and conversion are sound.
- The conversion table itself (which drugs count as antipsychotics, and
  their CPZ-equivalence factors) is documented in comments at the top of
  `18_clinical_prep.py`, citing Gardner et al. (2010) and standard updated
  equivalence tables for newer agents. Users can swap in a different reference table
  if they prefer a different source.

## Step 2 - correlate against network structure (run after the full-resolution pipeline)

Requires `results/strength_tabs_full.rds` (from `06_aggregate_strength_full.R`) and
`analysis/prepped_data.rds` (from `00_prep.R`) to already exist, plus
`clinical_data/clinical_variables_pds.csv` from Step 1 in the same directory.

```bash
Rscript 19_clinical_correlation.R 2>&1 | tee clinical_correlation.log
```

**What this computes**, and why it isn't a new DBN fit per participant: a
single participant's ~170-215 scans are nowhere near enough to fit a
reliable 118-node DBN on their own. Instead, for every edge retained in the
PDS group's bootstrap-thresholded network, each participant's own
correlation between the two involved regions (using only their own scans)
is computed, then averaged across all such edges. A "subject-level
expression of the group network template," computed separately for encoding
and retrieval. This is a standard individual-differences approach in
network neuroscience.

Output: `results/network_scores_with_clinical.csv` (per-participant scores + all
clinical variables) and `results/clinical_correlation_results.csv` (Spearman
correlations between network score and each of the 4 BPRS subscales + CPZ
dose, separately for encoding/retrieval, BH-FDR corrected across all 10
tests).

## Results (n = 45 PDS participants, Spearman, BH-FDR over 10 tests)

| Phase | Variable | rho | p raw | p FDR |
|---|---|---|---|---|
| Encoding | BPRS negative | -0.19 | 0.22 | 0.56 |
| Encoding | BPRS depression/anxiety | -0.26 | 0.088 | 0.56 |
| Encoding | BPRS mania | 0.12 | 0.43 | 0.66 |
| Encoding | BPRS positive | 0.01 | 0.94 | 0.94 |
| Encoding | CPZ-equivalent dose | -0.10 | 0.51 | 0.66 |
| Retrieval | BPRS negative | 0.19 | 0.20 | 0.56 |
| Retrieval | BPRS depression/anxiety | 0.11 | 0.49 | 0.66 |
| Retrieval | BPRS mania | 0.21 | 0.17 | 0.56 |
| Retrieval | BPRS positive | -0.10 | 0.53 | 0.66 |
| Retrieval | CPZ-equivalent dose | -0.02 | 0.88 | 0.94 |

Nothing passes FDR correction (smallest FDR p = 0.56).

**Power caveat:** with n = 45, only |rho| of about 0.4 or larger is reliably detectable, 
so these are null results of limited sensitivity, not evidence of absence. 
The correlation with CPZ dose is small, which is mildly reassuring regarding medication, 
but it is not a medicated vs unmedicated group comparison of network structure (analysis not done). 
The CPZ conversion table should be checked by the user(s).
