# Clinical Correlation Analysis (R1 Comment 4 + medication confound for R2 Major Comment 4)

Two things happen here:
1. **CPZ-equivalent antipsychotic dose** is computed from `medication.tsv`
   (completing the medication half of R2 Major Comment 4, alongside the
   smoking sensitivity analysis already built).
2. **BPRS subscales are correlated against a per-participant network
   score**, addressing R1 Comment 4.

## Step 1 — build the clinical variables file (run once, on any machine with Python)

```bash
python3 18_clinical_prep.py
```

This reads `medication.tsv` and `bprs.tsv`, and produces
`clinical_variables_pds.csv` (45 PDS participants × BPRS subscales + total
CPZ-equivalent dose). Read the printed output carefully:

- **Flagged dose adjustments** are printed explicitly, not applied silently.
  Two categories occur in this dataset:
  - A Risperdal **Consta** (long-acting injectable) dose was converted from
    its biweekly injection dose to an oral daily-equivalent before applying
    the standard oral CPZ factor — applying the oral factor directly to an
    injectable dose would have overstated exposure roughly 6-fold.
  - Two Olanzapine doses (60mg, 200mg/day) exceeded 1.5× the recognised
    maximum oral daily dose and were winsorised to that maximum (30mg/day) —
    almost certainly data-entry artifacts, but capped rather than silently
    dropped, and reported in `cpz_dose_flags.csv` for transparency.
- **Validation against the manuscript**: the resulting mean CPZ-equivalent
  dose (398.9mg, SD 445.3) is reasonably close to the manuscript's reported
  370.4mg (SD 299.3) — full agreement isn't expected, since published
  chlorpromazine-equivalence tables commonly differ by 10-30% depending on
  source, but this closeness is good corroborating evidence the
  classification and conversion are basically sound, not a methodology
  producing wildly different numbers from an independent construction.
- The conversion table itself (which drugs count as antipsychotics, and
  their CPZ-equivalence factors) is documented in comments at the top of
  `18_clinical_prep.py`, citing Gardner et al. (2010) and standard updated
  equivalence tables for newer agents — swap in a different reference table
  there if you or a coauthor prefers a different source.

## Step 2 — correlate against network structure (run on your machine, after the full-resolution pipeline)

Requires `strength_tabs_full.rds` (from `06_aggregate_strength_full.R`) and
`prepped_data.rds` (from `00_prep.R`) to already exist, plus
`clinical_variables_pds.csv` from Step 1 in the same directory.

```bash
Rscript 19_clinical_correlation.R 2>&1 | tee clinical_correlation.log
```

**What this computes, and why it isn't a new DBN fit per participant:** a
single participant's ~170-215 scans are nowhere near enough to fit a
reliable 118-node DBN on their own. Instead, for every edge retained in the
PDS group's bootstrap-thresholded network, each participant's own
correlation between the two involved regions (using only their own scans)
is computed, then averaged across all such edges — a "subject-level
expression of the group network template," computed separately for encoding
and retrieval. This is a standard individual-differences approach in
network neuroscience and needs no new model fitting.

Output: `network_scores_with_clinical.csv` (per-participant scores + all
clinical variables) and `clinical_correlation_results.csv` (Spearman
correlations between network score and each of the 4 BPRS subscales + CPZ
dose, separately for encoding/retrieval, BH-FDR corrected across all 10
tests).

## What to send back

The printed correlation table (or `clinical_correlation_results.csv`) —
I'll write up whatever survives FDR correction (or note honestly if nothing
does) for R1 Comment 4, and fold the CPZ-equivalent numbers into R2 Major
Comment 4 alongside the smoking sensitivity result.
