import pandas as pd
import numpy as np
import re

MED_PATH = "/mnt/user-data/uploads/medication.tsv"
BPRS_PATH = "/mnt/user-data/uploads/bprs.tsv"
FMRI_PATH = "/mnt/user-data/uploads/240621_Final_Data_for_Sparsity_Manuscript.csv"
OUT = "/home/claude/analysis"

# ---------------------------------------------------------------------------
# Antipsychotic classification + chlorpromazine-equivalent (CPZ-eq) conversion
#
# Factors below are approximate daily-dose equivalents to 100mg chlorpromazine,
# drawn from the widely-cited international consensus dosing study (Gardner
# et al., Am J Psychiatry 2010) and standard updated equivalence tables
# (Leucht et al., 2020, Schizophrenia Bulletin, for newer agents). These are
# necessarily approximations -- equivalence tables vary somewhat across
# sources -- and are stated explicitly here so the conversion is auditable
# and can be swapped for a different reference table if reviewers prefer one.
#
# Non-antipsychotic psychotropics (antidepressants, benzodiazepines, mood
# stabilizers/anticonvulsants, anticholinergics prescribed for extrapyramidal
# side effects, sleep aids, thyroid medication, ECT) are explicitly excluded
# from the antipsychotic dose sum -- they answer a different question
# (polypharmacy burden) than the "antipsychotic exposure" confound the
# reviewer raised, so folding them in would answer something other than what
# was asked.
# ---------------------------------------------------------------------------

# mg of drug equivalent to 100mg chlorpromazine (i.e. CPZ-eq = dose_mg * 100 / factor)
CPZ_FACTORS = {
    "chlorpromazine": 100, "haloperidol": 2, "haldol": 2,
    "fluphenazine": 2, "prolixin": 2,
    "risperidone": 2, "risperdal": 2,
    "olanzapine": 5, "zyprexa": 5,
    "quetiapine": 100, "seroquel": 100,
    "ziprasidone": 60, "geodon": 60,
    "aripiprazole": 7.5, "abilify": 7.5,
    "clozapine": 50, "clozaril": 50,
    "paliperidone": 1.5, "invega": 1.5,
    "iloperidone": 4, "fanapt": 4,
    "loxapine": 10, "loxitane": 10,
    "thiothixene": 5, "navane": 5,
    "asenapine": 5, "saphris": 5,
}

NON_ANTIPSYCHOTIC_KEYWORDS = [
    "zolpidem","ambien","trihexyphenidyl","artane","lorazepam","ativan",
    "diphenhydramine","benadryl","citalopram","celexa","clonazepam","klonopin",
    "benztropine","cogentin","divalproex","depakote","doxepin","sinequan",
    "adapin","ect","fluoxetine","prozac","hydroxyzine","atarax","vistaril",
    "lamotrigine","lamictal","levothyroxine","levotabs","levothroid",
    "escitalopram","lexapro","lithium","eskalith","lithonate","loxapine",
    "paroxetine","paxil","phentermine","propranolol","inderal","remeron",
    "mirtazapine","temazepam","restoril","rozerem","topiramate","topamax",
    "trazodone","desyrel","bupropion","wellbutrin","synthroid",
]

# Maximum clinically-recognised oral daily doses (mg/day), used only to flag
# values that look like data-entry errors or a non-oral dosing interval
# recorded as if it were a daily oral dose -- NOT to silently "correct" the
# data. Values exceeding ~1.5x this ceiling are flagged and winsorised to the
# ceiling, with the flag reported so reviewers can see exactly what was done.
MAX_ORAL_DAILY_MG = {
    "haloperidol": 40, "haldol": 40, "fluphenazine": 40, "prolixin": 40,
    "risperidone": 16, "risperdal": 16, "olanzapine": 30, "zyprexa": 30,
    "quetiapine": 1600, "seroquel": 1600, "ziprasidone": 200, "geodon": 200,
    "aripiprazole": 30, "abilify": 30, "clozapine": 900, "clozaril": 900,
    "paliperidone": 12, "invega": 12, "iloperidone": 24, "fanapt": 24,
    "loxapine": 250, "loxitane": 250, "thiothixene": 60, "navane": 60,
    "asenapine": 20, "saphris": 20, "chlorpromazine": 2000,
}

# Long-acting injectable (LAI) formulations are dosed per 2-4 weeks, not
# daily -- applying an oral CPZ factor directly to the injectable dose
# massively overstates exposure. Risperdal Consta is the only LAI appearing
# in this dataset; converted here via the standard biweekly-dose -> oral
# daily-equivalent approximation (dose_mg / 12.5 -> oral mg/day), then the
# ordinary oral risperidone CPZ factor is applied to that daily equivalent.
LAI_BIWEEKLY_DIVISOR = {"consta": 12.5}

flags = []

def classify_dose(name, dose, participant_id, slot):
    if pd.isna(name) or pd.isna(dose) or dose < 0:
        return None
    n = name.lower()

    if "consta" in n:
        oral_equiv = dose / LAI_BIWEEKLY_DIVISOR["consta"]
        flags.append({"participant_id": participant_id, "slot": slot, "name": name,
                       "issue": "LAI dose converted to oral-equivalent",
                       "raw_dose": dose, "adjusted_oral_equiv_mg": oral_equiv})
        dose = oral_equiv
        n = n.replace("consta", "")

    for key, factor in CPZ_FACTORS.items():
        if key in n:
            max_dose = MAX_ORAL_DAILY_MG.get(key)
            if max_dose is not None and dose > 1.5 * max_dose:
                flags.append({"participant_id": participant_id, "slot": slot, "name": name,
                               "issue": f"exceeds 1.5x recognised max oral daily dose ({max_dose}mg) -- winsorised",
                               "raw_dose": dose, "adjusted_oral_equiv_mg": max_dose})
                dose = max_dose
            return dose * 100.0 / factor
    return 0.0  # matched a known non-antipsychotic, or unrecognised -> contributes 0

def is_recognised(name):
    if pd.isna(name):
        return True
    n = name.lower()
    if any(key in n for key in CPZ_FACTORS):
        return True
    if any(key in n for key in NON_ANTIPSYCHOTIC_KEYWORDS):
        return True
    return False

# ---------------------------------------------------------------------------
fmri = pd.read_csv(FMRI_PATH, usecols=["Participant_id","conditiongroup"]).drop_duplicates("Participant_id")
fmri["dx"] = fmri["conditiongroup"].str[0]
pds_ids = fmri.loc[fmri["dx"] == "S", "Participant_id"].unique()

med = pd.read_csv(MED_PATH, sep="\t")
med = med[med["participant_id"].isin(pds_ids)].copy()

unrecognised = set()
rows = []
for _, r in med.iterrows():
    total_cpz = 0.0
    any_med = False
    for i in range(1, 21):
        name = r.get(f"med_name{i}")
        use = r.get(f"med_use{i}")
        dos = r.get(f"med_dos{i}")
        if pd.isna(name):
            continue
        any_med = True
        if not is_recognised(name):
            unrecognised.add(name)
        if use == 1:
            cpz = classify_dose(name, dos, r["participant_id"], i)
            if cpz is not None:
                total_cpz += cpz
    rows.append({"participant_id": r["participant_id"], "total_cpz_eq": total_cpz, "any_medication_listed": any_med})

cpz_df = pd.DataFrame(rows)
print("Unrecognised medication names (not classified either way):", unrecognised)
print()
if flags:
    print(f"Flagged dose adjustments ({len(flags)} total) -- reviewable, not silent:")
    for f in flags:
        print(f"  {f['participant_id']} slot{f['slot']} {f['name']}: {f['issue']} "
              f"(raw={f['raw_dose']} -> adjusted={f['adjusted_oral_equiv_mg']:.1f})")
    pd.DataFrame(flags).to_csv(f"{OUT}/cpz_dose_flags.csv", index=False)
print()
print("N PDS with medication.tsv row:", len(cpz_df))
print("N with zero medications listed at all:", (~cpz_df["any_medication_listed"]).sum())
print()
print("CPZ-eq summary (all 45 PDS, zero-filled for those with no antipsychotic):")
print(cpz_df["total_cpz_eq"].describe())
print(f"\nMean={cpz_df['total_cpz_eq'].mean():.1f}, SD={cpz_df['total_cpz_eq'].std():.1f}")
print("(manuscript Table 2 reports 370.4mg, SD 299.3mg CPZ-equivalent for PDS)")

cpz_df.to_csv(f"{OUT}/cpz_equivalent_dose.csv", index=False)

# ---------------------------------------------------------------------------
bprs = pd.read_csv(BPRS_PATH, sep="\t")
bprs_pds = bprs[bprs["participant_id"].isin(pds_ids)][
    ["participant_id","bprs_positive","bprs_negative","bprs_mania","bprs_depanx"]
]
print("\nBPRS coverage:", bprs_pds.notna().all(axis=1).sum(), "of", len(pds_ids), "PDS participants complete")
bprs_pds.to_csv(f"{OUT}/bprs_subscales_pds.csv", index=False)

clinical = bprs_pds.merge(cpz_df, on="participant_id", how="outer")
clinical.to_csv(f"{OUT}/clinical_variables_pds.csv", index=False)
print("\nSaved clinical_variables_pds.csv:", clinical.shape)
