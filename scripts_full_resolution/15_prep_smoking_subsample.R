suppressMessages({library(data.table)})

# Sensitivity analysis for R2 Major Comment 4 (smoking confound). Restricts
# the full-resolution sample to NON-SMOKERS ONLY in both groups, so that any
# group density difference found here can no longer be attributed to the
# smoking-status imbalance in the full sample (PDS were predominantly
# smokers, HC predominantly non-smokers in the full sample; see Table 1).
#
# NOTE: this addresses the SMOKING half of Reviewer 2's Major Comment 4 only.
# The manuscript also reports antipsychotic medication exposure as a
# potential confound; no per-participant medication/dose data has been
# provided, so that half of the comment remains a stated limitation, not an
# analysis, until such a file is available.

raw <- readRDS("analysis/prepped_data.rds")  # output of 00_prep.R -- must be run first
demo <- fread("data/Demographics.csv")  # <-- edit path if needed

# One row per participant is enough; smoking status doesn't vary by task/scan.
demo_small <- unique(demo[, .(Participant_id = participant_id, cigs)])
demo_small[, smoker := as.integer(cigs > 0)]

cat("Smoking status by group (unique participants):\n")
ids <- unique(raw[, .(Participant_id, conditiongroup)])
ids[, group_letter := substr(conditiongroup, 1, 1)]
ids <- unique(ids[, .(Participant_id, group_letter)])
merged_check <- merge(ids, demo_small, by = "Participant_id")
print(table(merged_check$group_letter, merged_check$smoker))

# Filter the full-resolution prepped data down to non-smokers only.
nonsmoker_ids <- demo_small[smoker == 0, Participant_id]
raw_ns <- raw[Participant_id %in% nonsmoker_ids]

cat("\nParticipants retained per group (non-smokers only):\n")
print(raw_ns[, .(n_participants = uniqueN(Participant_id)), by = conditiongroup])

saveRDS(raw_ns, "analysis/prepped_data_nonsmoker.rds")
cat("\nSaved prepped_data_nonsmoker.rds:", nrow(raw_ns), "rows\n")
