suppressMessages({
  library(data.table)
  library(bnlearn)
  library(dbnR)
})
set.seed(1)
raw <- fread("data/240621_Final_Data_for_Sparsity_Manuscript.csv")
raw[, V1 := NULL]

region_cols <- setdiff(names(raw), c("Participant_id","DVARS","Scan_no","trialtype2","conditiongroup"))

raw <- raw[!(trialtype2 %in% c("CONTROL","Firstob"))]
raw <- raw[!(conditiongroup %in% c("H.E","S.E") & trialtype2 != "TASK")]
raw <- raw[!(conditiongroup %in% c("H.R","S.R") & trialtype2 == "TASK")]

raw[, (region_cols) := lapply(.SD, function(x) sign(x) * sqrt(abs(x))), .SDcols = region_cols]
raw[, Scan_no := as.double(Scan_no)]
setorder(raw, conditiongroup, Participant_id, Scan_no)

saveRDS(raw, "analysis/prepped_data.rds")
cat("Saved:", nrow(raw), "x", ncol(raw), "\n")
