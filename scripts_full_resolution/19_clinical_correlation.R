# ---------------------------------------------------------------------------
# R1 Comment 4: link network metrics to symptom severity (BPRS) and, as a
# natural extension, to antipsychotic dose (CPZ-equivalent).
#
# The group-level pipeline (05/14 -> 06) produces ONE bootstrap-averaged
# network PER GROUP, not a per-participant network -- a single participant's
# ~170-215 scans is nowhere near enough data to fit a reliable 118-node DBN
# on their own. To get a per-participant number to correlate against BPRS/
# medication, we instead compute each PDS participant's own "expression" of
# the group-level network: for every edge retained in the PDS group's
# bootstrap-thresholded network, take that participant's own correlation
# between the two involved regions (their own scans only), then average
# across all such edges. This is a standard "subject-level projection onto a
# group template" approach in individual-differences neuroimaging, and it is
# computationally cheap (no new DBN fitting) since bnlearn already fit the
# group template.
#
# Intra-slice edges use the contemporaneous (same-scan) correlation between
# the two regions. Transition edges use the same participant's own lag-1
# (scan t vs scan t+1) correlation between the two regions, mirroring what
# the DBN's transition arcs represent.
# ---------------------------------------------------------------------------

GROUPS <- c("S.E","S.R")
strength_tabs <- readRDS("results/strength_tabs_full.rds")  # from 06_aggregate_strength_full.R
raw <- readRDS("analysis/prepped_data.rds")                   # from 00_prep.R
clinical <- read.csv("clinical_data/clinical_variables_pds.csv", stringsAsFactors = FALSE)

score_participant <- function(pid, dt, edges) {
  sub <- dt[Participant_id == pid]
  setorder(sub, Scan_no)
  vals <- numeric(0)
  for (i in seq_len(nrow(edges))) {
    from_base <- edges$from_base[i]; to_base <- edges$to_base[i]
    is_trans <- edges$is_transition[i]
    if (!(from_base %in% names(sub)) || !(to_base %in% names(sub))) next
    x <- sub[[from_base]]; y <- sub[[to_base]]
    if (is_trans) {
      n <- length(x)
      if (n < 3) next
      r <- suppressWarnings(cor(x[1:(n-1)], y[2:n], use = "complete.obs"))
    } else {
      r <- suppressWarnings(cor(x, y, use = "complete.obs"))
    }
    if (!is.na(r)) vals <- c(vals, abs(r))
  }
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

results <- list()
for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  thr <- bnlearn:::threshold(s)
  avg <- s[s$strength >= thr, ]

  from_is_past <- grepl("_t_1$", avg$from); to_is_present <- grepl("_t_0$", avg$to)
  is_transition <- from_is_past & to_is_present
  from_base <- sub("_t_[01]$", "", avg$from); to_base <- sub("_t_[01]$", "", avg$to)
  # keep only intra-slice-t0 and genuine transition edges once (drop the
  # duplicated intra-slice-t1 copy so each undirected/edge pair is scored once)
  keep <- is_transition | (grepl("_t_0$", avg$from) & grepl("_t_0$", avg$to))
  edges <- data.frame(from_base = from_base[keep], to_base = to_base[keep], is_transition = is_transition[keep])
  edges <- unique(edges)
  cat(sprintf("%s: %d edges used for participant scoring (%d transition, %d intra-slice)\n",
              grp, nrow(edges), sum(edges$is_transition), sum(!edges$is_transition)))

  grp_dt <- raw[conditiongroup == grp]
  pids <- clinical$participant_id
  scores <- sapply(pids, function(p) score_participant(p, grp_dt, edges))
  results[[grp]] <- data.frame(participant_id = pids, network_score = scores)
  names(results[[grp]])[2] <- paste0("network_score_", grp)
}

merged <- clinical
for (grp in GROUPS) merged <- merge(merged, results[[grp]], by = "participant_id", all.x = TRUE)
write.csv(merged, "results/network_scores_with_clinical.csv", row.names = FALSE)

cat("\n--- Correlations: network expression score vs. clinical variables (PDS only) ---\n")
clinical_vars <- c("bprs_positive","bprs_negative","bprs_mania","bprs_depanx","total_cpz_eq")
test_rows <- list()
for (grp in GROUPS) {
  score_col <- paste0("network_score_", grp)
  for (cv in clinical_vars) {
    ok <- complete.cases(merged[[score_col]], merged[[cv]])
    if (sum(ok) < 5) next
    ct <- suppressWarnings(cor.test(merged[[score_col]][ok], merged[[cv]][ok], method = "spearman"))
    test_rows[[paste(grp, cv)]] <- data.frame(task_phase = grp, clinical_variable = cv,
                                                n = sum(ok), rho = unname(ct$estimate), p_raw = ct$p.value)
  }
}
test_df <- do.call(rbind, test_rows)
test_df$p_fdr <- p.adjust(test_df$p_raw, method = "BH")
test_df <- test_df[order(test_df$p_fdr), ]
print(test_df, row.names = FALSE)
write.csv(test_df, "results/clinical_correlation_results.csv", row.names = FALSE)
