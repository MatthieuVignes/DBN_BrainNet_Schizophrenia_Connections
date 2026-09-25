suppressMessages({library(data.table)})
GROUPS <- c("H.E","H.R","S.E","S.R")
strength_tabs <- readRDS("results/strength_tabs_full.rds")
summary_rows <- list()
for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  thr <- bnlearn:::threshold(s)
  avg <- s[s$strength >= thr, ]
  # dbnR convention: "_t_0" = present, "_t_1" = one lag in the past;
  # allowed causal direction is past -> present, i.e. "_t_1" -> "_t_0".
  is_transition <- grepl("_t_1$", avg$from) & grepl("_t_0$", avg$to)
  n_total <- nrow(avg); n_trans <- sum(is_transition); n_intra <- n_total - n_trans
  summary_rows[[grp]] <- data.frame(group=grp, n_total=n_total, n_intra=n_intra,
                                     n_transition=n_trans, pct_transition=100*n_trans/n_total)
  cat(sprintf("%s: total=%d | intra=%d | transition=%d (%.1f%%)\n", grp, n_total, n_intra, n_trans, 100*n_trans/n_total))
}
write.csv(do.call(rbind, summary_rows), "results/directed_edge_summary_full.csv", row.names=FALSE)
