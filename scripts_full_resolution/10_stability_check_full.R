suppressMessages({source("scripts_full_resolution/01_helpers.R")})
mg <- readRDS("analysis/prepped_data.rds")
GROUPS <- c("H.E","H.R","S.E","S.R")
MAXP <- 8  # match your bootstrap run
set.seed(99)
stab_results <- list()
for (grp in GROUPS) {
  grp_dt <- mg[conditiongroup == grp]
  parts <- unique(grp_dt$Participant_id)
  shuf <- sample(parts)
  half1 <- shuf[1:(length(shuf)/2)]; half2 <- shuf[(length(shuf)/2+1):length(shuf)]
  d1 <- grp_dt[Participant_id %in% half1]; d2 <- grp_dt[Participant_id %in% half2]
  net1 <- fit_group_dbn_hc(d1, size=2, score="bic-g", maxp=MAXP)
  net2 <- fit_group_dbn_hc(d2, size=2, score="bic-g", maxp=MAXP)
  a1 <- extract_arcs(net1); a1[, arc:=paste0(from,"->",to)]
  a2 <- extract_arcs(net2); a2[, arc:=paste0(from,"->",to)]
  ov <- length(intersect(a1$arc,a2$arc)); jac <- ov/length(union(a1$arc,a2$arc))
  stab_results[[grp]] <- data.frame(group=grp, arcs_half1=nrow(a1), arcs_half2=nrow(a2), overlap=ov, jaccard=round(jac,3))
  cat(sprintf("%s: half1 arcs=%d | half2 arcs=%d | Jaccard=%.3f\n", grp, nrow(a1), nrow(a2), jac))
}
write.csv(do.call(rbind, stab_results), "results/split_half_stability_full.csv", row.names=FALSE)
