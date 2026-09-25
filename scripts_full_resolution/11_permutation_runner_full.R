suppressMessages({source("scripts_full_resolution/01_helpers.R"); library(parallel)})
mg <- readRDS("analysis/prepped_data.rds")
N_PERM <- 50
SCORE <- "bic-g"
MAXP <- 8
N_CORES <- 7
dir.create("perm_results_full", showWarnings = FALSE)
contrasts <- list(encoding = c("H.E","S.E"), retrieval = c("H.R","S.R"))

for (grp in unique(unlist(contrasts))) {
  outfile <- sprintf("perm_results_full/observed_%s.rds", grp)
  if (!file.exists(outfile)) {
    d <- mg[conditiongroup == grp]
    net <- fit_group_dbn_hc(d, size=2, score=SCORE, maxp=MAXP)
    saveRDS(extract_arcs(net), outfile)
  }
}

run_one_perm <- function(contrast_name, p, mg, groups2) {
  outfile <- sprintf("perm_results_full/%s_perm%03d.rds", contrast_name, p)
  if (file.exists(outfile)) return(invisible(NULL))
  set.seed(20000 + 1000*which(names(contrasts)==contrast_name) + p)
  pooled <- mg[conditiongroup %in% groups2]
  parts <- unique(pooled[, .(Participant_id, conditiongroup)])
  shuffled_labels <- sample(parts$conditiongroup)
  parts[, perm_label := shuffled_labels]
  pooled2 <- merge(pooled, parts[, .(Participant_id, perm_label)], by="Participant_id")
  dA <- pooled2[perm_label==groups2[1]]; dB <- pooled2[perm_label==groups2[2]]
  dA[, perm_label:=NULL]; dB[, perm_label:=NULL]
  netA <- fit_group_dbn_hc(dA, size=2, score=SCORE, maxp=MAXP)
  netB <- fit_group_dbn_hc(dB, size=2, score=SCORE, maxp=MAXP)
  aA <- extract_arcs(netA); aA[,arc:=paste0(from,"->",to)]
  aB <- extract_arcs(netB); aB[,arc:=paste0(from,"->",to)]
  edge_diff <- abs(nrow(aA)-nrow(aB))
  jacc <- length(intersect(aA$arc,aB$arc))/length(union(aA$arc,aB$arc))
  saveRDS(list(edge_diff=edge_diff, jaccard=jacc, nA=nrow(aA), nB=nrow(aB)), outfile)
  invisible(NULL)
}

jobs <- do.call(rbind, lapply(names(contrasts), function(cn) data.frame(cname=cn, p=1:N_PERM)))
jobs <- jobs[!file.exists(sprintf("perm_results_full/%s_perm%03d.rds", jobs$cname, jobs$p)), ]
cat("Permutation jobs remaining:", nrow(jobs), "\n")
if (nrow(jobs) > 0) {
  invisible(mclapply(seq_len(nrow(jobs)), function(i) {
    run_one_perm(jobs$cname[i], jobs$p[i], mg, contrasts[[jobs$cname[i]]])
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
}
cat("Done.\n")
