suppressMessages({source("scripts_full_resolution/01_helpers.R")})
mg <- readRDS("analysis/prepped_data.rds")
GROUPS <- c("H.E","H.R","S.E","S.R")
strength_tabs <- readRDS("results/strength_tabs_full.rds")
results <- list()
for (grp in GROUPS) {
  d <- mg[conditiongroup == grp]
  vars <- vars_for_model(d)
  X <- as.matrix(d[, ..vars])
  R <- cor(X); n <- nrow(X)
  pairs <- combn(vars, 2)
  z <- atanh(R[t(pairs)]); se <- 1/sqrt(n-3)
  pval <- 2*pnorm(-abs(z/se)); padj <- p.adjust(pval, method="BH"); sig <- padj < 0.05
  corr_edges <- sum(sig)
  s <- strength_tabs[[grp]]; thr <- bnlearn:::threshold(s); avg <- s[s$strength>=thr,]
  intra0 <- avg[grepl("_t_0$",avg$from) & grepl("_t_0$",avg$to),]
  base_from <- sub("_t_0$","",intra0$from); base_to <- sub("_t_0$","",intra0$to)
  skeleton <- unique(t(apply(cbind(base_from,base_to),1,sort)))
  dbn_edges <- nrow(skeleton)
  corr_pairs_sig <- t(pairs[,sig,drop=FALSE]); corr_pairs_sig <- unique(t(apply(corr_pairs_sig,1,sort)))
  corr_set <- apply(corr_pairs_sig,1,paste,collapse="|"); dbn_set <- apply(skeleton,1,paste,collapse="|")
  overlap <- length(intersect(corr_set,dbn_set)); jaccard <- overlap/length(union(corr_set,dbn_set))
  results[[grp]] <- data.frame(group=grp, correlation_edges=corr_edges, dbn_static_edges=dbn_edges,
                                overlap=overlap, jaccard=round(jaccard,3))
  cat(sprintf("%s: Pearson-FDR=%d | DBN static=%d | overlap=%d | Jaccard=%.3f\n", grp, corr_edges, dbn_edges, overlap, jaccard))
}
write.csv(do.call(rbind, results), "results/correlation_vs_dbn_comparison_full.csv", row.names=FALSE)
