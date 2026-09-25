suppressMessages({library(data.table)})
N_PERM <- 50
contrasts <- list(encoding=c("H.E","S.E"), retrieval=c("H.R","S.R"))
for (cname in names(contrasts)) {
  g <- contrasts[[cname]]
  oA <- readRDS(sprintf("perm_results_full/observed_%s.rds", g[1])); oA$arc <- paste0(oA$from,"->",oA$to)
  oB <- readRDS(sprintf("perm_results_full/observed_%s.rds", g[2])); oB$arc <- paste0(oB$from,"->",oB$to)
  obs_edge_diff <- abs(nrow(oA)-nrow(oB))
  obs_jaccard <- length(intersect(oA$arc,oB$arc))/length(union(oA$arc,oB$arc))
  null_edge_diff <- sapply(1:N_PERM, function(p) readRDS(sprintf("perm_results_full/%s_perm%03d.rds",cname,p))$edge_diff)
  null_jaccard <- sapply(1:N_PERM, function(p) readRDS(sprintf("perm_results_full/%s_perm%03d.rds",cname,p))$jaccard)
  cat(sprintf("[%s] observed |edge diff|=%d, p=%.3f | observed Jaccard=%.3f, p=%.3f\n",
      cname, obs_edge_diff, mean(null_edge_diff>=obs_edge_diff), obs_jaccard, mean(null_jaccard<=obs_jaccard)))
}
