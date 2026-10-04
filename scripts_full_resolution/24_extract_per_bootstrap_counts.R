# Extracts, from the bootstrap fits ALREADY ON DISK (from
# 14_bootstrap_parallel_runner.R -- no new model fitting needed), the total
# edge count produced by each individual bootstrap replicate.
#
# This is the input needed for a properly-calibrated group-density
# comparison (25_corrected_density_test.R), which treats each of the 100
# bootstrap replicates as ONE independent observation (real n=100 per
# group), rather than pooling every (bootstrap x possible-arc) pair into a
# single giant contingency table as the original Beta-Binomial approach
# (07_bayesian_density_full.R) did -- which Tier 2 validation showed
# produces a ~30% false-positive rate under the null instead of the nominal
# ~5%, almost certainly because it treats correlated bootstrap-arc pairs as
# if they were independent, understating uncertainty.

GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 100
SCORE <- "bic-g"
MAXP <- 8  # must match what 14_bootstrap_parallel_runner.R was run with

counts <- list()
for (grp in GROUPS) {
  files <- sprintf("boot_results_full/%s_b%03d_%s_maxp%d.rds", grp, 1:B, SCORE, MAXP)
  files <- files[file.exists(files)]
  cat(sprintf("%s: %d/%d bootstrap files found\n", grp, length(files), B))
  n_arcs <- sapply(files, function(f) nrow(readRDS(f)))
  counts[[grp]] <- n_arcs
}

max_len <- max(sapply(counts, length))
counts_df <- data.frame(bootstrap = seq_len(max_len))
for (grp in GROUPS) {
  v <- counts[[grp]]
  length(v) <- max_len  # pad with NA if any group has fewer completed files
  counts_df[[grp]] <- v
}

cat("\n--- Per-bootstrap edge count summary ---\n")
for (grp in GROUPS) cat(sprintf("%s: mean=%.1f sd=%.1f min=%d max=%d\n",
                                  grp, mean(counts[[grp]]), sd(counts[[grp]]),
                                  min(counts[[grp]]), max(counts[[grp]])))

write.csv(counts_df, "results/per_bootstrap_edge_counts.csv", row.names = FALSE)
cat("\nSaved per_bootstrap_edge_counts.csv -- send this file back for the corrected test.\n")
