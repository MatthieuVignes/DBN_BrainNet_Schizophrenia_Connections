# Extract the edge count of every bootstrap fit
# Usage: Rscript 24_extract_per_bootstrap_counts.R [BOOT_DIR] [B] [OUT_CSV]
#  main analysis : Rscript 24_extract_per_bootstrap_counts.R
#  non-smokers   : Rscript 24_extract_per_bootstrap_counts.R boot_results_nonsmoker 50 results/per_bootstrap_edge_counts_nonsmoker.csv
args <- commandArgs(trailingOnly = TRUE)
BOOT_DIR <- if (length(args) >= 1) args[1] else "boot_results_full"
B        <- if (length(args) >= 2) as.integer(args[2]) else 100L
OUT_CSV  <- if (length(args) >= 3) args[3] else "results/per_bootstrap_edge_counts.csv"
GROUPS <- c("H.E","H.R","S.E","S.R"); SCORE <- "bic-g"; MAXP <- 8

counts <- lapply(GROUPS, function(g) {
  f <- sprintf("%s/%s_b%03d_%s_maxp%d.rds", BOOT_DIR, g, seq_len(B), SCORE, MAXP)
  f <- f[file.exists(f)]
  cat(sprintf("%s: %d/%d files\n", g, length(f), B))
  vapply(f, function(x) nrow(readRDS(x)), numeric(1))
})
names(counts) <- GROUPS
n <- max(lengths(counts))
df <- data.frame(bootstrap = seq_len(n))
for (g in GROUPS) { v <- unname(counts[[g]]); length(v) <- n; df[[g]] <- v }
for (g in GROUPS) cat(sprintf("%s: mean=%.1f sd=%.1f\n", g, mean(df[[g]], na.rm=TRUE), sd(df[[g]], na.rm=TRUE)))
write.csv(df, OUT_CSV, row.names = FALSE)
cat("Saved", OUT_CSV, "\n")
