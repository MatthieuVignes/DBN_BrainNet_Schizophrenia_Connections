# Group comparison of network density from the per-bootstrap edge counts.
# Each group's bootstrap distribution of the edge count approximates the
# sampling distribution of that statistic. The difference of two independent
# groups is therefore summarised by all B x B pairwise differences:
# percentile 95% interval, plus a normal-approximation z/p with
# SE = sqrt(var1 + var2). 
# Usage: Rscript 25_bootstrap_difference_test.R [COUNTS_CSV] [OUT_CSV] [ADDR_SPACE]
args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "results/per_bootstrap_edge_counts.csv"
OUT <- if (length(args) >= 2) args[2] else "results/bootstrap_difference_test_results.csv"
ADDR <- if (length(args) >= 3) as.numeric(args[3]) else 13688  # addressable arc space (06_aggregate output)
d <- read.csv(IN)
contrasts <- list(
  "Encoding: PDS - HC (S.E - H.E)"        = c("S.E","H.E"),
  "Retrieval: PDS - HC (S.R - H.R)"       = c("S.R","H.R"),
  "HC: Retrieval - Encoding (H.R - H.E)"  = c("H.R","H.E"),
  "PDS: Retrieval - Encoding (S.R - S.E)" = c("S.R","S.E"))
res <- do.call(rbind, lapply(names(contrasts), function(nm) {
  x <- na.omit(d[[contrasts[[nm]][1]]]); y <- na.omit(d[[contrasts[[nm]][2]]])
  dif <- as.vector(outer(x, y, "-"))
  ci <- quantile(dif, c(.025, .975)); se <- sqrt(var(x) + var(y))
  z <- (mean(x) - mean(y)) / se
  data.frame(contrast = nm, n1 = length(x), n2 = length(y),
    mean_diff_edges = mean(dif), ci_lower = ci[1], ci_upper = ci[2],
    se = se, z = z, p_normal = 2 * pnorm(-abs(z)),
    excludes_zero = ci[1] > 0 | ci[2] < 0,
    mean_diff_density = mean(dif) / ADDR, row.names = NULL)
}))
print(res, digits = 4)
write.csv(res, OUT, row.names = FALSE)
cat("Saved", OUT, "\n")
