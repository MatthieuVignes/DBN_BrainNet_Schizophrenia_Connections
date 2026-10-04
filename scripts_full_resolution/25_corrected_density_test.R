# Corrected, properly-calibrated group-density comparison, addressing the
# miscalibration Tier 2 validation revealed in the original Beta-Binomial
# approach (07_bayesian_density_full.R). Uses per_bootstrap_edge_counts.csv
# from 24_extract_per_bootstrap_counts.R -- no new model fitting required,
# this re-analyzes bootstrap fits you already have.
#
# Each of the ~100 bootstrap replicates' own total edge count is treated as
# ONE independent observation (real n ~ 100 per group), rather than pooling
# every (bootstrap x possible-arc) pair into a single contingency table.
# This gives a real, defensible sample size instead of an inflated one, at
# the cost of a coarser statistic (total edge count per replicate, not
# individual-arc presence/absence) -- exactly the same statistic the
# permutation test (11/12) already used, but now with ~100 replicate draws
# per group instead of just one fit per group.
#
# Reports both a Wilcoxon rank-sum test (distribution-free, no normality
# assumption) and a simple percentile-based bootstrap CI on the difference
# in means, so the result doesn't depend on one test's assumptions alone.

counts <- read.csv("results/per_bootstrap_edge_counts.csv")

contrasts <- list(
  "Encoding: PDS - HC (S.E - H.E)" = c("S.E", "H.E"),
  "Retrieval: PDS - HC (S.R - H.R)" = c("S.R", "H.R"),
  "HC: Retrieval - Encoding (H.R - H.E)" = c("H.R", "H.E"),
  "PDS: Retrieval - Encoding (S.R - S.E)" = c("S.R", "S.E")
)

cat("--- Corrected group-density comparison (per-bootstrap-replicate as the unit) ---\n\n")

set.seed(1)
results <- list()
for (nm in names(contrasts)) {
  g1 <- contrasts[[nm]][1]; g2 <- contrasts[[nm]][2]
  x <- counts[[g1]][!is.na(counts[[g1]])]
  y <- counts[[g2]][!is.na(counts[[g2]])]

  wt <- suppressWarnings(wilcox.test(x, y, conf.int = TRUE))

  # percentile bootstrap CI on the difference in means, resampling
  # replicates WITH replacement within each group (a bootstrap-of-the-
  # bootstrap, appropriate since x and y are themselves already bootstrap
  # replicate statistics -- this resamples which of those replicates we
  # "believe", not the underlying participants again)
  diffs <- replicate(5000, {
    mean(sample(x, length(x), replace = TRUE)) - mean(sample(y, length(y), replace = TRUE))
  })
  ci <- quantile(diffs, c(0.025, 0.975))

  cat(sprintf("%s:\n  n=%d/%d | mean edge count: %.1f vs %.1f\n", nm, length(x), length(y), mean(x), mean(y)))
  cat(sprintf("  Wilcoxon p=%.4f | location-shift 95%% CI [%.2f, %.2f]\n", wt$p.value, wt$conf.int[1], wt$conf.int[2]))
  cat(sprintf("  Bootstrap-of-bootstrap: mean diff=%.2f | 95%% CI [%.2f, %.2f] | credible: %s\n\n",
              mean(diffs), ci[1], ci[2], ifelse(ci[1] > 0 || ci[2] < 0, "YES", "no")))

  results[[nm]] <- data.frame(contrast = nm, n1 = length(x), n2 = length(y),
                                mean1 = mean(x), mean2 = mean(y), wilcoxon_p = wt$p.value,
                                boot_mean_diff = mean(diffs), boot_ci_lower = ci[1], boot_ci_upper = ci[2])
}

write.csv(do.call(rbind, results), "results/corrected_density_test_results.csv", row.names = FALSE)
cat("Saved corrected_density_test_results.csv\n")
