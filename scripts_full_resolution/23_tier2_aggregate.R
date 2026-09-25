N_NODES <- 26
B_BOOT <- 25
N_REPLICATES <- 30

# Addressable arc space for this node count -- same computation as
# 06_aggregate_strength_full.R, needed for the Beta-Binomial density model.
f_dt_sample <- {
  node_names <- paste0("Region_", sprintf("%02d", 1:N_NODES))
  dummy <- data.table::data.table(Participant_id = rep("p1", 20), Scan_no = as.double(1:20))
  for (nm in node_names) dummy[[nm]] <- rnorm(20)
  fold_group(dummy, size = 2)
}
bl <- dbnR:::create_blacklist(names(f_dt_sample), 2)
n_nodes_total <- length(names(f_dt_sample))
ADDR_SPACE <- n_nodes_total * (n_nodes_total - 1) - nrow(bl)
cat("Addressable arc space:", ADDR_SPACE, "\n\n")

density_posterior <- function(scenario, replicate, group_label) {
  files <- sprintf("tier2_boot/%s_r%03d_%s_b%03d.rds", scenario, replicate, group_label, 1:B_BOOT)
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  arcs_list <- lapply(files, function(f) { a <- readRDS(f); paste0(a$from, "->", a$to) })
  tab <- table(unlist(arcs_list))
  total_present <- sum(as.numeric(tab))
  total_trials <- length(files) * ADDR_SPACE
  a_post <- 1 + total_present; b_post <- 1 + total_trials - total_present
  list(a = a_post, b = b_post, n_boot = length(files))
}

evaluate_replicate <- function(scenario, replicate) {
  postA <- density_posterior(scenario, replicate, "A")
  postB <- density_posterior(scenario, replicate, "B")
  if (is.null(postA) || is.null(postB)) return(NULL)
  set.seed(replicate)
  drawsA <- rbeta(20000, postA$a, postA$b)
  drawsB <- rbeta(20000, postB$a, postB$b)
  diff <- drawsB - drawsA  # B - A: in the effect scenario, B has MORE true edges, so a
                            # correctly-powered test should find diff credibly > 0
  ci <- quantile(diff, c(0.025, 0.975))
  credible <- (ci[1] > 0) || (ci[2] < 0)
  correct_direction <- credible && (ci[1] > 0)  # for effect scenario: B > A is the true direction
  data.frame(scenario = scenario, replicate = replicate,
             mean_diff = mean(diff), ci_lower = ci[1], ci_upper = ci[2],
             credible = credible, correct_direction = correct_direction)
}

results <- list()
for (scenario in c("null","effect")) {
  for (r in 1:N_REPLICATES) {
    res <- evaluate_replicate(scenario, r)
    if (!is.null(res)) results[[paste(scenario, r)]] <- res
  }
}
results_df <- do.call(rbind, results)
write.csv(results_df, "results/tier2_results.csv", row.names = FALSE)

cat("--- Tier 2 results ---\n\n")
null_res <- results_df[results_df$scenario == "null", ]
effect_res <- results_df[results_df$scenario == "effect", ]

fpr <- mean(null_res$credible)
cat(sprintf("NULL scenario (no true difference): %d/%d replicates completed\n", nrow(null_res), N_REPLICATES))
cat(sprintf("  False-positive rate (credible difference found when none exists): %.3f\n", fpr))
cat(sprintf("  (binomial 95%% CI on this rate: [%.3f, %.3f])\n",
            fpr - 1.96*sqrt(fpr*(1-fpr)/nrow(null_res)), fpr + 1.96*sqrt(fpr*(1-fpr)/nrow(null_res))))
cat("  A well-calibrated test should show this near the nominal ~0.05.\n\n")

power <- mean(effect_res$correct_direction)
cat(sprintf("EFFECT scenario (known true difference): %d/%d replicates completed\n", nrow(effect_res), N_REPLICATES))
cat(sprintf("  Power (correctly detected a credible difference in the right direction): %.3f\n", power))
cat(sprintf("  (binomial 95%% CI on this rate: [%.3f, %.3f])\n",
            power - 1.96*sqrt(power*(1-power)/nrow(effect_res)), power + 1.96*sqrt(power*(1-power)/nrow(effect_res))))
cat(sprintf("  Any credible difference found in the WRONG direction: %d/%d replicates\n",
            sum(effect_res$credible & !effect_res$correct_direction), nrow(effect_res)))
