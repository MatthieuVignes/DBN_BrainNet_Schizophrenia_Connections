suppressMessages({library(data.table)})
GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 100
ADDR_SPACE <- addr_space  # inherits the number printed by 06_aggregate_strength_full.R

if (is.null(ADDR_SPACE)) stop("Set ADDR_SPACE from the output of 06_aggregate_strength_full.R first.")

strength_tabs <- readRDS("results/strength_tabs_full.rds")
post_draws <- list()
set.seed(42)
for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  total_present <- sum(round(s$strength * B))
  total_trials <- B * ADDR_SPACE
  a_post <- 1 + total_present; b_post <- 1 + total_trials - total_present
  post_draws[[grp]] <- rbeta(200000, a_post, b_post)
  cat(sprintf("%s: posterior mean density=%.4f (95%% CrI %.4f-%.4f)\n",
              grp, a_post/(a_post+b_post), qbeta(0.025,a_post,b_post), qbeta(0.975,a_post,b_post)))
}
cat("\n--- Group-difference posterior contrasts ---\n")
contrasts <- list(
  "Encoding: PDS - HC (S.E - H.E)" = c("S.E","H.E"),
  "Retrieval: PDS - HC (S.R - H.R)" = c("S.R","H.R"),
  "HC: Retrieval - Encoding (H.R - H.E)" = c("H.R","H.E"),
  "PDS: Retrieval - Encoding (S.R - S.E)" = c("S.R","S.E")
)
for (nm in names(contrasts)) {
  g1 <- contrasts[[nm]][1]; g2 <- contrasts[[nm]][2]
  diff <- post_draws[[g1]] - post_draws[[g2]]
  ci <- quantile(diff, c(0.025, 0.975))
  cat(sprintf("%s: mean diff=%.4f, 95%% CrI [%.4f, %.4f], P(diff>0)=%.3f\n",
              nm, mean(diff), ci[1], ci[2], mean(diff > 0)))
}
saveRDS(post_draws, "results/bayesian_density_posteriors_full.rds")
