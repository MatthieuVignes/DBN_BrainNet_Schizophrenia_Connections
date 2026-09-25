GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 50
SCORE <- "bic-g"
MAXP <- 8  # must match 16_bootstrap_parallel_nonsmoker.R

strength_tabs <- list()
for (grp in GROUPS) {
  files <- sprintf("boot_results_nonsmoker/%s_b%03d_%s_maxp%d.rds", grp, 1:B, SCORE, MAXP)
  files <- files[file.exists(files)]
  cat(sprintf("%s: %d/%d bootstrap fits found\n", grp, length(files), B))
  arcs_list <- lapply(files, function(f) { a <- readRDS(f); paste0(a$from, "->", a$to) })
  all_arcs <- unlist(arcs_list)
  tab <- table(all_arcs)
  strength_tabs[[grp]] <- data.frame(
    from = sub("->.*$", "", names(tab)), to = sub("^.*->", "", names(tab)),
    strength = as.numeric(tab) / length(files)
  )
}
saveRDS(strength_tabs, "results/strength_tabs_nonsmoker.rds")

for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  thr <- bnlearn:::threshold(s)
  cat(sprintf("%s: %d unique arcs | optimal threshold=%.3f | arcs retained=%d\n",
              grp, nrow(s), thr, sum(s$strength >= thr)))
}

# Addressable arc space -- same node set as the primary full-resolution
# analysis, so this should match the number from 06_aggregate_strength_full.R
# (paste that number in below to skip recomputing it).
raw <- readRDS("analysis/prepped_data_nonsmoker.rds")
f_dt_sample <- fold_group(raw[conditiongroup == "H.E"][1:500], size = 2)
bl <- dbnR:::create_blacklist(names(f_dt_sample), 2)
n_nodes_total <- length(names(f_dt_sample))
addr_space <- n_nodes_total * (n_nodes_total - 1) - nrow(bl)
cat("\nAddressable arc space:", addr_space, "\n\n")

# --- Bayesian density comparison, restricted to non-smokers only ---
post_draws <- list()
set.seed(42)
for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  total_present <- sum(round(s$strength * B))
  total_trials <- B * addr_space
  a_post <- 1 + total_present; b_post <- 1 + total_trials - total_present
  post_draws[[grp]] <- rbeta(200000, a_post, b_post)
  cat(sprintf("%s: posterior mean density=%.4f (95%% CrI %.4f-%.4f)\n",
              grp, a_post/(a_post+b_post), qbeta(0.025,a_post,b_post), qbeta(0.975,a_post,b_post)))
}

cat("\n--- Non-smoker-only group-difference posterior contrasts ---\n")
contrasts <- list(
  "Encoding: PDS - HC (non-smokers only)" = c("S.E","H.E"),
  "Retrieval: PDS - HC (non-smokers only)" = c("S.R","H.R")
)
for (nm in names(contrasts)) {
  g1 <- contrasts[[nm]][1]; g2 <- contrasts[[nm]][2]
  diff <- post_draws[[g1]] - post_draws[[g2]]
  ci <- quantile(diff, c(0.025, 0.975))
  cat(sprintf("%s: mean diff=%.4f, 95%% CrI [%.4f, %.4f], P(diff>0)=%.3f\n",
              nm, mean(diff), ci[1], ci[2], mean(diff > 0)))
}
saveRDS(post_draws, "results/bayesian_density_posteriors_nonsmoker.rds")
