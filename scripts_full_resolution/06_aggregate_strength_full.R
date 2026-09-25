GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 100
SCORE <- "bic-g"
MAXP <- 8  # must match what you used in 14_bootstrap_parallel_runner.R

strength_tabs <- list()
for (grp in GROUPS) {
  files <- sprintf("boot_results_full/%s_b%03d_%s_maxp%d.rds", grp, 1:B, SCORE, MAXP)
  files <- files[file.exists(files)]
  cat(sprintf("%s: %d/%d bootstrap fits found\n", grp, length(files), B))
  arcs_list <- lapply(files, function(f) {
    a <- readRDS(f); paste0(a$from, "->", a$to)
  })
  all_arcs <- unlist(arcs_list)
  tab <- table(all_arcs)
  strength_tabs[[grp]] <- data.frame(
    from = sub("->.*$", "", names(tab)),
    to   = sub("^.*->", "", names(tab)),
    strength = as.numeric(tab) / length(files)
  )
}
saveRDS(strength_tabs, "results/strength_tabs_full.rds")

for (grp in GROUPS) {
  s <- strength_tabs[[grp]]
  thr <- bnlearn:::threshold(s)
  n_retained <- sum(s$strength >= thr)
  cat(sprintf("%s: %d unique arcs observed | optimal threshold=%.3f | arcs retained=%d\n",
              grp, nrow(s), thr, n_retained))
}

# Addressable arc space size for this node set (needed for the Bayesian
# density model in 07) -- recompute since full resolution has ~118 nodes,
# not the 26 used in the macro-ROI pilot.
raw <- readRDS("analysis/prepped_data.rds")
vars <- vars_for_model(raw[conditiongroup == "H.E"][1:500])
f_dt_sample <- fold_group(raw[conditiongroup == "H.E"][1:500], size = 2)
bl <- dbnR:::create_blacklist(names(f_dt_sample), 2)
n_nodes_total <- length(names(f_dt_sample))
addr_space <- n_nodes_total * (n_nodes_total - 1) - nrow(bl)
cat("\nAddressable arc space (for 07_bayesian_density.R's ADDR_SPACE):", addr_space, "\n")
