# Run this FIRST on your machine. It fits a single (non-bootstrap) DBN on
# one group at full 116-region resolution, at a few different `maxp` values,
# so you can see real timing on your hardware before committing to the full
# 400-fit bootstrap run. Takes the results as a guide to choose `maxp` and
# estimate total runtime (see printed extrapolation at the end).

raw <- readRDS("analysis/prepped_data.rds")
grp <- raw[conditiongroup == "H.E"]
cat("Group H.E:", nrow(grp), "rows,", length(vars_for_model(grp)), "nodes (unrolled:",
    2 * length(vars_for_model(grp)), ")\n\n")

f_dt <- fold_group(grp, size = 2)
cat("Folding done. Benchmarking fits at different maxp values...\n\n")

for (mp in c(4, 6, 8, 12)) {
  t0 <- Sys.time()
  net <- fit_group_dbn_hc(grp, size = 2, score = "bic-g", f_dt = f_dt, maxp = mp)
  t1 <- Sys.time()
  secs <- as.numeric(t1 - t0, units = "secs")
  cat(sprintf("maxp=%2d | fit time=%6.1fs | n_arcs=%d\n", mp, secs, nrow(bnlearn::arcs(net))))
  cat(sprintf("  -> extrapolated full run (400 fits, %d cores parallel): ~%.1f min\n\n",
              8L, secs * 400 / 8 / 60))
}

cat("If you want to also try unbounded (maxp = Inf, the manuscript's nominal\n")
cat("method with no cap), uncomment the block below -- but only after seeing\n")
cat("the capped timings above; it may take substantially longer or not finish\n")
cat("in a reasonable time at 116 regions.\n\n")
cat("# t0 <- Sys.time()\n")
cat("# net_inf <- fit_group_dbn_hc(grp, size=2, score='bic-g', f_dt=f_dt, maxp=Inf)\n")
cat("# cat('Unbounded fit time:', as.numeric(Sys.time()-t0, units='secs'), 's\\n')\n")
