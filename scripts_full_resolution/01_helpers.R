suppressMessages({
  library(data.table)
  library(bnlearn)
  library(dbnR)
})

# vars_for_model: which columns become DBN nodes. Scan_no IS included here
# (unlike the reduced pilot) to match the manuscript's stated use of scan
# number as a control node for row/time dependence.
vars_for_model <- function(dt) {
  setdiff(names(dt), c("Participant_id","inst_id","trialtype2","conditiongroup"))
}

fold_group <- function(dt, size = 2, id_col = "Participant_id") {
  vars <- vars_for_model(dt)
  ids <- unique(dt[[id_col]])
  out <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    sub <- dt[dt[[id_col]] == ids[i], ..vars]
    out[[i]] <- fold_dt(sub, size = size)
  }
  rbindlist(out)
}

extract_arcs <- function(net) {
  a <- bnlearn::arcs(net)
  data.table(from = a[,1], to = a[,2])
}

# Two-slice DBN fit, faithful to the manuscript's stated hill-climbing
# procedure (static network -> transition network -> merge), with an added
# `maxp` cap on parents-per-node. maxp bounds the combinatorial search space,
# which is what makes 116-region / 236-unrolled-node fitting tractable at
# all. maxp = Inf reproduces unbounded hill-climbing exactly (manuscript's
# nominal method) but may be very slow or fail to finish at this scale --
# benchmark on your machine (see 13_benchmark_full_res.R) before committing
# to a full bootstrap run at maxp = Inf.
fit_group_dbn_hc <- function(dt, size = 2, score = "bic-g", f_dt = NULL,
                              restart = 0, perturb = 1, id_col = "Participant_id",
                              maxp = Inf) {
  vars <- vars_for_model(dt)
  base <- dt[, ..vars]
  base_r <- dbnR:::time_rename(data.table::copy(base))
  net0 <- bnlearn::hc(base_r, score = score, restart = restart, perturb = perturb, maxp = maxp)

  if (is.null(f_dt)) f_dt <- fold_group(dt, size = size, id_col = id_col)
  blacklist <- dbnR:::create_blacklist(names(f_dt), size)
  net <- bnlearn::hc(f_dt, blacklist = blacklist, score = score,
                      restart = restart, perturb = perturb, maxp = maxp)

  bnlearn::arcs(net) <- dbnR:::merge_nets(net0, net, size)
  class(net) <- c("dbn", class(net))
  net
}
