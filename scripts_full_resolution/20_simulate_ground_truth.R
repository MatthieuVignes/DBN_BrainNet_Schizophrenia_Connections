suppressMessages({library(data.table)})
# Simulates data from a KNOWN two-slice DBN structure, using the same node
# conventions as the real pipeline (a Scan_no column, Participant_id-grouped
# time series), so the exact same fit_group_dbn_hc()/fold_group() functions
# in 01_helpers.R can be applied unmodified to simulated data.
#
# Ground-truth topology: a stochastic block model (modular, sparse) rather
# than a uniformly random graph -- brain networks are modular, and a random
# graph would be an easy, unrealistic target that invites the objection that
# the simulation doesn't resemble what the method is actually used on.
#
# Ground-truth dynamics: contemporaneous (intra-slice) structural equations,
# solved in topological order, plus lagged (transition, t-1 -> t) effects,
# plus AR(1)-autocorrelated innovations (mimicking hemodynamic-response
# smoothing in real BOLD signal) and between-participant random-effect
# jitter on every true coefficient (so no two simulated participants have
# identical dynamics, matching real between-subject heterogeneity).
#
# Scan_no is included as a genuine node with NO true edges to or from any
# other node -- exactly mirroring its role in the real analysis (a control
# variable, not part of the network) -- so any edges the recovery pipeline
# finds involving Scan_no are, by construction, false positives. This is a
# built-in specificity check, not just a convenience column.

#' Try to estimate a realistic AR(1) coefficient and innovation SD from the
#' user's own real full-resolution data, if available. Falls back to a
#' documented literature-typical default (fMRI frame-to-frame autocorrelation
#' is commonly reported in the ~0.2-0.5 range at TR ~ 1-2s) if the real data
#' isn't found -- so this script works standalone but calibrates itself
#' automatically when run alongside the real pipeline's prepped_data.rds.
estimate_ar1_params <- function(path = "analysis/prepped_data.rds", n_regions_sample = 10) {
  default <- list(phi = 0.35, innovation_sd = 1.0, source = "literature default (no real data found)")
  if (!file.exists(path)) return(default)
  tryCatch({
    raw <- readRDS(path)
    region_cols <- setdiff(names(raw), c("Participant_id","DVARS","Scan_no","trialtype2","conditiongroup"))
    sample_cols <- sample(region_cols, min(n_regions_sample, length(region_cols)))
    one_participant <- raw[Participant_id == raw$Participant_id[1]]
    setorder(one_participant, Scan_no)
    phis <- sapply(sample_cols, function(cl) {
      x <- one_participant[[cl]]
      if (length(x) < 10) return(NA_real_)
      acf_val <- cor(x[-length(x)], x[-1], use = "complete.obs")
      acf_val
    })
    phi_est <- mean(phis, na.rm = TRUE)
    if (is.na(phi_est) || phi_est <= 0 || phi_est >= 0.95) return(default)
    list(phi = phi_est, innovation_sd = 1.0,
         source = sprintf("estimated from real data (mean lag-1 autocorrelation over %d sampled regions)",
                           length(sample_cols)))
  }, error = function(e) default)
}

#' Build a ground-truth two-slice DBN structure.
#' @param node_names character vector of node names EXCLUDING Scan_no
#'   (Scan_no is added automatically as an edge-free node).
#' @param n_modules number of blocks for the stochastic block model.
#' @param p_intra_within / p_intra_between: contemporaneous edge probability
#'   within vs. between modules.
#' @param p_trans_within / p_trans_between: transition (lag-1) edge
#'   probability within vs. between modules.
#' @param coef_range absolute value range for true structural coefficients.
build_ground_truth <- function(node_names, n_modules = 4,
                                p_intra_within = 0.30, p_intra_between = 0.03,
                                p_trans_within = 0.12, p_trans_between = 0.01,
                                coef_range = c(0.15, 0.40), seed = 1) {
  set.seed(seed)
  n <- length(node_names)
  module_id <- sample(rep(1:n_modules, length.out = n))
  topo_order <- sample(node_names)  # random topological order for intra-slice DAG

  B0 <- matrix(0, n, n, dimnames = list(node_names, node_names))  # intra-slice: B0[from, to]
  B1 <- matrix(0, n, n, dimnames = list(node_names, node_names))  # transition: B1[from(t-1), to(t)]
  names(module_id) <- node_names

  rand_coef <- function() sample(c(-1, 1), 1) * runif(1, coef_range[1], coef_range[2])

  # Intra-slice: only allow edges forward in topo_order (guarantees acyclic)
  pos <- setNames(seq_along(topo_order), topo_order)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      from <- topo_order[i]; to <- topo_order[j]
      if (pos[from] >= pos[to]) next
      same_mod <- module_id[from] == module_id[to]
      p <- if (same_mod) p_intra_within else p_intra_between
      if (runif(1) < p) B0[from, to] <- rand_coef()
    }
  }

  # Transition edges: any ordered pair (including self-loops, i.e. AR
  # continuity of a region onto itself -- realistic for BOLD signal), no
  # acyclicity constraint needed since time itself enforces direction.
  for (from in node_names) {
    for (to in node_names) {
      same_mod <- module_id[from] == module_id[to]
      p <- if (from == to) 0.6 else if (same_mod) p_trans_within else p_trans_between
      if (runif(1) < p) B1[from, to] <- if (from == to) runif(1, 0.2, 0.5) else rand_coef()
    }
  }

  list(node_names = node_names, module_id = module_id, topo_order = topo_order,
       B0 = B0, B1 = B1,
       n_intra_edges = sum(B0 != 0), n_trans_edges = sum(B1 != 0))
}

#' Simulate one participant's time series from a ground-truth structure.
simulate_participant <- function(gt, n_timepoints, participant_id,
                                  ar_params = list(phi = 0.35, innovation_sd = 1.0),
                                  jitter_sd = 0.05) {
  node_names <- gt$node_names; n <- length(node_names)
  topo_order <- gt$topo_order

  # per-participant random-effect jitter on every true (nonzero) coefficient
  B0_p <- gt$B0; B1_p <- gt$B1
  B0_p[B0_p != 0] <- B0_p[B0_p != 0] + rnorm(sum(B0_p != 0), 0, jitter_sd)
  B1_p[B1_p != 0] <- B1_p[B1_p != 0] + rnorm(sum(B1_p != 0), 0, jitter_sd)

  phi <- ar_params$phi; sd_innov <- ar_params$innovation_sd
  X <- matrix(0, n_timepoints, n, dimnames = list(NULL, node_names))
  E <- matrix(0, n_timepoints, n, dimnames = list(NULL, node_names))  # AR(1) innovations

  E[1, ] <- rnorm(n, 0, sd_innov)
  for (t in 2:n_timepoints) E[t, ] <- phi * E[t - 1, ] + sqrt(1 - phi^2) * rnorm(n, 0, sd_innov)

  for (t in seq_len(n_timepoints)) {
    x_t <- setNames(numeric(n), node_names)
    lag_contrib <- if (t == 1) setNames(numeric(n), node_names) else
      setNames(as.numeric(t(B1_p) %*% X[t - 1, ]), node_names)
    for (nd in topo_order) {
      parents <- names(which(B0_p[, nd] != 0))
      intra_contrib <- if (length(parents) == 0) 0 else sum(B0_p[parents, nd] * x_t[parents])
      x_t[nd] <- intra_contrib + lag_contrib[nd] + E[t, nd]
    }
    X[t, ] <- x_t
  }

  dt <- as.data.table(X)
  dt[, Participant_id := participant_id]
  dt[, Scan_no := as.double(seq_len(n_timepoints))]
  dt
}

#' Simulate a full group (multiple participants).
simulate_group <- function(gt, n_participants, n_timepoints, group_label,
                            ar_params = list(phi = 0.35, innovation_sd = 1.0),
                            jitter_sd = 0.05, seed = 1) {
  set.seed(seed)
  parts <- lapply(seq_len(n_participants), function(i) {
    simulate_participant(gt, n_timepoints, participant_id = sprintf("sim-%s-%03d", group_label, i),
                          ar_params = ar_params, jitter_sd = jitter_sd)
  })
  rbindlist(parts)
}

#' Ground-truth edge tables in the same "arc" format extract_arcs() produces,
#' for direct comparison against a recovered network's arcs.
ground_truth_arcs <- function(gt) {
  intra <- which(gt$B0 != 0, arr.ind = TRUE)
  intra_df <- data.table(
    from = paste0(gt$node_names[intra[, 1]], "_t_0"),
    to   = paste0(gt$node_names[intra[, 2]], "_t_0"),
    type = "intra"
  )
  trans <- which(gt$B1 != 0, arr.ind = TRUE)
  trans_df <- data.table(
    from = paste0(gt$node_names[trans[, 1]], "_t_1"),
    to   = paste0(gt$node_names[trans[, 2]], "_t_0"),
    type = "transition"
  )
  rbind(intra_df, trans_df)
}
