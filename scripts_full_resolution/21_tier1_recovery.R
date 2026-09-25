# Tier 1: does the method recover a KNOWN network structure?
# Simulates one group's worth of data (default: matching H.E's dimensions --
# 45 participants x 171 scans) from a known ground-truth two-slice DBN,
# then runs the EXACT SAME bootstrap + optimal-threshold pipeline used on
# the real data (fit_group_dbn_hc, bnlearn:::threshold), and reports:
#   - precision/recall/F1 for intra-slice edges, evaluated on the UNDIRECTED
#     skeleton (see note below on why)
#   - precision/recall/F1 for transition edges, evaluated on the DIRECTED
#     (ordered) pair, since time itself provides identification here
#   - a full precision-recall curve across all candidate thresholds, so the
#     optimal-threshold procedure's chosen point can be shown against the
#     best achievable trade-off
#   - specificity check: are any edges spuriously found to/from Scan_no,
#     which has NO true edges by construction?
#
# NOTE on intra-slice direction: purely contemporaneous linear-Gaussian
# structural relationships are not always identifiable from observational
# data alone (Markov equivalence -- see e.g. Spirtes, Glymour & Scheines,
# 2000) -- a correctly-specified method can recover the right UNDIRECTED
# relationship while still sometimes getting the direction "wrong" simply
# because both directions are statistically indistinguishable given the
# data. We therefore score intra-slice recovery on the undirected skeleton,
# which is the fair and standard way to evaluate this, and reserve directed
# scoring for transition edges, where the temporal order breaks the
# symmetry and both directions genuinely are distinguishable in principle.

N_NODES <- 26          # matches the macro-ROI pilot's node count for direct comparability
N_MODULES <- 4
N_PARTICIPANTS <- 45
N_TIMEPOINTS <- 171    # matches H.E's per-participant scan count
B <- 50                # bootstrap count -- reduced from the primary analysis's 100,
                        # documented choice for a validation exercise, not a shortcut
SCORE <- "bic-g"
MAXP <- 8               # match whatever you used for the primary full-resolution run
SEED <- 1

node_names <- paste0("Region_", sprintf("%02d", 1:N_NODES))
gt <- build_ground_truth(node_names, n_modules = N_MODULES, seed = SEED)
cat(sprintf("Ground truth: %d intra-slice edges (%.1f%% of %d pairs), %d transition edges (%.1f%% of %d ordered pairs)\n",
            gt$n_intra_edges, 100 * gt$n_intra_edges / choose(N_NODES, 2), choose(N_NODES, 2),
            gt$n_trans_edges, 100 * gt$n_trans_edges / (N_NODES * N_NODES), N_NODES * N_NODES))

ar_params <- estimate_ar1_params("analysis/prepped_data.rds")
cat("AR(1) params:", ar_params$source, "- phi =", round(ar_params$phi, 3), "\n\n")

sim_data <- simulate_group(gt, N_PARTICIPANTS, N_TIMEPOINTS, group_label = "tier1",
                            ar_params = ar_params, seed = SEED)

cat("Fitting bootstrap DBNs on simulated data (", B, "resamples )...\n")
dir.create("tier1_boot", showWarnings = FALSE)
parts <- unique(sim_data$Participant_id)
for (b in seq_len(B)) {
  outfile <- sprintf("tier1_boot/b%03d.rds", b)
  if (file.exists(outfile)) next
  set.seed(1000 + b)
  samp <- sample(parts, length(parts), replace = TRUE)
  pieces <- lapply(seq_along(samp), function(i) {
    s <- sim_data[Participant_id == samp[i]]; s[, inst_id := paste0(samp[i], "_", i)]; s
  })
  bdt <- rbindlist(pieces)
  net <- fit_group_dbn_hc(bdt, size = 2, score = SCORE, id_col = "inst_id", maxp = MAXP)
  saveRDS(extract_arcs(net), outfile)
  if (b %% 10 == 0) cat("  ", b, "/", B, "\n")
}

# --- aggregate strength across bootstraps ---
arcs_list <- lapply(list.files("tier1_boot", full.names = TRUE), function(f) {
  a <- readRDS(f); paste0(a$from, "->", a$to)
})
tab <- table(unlist(arcs_list))
strength <- data.frame(from = sub("->.*$", "", names(tab)), to = sub("^.*->", "", names(tab)),
                        strength = as.numeric(tab) / B)
thr <- bnlearn:::threshold(strength)
cat("\nOptimal threshold:", round(thr, 3), "\n")

recovered <- strength[strength$strength >= thr, ]

# --- build recovered edge sets in comparable form to ground_truth_arcs() ---
rec_intra <- recovered[grepl("_t_0$", recovered$from) & grepl("_t_0$", recovered$to), ]
rec_trans <- recovered[grepl("_t_1$", recovered$from) & grepl("_t_0$", recovered$to), ]

gt_arcs <- ground_truth_arcs(gt)
gt_intra <- gt_arcs[type == "intra"]
gt_trans <- gt_arcs[type == "transition"]

# --- intra-slice: undirected skeleton comparison ---
undirected_pair <- function(from, to) apply(cbind(from, to), 1, function(x) paste(sort(x), collapse = "|"))
true_intra_pairs <- unique(undirected_pair(gt_intra$from, gt_intra$to))
rec_intra_pairs <- unique(undirected_pair(rec_intra$from, rec_intra$to))
all_possible_intra_pairs <- combn(paste0(node_names, "_t_0"), 2, function(x) paste(sort(x), collapse = "|"))

tp_intra <- length(intersect(true_intra_pairs, rec_intra_pairs))
fp_intra <- length(setdiff(rec_intra_pairs, true_intra_pairs))
fn_intra <- length(setdiff(true_intra_pairs, rec_intra_pairs))
prec_intra <- tp_intra / (tp_intra + fp_intra)
rec_intra_rate <- tp_intra / (tp_intra + fn_intra)
f1_intra <- 2 * prec_intra * rec_intra_rate / (prec_intra + rec_intra_rate)

cat(sprintf("\n--- Intra-slice edges (undirected skeleton) ---\nTP=%d FP=%d FN=%d | Precision=%.3f Recall=%.3f F1=%.3f\n",
            tp_intra, fp_intra, fn_intra, prec_intra, rec_intra_rate, f1_intra))

# --- transition: directed (ordered) pair comparison ---
directed_pair <- function(from, to) paste(from, to, sep = "|")
true_trans_pairs <- unique(directed_pair(gt_trans$from, gt_trans$to))
rec_trans_pairs <- unique(directed_pair(rec_trans$from, rec_trans$to))

tp_trans <- length(intersect(true_trans_pairs, rec_trans_pairs))
fp_trans <- length(setdiff(rec_trans_pairs, true_trans_pairs))
fn_trans <- length(setdiff(true_trans_pairs, rec_trans_pairs))
prec_trans <- tp_trans / (tp_trans + fp_trans)
rec_trans_rate <- tp_trans / (tp_trans + fn_trans)
f1_trans <- 2 * prec_trans * rec_trans_rate / (prec_trans + rec_trans_rate)

cat(sprintf("\n--- Transition edges (directed pair) ---\nTP=%d FP=%d FN=%d | Precision=%.3f Recall=%.3f F1=%.3f\n",
            tp_trans, fp_trans, fn_trans, prec_trans, rec_trans_rate, f1_trans))

# --- direction accuracy, conditional on skeleton being detected in some direction ---
true_trans_skeleton <- unique(undirected_pair(gt_trans$from, gt_trans$to))
detected_skeleton <- rec_trans[undirected_pair(rec_trans$from, rec_trans$to) %in% true_trans_skeleton, ]
if (nrow(detected_skeleton) > 0) {
  correct_dir <- sum(directed_pair(detected_skeleton$from, detected_skeleton$to) %in% true_trans_pairs)
  cat(sprintf("\nDirection accuracy (given the pair was detected as connected): %d/%d = %.3f\n",
              correct_dir, nrow(detected_skeleton), correct_dir / nrow(detected_skeleton)))
}

# --- specificity check: any spurious edges involving Scan_no? ---
scan_no_edges <- recovered[grepl("^Scan_no", recovered$from) | grepl("^Scan_no", recovered$to), ]
cat(sprintf("\nSpurious edges involving Scan_no (should have ~0 true edges): %d found\n", nrow(scan_no_edges)))

# --- full precision-recall curve across thresholds (intra-slice) ---
thresholds_to_try <- sort(unique(strength$strength))
pr_curve <- data.frame(threshold = numeric(0), precision = numeric(0), recall = numeric(0))
for (th in thresholds_to_try) {
  sel <- strength[strength$strength >= th, ]
  sel_intra <- sel[grepl("_t_0$", sel$from) & grepl("_t_0$", sel$to), ]
  sel_pairs <- unique(undirected_pair(sel_intra$from, sel_intra$to))
  tp <- length(intersect(true_intra_pairs, sel_pairs))
  fp <- length(setdiff(sel_pairs, true_intra_pairs))
  fn <- length(setdiff(true_intra_pairs, sel_pairs))
  p <- if (tp + fp == 0) NA else tp / (tp + fp)
  r <- if (tp + fn == 0) NA else tp / (tp + fn)
  pr_curve <- rbind(pr_curve, data.frame(threshold = th, precision = p, recall = r))
}
write.csv(pr_curve, "results/tier1_pr_curve.csv", row.names = FALSE)
cat("\nFull precision-recall curve saved to tier1_pr_curve.csv (", nrow(pr_curve), "points )\n")

saveRDS(list(gt = gt, strength = strength, threshold = thr), "results/tier1_result.rds")
cat("\nTier 1 complete. Summary saved to results/tier1_result.rds\n")
