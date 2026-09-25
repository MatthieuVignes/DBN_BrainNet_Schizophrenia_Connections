library(parallel)
data.table::setDTthreads(1)
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) install.packages("RhpcBLASctl")
RhpcBLASctl::blas_set_num_threads(1)

# Tier 2: are the formal statistical tests well-calibrated, and do they have
# power to detect an effect of the size actually claimed on the real data?
#
# Two scenarios, each repeated N_REPLICATES times with fresh random data:
#   "null"   -- group A and group B are simulated from the IDENTICAL ground
#               truth (only participant-level noise/jitter differs). Any
#               replicate where the Bayesian test finds a credible density
#               difference here is a FALSE POSITIVE. Across replicates, the
#               false-positive rate should sit near the nominal ~5% for a
#               well-calibrated test.
#   "effect" -- group B's ground truth has additional true edges layered on
#               top of group A's, producing a KNOWN, designed density
#               difference. The fraction of replicates where the test
#               correctly detects a credible difference IN THE CORRECT
#               DIRECTION is the method's empirical power at this effect size
#               and this sample size (45 participants/group, matching the
#               real data).
#
# B_BOOT and N_REPLICATES are both reduced from the primary/Tier-1 analyses,
# since the total cost multiplies across replicates x scenarios x groups x
# bootstraps. This is a deliberate, documented trade-off for a large-scale
# validation exercise -- see the printed cost estimate below before running.
# Increase N_REPLICATES for tighter calibration/power estimates if you have
# the time to let this run longer (each replicate is independent, so this
# scales linearly, not combinatorially).

N_NODES <- 26
N_MODULES <- 4
N_PARTICIPANTS <- 45
N_TIMEPOINTS <- 171
B_BOOT <- 25            # bootstraps per group per replicate (reduced -- see note above)
N_REPLICATES <- 30      # replicates per scenario (reduced -- see note above)
SCORE <- "bic-g"
MAXP <- 8
N_CORES <- 7
N_EXTRA_EFFECT_EDGES <- 15  # how many extra true edges group B gets in the "effect" scenario

est_total_fits <- 2 * N_REPLICATES * 2 * B_BOOT
cat(sprintf("Planned: %d scenarios x %d replicates x 2 groups x %d bootstraps = %d total fits.\n",
            2, N_REPLICATES, B_BOOT, est_total_fits))
cat("At this node count (26), your earlier macro-ROI pilot fits ran in ~6s each --\n")
cat(sprintf("at that rate this is roughly %.0f-%.0f minutes wall time across %d cores,\n",
            est_total_fits * 5 / N_CORES / 60, est_total_fits * 8 / N_CORES / 60, N_CORES))
cat("though actual per-fit time here may differ slightly (different maxp, different data).\n")
cat("Safe to interrupt after the first few dozen jobs and check the pace before committing\n")
cat("to the full run -- this script is resumable, nothing is lost by stopping early.\n\n")

node_names <- paste0("Region_", sprintf("%02d", 1:N_NODES))
ar_params <- estimate_ar1_params("analysis/prepped_data.rds")
cat("AR(1) params:", ar_params$source, "- phi =", round(ar_params$phi, 3), "\n\n")

dir.create("tier2_boot", showWarnings = FALSE)

# --- build the two scenarios' ground truths, once, deterministically ---
gt_null_A <- build_ground_truth(node_names, n_modules = N_MODULES, seed = 100)
gt_null_B <- gt_null_A  # identical structure -> no true difference

gt_effect_A <- build_ground_truth(node_names, n_modules = N_MODULES, seed = 200)
gt_effect_B <- gt_effect_A
# add N_EXTRA_EFFECT_EDGES true intra-slice edges to group B only, at
# currently-zero positions, giving a known, designed density difference
set.seed(201)
zero_positions <- which(gt_effect_B$B0 == 0 & row(gt_effect_B$B0) != col(gt_effect_B$B0), arr.ind = TRUE)
extra_idx <- zero_positions[sample(nrow(zero_positions), N_EXTRA_EFFECT_EDGES), , drop = FALSE]
for (i in seq_len(nrow(extra_idx))) {
  gt_effect_B$B0[extra_idx[i, 1], extra_idx[i, 2]] <- sample(c(-1, 1), 1) * runif(1, 0.15, 0.40)
}
gt_effect_B$n_intra_edges <- sum(gt_effect_B$B0 != 0)

cat(sprintf("Null scenario: group A and B both have %d intra-slice edges (identical truth)\n",
            gt_null_A$n_intra_edges))
cat(sprintf("Effect scenario: group A has %d intra-slice edges, group B has %d (%d extra, +%.1f%% of pair-space)\n\n",
            gt_effect_A$n_intra_edges, gt_effect_B$n_intra_edges, N_EXTRA_EFFECT_EDGES,
            100 * N_EXTRA_EFFECT_EDGES / choose(N_NODES, 2)))

saveRDS(list(gt_null_A = gt_null_A, gt_null_B = gt_null_B,
             gt_effect_A = gt_effect_A, gt_effect_B = gt_effect_B),
        "results/tier2_ground_truths.rds")

# --- simulate all replicate datasets up front (cheap, deterministic given seed) ---
simulate_replicate_data <- function(scenario, replicate) {
  seed_base <- 10000 * (if (scenario == "null") 1 else 2) + replicate * 17
  gtA <- if (scenario == "null") gt_null_A else gt_effect_A
  gtB <- if (scenario == "null") gt_null_B else gt_effect_B
  dA <- simulate_group(gtA, N_PARTICIPANTS, N_TIMEPOINTS, group_label = paste0(scenario, replicate, "A"),
                        ar_params = ar_params, seed = seed_base + 1)
  dB <- simulate_group(gtB, N_PARTICIPANTS, N_TIMEPOINTS, group_label = paste0(scenario, replicate, "B"),
                        ar_params = ar_params, seed = seed_base + 2)
  list(A = dA, B = dB)
}

# --- pre-simulate and cache ALL replicate datasets to disk, SERIALLY, before
# any parallel fitting starts. This is deliberate: doing this inside the
# parallel bootstrap loop would let multiple worker processes race to
# generate/write the same cache file simultaneously (harmless in principle
# since the seed is deterministic, but file writes from concurrent processes
# can still corrupt the file on disk) -- simulating once, up front, avoids
# this entirely and is cheap relative to the DBN fits themselves.
dir.create("tier2_data", showWarnings = FALSE)
cat("Pre-simulating all replicate datasets (serial, cheap)...\n")
for (scenario in c("null","effect")) {
  for (replicate in 1:N_REPLICATES) {
    fA <- sprintf("tier2_data/%s_r%03d_A.rds", scenario, replicate)
    fB <- sprintf("tier2_data/%s_r%03d_B.rds", scenario, replicate)
    if (file.exists(fA) && file.exists(fB)) next
    sim <- simulate_replicate_data(scenario, replicate)
    saveRDS(sim$A, fA); saveRDS(sim$B, fB)
  }
}
cat("Done.\n\n")

get_replicate_data <- function(scenario, replicate, group_label) {
  readRDS(sprintf("tier2_data/%s_r%03d_%s.rds", scenario, replicate, group_label))
}

# one bootstrap fit job
run_one_fit <- function(scenario, replicate, group_label, b) {
  outfile <- sprintf("tier2_boot/%s_r%03d_%s_b%03d.rds", scenario, replicate, group_label, b)
  if (file.exists(outfile)) return(invisible(NULL))

  dt <- get_replicate_data(scenario, replicate, group_label)

  set.seed(50000 + replicate * 97 + b + (if (group_label == "B") 500 else 0))
  parts <- unique(dt$Participant_id)
  samp <- sample(parts, length(parts), replace = TRUE)
  pieces <- lapply(seq_along(samp), function(i) {
    s <- dt[Participant_id == samp[i]]; s[, inst_id := paste0(samp[i], "_", i)]; s
  })
  bdt <- rbindlist(pieces)

  net <- tryCatch(fit_group_dbn_hc(bdt, size = 2, score = SCORE, id_col = "inst_id", maxp = MAXP),
                   error = function(e) { message(sprintf("FAILED %s r%d %s b%d: %s",
                                                            scenario, replicate, group_label, b, e$message)); NULL })
  if (is.null(net)) return(invisible(NULL))
  saveRDS(extract_arcs(net), outfile)
  invisible(NULL)
}

jobs <- expand.grid(scenario = c("null","effect"), replicate = 1:N_REPLICATES,
                     group_label = c("A","B"), b = 1:B_BOOT, stringsAsFactors = FALSE)
jobs$outfile <- sprintf("tier2_boot/%s_r%03d_%s_b%03d.rds", jobs$scenario, jobs$replicate, jobs$group_label, jobs$b)
jobs <- jobs[!file.exists(jobs$outfile), ]
cat("Jobs remaining:", nrow(jobs), "of", est_total_fits, "\n")

if (nrow(jobs) > 0) {
  t0 <- Sys.time()
  invisible(mclapply(seq_len(nrow(jobs)), function(i) {
    run_one_fit(jobs$scenario[i], jobs$replicate[i], jobs$group_label[i], jobs$b[i])
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  cat("Batch wall time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")
}

n_done <- length(list.files("tier2_boot", pattern = "\\.rds$"))
cat("Total fits completed so far:", n_done, "of", est_total_fits, "\n")
if (n_done < est_total_fits) cat("Some jobs remain (or failed) -- just re-run this script; completed fits are skipped.\n")
