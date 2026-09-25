# Bootstrap DBN runner for the smoking-restricted sensitivity analysis
# (R2 Major Comment 4). Same method as 14_bootstrap_parallel_runner.R, run on
# the non-smoker-only subsample instead of the full sample.
#
# B is set lower than the primary analysis (50 vs 100) since this is a
# secondary/sensitivity check, not the primary result -- this is a
# documented, deliberate choice to keep runtime reasonable, not a shortcut
# hiding anything: state it as such in the response letter.

data.table::setDTthreads(1)  # required to avoid fork+thread deadlock with mclapply
                              # -- you must ALSO export OPENBLAS_NUM_THREADS=1 and
                              # OMP_NUM_THREADS=1 in the shell before running this,
                              # exactly as for 14_bootstrap_parallel_runner.R.

mg <- readRDS("analysis/prepped_data_nonsmoker.rds")

GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 50
SCORE <- "bic-g"
MAXP <- 8          # match whatever you used for the primary full-resolution run
N_CORES <- 7

dir.create("boot_results_nonsmoker", showWarnings = FALSE)

run_one_boot <- function(grp_name, b, mg) {
  outfile <- sprintf("boot_results_nonsmoker/%s_b%03d_%s_maxp%d.rds", grp_name, b, SCORE, MAXP)
  if (file.exists(outfile)) return(invisible(NULL))

  set.seed(30000 * which(GROUPS == grp_name) + b)
  grp <- mg[conditiongroup == grp_name]
  parts <- unique(grp$Participant_id)
  samp <- sample(parts, length(parts), replace = TRUE)

  pieces <- vector("list", length(samp))
  for (i in seq_along(samp)) {
    sub <- grp[Participant_id == samp[i]]
    sub[, inst_id := paste0(samp[i], "_", i)]
    pieces[[i]] <- sub
  }
  bdt <- rbindlist(pieces)

  net <- tryCatch(
    fit_group_dbn_hc(bdt, size = 2, score = SCORE, id_col = "inst_id", maxp = MAXP),
    error = function(e) { message(sprintf("FAILED %s b=%d: %s", grp_name, b, e$message)); NULL }
  )
  if (is.null(net)) return(invisible(NULL))
  saveRDS(extract_arcs(net), outfile)
  invisible(NULL)
}

jobs <- expand.grid(grp_name = GROUPS, b = 1:B, stringsAsFactors = FALSE)
jobs <- jobs[!file.exists(sprintf("boot_results_nonsmoker/%s_b%03d_%s_maxp%d.rds",
                                    jobs$grp_name, jobs$b, SCORE, MAXP)), ]
cat("Jobs remaining:", nrow(jobs), "of", 4 * B, "\n")

if (nrow(jobs) > 0) {
  t0 <- Sys.time()
  invisible(mclapply(seq_len(nrow(jobs)), function(i) {
    run_one_boot(jobs$grp_name[i], jobs$b[i], mg)
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  cat("Batch wall time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")
}

n_done <- length(list.files("boot_results_nonsmoker", pattern = sprintf("_%s_maxp%d\\.rds$", SCORE, MAXP)))
cat("Total fits completed so far:", n_done, "of", 4 * B, "\n")
if (n_done < 4 * B) cat("Some jobs remain (or failed) -- just re-run this script; completed fits are skipped.\n")
