# Full-resolution (116-region) bootstrap DBN runner, parallelised across
# cores. Run 13_benchmark_full_res.R FIRST to pick a `maxp` value and
# estimate total runtime before launching this.
#
# Uses parallel::mclapply (fork-based parallelism and works on Linux/Mac,
# NOT on Windows)

library(parallel)

data.table::setDTthreads(1)
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) install.packages("RhpcBLASctl")
RhpcBLASctl::blas_set_num_threads(1)

mg <- readRDS("analysis/prepped_data.rds")

GROUPS <- c("H.E","H.R","S.E","S.R")
B <- 100
SCORE <- "bic-g"
MAXP <- 8          # set this from your benchmark results (13_benchmark_full_res.R)
N_CORES <- 7       # leaves 1 core free on an 8-core machine for the OS/other work

dir.create("boot_results_full", showWarnings = FALSE)

run_one_boot <- function(grp_name, b, mg) {
  outfile <- sprintf("boot_results_full/%s_b%03d_%s_maxp%d.rds", grp_name, b, SCORE, MAXP)
  if (file.exists(outfile)) return(invisible(NULL))

  set.seed(10000 * which(GROUPS == grp_name) + b)
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

# Build the full job list (group, b) pairs not yet done, run in parallel.
jobs <- expand.grid(grp_name = GROUPS, b = 1:B, stringsAsFactors = FALSE)
jobs <- jobs[!file.exists(sprintf("boot_results_full/%s_b%03d_%s_maxp%d.rds",
                                    jobs$grp_name, jobs$b, SCORE, MAXP)), ]
cat("Jobs remaining:", nrow(jobs), "of", 4 * B, "\n")

if (nrow(jobs) > 0) {
  t0 <- Sys.time()
  invisible(mclapply(seq_len(nrow(jobs)), function(i) {
    run_one_boot(jobs$grp_name[i], jobs$b[i], mg)
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  cat("Batch wall time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")
}

n_done <- length(list.files("boot_results_full", pattern = sprintf("_%s_maxp%d\\.rds$", SCORE, MAXP)))
cat("Total fits completed so far:", n_done, "of", 4 * B, "\n")
if (n_done < 4 * B) {
  cat("Some jobs remain (or failed) -- just re-run this script; completed\n")
  cat("fits are skipped automatically, only missing ones are (re)run.\n")
}
