## Master R script launcher for ...

suppressMessages({
  library(rstudioapi)
})
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

# Data preparation
source("scripts_full_resolution/00_prep.R")
# Shared functions
source("scripts_full_resolution/01_helpers.R")
## Build the network ##
# Hardware benchmarking
source("scripts_full_resolution/13_benchmark_full_res.R")
# Primary bootstrap beconstruction & thresholding
source("scripts_full_resolution/14_bootstrap_parallel_runner.R")
source("scripts_full_resolution/06_aggregate_strength_full.R")
## Test and characterise the network ##
# Formal Bayesian inference: Beta-Binomial network-density comparison between groups
source("scripts_full_resolution/07_bayesian_density_full.R")
# Directed/temporal edge analysis - intra-slice vs. transition edge breakdown
source("scripts_full_resolution/08_directed_edges_full.R")
# Conventional-method comparison: DBN network vs. FDR-corrected Pearson correlation network
source("scripts_full_resolution/09_correlation_comparison_full.R")
# Stability/overfitting check: split-half Jaccard overlap
source("scripts_full_resolution/10_stability_check_full.R")
# Formal inference: frequentist (permutations)
source("scripts_full_resolution/11_permutation_runner_full.R")
source("scripts_full_resolution/12_permutation_pvalues_full.R")
## Confound and clinical follow-ups ##
# Smoking sensitivity
source("scripts_full_resolution/15_prep_smoking_subsample.R")
source("scripts_full_resolution/16_bootstrap_parallel_nonsmoker.R")
source("scripts_full_resolution/17_aggregate_bayesian_nonsmoker.R")
# Clinical correlation question
# from scratch, one would need to run 18_clinical_prep.py, but needed material in the clinical data folder
source("scripts_full_resolution/19_clinical_correlation.R")
## Validate the network inference method ##
# "Realistic" ground-truth simulation
source("scripts_full_resolution/20_simulate_ground_truth.R")
# Tier 1: one ground-truth network, no group difference, report precision/recall/F1 + direction accuracy
source("scripts_full_resolution/21_tier1_recovery.R")
# Tier 2: adds to Tier 1 the two-condition design with/without an embedded difference, replicated enough times for null-calibration and power estimates
source("scripts_full_resolution/22_tier2_runner.R")
source("scripts_full_resolution/23_tier2_aggregate.R")
