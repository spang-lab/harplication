# SIMULATE ARTIFICIAL BULKS from real single cell data
# DEPENDENCIES
library(harplication)
library(tidyverse)
library(SummarizedExperiment)
library(DTD)
# set parameters
config <- config::get(file = system.file("R/inst/config", "config.yml", package = "harplication"), use_parent = FALSE)
if (config$standalone) {
    source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
}
con <- file(file.path(output_dir, "distort_input_simulation.log"))
sink(con)
sink(con, type = "message")

# REPRODUCABILITY
set.seed(1)

percentage_genes_to_distort <- 0.4

proportions_all <- tibble()

for (irun in 1:number_simulation_runs) {
    print(paste0("run ", irun, " / ", number_simulation_runs))

    print(paste0("run ", irun, " / ", number_simulation_runs))
    sc_summarized_experiment <- readRDS(file.path(output_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
    bulk_train_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_train_summarized_experiment_run_", irun, ".rds"))
    )
    bulk_pheno_train_true_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_pheno_train_true_summarized_experiment_run_", irun, ".rds"))
    )
    bulk_test_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_test_summarized_experiment_run_", irun, ".rds"))
    )
    bulk_pheno_train_cdeath_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_pheno_train_cdeath_summarized_experiment_run_", irun, ".rds"))
    )
    proportions_true <- readRDS(
        file.path(output_dir, paste0("proportions_true_run_", irun, ".rds"))
    )

    original_counts <- assay(bulk_train_summarized_experiment, "counts")
    dist_bulks <- distort_bulk_counts(original_counts, percentage_genes_to_distort, distortion_factor_mean, distortion_factor_std)
}

sink()
sink(type = "message")