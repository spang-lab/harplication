# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(harp)

# REPRODUCABILITY
set.seed(42)
input_dir <- output_dir
output_dir <- file.path(output_dir, 'reference_varying_lambda')
print(paste("Creating", output_dir))
dir.create(file.path(output_dir))
con <- file(file.path(output_dir, "infer_reference_simulation.log"))
sink(con)
sink(con, type = "message")

for (irun in 1:number_simulation_runs) {
    # load input data
    print(paste0("run ", irun, " / ", number_simulation_runs))
    sc_summarized_experiment <- readRDS(file.path(input_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
    bulk_train_summarized_experiment <- readRDS(file.path(input_dir, paste0("bulk_train_summarized_experiment_run_", irun, ".rds")))
    bulk_test_summarized_experiment <- readRDS(file.path(input_dir, paste0("bulk_test_summarized_experiment_run_", irun, ".rds")))
    bulk_pheno_train_true_summarized_experiment <- readRDS(file.path(input_dir, paste0("bulk_pheno_train_true_summarized_experiment_run_", irun, ".rds")))

    # infer proportions of HARP
    for (lambda in lambda_seq) {
        print(paste("InferringHarp proportions with lambda being", lambda))
        harp_true <- benchmark_harp(
            sc_summarized_experiment = sc_summarized_experiment,
            bulk_train_summarized_experiment = bulk_train_summarized_experiment,
            bulk_test_summarized_experiment = bulk_test_summarized_experiment,
            bulk_pheno_train_summarized_experiment = bulk_pheno_train_true_summarized_experiment,
            lambda_seq = lambda,
            n_folds = 5
        )
        print("Done with Harp")
        saveRDS(harp_true$output_harp, file.path(output_dir, paste0("output_harp_true_lambda_", lambda,"_run_", irun, ".rds")))

        harp_p_true <- harp_true$proportions %>% add_column(run = irun, algo = "harp_true")
        saveRDS(harp_p_true, file.path(output_dir, paste0("proportions_harp_true_lambda_run_", lambda, "_run_", irun, ".rds")))
    }

}


sink()
sink(type = "message")