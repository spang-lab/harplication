# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(gridExtra)
library(tidyverse)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)
library(uwot)


input_dir <- output_dir
output_dir <- file.path(output_dir, "reference_varying_lambda")

# Load data to compose the signature data list for input to plot_umap_signature
for (irun in 1:number_simulation_runs){
        # choose optimal lambda
        output_harp_cross_val <- readRDS(file.path(input_dir, paste0("output_harp_true_run_", irun, ".rds")))
        lambda_optimal <- output_harp_cross_val$lambda$optimal_lambda_first

        # get sc data
        sc_summarized_experiment <- readRDS(file.path(input_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
        sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
        sc_mat <- sc_library %>%
                dplyr::select(-celltype, -cell_id) %>%
                as.matrix() %>%
                t()
        colnames(sc_mat) <- sc_library %>% pull(celltype)
        signature_data <- list(
                Single_cell_Signature = sc_mat
        )
        # first reference
        # get all references (for varying lambdas)
        for (i in seq_along(lambda_seq)) {
                output_harp_lambda <- readRDS(
                        file.path(output_dir, paste0("output_harp_true_lambda_", lambda_seq[i], "_run_", irun, ".rds")))
                harp <- output_harp_lambda$reference_profiles$estimated_reference_first
                signature_data[[paste0("lambda_", i)]] <- harp
        }
        # make a single UMAP plot
        data_name <- paste0("_simulation_run_", irun, "_reference_first_")

        # second reference
        plot_umap_signature(signature_data, data_name, output_dir, lambda_seq, lambda_optimal, min_dist = 0.3, spread = 0.3, n_neighbors = 3)
        # get all references (for varying lambdas)
        for (i in seq_along(lambda_seq)) {
                output_harp_lambda <- readRDS(
                        file.path(output_dir, paste0("output_harp_true_lambda_", lambda_seq[i], "_run_", irun, ".rds"))
                )
                harp <- output_harp_lambda$reference_profiles$estimated_reference_second
                signature_data[[paste0("lambda_", i)]] <- harp
        }
        data_name <- paste0("_simulation_run_", irun, "_reference_second_")
        plot_umap_signature(signature_data, data_name, output_dir, lambda_seq, lambda_optimal, min_dist = 0.3, spread = 0.3, n_neighbors = 3)
}
