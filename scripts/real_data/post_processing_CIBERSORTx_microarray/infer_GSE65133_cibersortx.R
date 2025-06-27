# infer proportion from cibersortx
devtools::document()
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(patchwork)
library(uwot)
library(harp)

train_size <- 12 # update this for the number of train samples
test_size <- 8 # update this for the number of test samples
num_runs <- 5 # Update this to match the number of runs you have #5

# Define output directory and data name
output_dir <- file.path(
    system.file("data/generated", package = "harplication"),
    "real_data", paste0("GSE65133_", train_size, "train_", test_size, "test_", num_runs, "run")
)
if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first to get the input for cibersortx;
    after using the inputs on CIBESORTx, postprocess the outputs.")
}

cibersort_bulk_correlations <- tibble()
cibersort_prp <- tibble()
for (irun in 1:num_runs) {
    cibersort_lm22 <- paste0("CIBERSORTx_Adjusted_lm22_run_", irun, ".txt") # add the name of file of cibersortx output when using lm22 as rerefence here for example "cibersortx_job2.txt"
    cibersort_harp <- paste0("CIBERSORTx_Adjusted_harp_run_", irun, ".txt") # add the name of file of cibersortx output when using the reference from harp (when using lm22 as anchor) for example "cibersortx_job3.txt"
    test_bulks <- readRDS(file.path(output_dir, paste0("bulk_counts_test_run_", irun, ".rds")))
    test_samples <- colnames(test_bulks)
    # this microarray with lm22
    cibersort_prp_lm22 <- read_cibersort_proportions_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_lm22,
        algoname = "lm22",
        irun = irun
    )
    cibersort_prp_lm22 <- cibersort_prp_lm22 %>% filter(bulk_id %in% test_samples)

    cibersort_prp_mapped_lm22 <- mapping_proportions_data_celltypes(proportions = cibersort_prp_lm22, map_name = "GSE65133_lm22")
    cibersort_prp <- rbind(cibersort_prp, cibersort_prp_mapped_lm22)

    # this microarray with reference from harp
    cibersort_prp_harp <- read_cibersort_proportions_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_harp,
        algoname = "harp",
        irun = irun
    )

    cibersort_prp <- rbind(cibersort_prp, cibersort_prp_harp)


    # cibersort_bulk_correlations
    cibersort_lm22_bulk_correlations <- read_cibersort_bulk_correlations_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_lm22,
        algoname = "lm22",
        irun = irun
    )
    cibersort_lm22_bulk_correlations <- cibersort_lm22_bulk_correlations %>% filter(bulk_id %in% test_samples)

    cibersort_bulk_correlations <- rbind(cibersort_bulk_correlations, cibersort_lm22_bulk_correlations)

    cibersort_harp_bulk_correlations <- read_cibersort_bulk_correlations_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_harp,
        algoname = "harp",
        irun = irun
    )

    cibersort_bulk_correlations <- rbind(cibersort_bulk_correlations, cibersort_harp_bulk_correlations)
}

saveRDS(cibersort_prp, file.path(output_dir, "proportions_cibersort_GSE65133.rds"))
saveRDS(cibersort_bulk_correlations, file.path(output_dir, "bulk_correlation_cibersort_GSE65133.rds"))
