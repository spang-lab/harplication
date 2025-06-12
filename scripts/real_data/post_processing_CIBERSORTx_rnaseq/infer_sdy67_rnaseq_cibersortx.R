# infer proportion from cibersortx for zimmermann data with lm22 as referecne
devtools::document()
library(tidyverse)
library(MuSiC)
library(ggplot2)
library(gridExtra)
library(stringr)
library(patchwork)
library(uwot)
library(harp)

train_size <- 150 # update this for the number of train samples #20
test_size <- 100 # update this for the number of test samples #100
num_runs <- 5 # Update this to match the number of runs you have #5

# Define output directory and data name
output_dir <- file.path(
    system.file("data/generated", package = "harplication"),
    "real_data", paste0("sdy67_rnaseq_reference_", train_size, "train_", test_size, "test_", num_runs, "run")
)
if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first to get the input for cibersortx;
    after using the inputs on CIBESORTx, postprocess the outputs.")
}

cibersort_bulk_correlations <- tibble()
cibersort_prp <- tibble()
for (irun in 1:num_runs) {
    cibersort_rnasig <- paste0("CIBERSORTx_Adjusted_rna_run_", irun, ".txt") # add the name of file of cibersortx output when using the rna-seq signature with cibersortx for example "cibersortx_job1.txt"
    cibersort_lm22 <- paste0("CIBERSORTx_Adjusted_lm22_run_", irun, ".txt") # add the name of file of cibersortx output when using lm22 as rerefence here for example "cibersortx_job2.txt"
    cibersort_harp <- paste0("CIBERSORTx_Adjusted_harp_rna_run_", irun, ".txt") # add the name of file of cibersortx output when using the reference from harp (when using rna-seq as anchor) for example "cibersortx_job3.txt"

    cibersort_prp_rna <- read_cibersort_proportions_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_rnasig,
        algoname = "rna",
        irun = irun
    )

    cibersort_prp_mapped_rna <- mapping_proportions_data_celltypes(proportions = cibersort_prp_rna, map_name = "zimmermann_monaco")

    cibersort_prp <- rbind(cibersort_prp, cibersort_prp_mapped_rna)
    # this zimmermann with lm22 add all genes to the bulks for lm22
    cibersort_prp_lm22 <- read_cibersort_proportions_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_lm22,
        algoname = "lm22",
        irun = irun
    )
    cibersort_prp_mapped_lm22 <- mapping_proportions_data_celltypes(proportions = cibersort_prp_lm22, map_name = "zimmermann_lm22")
    cibersort_prp <- rbind(cibersort_prp, cibersort_prp_mapped_lm22)


    cibersort_prp_harp <- read_cibersort_proportions_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_harp,
        algoname = "harp",
        irun = irun
    )
    cibersort_prp <- rbind(cibersort_prp, cibersort_prp_harp)


    ##########

    # cibersort_bulk_correlations_result
    cibersort_lm22_bulk_correlations <- read_cibersort_bulk_correlations_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_lm22,
        algoname = "lm22",
        irun = irun
    )

    cibersort_harp_bulk_correlations <- read_cibersort_bulk_correlations_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_harp,
        algoname = "harp",
        irun = irun
    )

    cibersort_bulk_correlations_rna <- read_cibersort_bulk_correlations_result(
        output_dir = output_dir,
        cibersort_filename = cibersort_rnasig,
        algoname = "rna",
        irun = irun
    )

    cibersort_bulk_correlations <- rbind(
        cibersort_bulk_correlations,
        cibersort_lm22_bulk_correlations,
        cibersort_bulk_correlations_rna,
        cibersort_harp_bulk_correlations
    )
}


saveRDS(cibersort_prp, file.path(output_dir, "proportions_cibersort_sdy67_rna.rds"))
saveRDS(cibersort_bulk_correlations, file.path(output_dir, "bulk_correlation_cibersort_sdy67_rna.rds"))
