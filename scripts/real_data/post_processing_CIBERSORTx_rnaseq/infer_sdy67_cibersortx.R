# infer proportion from cibersortx for zimmermann data
devtools::document()
library(tidyverse)
library(MuSiC)
library(ggplot2)
library(gridExtra)
library(stringr)
library(patchwork)
library(uwot)
library(harp)

output_dir <- "data/generated/real_data/sdy67_rnaseq_reference"
if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first to get the input for cibersortx;
    after using the inputs on CIBESORTx, postprocess the outputs.")
}

irun <- 1

cibersort_rnasig <- "" # add the name of file of cibersortx output when using the rna-seq signature with cibersortx for example "cibersortx_job1.txt"
cibersort_lm22 <- "" # add the name of file of cibersortx output when using lm22 as rerefence here for example "cibersortx_job2.txt"
cibersort_harp <- "" # add the name of file of cibersortx output when using the reference from harp (when using rna-seq as anchor) for example "cibersortx_job3.txt"
cibersort_harp_lm22 <- "" # add the name of file of cibersortx output when using the reference from harp (when using lm22 as anchor) for example "cibersortx_job3.txt"


cibersort_prp <- tibble()
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


cibersort_prp_harp_lm22 <- read_cibersort_proportions_result(
    output_dir = output_dir,
    cibersort_filename = cibersort_harp_lm22,
    algoname = "harp_lm22",
    irun = irun
)

cibersort_prp <- rbind(cibersort_prp, cibersort_prp_harp_lm22)

saveRDS(cibersort_prp, file.path(output_dir, "proportions_cibersort.rds"))
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

cibersort_bulk_correlations_harp_lm22 <- read_cibersort_bulk_correlations_result(
    output_dir = output_dir,
    cibersort_filename = cibersort_harp_lm22,
    algoname = "harp_lm22",
    irun = irun
)
cibersort_bulk_correlations <- rbind(
    cibersort_lm22_bulk_correlations,
    cibersort_harp_bulk_correlations,
    cibersort_bulk_correlations_rna,
    cibersort_bulk_correlations_harp_lm22
)

saveRDS(cibersort_bulk_correlations, file.path(output_dir, "bulk_correlation_cibersort.rds"))
