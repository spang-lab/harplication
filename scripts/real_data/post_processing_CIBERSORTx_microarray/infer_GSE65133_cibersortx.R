# infer proportion from cibersortx
devtools::document()
library(tidyverse)
library(dplyr)
library(harp)

output_dir <- "data/generated/real_data/GSE65133"

if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first to get the harp reference as input for
    cibersortx; after using the inputs in CIBESORTx, postprocess the outputs here.")
}

irun <- 1

cibersort_harp <- "" # add the name of file of cibersortx output when using the reference from harp for example "cibersortx_job1.txt"
cibersort_lm22 <- "" # add the name of file of cibersortx output when using lm22 as rerefence here for example "cibersortx_job2.txt"

cibersort_prp <- tibble()

# this microarray with lm22
cibersort_prp_lm22 <- read_cibersort_proportions_result(
    output_dir = output_dir,
    cibersort_filename = cibersort_lm22,
    algoname = "lm22",
    irun = irun
)


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

saveRDS(cibersort_prp, file.path(output_dir, "proportions_cibersort.rds"))
##############

# cibersort_bulk_correlations
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

cibersort_bulk_correlations <- rbind(cibersort_lm22_bulk_correlations, cibersort_harp_bulk_correlations)

saveRDS(cibersort_bulk_correlations, file.path(output_dir, "bulk_correlation_cibersort.rds"))
