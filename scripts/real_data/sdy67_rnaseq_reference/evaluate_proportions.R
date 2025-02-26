# this script is to generate  cell proportions related quality scores for the rna-seq data using a sorted rna-seq reference
# DEPENDENCIES
devtools::document()
library(tidyverse)
library(MuSiC)
library(ggplot2)
library(gridExtra)
library(stringr)
library(patchwork)
library(uwot)
library(harp)
library(DTD)

output_dir <- "data/generated/real_data/sdy67_rnaseq_reference"

if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}

proportions <- readRDS(file.path(output_dir, "proportions_all.rds")) # proportions from methods except cibersort

# cibersortx precessed results
if (!file.exists(file.path(output_dir, "proportions_cibersort.rds"))) {
    stop("Error: The CIBERSORTx results file 'proportions_cibersort.rds' not found.")
}

proportions_cibersort <- readRDS(file.path(output_dir, "proportions_cibersort.rds")) # processed results from cibersortx

proportions_all <- rbind(proportions, proportions_cibersort) # proportions from all algorithms

# plot quality scores for the main benchmark
main_plots <- c("true", "harp", "cibersort_lm22", "bp_subtypes")
proportions_all <- proportions_all %>% filter(celltype != "extra") # remove the unidentified cells for cell composition quality scores

proportions_all_main <- proportions_all %>% filter(algo %in% main_plots)
plot_list_main <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_main,
    output_dir = output_dir,
    data_name = "zimmermann_main",
    file_type = "pdf",
    combined_width = 10, combined_height = 6,
    individual_width = 12, individual_height = 6
)
# plot quality scores for hybrid deconvolution
harp_ref <- c("true", "cibersort_rna", "cibersort_harp", "bp_subtypes", "bp_harp", "cibersort_lm22")
proportions_all_harp_ref <- proportions_all %>% filter(algo %in% harp_ref)
plot_list_harp_ref <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_harp_ref,
    output_dir = output_dir,
    data_name = "zimmermann_harp_ref",
    file_type = "pdf",
    combined_width = 10, combined_height = 6,
    individual_width = 12, individual_height = 7,
    y_min_r = 0.6
)
