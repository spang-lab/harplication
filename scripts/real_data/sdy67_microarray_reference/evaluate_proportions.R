# the script is to generate plots for cell proportions quality scores for the zimmermann (rna-seq) data with lm22 as reference

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


output_dir <- "data/generated/real_data/sdy67_microarray_reference"

if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}
result_dir <- "data/generated/real_data/sdy67_rnaseq_reference" # dir for the results of cibersort
proportions <- readRDS(file.path(output_dir, "proportions_all_lm22.rds")) # proportions from methods except cibersortx

if (!file.exists(file.path(result_dir, "proportions_cibersort.rds"))) {
    stop("Error: The CIBERSORTx results file 'proportions_cibersort.rds' not found.")
}

proportions_cibersort <- readRDS(file.path(result_dir, "proportions_cibersort.rds")) # result from cibersortx
proportions_all <- rbind(proportions, proportions_cibersort) # proportions from all algorithms

# plot quality scores
main_plots <- c("true", "harp", "cibersort_lm22", "bp_subtypes")
proportions_all <- proportions_all %>% filter(celltype != "extra") # remove the unidentified cells for cell composition quality scores
proportions_all_main <- proportions_all %>% filter(algo %in% main_plots)

# plot quality scores for main benchmark
plot_list_main <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_main,
    output_dir = output_dir,
    data_name = "zimmermann_lm22_main",
    file_type = "pdf",
    combined_width = 16, combined_height = 6,
    individual_width = 10, individual_height = 8,
    y_min_r = 0.6
)
# plot quality scores for hybrid deconvolution
harp_ref <- c("true", "cibersort_harp_lm22", "bp_subtypes", "bp_harp", "cibersort_lm22")
proportions_all_harp_ref <- proportions_all %>% filter(algo %in% harp_ref)
plot_list_ref_har <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_harp_ref,
    output_dir = output_dir,
    data_name = "zimmermann_lm22_harp_ref",
    file_type = "pdf",
    combined_width = 16, combined_height = 6,
    individual_width = 10, individual_height = 8,
    y_min = 0.03,
    y_min_r = 0.6,
    angle_celltype_names = 0,
    hjust_celltype_names = 0.5
)
