# the script is to generate plots for cell proportions quality scores for the microarray data( GSE65133) using lm22 as reference

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

output_dir <- "data/generated/real_data/GSE65133"

if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}
proportions <- readRDS(file.path(output_dir, "proportions_all.rds")) # proportions from methods except cibersort
test.bulks <- unique(proportions$bulk_id) # filter for test samples

if (!file.exists(file.path(output_dir, "proportions_cibersort.rds"))) {
    stop("Error: The CIBERSORTx results file 'proportions_cibersort.rds' not found.")
}

proportions_cibersort <- readRDS(file.path(output_dir, "proportions_cibersort.rds"))
# result from cibersort

proportions_cibersort <- proportions_cibersort %>% filter(bulk_id %in% test.bulks) # since we deconvoled all samples for CIBERSORT(LM22) but we only want to calculate the scores for test samples
proportions_all <- rbind(proportions, proportions_cibersort) # proportions from all algorithms

# removing gamma delta tcell from the evaluation of cell proportions quality scores, it has negative correlation in all algorithms.
proportions_all <- proportions_all %>% filter(celltype != "gamma_delta_T_cell")
proportions_all <- proportions_all %>% filter(celltype != "extra") # remove the unidentified cells for cell composition quality scores

# quality score plots for  main algoritms benchmark
main_plots <- c("true", "harp", "cibersort_lm22", "bp_subtypes")
proportions_all_main <- proportions_all %>% filter(algo %in% main_plots)
plot_list_main <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_main,
    output_dir = output_dir,
    data_name = "GSE65133_main",
    file_type = "pdf",
    combined_width = 12, combined_height = 8,
    individual_width = 9, individual_height = 8,
    angle_celltype_names = 45, hjust_celltype_names = 1,
    y_min_r = 0.6
)

# filter the algorithms for impact of harp's reference
harp_ref <- c("true", "cibersort_lm22", "cibersort_harp", "bp_subtypes", "bp_harp")
proportions_all_harp_ref <- proportions_all %>% filter(algo %in% harp_ref)
plot_list_harp_ref <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_harp_ref,
    output_dir = output_dir,
    data_name = "GSE65133_harp_ref",
    file_type = "pdf",
    combined_width = 12, combined_height = 8,
    individual_width = 9, individual_height = 8,
    angle_celltype_names = 45, hjust_celltype_names = 1,
    y_min_r = 0.6
)
