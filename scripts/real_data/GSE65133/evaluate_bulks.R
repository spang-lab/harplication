# this script is to benchmark recounstructed bulk samples from Harp against other algorithms using microarray data and lm22 as reference

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

GSE65133_data <- readRDS("data/source/processed_GSE65133.rds")
# loading the data

harp <- readRDS(file.path(output_dir, "output_harp_run_1.rds"))

# caclulate correlation of bulks for bayesprim
bp_res <- readRDS(file.path(output_dir, "bp_res_bayesPrism_run_1.rds"))
bp_res_harp <- readRDS(file.path(output_dir, "bp_res_bayesPrism_harp_run_1.rds"))
true_bulk_expression <- GSE65133_data$bulk_counts_test

# load the bulk correlations from cibersort
if (!file.exists(file.path(output_dir, "bulk_correlation_cibersort.rds"))) {
    stop("Error: The CIBERSORTx results file 'bulk_correlation_cibersort.rds' not found.")
}
correlations_bulk_cibersort <- readRDS(file.path(output_dir, "bulk_correlation_cibersort.rds"))
# bayesprism
correlations_bulk_bp <- calculate_bulk_correlation_bp(
    bp_res = bp_res,
    true_bulk_expression = true_bulk_expression,
    exclude_celltypes = NULL,
    algo = "bp_subtypes"
)
# bayesprism with harp reference
correlations_bulk_bp_harp <- calculate_bulk_correlation_bp(
    bp_res = bp_res_harp,
    true_bulk_expression = true_bulk_expression,
    exclude_celltypes = NULL,
    algo = "bp_harp"
)

# calculated bulk correlations for harp
correlations_bulk_harp <- calculate_bulk_correlation_harp(
    harp_output = harp,
    true_bulk_expression = true_bulk_expression
)


# combine all correlation for all algorithms
bulk_correlations_all <- rbind(
    correlations_bulk_harp,
    correlations_bulk_cibersort,
    correlations_bulk_bp,
    correlations_bulk_bp_harp
)

# bulk recoustruction plots for main algorithms
main_plots <- c("harp", "cibersort_lm22", "bp_subtypes")
bulk_correlations_all_main <- bulk_correlations_all %>% filter(algo %in% main_plots)

box_plot_bulk_correlation(
    bulk_correlations_all = bulk_correlations_all_main,
    data_name = "GSE65133_main",
    output_dir = output_dir,
    individual_width = 9,
    individual_height = 8,
    file_type = "pdf"
)
# bulk recounstruction for algorithms with harp refernce
harp_ref <- c("cibersort_harp", "bp_subtypes", "bp_harp", "cibersort_lm22")
bulk_correlations_harp_ref <- bulk_correlations_all %>% filter(algo %in% harp_ref)

box_plot_bulk_correlation(
    bulk_correlations_all = bulk_correlations_harp_ref,
    data_name = "GSE65133_harp_ref",
    output_dir = output_dir,
    individual_width = 9,
    individual_height = 8,
    file_type = "pdf"
)

# box plots of recoustructed bulk samples with harp and real data
predicted_samples <- predict_bulk_expression(harp_output = harp, data = GSE65133_data)
# box plot and historam of the predicted expressions
hist_plot <- plot_box_predicted_bulks(
    predicted_bulk_expression = predicted_samples$predicted_bulk_expression_harp,
    true_bulk_expression = predicted_samples$true_bulk_expression,
    predicted_bulk_from_experiment = predicted_samples$predicted_bulk_expression_experiment,
    file_type = "pdf",
    data_name = "GSE65133",
    output_dir = output_dir,
    reference_name = "LM22 reference",
    height = 8,
    width = 8
)


ref_xsc <- compute_reference_harp(sc_library = GSE65133_data$sc_library) # the anchor
df_bulk <- bulk_expression_correlations_data_frame(
    harp_model = harp,
    true_bulk_expression = predicted_samples$true_bulk_expression,
    bulk_pheno_test = GSE65133_data$bulk_pheno_test_true,
    ref_xsc = ref_xsc
)
