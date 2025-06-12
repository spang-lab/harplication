# Define the number of runs to process
devtools::document()
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(MuSiC)
library(gridExtra)
library(ggplot2)
library(BayesPrism)
library(DTD)
train_size <- 150 # update this for the number of train samples # 20
test_size <- 100 # update this for the number of test samples #100
num_runs <- 5 # Update this to match the number of runs you have #5

# Define output directory and data name
output_dir <- file.path(
    system.file("data/generated", package = "harplication"),
    "real_data", paste0("sdy67_rnaseq_reference_", train_size, "train_", test_size, "test_", num_runs, "run")
)
if (!dir.exists(output_dir)) {
    stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}

# Create a directory for statistical analysis
stats_dir <- file.path(output_dir, "statistical_analysis")
dir.create(stats_dir)

##############################################################################################
#
# proportions_harp_bp contains data from all runs for bayesPrism and harp
#
##############################################################################################
proportions_harp_bp <- NULL

for (irun in 1:num_runs) {
    # Read the proportions file for this run
    proportion_file <- file.path(output_dir, paste0("proportions_all_run_", irun, ".rds"))

    if (file.exists(proportion_file)) {
        cat(paste0("Processing run ", irun, "\n"))

        # Read the proportions
        proportions_run <- readRDS(proportion_file)

        # Add a column to identify the run number (optional)
        proportions_run$run <- irun

        # Combine with the main data frame
        if (is.null(proportions_harp_bp)) {
            proportions_harp_bp <- proportions_run
        } else {
            proportions_harp_bp <- rbind(proportions_harp_bp, proportions_run)
        }
    }
}


if (!file.exists(file.path(output_dir, "proportions_cibersort_sdy67_rna.rds"))) {
    stop("Error: The CIBERSORTx results file 'proportions_cibersort_sdy67_rna.rds' not found.")
}
proportions_cibersort <- readRDS(file.path(output_dir, "proportions_cibersort_sdy67_rna.rds"))

proportions_all <- rbind(proportions_harp_bp, proportions_cibersort)

proportions_all <- proportions_all %>% filter(celltype != "extra")

saveRDS(proportions_all, file.path(stats_dir, "proportions_all_combined.rds"))

######################################################################################################
#
# COMPUTING METRICS
#
######################################################################################################

music_metrics <- tibble()
celltypewise_correlation <- tibble()
sample_correlation <- tibble()

# For the plots we generate lists in order to plot all runs
table_music_metrics <- list()
proportions_heatmap <- list()
table_celltypewise_correlation <- list()
correlations_heatmap <- list()

for (irun in 1:num_runs) {
    print("Computing quality metrics for comparison")
    proportions_run <- proportions_all %>% filter(run == irun)
    data_name <- "main_benchmark"
    metric_list <- get_all_individual_eval_statistics_bars(proportions = proportions_run)
    music_metrics <- rbind(music_metrics, metric_list$music_metrics %>% cbind(run = irun))
    celltypewise_correlation <- rbind(celltypewise_correlation, metric_list$correlations %>% cbind(run = irun))
    sample_correlation <- rbind(sample_correlation, metric_list$sample_correlations %>% cbind(run = irun))
}
#-0.02 for harp_bp Bcell for 150 training
######################################################################################################
#
# SUMMARIZING STATISTICS ACROSS RUNS
#
######################################################################################################
# Summarize statistics for the R,RMSD and AAD plot
eval_statistics_summary <- music_metrics %>%
    group_by(algo) %>%
    mutate(
        mean_R = mean(R),
        mean_RMSD = mean(RMSD),
        mean_mAD = mean(mAD),
        sd_R = sd(R),
        sd_RMSD = sd(RMSD),
        sd_mAD = sd(mAD)
    ) %>%
    ungroup() %>%
    dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct() %>%
    arrange(desc(mean_R))

# celltype wise correlation
eval_statistics_celltype_summary <- celltypewise_correlation %>%
    group_by(algo, celltype) %>%
    mutate(
        mean_R = mean(correlation),
        sd_R = sd(correlation)
    ) %>%
    ungroup() %>%
    dplyr::select(algo, celltype, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct()

# sample wise correlation
eval_statistics_sample_summary <- sample_correlation %>%
    group_by(algo, run) %>%
    mutate(
        mean_R = mean(correlation),
    ) %>%
    ungroup() %>%
    dplyr::select(algo, run, starts_with("mean_")) %>%
    dplyr::distinct() %>%
    group_by(algo) %>%
    mutate(sd_R = sd(mean_R)) %>%
    mutate(mean_R = mean(mean_R)) %>%
    dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct()

print(eval_statistics_celltype_summary)
print(eval_statistics_summary)
print(eval_statistics_sample_summary)

# save the overall quality score summaries
saveRDS(eval_statistics_summary,
    file = file.path(
        stats_dir,
        paste0(num_runs, "_runs_", train_size, "_training_", test_size, "_test_", "eval_statistics_summary.rds")
    )
)
# save cell type correaltions summary
saveRDS(eval_statistics_celltype_summary,
    file = file.path(
        stats_dir,
        paste0(num_runs, "_runs_", train_size, "_training_", test_size, "_test_", "eval_statistics_celltype_summary.rds")
    )
)
# save sample wise correlation summary
saveRDS(eval_statistics_sample_summary,
    file = file.path(
        stats_dir,
        paste0(num_runs, "_runs_", train_size, "_training_", test_size, "_test_", "sample_correlation_runs_summary.rds")
    )
)

###################################################################################################
#
# Z-TEST FOR CELL TYPE CORRELATION AND SAMPLE WISE CORRELATION
#
###################################################################################################
irun <- 5
exclude <- c("true", "harp")
algos <- proportions_all %>%
    filter(!(algo %in% exclude)) %>%
    pull(algo) %>%
    unique()

for (algo in algos) {
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("harp", algo), irun = irun))
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("bp_harp", algo)), irun = irun)
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("cibersort_harp", algo)), irun = irun)
}
#################################################################################################
#
# CALCULATE BULK GENE EXPRESSION CORRELATIONS FOR BENCHMARKING
#
#################################################################################################


bulk_correlations_all <- tibble()

# load the bulk correlations from cibersort
if (!file.exists(file.path(output_dir, "bulk_correlation_cibersort_sdy67_rna.rds"))) {
    stop("Error: The CIBERSORTx results file 'bulk_correlation_cibersort.rds' not found.")
}
correlations_bulk_cibersort <- readRDS(file.path(output_dir, "bulk_correlation_cibersort_sdy67_rna.rds"))
bulk_correlations_all <- rbind(
    bulk_correlations_all,
    correlations_bulk_cibersort
)
for (irun in 1:num_runs) {
    harp <- readRDS(file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
    # caclulate correlation of bulks for bayesprism

    bp_res <- readRDS(file.path(output_dir, paste0("bp_res_bayesPrism_run_", irun, ".rds")))
    bp_res_harp <- readRDS(file.path(output_dir, paste0("bp_res_bayesPrism_harp_run_", irun, ".rds")))
    true_bulk_expression <- readRDS(file.path(output_dir, paste0("bulk_counts_test_run_", irun, ".rds")))
    # bayesprism
    correlations_bulk_bp <- calculate_bulk_correlation_bp(
        bp_res = bp_res,
        true_bulk_expression = true_bulk_expression,
        exclude_celltypes = NULL,
        algo = "bp_subtypes",
        irun = irun
    )

    # bayesprism with harp reference
    correlations_bulk_bp_harp <- calculate_bulk_correlation_bp(
        bp_res = bp_res_harp,
        true_bulk_expression = true_bulk_expression,
        exclude_celltypes = NULL,
        algo = "bp_harp",
        irun = irun
    )

    # calculate bulk correlations for harp
    correlations_bulk_harp <- calculate_bulk_correlation_harp(
        harp_output = harp,
        true_bulk_expression = true_bulk_expression,
        irun = irun
    )

    # combine all correlations
    bulk_correlations_all <- rbind(
        bulk_correlations_all,
        correlations_bulk_harp,
        correlations_bulk_bp,
        correlations_bulk_bp_harp
    )
}

saveRDS(bulk_correlations_all, file = file.path(stats_dir, "bulk_correlations_all_combined.rds"))


eval_statistics_bulk_experssion_summary <- bulk_correlations_all %>%
    group_by(algo, run) %>%
    mutate(
        mean_R = mean(correlation),
    ) %>%
    ungroup() %>%
    dplyr::select(algo, run, starts_with("mean_")) %>%
    dplyr::distinct() %>%
    group_by(algo) %>%
    mutate(sd_R = sd(mean_R)) %>%
    mutate(mean_R = mean(mean_R)) %>%
    dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct()

print(eval_statistics_bulk_experssion_summary)


saveRDS(eval_statistics_bulk_experssion_summary,
    file = file.path(
        stats_dir,
        paste0(num_runs, "_runs_", train_size, "_training_", test_size, "_test_", "correlation_bulk_experssion_summary.rds")
    )
)
###############################################################################################################
#
# GERNERATE PLOTS FOR THE PAPER FOR 1 RUN
#
###############################################################################################################
irun <- 5 # this run is the closest to the mean in this case using 150 training samples and 100 test samples

main_plots <- c("true", "harp", "cibersort_lm22", "bp_subtypes")
proportions_all_main <- proportions_all %>%
    filter(algo %in% main_plots) %>%
    filter(run == irun)

plot_list_main <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_main,
    output_dir = stats_dir,
    data_name = "zimmermann_main",
    file_type = "pdf",
    combined_width = 10, combined_height = 6,
    individual_width = 12, individual_height = 6
)
# plot quality scores for hybrid deconvolution
harp_ref <- c("true", "cibersort_rna", "cibersort_harp", "bp_subtypes", "bp_harp", "cibersort_lm22")
proportions_all_harp_ref <- proportions_all %>%
    filter(algo %in% harp_ref) %>%
    filter(run == irun)
plot_list_harp_ref <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_harp_ref,
    output_dir = stats_dir,
    data_name = "zimmermann_harp_ref",
    file_type = "pdf",
    combined_width = 10, combined_height = 6,
    individual_width = 12, individual_height = 7,
    y_min_r = 0.6
)
############################################################################
#
# CORRELATION PLOTS FOR BULK GENE EXPRESSION PROFILES FROM ALL ALGORITHMS
#
############################################################################
main_plots <- c("harp", "cibersort_lm22", "bp_subtypes")
bulk_correlations_all_main <- bulk_correlations_all %>%
    filter(algo %in% main_plots) %>%
    filter(run == irun)

bulk_correlations_plot_main <- box_plot_bulk_correlation(
    bulk_correlations_all = bulk_correlations_all_main,
    data_name = "zimmermann_main",
    output_dir = stats_dir,
    individual_width = 8,
    individual_height = 8,
    file_type = "pdf"
)

# bulk correlation plota for the algorithms using harp referecne
harp_ref <- c("cibersort_rna", "cibersort_harp", "bp_subtypes", "bp_harp", "cibersort_lm22")
bulk_correlations_harp_ref <- bulk_correlations_all %>%
    filter(algo %in% harp_ref) %>%
    filter(run == irun)

bulk_correlations_plot_harp_ref <- box_plot_bulk_correlation(
    bulk_correlations_all = bulk_correlations_harp_ref,
    data_name = "zimmermann_harp_ref",
    output_dir = stats_dir,
    individual_width = 8,
    individual_height = 8,
    file_type = "pdf"
)
#################################################################################
# CORRELATION PLOTS FOR BULK GENE EXPRESSION PROFILES FROM HARP AND REAL DATA
#################################################################################
irun <- 5
harp <- readRDS(file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
bulk_count_test <- readRDS(file.path(output_dir, paste0("bulk_counts_test_run_", irun, ".rds")))
test.samples <- colnames(bulk_count_test)
zimmermann_data <- readRDS("data/source/sdy67_rnaseq_ref.rds")
bulk_pheno_test <- zimmermann_data$bulk_pheno[, test.samples]

predicted_samples <- predict_bulk_expression(
    harp_output = harp,
    bulk_pheno_test = bulk_pheno_test,
    bulk_count_test = bulk_count_test,
    sc_library = zimmermann_data$sc_library
)

box_plot <- plot_box_predicted_bulks(
    predicted_bulk_expression = predicted_samples$predicted_bulk_expression_harp,
    true_bulk_expression = predicted_samples$true_bulk_expression,
    predicted_bulk_from_experiment = predicted_samples$predicted_bulk_expression_experiment,
    file_type = "pdf",
    data_name = "zimmermann",
    output_dir = stats_dir,
    width = 12,
    heigh = 8
)


# # table in the appendix
ref_xsc <- compute_reference_harp(sc_library = zimmermann_data$sc_library) # the anchor
df_bulk <- bulk_expression_correlations_data_frame(
    harp_model = harp,
    true_bulk_expression = predicted_samples$true_bulk_expression,
    bulk_pheno_test = bulk_pheno_test,
    ref_xsc = ref_xsc
)
df_bulk

# plot the combined figure for the paper-horizontal format
combined_plot_vertical_reduced <- (plot_list_main$individual_plots$R + plot_list_main$sample_wise) /
    (plot_list_main$celltype_wise + bulk_correlations_plot_main) /
    (box_plot + umap_bulk) /
    plot_layout(heights = c(1, 1, 1, 1, 1, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(3.5),
            color = "black",
            face = "bold"
        ),
        plot.tag.position = c(0.01, 0.95)
    )

ggsave(file.path(stats_dir, "zimmerman_monaco_vertical_reduced.pdf"),
    combined_plot_vertical_reduced,
    device = cairo_pdf,
    width = 700,
    height = 1200,
    units = "mm"
)
# plot the combined figure for the paper-horizontal format
combined_plot_horizontal <- (plot_list_main$combined_plot) /
    (plot_list_main$celltype_wise + plot_list_main$sample_wise) /
    (bulk_correlations_plot_main + box_plot) /
    (plot_list_harp_ref$combined_plot) /
    (plot_list_harp_ref$celltype_wise + plot_list_harp_ref$sample_wise) /
    (bulk_correlations_plot_harp_ref + umap_bulk) +
    plot_layout(heights = c(1, 1, 1, 1, 1, 1), guides = "collect", ncol = 2) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(3.5),
            color = "black",
            face = "bold"
        ),
        plot.tag.position = c(0.01, 0.95)
    )

ggsave(file.path(stats_dir, "zimmerman_monaco_horizontal.pdf"),
    combined_plot_horizontal,
    device = cairo_pdf,
    width = 1200,
    height = 1200,
    units = "mm"
)
