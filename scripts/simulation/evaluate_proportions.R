# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(gridExtra)
library(tidyverse)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)
library(patchwork)

######################################################################################################
#
# FETCHING PROPORTIONS
#
######################################################################################################

proportions_all <- readRDS(file.path(output_dir, "proportions_all.rds"))
# subset tibble to contain only test ids (becuase true contains all)
test_ids <- proportions_all %>%
    group_by(algo) %>%
    summarise(bulk_ids = list(unique(bulk_id))) %>%
    pull(bulk_ids) %>%
    Reduce(intersect, .)
proportions_all <- proportions_all %>% filter(bulk_id %in% test_ids)

downstream_output_dir <- file.path(output_dir, "downstream")
dir.create(downstream_output_dir)

if (!is.null(cibersort_algorithms)) {
    proportions_all <- add_cibersort_proportions(
        output_dir,
        proportions_all,
        cibersort_filename,
        cibersort_algorithms,
        number_simulation_runs,
        simulation_name
    )
}
print(proportions_all %>% pull(algo) %>% unique())
saveRDS(proportions_all, file.path(downstream_output_dir, "proportions.rds"))

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

for (irun in 1:number_simulation_runs) {
    print("Computing quality metrics for comparison")
    proportions_run <- proportions_all %>% filter(run == irun)
    # Now evaluate
    # MAIN benchmark
    data_name <- "main_benchmark"
    metric_list <- get_all_individual_eval_statistics_bars(proportions = proportions_run)
    music_metrics <- rbind(music_metrics, metric_list$music_metrics %>% cbind(run = irun))
    celltypewise_correlation <- rbind(celltypewise_correlation, metric_list$correlations %>% cbind(run = irun))
    sample_correlation <- rbind(sample_correlation, metric_list$sample_correlations %>% cbind(run = irun))
}

######################################################################################################
#
# SUMMARIZING STATISTICS ACROSS RUNS
#
######################################################################################################
# Summarize statistics for the R,RMSD and AAD plot
eval_statistics_summary <- music_metrics %>%
    group_by(algo) %>%
    mutate(
        mean_RMSD = mean(RMSD),
        mean_mAD = mean(mAD),
        mean_R = mean(R),
        sd_RMSD = sd(RMSD),
        sd_mAD = sd(mAD),
        sd_R = sd(R)
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

######################################################################################################
#
# SIGNIFICANCE TESTING
#
######################################################################################################
exclude <- c("true", "harp_true", "harp_cdeath", "harp_noise")
algos <- proportions_all %>%
    filter(!(algo %in% exclude)) %>%
    pull(algo) %>%
    unique()
for (algo in algos) {
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("harp_true", algo)))
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("harp_noise", algo)))
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("harp_cdeath", algo)))
    print(one_sided_z_test(proportions_all = proportions_all, algos = c("cibersort_harp", algo)))
}

######################################################################################################
#
# PLOTTING (only for first run)
#
######################################################################################################
irun <- 1
proportions_run <- proportions_all %>% filter(run == irun)

lm_results <- lm_celltypes(proportions = proportions_run)
saveRDS(lm_results, file.path(output_dir, paste0("lm_results_main_", irun, ".rds")))
lm_algos <- c("harp_true", "bp_subtypes", "music", "cibersort_lm22")
all_alg_list <- list()
for (alg in lm_algos) {
    plot_list <- list()
    alg_results <- lm_results[[paste(alg)]]
    celltypes <- sort(names(alg_results))
    for (ct in celltypes) {
        plot_list[[paste(ct)]] <- alg_results[[paste(ct)]]$plot
    }
    combined_plot <- wrap_plots(plot_list, ncol = 4)
    all_alg_list[[paste(alg)]] <- combined_plot
    ggsave(file.path(output_dir, paste0("lm_", alg, ".pdf")),
        combined_plot,
        device = cairo_pdf, width = 400, height = 350, units = "mm"
    )
}
all_alg <- wrap_plots(all_alg_list, ncol = 1)
ggsave(file.path(output_dir, paste0("lm_main.pdf")),
    all_alg,
    device = cairo_pdf, width = 400, height = 300, units = "mm"
)
print(names(all_alg_list))

main_plots <- c("true", "harp_true", "music", "cibersort_lm22", "bp_subtypes")
proportions_all_main <- proportions_run %>% filter(algo %in% main_plots)
plot_list_main <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_main,
    output_dir = output_dir,
    data_name = data_name,
    file_type = "pdf",
    combined_width = 12, combined_height = 6,
    individual_width = 6, individual_height = 6,
    y_min = 0,
    y_min_ct = 0.3,
    y_min_r = 0.8,
    simulation = TRUE
)

# Cibersort BayesPrism with Harp
cb_plots <- c("true", "cibersort_lm22", "cibersort_sc", "cibersort_harp", "bp_subtypes", "bp_harp")
proportions_all_cb <- proportions_run %>% filter(algo %in% cb_plots)
data_name <- "cibersort_bayesPrism_benchmark"
plot_list_ref <- plot_all_individual_eval_statistics_bars(
    proportions = proportions_all_cb,
    output_dir = output_dir,
    data_name = data_name,
    file_type = "pdf",
    combined_width = 12, combined_height = 6,
    individual_width = 6, individual_height = 6,
    y_min = 0,
    y_min_ct = 0.3,
    y_min_r = 0.8,
    simulation = TRUE
)

combined_plot <- (plot_list_main$combined_plot) /
    (plot_list_main$celltype_wise + plot_list_main$sample_wise) /
    (plot_list_ref$combined_plot) /
    (plot_list_ref$celltype_wise + plot_list_ref$sample_wise) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(2.5), # Increase size (adjust number as needed)
            color = "black", # Change color if desired
            face = "bold"
        )
    )
ggsave(file.path(output_dir, "simulation_appendix.pdf"),
    combined_plot,
    device = cairo_pdf, width = 460, height = 700, units = "mm"
)

combined_plot <- ((plot_list_main$individual_plots$R & theme(legend.position = "none")) +
    (plot_list_main$sample_wise & theme(legend.position = "none"))) /
    (plot_list_main$celltype_wise & theme(legend.position = "bottom")) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(2.5), # Increase size (adjust number as needed)
            color = "black", # Change color if desired
            face = "bold"
        )
    )
ggsave(file.path(output_dir, "simulation_main.pdf"),
    combined_plot,
    device = cairo_pdf, width = 300, height = 300, units = "mm"
)
