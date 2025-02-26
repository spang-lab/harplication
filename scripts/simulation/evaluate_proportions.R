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

music_metrics <- tibble()
celltypewise_correlation <- tibble()

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
    main_plots <- c("true", "harp_true", "music", "cibersort_lm22", "bp_subtypes")
    proportions_all_main <- proportions_run %>% filter(algo %in% main_plots)
    data_name <- "main_benchmark"
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
    ggsave(file.path(output_dir, "simulation_main.pdf"),
        combined_plot,
        device = cairo_pdf, width = 460, height = 700, units = "mm"
    )
}
saveRDS(proportions_all, file.path(downstream_output_dir, "proportions.rds"))