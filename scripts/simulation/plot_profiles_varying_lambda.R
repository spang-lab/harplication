# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(harp)
library(patchwork)

# REPRODUCABILITY
set.seed(42)
input_dir <- output_dir
output_dir <- file.path(output_dir, 'reference_varying_lambda')

plot_ct_correlation <- function(correlation_tibble, output_dir, file_name, lim = NULL, log_scale = NULL) {
    plot <- ggplot(data = correlation_tibble, aes(x = lambda, y = correlation, color = celltype)) +
        geom_line(linewidth = 1) +
        {
            if (!is.null(log_scale)) {
                scale_x_log10()
            }
        } +
        geom_point() +
        theme_minimal() +
        ggtitle(expression(bold(paste(R[c])))) +
        labs(
            x = "Regularization",
            y = "",
            color = "Celltype"
        ) +
        theme(text = element_text(size = 20)) +
        if (!is.null(lim)) {
            ylim(lim[1], lim[2])
        }
    output_path <- file.path(output_dir, file_name)
    ggsave(output_path, plot = plot, dpi = 300)
    return(plot)
}

plot_sample_correlations <- function(mean_sample_correlations, file_name, lim = NULL, log_scale = NULL) {
    box_plot <- ggplot(data = mean_sample_correlations, aes(x = lambda, y = correlation, group = lambda)) +
        {
            if (!is.null(log_scale)) {
                scale_x_log10()
            }
        } +
        ggtitle(expression(bold(paste(R[s])))) +
        geom_boxplot() +
        theme_minimal() +
        labs(
            x = "Regularization",
            y = ""
        ) +
        theme(text = element_text(size = 20)) +
        if (!is.null(lim)) {
            ylim(lim[1], lim[2])
        }
    print(paste("Saving plot to:", file.path(output_dir, file_name)))
    ggsave(file.path(output_dir, file_name),
        plot = box_plot, bg = "white", width = 12, height = 12
    )
    return(box_plot)
}


for (irun in 1:number_simulation_runs) {
    # load input data
    print(paste0("run ", irun, " / ", number_simulation_runs))
    proportions_all <- readRDS(file.path(input_dir, paste0("proportions_true_run_", irun, ".rds")))
    for (lambda in lambda_seq) {
        proportions_lambda <- readRDS(file.path(output_dir, paste0("proportions_harp_true_lambda_run_", lambda, "_run_", irun, ".rds")))
        proportions_all <- rbind(proportions_all, 
            proportions_lambda %>% mutate(algo = paste0(algo, "_", lambda)))
    }
    # Now evaluate
    test_ids <- proportions_all %>%
        group_by(algo) %>%
        summarise(bulk_ids = list(unique(bulk_id))) %>%
        pull(bulk_ids) %>%
        Reduce(intersect, .)
    proportions_all <- proportions_all %>% filter(bulk_id %in% test_ids)
    correlations <- correlate_celltype(proportions_all)
    correlations <- correlations %>% mutate(lambda = as.numeric(sub(".*_", "", algo)))
    correlations_small <- correlations %>% filter(lambda <= 1)
    correlations_large <- correlations %>% filter(lambda > 1)
    ct_plot_small <- plot_ct_correlation(correlation_tibble = correlations_small, output_dir = output_dir, file_name = "lambda_ct_correlations_small.pdf", lim = c(0.97, 1.0))
    ct_plot_large <- plot_ct_correlation(correlation_tibble = correlations_large, output_dir = output_dir, file_name = "lambda_ct_correlations_large.pdf", lim = c(0.4, 1.0), log_scale = TRUE)

    sample_correlations <- correlate_samples(proportions = proportions_all) %>% mutate(lambda = as.numeric(sub(".*_", "", algo)))
    sample_correlations_small <- sample_correlations %>% filter(lambda <= 1)
    sample_correlations_large <- sample_correlations %>% filter(lambda > 1)
    sample_plot_small <- plot_sample_correlations(mean_sample_correlations = sample_correlations_small, file_name = "lambda_size_box_small.pdf", lim = c(0.997, 1))
    sample_plot_large <- plot_sample_correlations(mean_sample_correlations = sample_correlations_large, file_name = "lambda_size_box_large.pdf", lim = c(0, 1), log_scale = TRUE)


    combined_plot <- (ct_plot_small + ct_plot_large) / sample_plot_small / sample_plot_large +
    plot_layout(
        heights = c(2, 1, 1)
    ) +
    plot_annotation(tag_levels = "a") &
        theme(
            plot.tag = element_text(
                size = rel(1),
                color = "black"
            )
        )
    ggsave(file.path(output_dir, "varying_lambda_main.pdf"),
        combined_plot,
        device = cairo_pdf, width = 300, height = 300, units = "mm"
    )

}




sink()
sink(type = "message")