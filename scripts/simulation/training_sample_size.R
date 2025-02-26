# DEPENDENCIES
library(harplication)
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(harp)
library(MuSiC)
library(gridExtra)
library(ggplot2)
library(patchwork)

plot_music_metrics <- function(metric, data, output_dir, file_name, y_axis = "Pearson Correlation") {
    plot <- ggplot(data = data, aes(x = sample, y = !!sym(metric))) +
        geom_line(linewidth = 1) +
        geom_point() +
        theme_minimal() +
        labs(
            x = "Amount Training Samples",
            y = ""
        ) +
        ggtitle(as.character(y_axis)) +
        theme(text = element_text(size = 25),
        plot.title = element_text(face = "bold"))
    output_path <- file.path(output_dir, file_name)
    ggsave(output_path, plot = plot, dpi = 300, device = cairo_pdf)
    return(plot)
}

plot_ct_correlation <- function(correlation_tibble, output_dir, file_name, lim = NULL) {
    plot <- ggplot(data = correlation_tibble, aes(x = sample, y = correlation, color = celltype)) +
        geom_line(linewidth = 1) +
        geom_point() +
        theme_minimal() +
        labs(
            x = "Amount Training Samples",
            y = "",
            color = "Celltype"
        ) + theme(text = element_text(size = 25)) +
        ggtitle(expression(bold(R[c]))) +
        if (!is.null(lim)) { ylim(lim[1], lim[2])}
    output_path <- file.path(output_dir, file_name)
    ggsave(output_path, plot = plot, dpi = 300, device = cairo_pdf)
    return(plot)
}

plot_sample_correlations <- function(mean_sample_correlations, file_name, lim = NULL) {
    box_plot <- ggplot(data = mean_sample_correlations, aes(x = sample_size, y = correlation, group = sample_size)) +
        geom_boxplot(width = 5) +
        theme_minimal() +
        labs(
            x = "Amount Training Samples",
            y = ""
        ) + theme(text = element_text(size = 25)) +
        ggtitle(expression(bold(R[s])))
        if (!is.null(lim)) {
            ylim(lim[1], lim[2])
        }
    ggsave(file.path(output_dir, file_name),
        plot = box_plot, bg = "white", dpi = 300, device = cairo_pdf
    )
    return(box_plot)
}

source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))


source(system.file("scripts/simulation", "generate_input.R", package = "harplication"))


subsets <- c(seq(6, 9, by = 1), seq(10, 160, by = 5))

# INFERENCE
for (subset in subsets) {
    subset_bulks <- subset
    source(system.file("scripts/simulation", "infer_proportions.R", package = "harplication"))
}

# EVALUATION
proportions_all <- readRDS(file.path(output_dir, paste0("proportions_true_run_1.rds")))
for (subset in subsets) {
    results_dir <- file.path(output_dir, paste0("sample_", subset))
    proportions_subset <- readRDS(file.path(results_dir, paste0("proportions_harp_true_run_1.rds")))
    proportions_all <- rbind(
        proportions_all,
        proportions_subset %>% mutate(algo = paste0("harp_sample_", subset))
    )
}
saveRDS(proportions_all, file.path(output_dir, paste0("proportions_sample.rds")))

# METRICS

# music metrics
proportions_all <- readRDS(file.path(output_dir, "proportions_sample.rds"))
eval_statistics <- calculate_metrics(proportions_all, output_dir, file_name = "sample_size_statistcs.pdf")
data <- eval_statistics %>% mutate(sample = as.numeric(sub(".*_", "", algo)))
saveRDS(data, file.path(output_dir, "sample_size_music_metrics.rds"))
r_plot <- plot_music_metrics("R", data, output_dir, file_name = "sample_size_R.pdf", y_axis = expression(R))
mad_plot <- plot_music_metrics("mAD", data, output_dir, file_name = "sample_size_AAD.pdf", y_axis = expression(mAD))
rmsd_plot <- plot_music_metrics("RMSD", data, output_dir, file_name = "sample_size_RMSD.pdf", y_axis = expression(RMSD))

# cell type wise
test_ids <- proportions_all %>%
    group_by(algo) %>%
    summarise(bulk_ids = list(unique(bulk_id))) %>%
    pull(bulk_ids) %>%
    Reduce(intersect, .)
proportions_all <- proportions_all %>% filter(bulk_id %in% test_ids)
correlations <- correlate_celltype(proportions_all) %>% mutate(sample = as.numeric(sub(".*_", "", algo)))
saveRDS(correlations, file.path(output_dir, "sample_size_ct_correlation.rds"))
ct_plot <- plot_ct_correlation(correlation_tibble = correlations, output_dir = output_dir, file_name = "sample_size_ct_correlations.pdf")

# sample wise
sample_correlations <- correlate_samples(proportions = proportions_all) %>% mutate(sample_size = as.numeric(sub(".*_", "", algo)))
mean_sample_correlations <- sample_correlations %>%
    mutate(sample_size = as.numeric(sub(".*_", "", algo))) %>%
    dplyr::select(-algo) %>%
    group_by(sample_size) %>%
    mutate(mean = mean(correlation)) %>%
    mutate(sd = sd(correlation)) %>%
    ungroup() %>%
    dplyr::select(sample_size, mean, sd) %>%
    dplyr::distinct()
s_plot <- plot_sample_correlations(mean_sample_correlations = sample_correlations, file_name = "sample_size_box.pdf", lim = c(0.9, 1))

combined_plot <- (r_plot + mad_plot) /
    (rmsd_plot + ct_plot) /
    s_plot +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(1.5), # Increase size (adjust number as needed)
            color = "black" # Change color if desired

        )
    )
ggsave(file.path(output_dir, "training_amount_main.pdf"),
    combined_plot,
    device = cairo_pdf, width = 400, height = 460, units = "mm"
)
