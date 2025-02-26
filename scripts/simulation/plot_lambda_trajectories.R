# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(patchwork)

plt_list <- list()

# lambda trajectories
for (irun in 1:number_simulation_runs) {
    for (harp_type in c("true")) {
        output_harp <- readRDS(file.path(output_dir, paste0("output_harp_", harp_type, "_run_", irun, ".rds")))
        for (cor_run in c("mean_cor_first", "mean_cor_second")) {
            mean_cor <- output_harp$cv[[cor_run]]
            lambda_seq <- output_harp$cv$lambda_seq
            if (harp_type == "cdeath" && cor_run == "mean_cor_first") {
                filter_lambda <- lambda_seq[lambda_seq < 10]
            } else {
                filter_lambda <- lambda_seq[lambda_seq < 10]
            }
            plt_list[[cor_run]] <- plot_lambda_trajectories(mean_cor = mean_cor,
                                        lambda_seq = lambda_seq,
                                        output_dir = output_dir,
                                        harp_type = harp_type,
                                        filter_lambda = filter_lambda,
                                        file_name = paste0("lambda_trajectories_", cor_run, "_run_", irun, ".pdf"),
                                        lim = c(0.945, 0.96))
        }
    }
}

combined_plot <- (plt_list$mean_cor_first + plt_list$mean_cor_second) +
    plot_annotation(tag_levels = "a") &
        theme(
            plot.tag = element_text(
                size = rel(1),
                color = "black"
            )
        )
    ggsave(file.path(output_dir, "lambda_trajectory.pdf"),
        combined_plot,
        device = cairo_pdf, width = 400, height = 200, units = "mm"
    )

