# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(harp)
library(MuSiC)

# REPRODUCABILITY
set.seed(42)

con <- file(file.path(output_dir, "compare_cdeath_simulation.log"))
sink(con)
sink(con, type = "message")


# HELPERS
convert_alpha_tibble <- function(alpha, mode) {
    celltypes <- alpha %>% names()
    alpha_tibble <- alpha %>%
        as_tibble() %>%
        mutate(celltype = celltypes) %>%
        dplyr::rename(alpha = value) %>%
        add_column(mode = mode)
    return(alpha_tibble)
}

plot_alpha <- function(alpha_tibble, output_dir, data_name) {
    ggplot(alpha_tibble, aes(x = celltype, y = alpha, fill = mode)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_minimal() +
        labs(
            x = "",
            y = ""
        ) +
        theme(
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white"),
            panel.border = element_blank(),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.line = element_line(colour = "black"),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold")
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(fill = paste("Correction/Distortion Value")) +
        scale_fill_manual(
            values = c("true" = "#eb9d76", "cdeath" = "#b95726", "delta_inv" = "#3182BD"),
            labels = c("true" = expression(alpha(paste(True))), "cdeath" = expression(alpha(paste(Distorted))), "delta_inv" = expression(delta^{
                -1
            })),
            breaks = c("true", "cdeath", "delta_inv")
        )
    ggsave(file.path(output_dir, data_name), width = 200, height = 150, units = "mm", device = cairo_pdf)
}

# INFER HARP PROPORTIONS
proportions_all <- tibble()

for (irun in 1:number_simulation_runs) {
    # load input data
    print(paste0("run ", irun, " / ", number_simulation_runs))
    sc_summarized_experiment <- readRDS(file.path(output_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
    bulk_train_summarized_experiment <- readRDS(file.path(output_dir, paste0("bulk_train_summarized_experiment_run_", irun, ".rds")))
    bulk_test_summarized_experiment <- readRDS(file.path(output_dir, paste0("bulk_test_summarized_experiment_run_", irun, ".rds")))
    bulk_pheno_train_true_summarized_experiment <- readRDS(file.path(output_dir, paste0("bulk_pheno_train_true_summarized_experiment_run_", irun, ".rds")))
    bulk_pheno_train_cdeath_summarized_experiment <- readRDS(file.path(output_dir, paste0("bulk_pheno_train_cdeath_summarized_experiment_run_", irun, ".rds")))
    proportions_true <- readRDS(file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        proportions_true
    )

    # if we ran the main benchmark, we already have the version with the true proportions
    if (file.exists(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))) {
        output_harp_true <- readRDS(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))
        harp_proportions_true <- readRDS(file.path(output_dir, paste0("proportions_harp_true_run_", irun, ".rds")))
    } else {
        # infer proportions of HARP
        print("Inferring Harp true proportions with lambda being")
        print(lambda_seq)
        # first with TRUE proportions, i.e., no cell death
        harp_true <- benchmark_harp(
            sc_summarized_experiment = sc_summarized_experiment,
            bulk_train_summarized_experiment = bulk_train_summarized_experiment,
            bulk_test_summarized_experiment = bulk_test_summarized_experiment,
            bulk_pheno_train_summarized_experiment = bulk_pheno_train_true_summarized_experiment,
            lambda_seq = lambda_seq,
            n_folds = 5
        )
        print("Done with Harp true proportions")
        output_harp_true <- harp_true$output_harp
        saveRDS(output_harp_true, file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))

        harp_proportions_true <- harp_true$proportions %>% add_column(run = irun, algo = "harp_true")
        saveRDS(harp_proportions_true, file.path(output_dir, paste0("proportions_harp_true_run_", irun, ".rds")))
    }
    proportions_all <- rbind(
        proportions_all,
        harp_proportions_true
    )


    # now with cell death
    print("Inferring Harp cdeath proportions with lambda being")
    print(lambda_seq)
    harp_cdeath <- benchmark_harp(
        sc_summarized_experiment = sc_summarized_experiment,
        bulk_train_summarized_experiment = bulk_train_summarized_experiment,
        bulk_test_summarized_experiment = bulk_test_summarized_experiment,
        bulk_pheno_train_summarized_experiment = bulk_pheno_train_cdeath_summarized_experiment,
        lambda_seq = lambda_seq,
        n_folds = 5
    )
    print("Done with Harp cell death proportions")
    harp_proportions_cdeath <- harp_cdeath$proportions %>% add_column(run = irun, algo = "harp_cdeath")
    saveRDS(harp_proportions_cdeath, file.path(output_dir, paste0("proportions_harp_cdeath_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        harp_proportions_cdeath
    )
    # store output of harp pipeline
    output_harp_cdeath <- harp_cdeath$output_harp
    saveRDS(output_harp_cdeath, file.path(output_dir, paste0("output_harp_cdeath_run_", irun, ".rds")))
}

print("Saving all proportions")
saveRDS(proportions_all, file.path(output_dir, "proportions_cdeath_comparison.rds"))

# Now evaluate
irun <- 1
proportions_all <- readRDS(file.path(output_dir, "proportions_cdeath_comparison.rds"))
output_harp_true <- readRDS(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))
output_harp_cdeath <- readRDS(file.path(output_dir, paste0("output_harp_cdeath_run_", irun, ".rds")))


data_name <- "cdeath_comparison"
test_ids <- proportions_all %>%
    group_by(algo) %>%
    summarise(bulk_ids = list(unique(bulk_id))) %>%
    pull(bulk_ids) %>%
    Reduce(intersect, .)
proportions_all <- proportions_all %>% filter(bulk_id %in% test_ids)
proportions_all <- proportions_all %>% filter(algo %in% c("true", "harp_true", "harp_cdeath"))
plot_list <- plot_all_individual_eval_statistics_bars(proportions = proportions_all,
                                            output_dir = output_dir,
                                            data_name = data_name,
                                            file_type = "pdf",
                                            combined_width = 12, combined_height = 6,
                                            individual_width = 6, individual_height = 6,
                                            y_min = 0.001,
                                            y_min_r = 0.99,
                                            y_min_ct = 0.95,
                                            simulation = TRUE
                                            )

combined_plot <- (plot_list$combined_plot) /
    (plot_list$celltype_wise + plot_list$sample_wise) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(2.5),
            color = "black",
            face = "bold"
        )
    )
ggsave(file.path(output_dir, "cell_death_main.pdf"),
    combined_plot,
    device = cairo_pdf, width = 400, height = 350, units = "mm"
)

# Load actual delta values (the amount of distorting bulks)
delta <- readRDS(file.path(output_dir, "bulk_pheno_cell_death_rates_1.rds"))

# Analyze alpha
alpha_true <- output_harp_true$alpha
alpha_true <- convert_alpha_tibble(alpha_true, mode = "true")

alpha_cdeath <- output_harp_cdeath$alpha
alpha_cdeath <- convert_alpha_tibble(alpha_cdeath, mode = "cdeath")

delta_inverse <- delta %>% mutate(alpha = 1/rate) %>% add_column(mode = "delta_inv") %>% dplyr::select(-rate)

alpha_tibble <- rbind(alpha_true, alpha_cdeath, delta_inverse)
saveRDS(alpha_tibble, file.path(output_dir, "alpha_tibble.rds"))

plot_alpha(alpha_tibble, output_dir, "alpha_celltype.pdf")

sink()
sink(type = "message")
