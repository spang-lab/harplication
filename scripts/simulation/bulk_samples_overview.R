# DEPENDENCIES
devtools::document()
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)


proportions <- readRDS(file.path(output_dir, "proportions_true_run_1.rds"))
train_proportions <- proportions %>% filter(bulk_id <= amount_train_bulks)
test_proportions <- proportions %>% filter(bulk_id > amount_train_bulks)

test_indices <- (amount_train_bulks + 1):(amount_train_bulks + amount_test_bulks)
test_bulks_plt <- ggplot(test_proportions, aes(x = celltype, y = prp)) +
    geom_boxplot() +
    labs(
        x = "",
        y = "Proportion"
    ) + theme(text = element_text(size = 20))
ggsave(file.path(output_dir, "test_bulks_overview.pdf"), last_plot(), width = 5, height = 5, device = cairo_pdf)

set.seed(1)
test_pie <- ggplot(test_proportions %>% filter(bulk_id %in% sample(x = test_indices, size = 10, replace = FALSE)), aes(x = "", y = prp, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(~bulk_id, ncol = 5) +
    theme_void() +
    labs(
        fill = "Celltype"
    ) +
    theme(text = element_text(size = 20))
ggsave(file.path(output_dir, "pie_bulks_test.pdf"), last_plot(), width = 8, height = 5, device = cairo_pdf)

train_indices <- 1:amount_train_bulks
train_bulks_plt <- ggplot(train_proportions, aes(x = celltype, y = prp)) +
    geom_boxplot() +
    labs(
        x = "",
        y = "Proportion"
    ) + theme(text = element_text(size = 20))
ggsave(file.path(output_dir, "train_bulks_overview.pdf"), last_plot(), width = 5, height = 5, device = cairo_pdf)

set.seed(1)
train_pie <- ggplot(train_proportions %>% filter(bulk_id %in% sample(x = train_indices, size = 10, replace = FALSE)), aes(x = "", y = prp, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(~bulk_id, ncol = 5) +
    theme_void() +
    labs(
        fill = "Celltype"
    ) +
    theme(text = element_text(size = 20))
ggsave(file.path(output_dir, "pie_bulks_train.pdf"), last_plot(), width = 8, height = 5, device = cairo_pdf)

combined_plot <- (test_bulks_plt +  train_bulks_plt) /
    (test_pie + train_pie) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(
            size = rel(1.5),
            color = "black" 
        )
    )
ggsave(file.path(output_dir, "overview_bulks_main.pdf"),
    combined_plot,
    device = cairo_pdf, width = 400, height = 250, units = "mm"
)
