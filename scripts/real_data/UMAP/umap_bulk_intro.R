# this is to generate the plot in the introduction
devtools::document()
library(tidyverse)
library(ggplot2)
library(SRSU)
library(DTD)
library(stringr)
zimmermann_data <- readRDS("data/source/sdy67_rnaseq_ref.rds")

base_dir <- system.file("data/generated", package = "harplication")
real_data_dir <- file.path(base_dir, "real_data")

# Create real_data directory if it doesn't exist
if (!dir.exists(real_data_dir)) {
    dir.create(real_data_dir)
}

output_dir <- file.path(real_data_dir, "sdy67_rnaseq_reference")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}
# bind the smaples in train and test
bulk.counts.all <- zimmermann_data$bulk_counts
bulk.pheno.all <- zimmermann_data$bulk_pheno

# create the average reference
ref_xsc <- compute_reference_harp(sc_library = zimmermann_data$sc_library)
# celltypes <- rownames(zimmermann_data$sc_library)
celltypes <- colnames(ref_xsc)
samples <- colnames(bulk.counts.all)
# genes <- rownames(zimmermann_data$sc_library)
genes <- rownames(ref_xsc)
reference.profile <- ref_xsc[genes, celltypes]
predicted_bulks <- reference.profile %*% bulk.pheno.all[celltypes, samples]
colnames(predicted_bulks) <- paste0("predicted_", colnames(predicted_bulks))

bulk.df <- cbind(bulk.counts.all[genes, ], predicted_bulks)
# UMAP Embedding
set.seed(12)
message("Starting UMAP embedding...")
umap.embedded.bulks <- uwot::umap(
    X = t(bulk.df),
    metric = "correlation",
    ret_model = TRUE,
    n_neighbors = 20,
    min_dist = 2,
    spread = 1
)
umap.frame <- data.frame(
    "umap1" = umap.embedded.bulks$embedding[, 1],
    "umap2" = umap.embedded.bulks$embedding[, 2],
    "samples" = colnames(bulk.df)
)

umap.frame$origin <- ifelse(grepl("predicted_", rownames(umap.frame)), "Predicted Bulk Epxression", "PBMC (Zimmermann et al.)")
umap_bulk <- ggplot(data = umap.frame, aes(x = umap1, y = umap2, col = origin)) +
    geom_point(size = 2) +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = c(0.89, 1.05),
        legend.justification = c(0.89, 1.05)
    ) +
    labs(
        title = "",
        x = "UMAP 1",
        y = "UMAP 2",
        color = ""
    ) +
    scale_color_manual(
        values = alpha(c("#D16103", "#293352"))
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))

# save the plot
ggsave(file.path(output_dir, "zimmerman_monaco_umap.pdf"), plot = umap_bulk)
