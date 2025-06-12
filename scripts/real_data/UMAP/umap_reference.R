# this is to generate the umap of the sorted rna-seq references and the averaged reference
devtools::document()
library(tidyverse)
library(ggplot2)
library(SRSU)
library(DTD)
library(stringr)

zimmermann_data <- readRDS("data/source/sdy67_rnaseq_ref.rds")

# create outputdir
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


ref_xsc <- compute_reference_harp(sc_library = zimmermann_data$sc_library)
# load the data with original cell types
sc_counts_unmapped <- t(zimmermann_data$sc_library_unmapped %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix())

sc_df <- sc_counts_unmapped
colnames(sc_df) <- gsub("^[^_]*_", "", colnames(sc_df))
colnames(ref_xsc) <- paste0("ref_", colnames(ref_xsc))
sc_df <- cbind(sc_df[rownames(ref_xsc), ], ref_xsc)
mapping_monaco_cells_to_coarse_celltypes <- function(X) {
    if (str_detect(X, "CD4|Th|TFH|Treg|CD8|VD|MAIT")) {
        return("T_cells")
    }
    if (str_detect(X, "mono|mDC|pDC")) {
        return("Monocytes")
    }
    if (str_detect(X, "B_")) {
        return("B_cells")
    }
    if (str_detect(X, "NK")) {
        return("NK_cells")
    }
    if (str_detect(X, "Plasmablasts")) {
        return("Plasmablasts")
    }
    if (str_detect(X, "Progenitor|Neut|Basophil")) {
        return("extra")
    }
    return(X)
}

set.seed(12)
message("Starting UMAP embedding...")
umap.embedded <- uwot::umap(
    X = t(sc_df),
    metric = "correlation",
    ret_model = TRUE,
    n_neighbors = 20,
    min_dist = 0.3,
    spread = 1
)
umap.frame <- data.frame(
    "umap1" = umap.embedded$embedding[, 1],
    "umap2" = umap.embedded$embedding[, 2],
    "celltype" = colnames(sc_df)
)
umap.frame$origin <- umap.frame$celltype

umap.frame$origin <- vapply(
    umap.frame$origin,
    mapping_monaco_cells_to_coarse_celltypes,
    character(1)
)
umap.frame$origin <- gsub("ref_", "", umap.frame$origin)

umap.frame$group <- ifelse(grepl("ref_", umap.frame$celltype), "average",
    "references"
)
umap_cells <- ggplot(data = umap.frame, aes(x = umap1, y = umap2, col = origin, shape = group)) +
    scale_color_discrete(labels = function(x) gsub("_", " ", gsub("extra", "Unidentified", x))) +
    geom_point(size = 3, stroke = 0) +
    scale_shape_manual(
        values = c(
            "references" = 16,
            "average" = 17
        )
    ) +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22)
    ) +
    labs(
        title = "",
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Cell type",
        group = "type"
    )

ggsave(file.path(output_dir, "umap_monaco_reference.pdf"), plot = umap_cells)
