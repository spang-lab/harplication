# umap for the predicted bulks from harp and compare it to real and predicted from sortedrna-seq
devtools::document()
library(tidyverse)
library(ggplot2)
library(SRSU)
library(DTD)
library(stringr)
GSE65133_data <- readRDS("data/source/GSE65133_microarray.rds")

train_size <- 12 # update this for the number of train samples #20
test_size <- 8 # update this for the number of test samples #100
num_runs <- 5 # Update this to match the number of runs you have #5
irun <- 4
# Define output directory and data name
output_dir <- file.path(
  system.file("data/generated", package = "harplication"),
  "real_data", paste0("GSE65133_", train_size, "train_", test_size, "test_", num_runs, "run")
)
# TODO:change the dir name later

if (!dir.exists(output_dir)) {
  stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}
stats_dir <- file.path(output_dir, "statistical_analysis")



harp <- readRDS(file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
harp_ref <- harp$reference_profiles$estimated_reference_second
bulk_count_test <- readRDS(file.path(output_dir, paste0("bulk_counts_test_run_", irun, ".rds")))

bulk_count_train <- readRDS(file.path(output_dir, paste0("bulk_counts_train_run_", irun, ".rds")))


celltype <- colnames(harp_ref)
gene_names <- rownames(harp_ref)
test.samples <- colnames(bulk_count_test)
train.samples <- colnames(bulk_count_train)
true_bulk_expression_test <- bulk_count_test[gene_names, test.samples]


colnames(true_bulk_expression_test) <- paste0("test_", colnames(true_bulk_expression_test))
true_bulk_expression_train <- bulk_count_train[gene_names, train.samples]

colnames(true_bulk_expression_train) <- paste0("train_", colnames(true_bulk_expression_train))

# true bulk_rna expresiion of test samples
predicted_bulk_expression_harp_test <- harp_ref %*% GSE65133_data$bulk_pheno[celltype, test.samples]

predicted_bulk_expression_harp_train <- harp_ref %*% GSE65133_data$bulk_pheno[celltype, train.samples]

colnames(predicted_bulk_expression_harp_test) <- paste0("predicted_test", colnames(predicted_bulk_expression_harp_test))

colnames(predicted_bulk_expression_harp_train) <- paste0("predicted_train", colnames(predicted_bulk_expression_harp_train))


bulk.pheno.all <- GSE65133_data$bulk_pheno
reference.profile <- compute_reference_harp(sc_library = GSE65133_data$sc_library)
reference.profile <- reference.profile[gene_names, celltype]
predicted_bulks_sc <- reference.profile %*% bulk.pheno.all[celltype, colnames(bulk.pheno.all)]
colnames(predicted_bulks_sc) <- paste0("sc_", colnames(predicted_bulks_sc))

bulk.df <- cbind(
  true_bulk_expression_test, true_bulk_expression_train,
  predicted_bulk_expression_harp_test, predicted_bulk_expression_harp_train, predicted_bulks_sc
)


# UMAP Embedding
set.seed(12)
message("Starting UMAP embedding...")
umap.embedded.bulks <- uwot::umap(
  X = t(bulk.df),
  metric = "correlation",
  ret_model = TRUE,
  n_neighbors = 15,
  min_dist = 2,
  spread = 2
)
umap.frame <- data.frame(
  "umap1" = umap.embedded.bulks$embedding[, 1],
  "umap2" = umap.embedded.bulks$embedding[, 2],
  "samples" = colnames(bulk.df)
)

umap.frame$origin <- ifelse(grepl("predicted_", rownames(umap.frame)), "Predicted Bulk Expression (Harp)",
  ifelse(grepl("sc", rownames(umap.frame)), "Predicted Bulk Expression (LM22)",
    "PBMC micorarray Bulks"
  )
)

umap.frame$group <- ifelse(grepl("test", rownames(umap.frame)), "test", "train")

umap_bulk <- ggplot(data = umap.frame, aes(
  x = umap1, y = umap2, col = origin,
  alpha = ifelse(origin == "Predicted Bulk Expression (LM22)", 1, group)
)) +
  geom_point(size = 2) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.position = c(0.6, 0.95)
  ) +
  labs(
    title = "",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Bulk samples",
    alpha = "Group samples"
  ) +
  scale_color_manual(
    values = alpha(c("#D16103", "#52854C", "#882255"))
  ) +
  scale_alpha_manual(
    values = c("test" = 1, "train" = 0.5),
    labels = c("Test Set", "Train Set")
  ) +
  guides(alpha = "none")
# guides(color = guide_legend(override.aes = list(size = 2)))

# save the plot
ggsave(file.path(stats_dir, "umap_GSE65133_lm22_harp_predicted.pdf"),
  plot = umap_bulk, width = 11, height = 11
)
