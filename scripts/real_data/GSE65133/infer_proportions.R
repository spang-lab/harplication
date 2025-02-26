# benchmarking the performance of HARP against BayesPrism on microarray data- input files for cibersortx are provided seprately

# DEPENDENCIES
devtools::document()

library(tidyverse)
# library(edgeR)
library(ggplot2)
library(SummarizedExperiment)
library(gridExtra)
library(MuSiC)
library(DTD)

# creating output dir
base_dir <- system.file("data/generated", package = "harplication")
dir.create(file.path(base_dir, "real_data"))
output_dir <- file.path(base_dir, "real_data", "GSE65133")
dir.create(output_dir)

# creal .log file
con <- file(file.path(output_dir, "infer_proportions_GSE65133.log"))
sink(con)
sink(con, type = "message")

# select algorithms to be benchmarked
algorithms <- c("Harp", "bp_subtypes", "cibersort")


if ("Harp" %in% algorithms) {
    library(harp)
}
if ("DTD" %in% algorithms) {
    library(DTD)
}
if ("bp_subtypes" %in% algorithms) {
    library(Seurat)
    # Note that in the bayesPrism we make use of
    # 50 cores for parallel computing specify the respective parts to your system as needed
    library(BayesPrism)
}
if ("music" %in% algorithms) {
    library(SingleCellExperiment)
    library(MuSiC)
}

# loading the data

GSE65133_data <- readRDS("data/source/processed_GSE65133.rds")

# loading the reference lm22
sc_library <- GSE65133_data$sc_library

# training data
bulk_counts_train <- GSE65133_data$bulk_counts_train
bulk_pheno_train <- GSE65133_data$bulk_pheno_train_true

# test data
bulk_counts_test <- GSE65133_data$bulk_counts_test
bulk_pheno_test <- GSE65133_data$bulk_pheno_test_true # this is only for evaluation



# preparing the input data
irun <- 1
sc_counts <- sc_library %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix()

print("Selecting gene intersect")
sc_genes <- colnames(sc_counts)
print(paste("Single cell library consists of ", length(sc_genes), "genes"))
celltype <- sc_library %>% pull(celltype)
cell_id <- sc_library %>% pull(cell_id)
bulk_genes <- rownames(bulk_counts_train)
print(paste("bulk consists of ", length(bulk_genes), "genes"))
isec <- intersect(sc_genes, bulk_genes)
print(paste(length(isec), "genes remain after filtering"))

# filtering genes
sc_counts <- sc_counts[, isec]
sc_counts_normalized <- normalize_to_count(sc_counts[, isec])
bulk_counts_train <- bulk_counts_train[isec, ]
bulk_counts_test <- bulk_counts_test[isec, ]


# prepraring inputs formats
sc_summarized_experiment <- SummarizedExperiment(list(counts = t(sc_counts)),
    rowData = DataFrame(gene.name = colnames(sc_counts)),
    colData = DataFrame(cell_id = cell_id, celltype = celltype)
)

sc_summarized_experiment_normalized <- SummarizedExperiment(list(counts = t(sc_counts_normalized)),
    rowData = DataFrame(gene.name = colnames(sc_counts)),
    colData = DataFrame(cell_id = cell_id, celltype = celltype)
)


bulk_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_train),
    rowData = DataFrame(gene.name = rownames(bulk_counts_train)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_train))
)
bulk_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_test),
    rowData = DataFrame(gene.name = rownames(bulk_counts_test)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_test))
)

bulk_pheno_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_train),
    rowData = DataFrame(celltype = rownames(bulk_pheno_train)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_train))
)

bulk_pheno_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_test),
    rowData = DataFrame(celltype = rownames(bulk_pheno_test)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_test))
)



# INFER PROPORTION FOR ALL ALGO
proportions_all <- tibble()
# adding the ground truth data for the evaluation
true_proportion <- t(bulk_pheno_test) %>%
    as_tibble(rownames = "bulk_id") %>%
    pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
    add_column(run = irun, algo = "true")

proportions_all <- rbind(proportions_all, true_proportion)

if ("Harp" %in% algorithms) {
    # infer proportions of HARP (and DTD)
    print("Inferring Harp proportions")
    harp <- benchmark_harp(
        sc_summarized_experiment = sc_summarized_experiment,
        bulk_train_summarized_experiment = bulk_train_summarized_experiment,
        bulk_test_summarized_experiment = bulk_test_summarized_experiment,
        bulk_pheno_train_summarized_experiment = bulk_pheno_train_summarized_experiment,
        n_folds = 4,
        lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1))
    )
    print("Done with Harp proportions")
    saveRDS(harp$output_harp, file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
    harp_ref <- harp$output_harp$reference_profiles$estimated_reference_second

    print("creating reference from harp file for cibersort")
    harp_ref_cibersort <- write_harp_reference_input_cibersort(reference = harp_ref)
    harp_p <- harp$proportions %>% add_column(run = irun, algo = "harp")
    saveRDS(harp_p, file.path(output_dir, paste0("proportions_harp_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        harp_p
    )
}

if ("DTD" %in% algorithms) {
    output_DTD <- benchmark_dtd(
        sc_summarized_experiment,
        bulk_test_summarized_experiment
    )
    proportions_DTD <- output_DTD$proportions
    model <- output_DTD$model
    DTD_p <- proportions_DTD %>% add_column(run = irun, algo = "DTD")
    saveRDS(DTD_p, file.path(output_dir, paste0("proportions_DTD_run_", irun, ".rds")))
    saveRDS(model, file.path(output_dir, paste0("DTD_model_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        DTD_p
    )
}
if ("cibersort" %in% algorithms) {
    print("Generate input for Cibersort")
    input_cibersort <- compute_input_cibersort(
        sc_summarized_experiment = NULL,
        bulk_summarized_experiment_test = bulk_test_summarized_experiment,
        output_harp = NULL
    )
    write_cibersort_input_txt(
        bulk_counts = input_cibersort$bulk_counts$bulk_counts_matrix,
        output_dir = output_dir,
        irun = irun
    )
}
# COMPETING ALGORITHMS on TEST SET
if ("music" %in% algorithms) {
    print("Inferring MuSiC proportions")
    music_p <- benchmark_music(sc_summarized_experiment, bulk_test_summarized_experiment) %>%
        add_column(run = irun, algo = "music")
    print("Done with music proportions")
    saveRDS(music_p, file.path(output_dir, paste0("proportions_music_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        music_p
    )
}

if ("bp_subtypes" %in% algorithms) {
    print("Inferring BayesPrism proportions")
    bp.res <- fit_bulks_bayesPrism(sc_summarized_experiment,
        bulk_test_summarized_experiment,
        clustering = NULL,
        input.type = "GEP"
    )

    print("Done with BayesPrism")
    saveRDS(bp.res, file.path(output_dir, paste0("bp_res_bayesPrism_run_", irun, ".rds")))
    bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_subtypes")
    saveRDS(bayesPrism_p, file.path(output_dir, paste0("proportions_bayesPrism_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        bayesPrism_p
    )
    print("Inferring BayesPrism proportions with reference from Harp")
    sc_summarized_experiment_harp <- SummarizedExperiment(list(counts = harp_ref),
        rowData = DataFrame(gene.name = rownames(harp_ref)),
        colData = DataFrame(cell_id = colnames(harp_ref), celltype = colnames(harp_ref))
    )
    bp.res_harp <- fit_bulks_bayesPrism(
        sc_summarized_experiment = sc_summarized_experiment_harp,
        bulk_test_summarized_experiment = bulk_test_summarized_experiment,
        clustering = NULL,
        input.type = "GEP"
    )
    saveRDS(bp.res_harp, file.path(output_dir, paste0("bp_res_bayesPrism_harp_run_", irun, ".rds")))
    bayesPrism_p_harp <- get_bayesPrism_proportions(bp.res_harp) %>% add_column(run = irun, algo = "bp_harp")
    proportions_all <- rbind(
        proportions_all,
        bayesPrism_p_harp
    )
}


print("Saving all proportions")
saveRDS(proportions_all, file.path(output_dir, "proportions_all.rds"))

sink()
sink(type = "message")
