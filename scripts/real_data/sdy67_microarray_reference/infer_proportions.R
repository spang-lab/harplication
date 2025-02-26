# benchmarking the performance of HARP against BayesPrism  on  bulk rna-seq data using a microarray reference (lm22)
# this also generate the .txt files as inputs for cibersortx
# DEPENDENCIES

devtools::document()
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(MuSiC)
library(gridExtra)
library(ggplot2)

algorithms <- c("Harp", "bp_subtypes", "cibersort")

# library(parallel)
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
# # creating output dir
base_dir <- system.file("data/generated", package = "harplication")

dir.create(file.path(base_dir, "real_data"))
output_dir <- file.path(base_dir, "real_data", "sdy67_microarray_reference")
dir.create(output_dir)

# create log file
con <- file(file.path(output_dir, "infer_proportions_sdy67_lm22.log"))
sink(con)
sink(con, type = "message")


# reading the data the data
zimmermann_data <- readRDS("data/source/processed_sdy67_microarray_ref.rds")

#  preprare the inputs for harplication
irun <- 1
# lm22 with mapped cell types
sc_library_lm22 <- zimmermann_data$sc_library

sc_counts_lm22 <- sc_library_lm22 %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix()

gene_names <- colnames(sc_counts_lm22)
celltype <- sc_library_lm22 %>% pull(celltype)
cell_id <- sc_library_lm22 %>% pull(cell_id)

# lm22 with the original cell type lables
sc_library_unmapped <- zimmermann_data$sc_library_unmapped # this for BayesPrism-  cibersort already have this reference(lm22) available

sc_counts_unmapped <- sc_library_unmapped %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix()
gene_names <- colnames(sc_counts_unmapped)
celltype_all <- sc_library_unmapped %>% pull(celltype)
cell_id_all <- sc_library_unmapped %>% pull(cell_id)


# restrict to common genes by filtering sc genes within bulk data
isec <- intersect(gene_names, rownames(zimmermann_data$bulk_counts_train))
bulk_counts_train <- zimmermann_data$bulk_counts_train[isec, ]
bulk_counts_test <- zimmermann_data$bulk_counts_test[isec, ]

bulk_counts_test_all_genes <- zimmermann_data$bulk_counts_test # we keep all of the genes here for cibersort

# cell proportions
bulk_pheno_train_lm22 <- zimmermann_data$bulk_pheno_train_true
bulk_pheno_test_lm22 <- zimmermann_data$bulk_pheno_test_true # this is only for evaluation


# this is for proper format of input data for harplication package
sc_summarized_experiment <- SummarizedExperiment(list(counts = t(sc_counts_lm22)),
    rowData = DataFrame(gene.name = gene_names),
    colData = DataFrame(cell_id = cell_id, celltype = celltype)
)

sc_summarized_experiment_unmapped <- SummarizedExperiment(list(counts = t(sc_counts_unmapped)),
    rowData = DataFrame(gene.name = gene_names),
    colData = DataFrame(cell_id = cell_id_all, celltype = celltype_all)
)

bulk_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_train),
    rowData = DataFrame(gene.name = rownames(bulk_counts_train)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_train))
)
bulk_pheno_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_train_lm22),
    rowData = DataFrame(celltype = rownames(bulk_pheno_train_lm22)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_train_lm22))
)
bulk_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_test),
    rowData = DataFrame(gene.name = rownames(bulk_counts_test)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_test))
)
bulk_pheno_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_test_lm22),
    rowData = DataFrame(celltype = rownames(bulk_pheno_test_lm22)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_test_lm22))
)

# this is used to provide to cibersortx with all genes for deconvolution with lm22
bulk_test_summarized_experiment_all_genes <- SummarizedExperiment(list(counts = bulk_counts_test_all_genes),
    rowData = DataFrame(gene.name = rownames(bulk_counts_test_all_genes)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_test_all_genes))
)



# INFER PROPORTION FOR ALL ALGOS
proportions_all_lm22 <- tibble()
# add the ground truth data so it can be used in the evaluations later on

true_proportion <- t(bulk_pheno_test_lm22) %>%
    as_tibble(rownames = "bulk_id") %>%
    pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
    add_column(run = irun, algo = "true")

proportions_all_lm22 <- rbind(proportions_all_lm22, true_proportion)


if ("Harp" %in% algorithms) {
    # infer proportions of HARP (and DTD)
    print("Inferring DTD/Harp proportions without cell death")
    harp <- benchmark_harp(sc_summarized_experiment,
        bulk_train_summarized_experiment,
        bulk_test_summarized_experiment,
        bulk_pheno_train_summarized_experiment,
        lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1))
    )


    print("Done with DTD/Harp")
    saveRDS(harp$output_harp, file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
    harp_ref <- harp$output_harp$reference_profiles$estimated_reference_second
    harp_p <- harp$proportions %>% add_column(run = irun, algo = "harp")
    saveRDS(harp_p, file.path(output_dir, paste0("proportions_harp_run_", irun, ".rds")))
    proportions_all_lm22 <- rbind(
        proportions_all_lm22,
        harp_p
    )
}
# Generate .txt files for CIBERSORTx to be uploaded on their website
if ("cibersort" %in% algorithms) {
    print("Generate input for Cibersort")
    input_cibersort <- compute_input_cibersort(
        sc_summarized_experiment = NULL,
        bulk_summarized_experiment_test = bulk_test_summarized_experiment_all_genes,
        output_harp = harp$output_harp
    )
    write_cibersort_input_txt(
        harp_reference = input_cibersort$harp_reference,
        bulk_counts = input_cibersort$bulk_counts$bulk_counts_matrix,
        output_dir = output_dir,
        irun = irun
    )
}

# COMPETING ALGORITHMS on TEST SET
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
    proportions_all_lm22 <- rbind(
        proportions_all_lm22,
        DTD_p
    )
}

if ("bp_subtypes" %in% algorithms) {
    print("Inferring BayesPrism proportions")
    bp.res <- fit_bulks_bayesPrism(sc_summarized_experiment_unmapped,
        bulk_test_summarized_experiment,
        clustering = NULL,
        input.type = "GEP"
    )
    saveRDS(bp.res, file.path(output_dir, paste0("bp_res_bayesPrism_run_", irun, ".rds")))
    bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_subtypes")
    bayesPrism_p <- mapping_proportions_data_celltypes(proportions = bayesPrism_p, map_name = "zimmermann_lm22") # this to map the cell types after deconvolution
    saveRDS(bayesPrism_p, file.path(output_dir, paste0("proportions_bayesPrism_run_", irun, ".rds")))
    proportions_all_lm22 <- rbind(
        proportions_all_lm22,
        bayesPrism_p
    )
    print("Inferring BayesPrism proportions with reference from harp")
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
    saveRDS(bayesPrism_p_harp, file.path(output_dir, paste0("proportions_bayesPrism_harp_ref_run_", irun, ".rds")))

    proportions_all_lm22 <- rbind(
        proportions_all_lm22,
        bayesPrism_p_harp
    )
}


print("Saving all proportions")
saveRDS(proportions_all_lm22, file.path(output_dir, "proportions_all_lm22.rds"))

sink()
sink(type = "message")
