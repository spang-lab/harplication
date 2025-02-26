# benchmarking the performance of HARP against BayesPrism on bulk rna-seq data using sorted rna-seq reference- also the .txt files for cibersortx input are generated
# using this script

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

# creating output dir
base_dir <- system.file("data/generated", package = "harplication")
dir.create(file.path(base_dir, "real_data"))

output_dir <- file.path(base_dir, "real_data", "sdy67_rnaseq_reference")
dir.create(output_dir)

# create log file

con <- file(file.path(output_dir, "infer_proportions_sdy67.log"))
sink(con)

sink(con, type = "message")

# reading the data
zimmermann_data <- readRDS("data/source/processed_sdy67_rnaseq_ref.rds")

sc_library <- zimmermann_data$sc_library

sc_counts <- sc_library %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix()
gene_names <- colnames(sc_counts)

celltype <- sc_library %>% pull(celltype)
cell_id <- sc_library %>% pull(cell_id)

sc_library_unmapped <- zimmermann_data$sc_library_unmapped # this for cibersort and BayesPrism, with original cell types

# training data
bulk_counts_train <- zimmermann_data$bulk_counts_train[gene_names, ] # restrict to common genes by filtering sc genes within bulk data
bulk_pheno_train <- zimmermann_data$bulk_pheno_train_true

# test data
bulk_pheno_test <- zimmermann_data$bulk_pheno_test_true # this only will be used to evaluate the results
bulk_counts_test <- zimmermann_data$bulk_counts_test[gene_names, ] # restrict to common genes by filtering sc genes within bulk data


bulk_counts_test_all_genes <- zimmermann_data$bulk_counts_test # we keep all of the genes here,
# since it will be used to generate .txt files for cibersort (for deconvolution with lm22)

# prepare the dataformat for input of harplication
irun <- 1
sc_summarized_experiment <- SummarizedExperiment(list(counts = t(sc_counts)),
    rowData = DataFrame(gene.name = gene_names),
    colData = DataFrame(cell_id = cell_id, celltype = celltype)
)

# single cell monaco with original cell types to use in bayesprism
sc_counts_unmapped <- sc_library_unmapped %>%
    dplyr::select(-c(celltype, cell_id)) %>%
    as.matrix()
celltype_all <- sc_library_unmapped %>% pull(celltype)
cell_id_all <- sc_library_unmapped %>% pull(cell_id)

# single cell monaco with original cell types to use in bayesprism
sc_summarized_experiment_unmapped <- SummarizedExperiment(list(counts = t(sc_counts_unmapped)),
    rowData = DataFrame(gene.name = gene_names),
    colData = DataFrame(cell_id = cell_id_all, celltype = celltype_all)
)


bulk_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_train),
    rowData = DataFrame(gene.name = rownames(bulk_counts_train)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_train))
)
bulk_pheno_train_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_train),
    rowData = DataFrame(celltype = rownames(bulk_pheno_train)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_train))
)

bulk_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_counts_test),
    rowData = DataFrame(gene.name = rownames(bulk_counts_test)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_test))
)
bulk_pheno_test_summarized_experiment <- SummarizedExperiment(list(counts = bulk_pheno_test),
    rowData = DataFrame(celltype = rownames(bulk_pheno_test)),
    colData = DataFrame(bulk_id = colnames(bulk_pheno_test))
)

# this is used to provide cibersort with all genes for deconvolution with lm22
bulk_test_summarized_experiment_all_genes <- SummarizedExperiment(list(counts = bulk_counts_test_all_genes),
    rowData = DataFrame(gene.name = rownames(bulk_counts_test_all_genes)),
    colData = DataFrame(bulk_id = colnames(bulk_counts_test_all_genes))
)

# INFER PROPORTION FOR ALL ALGOS
proportions_all <- tibble()
true_proportion <- t(bulk_pheno_test) %>%
    as_tibble(rownames = "bulk_id") %>%
    pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
    add_column(run = irun, algo = "true")

proportions_all <- rbind(proportions_all, true_proportion)

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
    proportions_all <- rbind(
        proportions_all,
        harp_p
    )
}
# Generate .txt files for CIBERSORTx
if ("cibersort" %in% algorithms) {
    print("Generate input for Cibersort")
    input_cibersort <- compute_input_cibersort(
        sc_summarized_experiment = sc_summarized_experiment_unmapped,
        bulk_summarized_experiment_test = bulk_test_summarized_experiment_all_genes,
        output_harp = harp$output_harp
    )
    write_cibersort_input_txt(
        mean_reference = input_cibersort$mean_reference,
        harp_reference = input_cibersort$harp_reference,
        refsample = input_cibersort$refsample,
        bulk_counts = input_cibersort$bulk_counts$bulk_counts_matrix,
        phenotype_classes = input_cibersort$phenotype_classes,
        bulk_counts_hgnc = NULL,
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
    proportions_all <- rbind(
        proportions_all,
        DTD_p
    )
}

if ("bp_subtypes" %in% algorithms) {
    print("Done computing subclustering")
    print("Inferring BayesPrism proportions")
    bp.res <- fit_bulks_bayesPrism(sc_summarized_experiment_unmapped,
        bulk_test_summarized_experiment,
        clustering = NULL,
        input.type = "count.matrix"
    )
    saveRDS(bp.res, file.path(output_dir, paste0("bp_res_bayesPrism_run_", irun, ".rds")))
    bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_subtypes")
    bayesPrism_p <- mapping_proportions_data_celltypes(proportions = bayesPrism_p, map_name = "zimmermann_monaco") # this to map the cell types after deconvolution

    saveRDS(bayesPrism_p, file.path(output_dir, paste0("proportions_bayesPrism_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
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

    proportions_all <- rbind(
        proportions_all,
        bayesPrism_p_harp
    )
}


print("Saving all proportions")
saveRDS(proportions_all, file.path(output_dir, "proportions_all.rds"))

sink()
sink(type = "message")
