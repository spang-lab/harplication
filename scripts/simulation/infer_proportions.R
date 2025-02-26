# DEPENDENCIES
library(harplication)
# set parameters
config <- config::get(file = system.file("R/inst/config", "config.yml", package = "harplication"), use_parent = FALSE)
if (config$standalone) {
    source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
}


library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
if ("Harp" %in% algorithms) {
    library(harp)
}
if ("DTD" %in% algorithms) {
    library(DTD)
}
if ("bp_subtypes" %in% algorithms || "bp_harp" %in% algorithms) {
    library(Seurat)
    # Note that in the bayesPrism we make use of
    # 50 cores for parallel computing specify the respective parts to your system as needed
    library(BayesPrism)
}
if ("music" %in% algorithms) {
    library(SingleCellExperiment)
    library(MuSiC)
}

# REPRODUCABILITY
set.seed(42)
if (!is.null(subset_bulks)) {
    print(paste0("sample_", subset_bulks))
    results_dir <- file.path(output_dir, paste0("sample_", subset_bulks))
    print(paste("Creating", results_dir))
    dir.create(file.path(results_dir))
} else {
    results_dir <- output_dir
}
con <- file(file.path(results_dir, "infer_proportions_simulation.log"))
sink(con)
sink(con, type = "message")

# INFER PROPORTION FOR ALL ALGOS
proportions_all <- tibble()

for (irun in 1:number_simulation_runs) {
    # load input data
    print(paste0("run ", irun, " / ", number_simulation_runs))
    sc_summarized_experiment <- readRDS(file.path(output_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
    bulk_train_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_train_summarized_experiment_run_", irun, ".rds"))
    )
    bulk_pheno_train_true_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_pheno_train_true_summarized_experiment_run_", irun, ".rds"))
    )
    # We provide the option to subset bulks, i.e. use only a subset of bulks in training.
    if (!is.null(subset_bulks)) {
        bulk_ids <- sample(bulk_train_summarized_experiment$bulk_id, subset_bulks, replace = FALSE)
        bulk_train_summarized_experiment <- bulk_train_summarized_experiment[,bulk_train_summarized_experiment$bulk_id %in% bulk_ids]
        bulk_pheno_train_true_summarized_experiment <- bulk_pheno_train_true_summarized_experiment[, bulk_pheno_train_true_summarized_experiment$bulk_id %in% bulk_ids]
    }
    # we take maximally 5 folds, but if the number of train bulks is to small we choose the number of folds,
    # such that at least 5 samples are in each fold
    number_bulks <- length(bulk_train_summarized_experiment$bulk_id)
    n_folds <- round(number_bulks / 5)
    if (n_folds < 2) {
        n_folds <- round(number_bulks / 3)
    }
    if (n_folds > 5) {
        n_folds <- 5
    }
    bulk_test_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_test_summarized_experiment_run_", irun, ".rds")))
    bulk_pheno_train_cdeath_summarized_experiment <- readRDS(
        file.path(output_dir, paste0("bulk_pheno_train_cdeath_summarized_experiment_run_", irun, ".rds")))
    proportions_true <- readRDS(
        file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        proportions_true
    )


    if ("Harp" %in% algorithms) {
        # infer proportions of HARP
        print("Inferring Harp proportions with lambda being")
        print(lambda_seq)
        harp_true <- benchmark_harp(
            sc_summarized_experiment = sc_summarized_experiment,
            bulk_train_summarized_experiment = bulk_train_summarized_experiment,
            bulk_test_summarized_experiment = bulk_test_summarized_experiment,
            bulk_pheno_train_summarized_experiment = bulk_pheno_train_true_summarized_experiment,
            n_folds = n_folds,
            lambda_seq = lambda_seq)
        print("Done with Harp")
        saveRDS(harp_true$output_harp, file.path(results_dir, paste0("output_harp_true_run_", irun, ".rds")))
        harp_p_true <- harp_true$proportions %>% add_column(run = irun, algo = "harp_true")
        saveRDS(harp_p_true, file.path(results_dir, paste0("proportions_harp_true_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            harp_p_true
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
        saveRDS(DTD_p, file.path(results_dir, paste0("proportions_DTD_run_", irun, ".rds")))
        saveRDS(model, file.path(results_dir, paste0("DTD_model_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            DTD_p)
    }

    # COMPETING ALGORITHMS on TEST SET
    # Preprocessing for CIBERSORT in order to be passed to Docker Container
    if ("cibersort" %in% algorithms && "Harp" %in% algorithms) {
        # Generate txt files for CIBERSORTx (the actual CIBERSORT algorithm is applied in an external docker container)
        output_harp <- readRDS(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))
        input_cibersort <- compute_input_cibersort(
            sc_summarized_experiment = sc_summarized_experiment,
            bulk_summarized_experiment_test = bulk_test_summarized_experiment,
            output_harp = output_harp
        )
        write_cibersort_input_txt(
            harp_reference = input_cibersort$harp_reference,
            mean_reference = input_cibersort$mean_reference,
            refsample = input_cibersort$refsample,
            bulk_counts = input_cibersort$bulk_counts$bulk_counts_matrix,
            bulk_counts_hgnc = input_cibersort$bulk_counts$bulk_counts_matrix_hgnc,
            output_dir = output_dir,
            irun = irun
        )
    }
    if ("music" %in% algorithms) {
        print("Inferring MuSiC proportions")
        music_p <- benchmark_music(sc_summarized_experiment, bulk_test_summarized_experiment) %>%
            add_column(run = irun, algo = "music")
        print("Done with music proportions")
        saveRDS(music_p, file.path(results_dir, paste0("proportions_music_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            music_p
        )
    }

    if ("bp_subtypes" %in% algorithms) {
        print("Computing subclustering of sc data for BayesPrism")
        sub_clustering <- compute_subclustering(sc_summarized_experiment)
        saveRDS(sub_clustering, file.path(results_dir, paste0("clustering_bp_run_", irun, ".rds")))
        print("Done computing subclustering")
        print("Inferring BayesPrism proportions using single cell data")
        bp.res <- fit_bulks_bayesPrism(sc_summarized_experiment, bulk_test_summarized_experiment, sub_clustering)
        bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_subtypes")
        saveRDS(bayesPrism_p, file.path(results_dir, paste0("proportions_bayesPrism_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            bayesPrism_p
        )
    }

    if ("bp_harp" %in% algorithms) {
        # in order to have a fair comparison we compute also BayesPrism with the "naive" X1
        # reference
        sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
        print("Computing X1 for use in BayesPrism")
        reference <- compute_reference_harp(sc_library)
        reference_summarized_experiment <- SummarizedExperiment(list(counts = reference),
            rowData = DataFrame(gene.name = rownames(reference)),
            colData = DataFrame(cell_id = colnames(reference), celltype = colnames(reference))
        )
        print("Inferring BayesPrism proportions using X1")
        bp.res <- fit_bulks_bayesPrism(reference_summarized_experiment,
            bulk_test_summarized_experiment,
            clustering = NULL,
            input.type = "GEP"
        )
        bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_X1")
        saveRDS(bayesPrism_p, file.path(results_dir, paste0("proportions_bayesPrism_X1_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            bayesPrism_p
        )
        # now use reference computed by Harp
        if (file.exists(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))) {
            output_harp <- readRDS(file.path(output_dir, paste0("output_harp_true_run_", irun, ".rds")))
            harp_ref <- output_harp$reference_profiles$estimated_reference_second
            print("Inferring BayesPrism proportions using Harp reference")
            sc_summarized_experiment_harp <- SummarizedExperiment(list(counts = harp_ref),
                rowData = DataFrame(gene.name = rownames(harp_ref)),
                colData = DataFrame(cell_id = colnames(harp_ref), celltype = colnames(harp_ref)))
            bp.res <- fit_bulks_bayesPrism(sc_summarized_experiment_harp,
                bulk_test_summarized_experiment,
                clustering = NULL,
                input.type = "GEP")
            bayesPrism_p <- get_bayesPrism_proportions(bp.res) %>% add_column(run = irun, algo = "bp_harp")
            saveRDS(bayesPrism_p, file.path(results_dir, paste0("proportions_bayesPrism_harp_run_", irun, ".rds")))
            proportions_all <- rbind(
                proportions_all,
                bayesPrism_p
            )
        } else {
            print("Please run inferrence for Harp first, otherwise BayesPrism with Harp reference can not be computed.")
        }
    }
}

print("Saving all proportions")
saveRDS(proportions_all, file.path(results_dir, "proportions_all.rds"))

sink()
sink(type = "message")