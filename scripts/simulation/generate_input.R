# SIMULATE ARTIFICIAL BULKS from real single cell data
# DEPENDENCIES
library(harplication)
library(tidyverse)
library(SummarizedExperiment)
library(DTD)
# set parameters
config <- config::get(file = system.file("R/inst/config", "config.yml", package = "harplication"), use_parent = FALSE)
if (config$standalone) {
    source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
}
con <- file(file.path(output_dir, "generate_input_simulation.log"))
sink(con)
sink(con, type = "message")

# REPRODUCABILITY
set.seed(1)

proportions_all <- tibble()

for (irun in 1:number_simulation_runs) {
    print(paste0("run ", irun, " / ", number_simulation_runs))

    print("Preparing single cell data")
    sclong <- preprocess_sc_bulk(single_cell_file, sc_library_patients, bulk_patients, simulation_name)
    sclong <- filter_sc_bulk_common_celltypes(sclong)
    saveRDS(sclong, file.path(output_dir, paste0("sclong_run_", irun, ".rds")))

    # the part of the sc that simulates the bulks
    sc_bulksim <- sclong %>%
        filter(testtrain == "test") %>%
        dplyr::select(cell_id, celltype, expression, gene, Patient) %>%
        pivot_wider(names_from = "gene", values_from = "expression")

    # the part of the sc that constitutes the sc library
    sc_library <- sclong %>%
        filter(testtrain == "train") %>%
        dplyr::select(cell_id, celltype, expression, gene) %>%
        pivot_wider(names_from = "gene", values_from = "expression")

    # generate single cell experiment out of that
    sc_counts <- sc_library %>%
        dplyr::select(-cell_id, -celltype) %>%
        as.matrix()
    gene_names <- colnames(sc_counts)
    celltype <- sc_library %>% pull(celltype)
    cell_id <- sc_library %>% pull(cell_id)
    sc_counts <- DTD::normalize_to_count(t(sc_counts))


    sc_summarized_experiment <- SummarizedExperiment(list(counts = sc_counts),
        rowData = DataFrame(gene.name = gene_names),
        colData = DataFrame(cell_id = cell_id, celltype = celltype)
    )
    saveRDS(sc_summarized_experiment, file.path(output_dir, paste0("sc_summarized_experiment_run_", irun, ".rds")))
    # in order to separate the training samples for HARP on the patient level
    bulk_patients <- sc_bulksim %>%
        pull(Patient) %>%
        unique()
    # random split of single cell patients
    if (!is.null(bulks_simulation_train_split)) {
        train_pos <- sample(seq_along(bulk_patients),
            round(bulks_simulation_train_split * length(bulk_patients)),
            replace = FALSE
        )
        train_patients <- bulk_patients[train_pos]
        test_patients <- bulk_patients[-train_pos]
        # customized split
    } else {
        if (length(train_sc_patients) > 0 && length(test_sc_patients) > 0) {
            train_patients <- train_sc_patients
            test_patients <- test_sc_patients
            if (length(intersect(train_patients, test_patients)) > 0) {
                stop("Training and test set must be disjoint. Please provide another split.")
            }
        } else {
            stop("Please either provide explicit train and test patients or a split in the config.yaml.")
        }
    }
    print(paste("Simulating", amount_train_bulks, "training bulks from the patients:"))
    print(train_patients)
    print(paste("Simulating", amount_test_bulks, "test bulks from the patients:"))
    print(test_patients)

    # generate now train/test artificial bulks
    train_bulks_sim <- simulate_bulks(sc_bulksim,
        amount_train_bulks,
        train_patients,
        fuzzyness = 0.3
    )
    test_bulks_sim <- simulate_bulks(sc_bulksim,
        amount_test_bulks,
        test_patients,
        fuzzyness = 0.3
    )

    # Combine train and test bulks to a single tibble with a training and test block
    test_bulks_sim$bulks <- test_bulks_sim$bulks %>% mutate(bulk_id = bulk_id + amount_train_bulks)
    test_bulks_sim$composition <- test_bulks_sim$composition %>% mutate(bulk_id = bulk_id + amount_train_bulks)

    # bind train and test bulks together
    bulks_sim <- c()
    bulks_sim$bulks <- rbind(train_bulks_sim$bulks, test_bulks_sim$bulks)
    bulks_sim$composition <- rbind(train_bulks_sim$composition, test_bulks_sim$composition)
    saveRDS(bulks_sim, file.path(output_dir, paste0("bulks_sim_run_", irun, ".rds")))

    # record the train/test bulk_ids in split
    train_samples <- 1:amount_train_bulks
    test_samples <- (amount_train_bulks + 1):(amount_train_bulks + amount_test_bulks)
    print("Done simulating bulks")

    # infer the true proportions (on actual bulks, NON-distorted)
    print("Inferring true proportions")
    true <- true_quantities(bulks_sim$bulks, bulks_sim$composition, sclong)
    true_prp <- true$proportions
    print("Done with true proportions")
    tp_save <- true_prp %>% add_column(run = irun, algo = "true")
    saveRDS(tp_save, file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))

    # now we distort the bulk expression profile by introducting random counts in
    # a fixed percentage of genes
    if (!is.null(percentage_genes_to_distort)) {
        print("Distorting bulks")
        dist_bulks <- distort_bulks(bulks_sim, percentage_genes_to_distort, distortion_factor_mean, distortion_factor_std)
        # overwriting true bulks, everything downstream works on distorted bulks!
        bulks_sim <- dist_bulks$bulks_sim
        modified_genes <- dist_bulks$modified_genes
        saveRDS(modified_genes, file.path(output_dir, paste0("modified_genes_run_", irun, ".rds")))
        print("Done distorting bulks")
    }
    train_bulks <- bulks_sim$bulks %>% filter(bulk_id %in% train_samples)
    bulk_counts <- train_bulks %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()

    gene_names <- colnames(bulk_counts)
    bulk_id <- train_bulks %>% pull(bulk_id)

    bulk_counts <- DTD::normalize_to_count(t(bulk_counts))
    bulk_summarized_experiment_train <- SummarizedExperiment(list(counts = bulk_counts),
        rowData = DataFrame(gene.name = gene_names),
        colData = DataFrame(bulk_id = bulk_id)
    )
    saveRDS(
        bulk_summarized_experiment_train,
        file.path(output_dir, paste0("bulk_train_summarized_experiment_run_", irun, ".rds"))
    )

    test_bulks <- bulks_sim$bulks %>% filter(bulk_id %in% test_samples)
    bulk_counts <- test_bulks %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()

    gene_names <- colnames(bulk_counts)
    bulk_id <- test_bulks %>% pull(bulk_id)

    bulk_counts <- DTD::normalize_to_count(t(bulk_counts))
    bulk_summarized_experiment_test <- SummarizedExperiment(list(counts = bulk_counts),
        rowData = DataFrame(gene.name = gene_names),
        colData = DataFrame(bulk_id = bulk_id)
    )
    saveRDS(
        bulk_summarized_experiment_test,
        file.path(output_dir, paste0("bulk_test_summarized_experiment_run_", irun, ".rds"))
    )

    print("Computing bulk pheno for Harp (true/death)")
    bulk_pheno <- compute_bulk_pheno_harp(true_prp %>% filter(bulk_id %in% train_samples), sc_library, kill_cells_mean, kill_cells_std, noise_mean, noise_std)
    celltype <- rownames(bulk_pheno$true)
    bulk_id <- colnames(bulk_pheno$true)
    bulk_pheno_experiment <- SummarizedExperiment(list(counts = bulk_pheno$true),
        rowData = DataFrame(celltype = celltype),
        colData = DataFrame(bulk_id = bulk_id)
    )
    saveRDS(bulk_pheno_experiment, file.path(output_dir, paste0("bulk_pheno_train_true_summarized_experiment_run_", irun, ".rds")))
    celltype <- rownames(bulk_pheno$cdeath)
    bulk_id <- colnames(bulk_pheno$cdeath)
    bulk_pheno_experiment <- SummarizedExperiment(list(counts = bulk_pheno$cdeath),
        rowData = DataFrame(celltype = celltype),
        colData = DataFrame(bulk_id = bulk_id)
    )
    saveRDS(bulk_pheno_experiment, file.path(output_dir, paste0("bulk_pheno_train_cdeath_summarized_experiment_run_", irun, ".rds")))
    saveRDS(bulk_pheno$cell_death_rates, file.path(output_dir, paste0("bulk_pheno_cell_death_rates_", irun, ".rds")))
    celltype <- rownames(bulk_pheno$noise)
    bulk_id <- colnames(bulk_pheno$noise)
    bulk_pheno_experiment <- SummarizedExperiment(list(counts = bulk_pheno$noise),
        rowData = DataFrame(celltype = celltype),
        colData = DataFrame(bulk_id = bulk_id)
    )
    saveRDS(bulk_pheno_experiment, file.path(output_dir, paste0("bulk_pheno_train_noise_summarized_experiment_run_", irun, ".rds")))
    saveRDS(bulk_pheno$noise_rates, file.path(output_dir, paste0("bulk_pheno_noise_rates_", irun, ".rds")))
}

sink()
sink(type = "message")
