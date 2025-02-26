#' @export
compute_harp_reference_cibersort <- function(output_harp) {
    reference <- output_harp$reference_profiles$estimated_reference_second
    gene_symbols <- rownames(reference)
    cell_names <- colnames(reference)
    reference_gene_names <- cbind("gene_name" = gene_symbols, reference)
    reference_matrix <- rbind(c("name", cell_names), reference_gene_names)
    return(reference_matrix)
}

#' @export
compute_mean_reference_cibersort <- function(sc_library) {
    # reference is the cell type specific expression, so each
    # celltype gets a column
    reference <- tibble()
    reference_tibble <- sc_library %>%
        pivot_longer(cols = -c(celltype, cell_id), names_to = "gene", values_to = "count") %>%
        group_by(celltype, gene) %>%
        mutate(cs_expr = mean(count)) %>%
        ungroup() %>%
        dplyr::select(celltype, gene, cs_expr) %>%
        unique()
    reference_tibble <- reference_tibble %>%
        pivot_wider(names_from = gene, values_from = cs_expr)
    reference <- reference_tibble %>%
        dplyr::select(-celltype) %>%
        as.matrix()
    cell_names <- reference_tibble %>% pull(celltype)
    rownames(reference) <- cell_names
    reference <- t(reference)
    gene_symbols <- rownames(reference)
    reference <- DTD::normalize_to_count(reference)
    reference_gene_names <- cbind("gene_name" = gene_symbols, reference)
    reference_matrix <- rbind(c("name", cell_names), reference_gene_names)
    return(reference_matrix)
}


#' @export
compute_refsample_cibersort <- function(sc_library) {
    # refsample is basically the WHOLE single cell data
    celltypes <- sc_library %>% pull(celltype)
    refsample <- sc_library %>%
        dplyr::select(-c(celltype, cell_id)) %>%
        as.matrix()
    refsample <- t(refsample)
    gene_symbols <- rownames(refsample)
    refsample_gene_names <- cbind("gene_name" = gene_symbols, refsample)
    refsample_matrix <- rbind(c("name", celltypes), refsample_gene_names)
    return(refsample_matrix)
}

#' @export
compute_phenotype_classes_cibersort <- function(sc_library) {
    celltypes <- sc_library %>% pull(celltype)
    unique_celltypes <- unique(celltypes)
    phenotype <- matrix(2,
        nrow = length(unique_celltypes), ncol = length(celltypes),
        dimnames = list(unique_celltypes, celltypes)
    )
    for (i in seq_along(unique_celltypes)) {
        row_name <- unique_celltypes[i]
        phenotype[row_name, celltypes == row_name] <- 1
    }
    unique_celltypes <- gsub("\\..", "_", unique_celltypes)
    phenotype_classes <- cbind("celltye" = unique_celltypes, phenotype)
}

#' @export
compute_bulk_counts_cibersort <- function(bulk_counts, output_dir = NULL) {
    bulk_matrix <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    sample_names <- bulk_counts %>% pull(bulk_id)
    bulk_matrix <- t(bulk_matrix)
    gene_symbols <- rownames(bulk_matrix)
    # Check if gene symbols start with "ENSG"
    if (all(startsWith(gene_symbols, "ENSG"))) {
    tryCatch({
        mapping <- map_feature_names(
            gene_symbols,
            input.type = "ensembl_gene_id",
            output.type = "hgnc_symbol",
            n.tries = 10,
            undo.safety.check = TRUE
        )
    }, error = function(e) {
        print("Ensembl mapping did not succed")
    })
    tryCatch({
        print(paste("Trying to use precomputed mapping"))
        mapping <- readRDS(mapping_file)
    }, error = function(e) {
        print("No precomputed mapping available")
    })
    tryCatch({
        # as LM22 works with HGNC symbols we need these also in our bulk data to deconvolute it
        # with LM22
        hgnc_gene_symbols <- as_tibble_col(gene_symbols, column_name = "ensembl_gene_id") %>%
            inner_join(tibble(mapping), by = "ensembl_gene_id") %>%
            pull(hgnc_symbol)
        bulk_counts_gene_names <- cbind("gene_name" = gene_symbols, bulk_matrix)
        bulk_counts_matrix <- rbind(c("name", sample_names), bulk_counts_gene_names)
        bulk_counts_gene_names_hgnc <- cbind("gene_name" = hgnc_gene_symbols, bulk_matrix)
        bulk_counts_matrix_hgnc <- rbind(c("name", sample_names), bulk_counts_gene_names_hgnc)
        return(list(bulk_counts_matrix = bulk_counts_matrix, bulk_counts_matrix_hgnc = bulk_counts_matrix_hgnc))
    }, error = function(e) {
            print("Mapping to HGNC did not succeed, so no HGNC bulks will be exported")
            bulk_counts_gene_names <- cbind("gene_name" = gene_symbols, bulk_matrix)
            bulk_counts_matrix <- rbind(c("name", sample_names), bulk_counts_gene_names)
            return(list(bulk_counts_matrix = bulk_counts_matrix))
    })
     }else {
          # If gene symbols do not start with "ENSG", skip mapping
        bulk_counts_gene_names <- cbind("gene_name" = gene_symbols, bulk_matrix)
        bulk_counts_matrix <- rbind(c("name", sample_names), bulk_counts_gene_names)
        return(list(bulk_counts_matrix = bulk_counts_matrix))
    }
}


#' @export
compute_input_cibersort <- function(sc_summarized_experiment = NULL, bulk_summarized_experiment_test = NULL, output_harp = NULL) {
    results <- list()
    if (!is.null(sc_summarized_experiment)) {
        sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
        results$mean_reference <- compute_mean_reference_cibersort(sc_library)
        results$refsample <- compute_refsample_cibersort(sc_library)
        results$phenotype_classes <- compute_phenotype_classes_cibersort(sc_library)
    }
    if (!is.null(bulk_summarized_experiment_test)) {
        bulk_counts <- convert_summarized_experiment_bulk_counts(bulk_summarized_experiment_test)

        results$bulk_counts <- compute_bulk_counts_cibersort(bulk_counts)
    }
    if (!is.null(output_harp)) {
        results$harp_reference <- compute_harp_reference_cibersort(output_harp)
    }

    results
}

#' @export
write_cibersort_input_txt <- function(
    harp_reference = NULL,
    mean_reference = NULL,
    refsample = NULL,
    bulk_counts = NULL,
    bulk_counts_hgnc = NULL,
    phenotype_classes = NULL,
    output_dir,
    irun) {
    if (!is.null(harp_reference)) {
        write.table(harp_reference,
            file = file.path(output_dir, paste0("CIBERSORTx_sigmatrix_harp_reference_run_", irun, ".txt")),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE
        )
    }
    if (!is.null(mean_reference)) {
        write.table(mean_reference,
            file = file.path(output_dir, paste0("CIBERSORTx_sigmatrix_mean_reference_run_", irun, ".txt")),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE
        )
    }
    if (!is.null(refsample)) {
        write.table(refsample,
            file = file.path(output_dir, paste0("CIBERSORTx_refsample_run_", irun, ".txt")),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE
        )
    }

    if (!is.null(bulk_counts)) {
        write.table(bulk_counts,
            file = file.path(output_dir, paste0("CIBERSORTx_bulks_run_", irun, ".txt")),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE
        )
    }
    if (!is.null(phenotype_classes)) {
        write.table(phenotype_classes,
            file = file.path(output_dir, paste0("CIBERSORTx_phenotype_classes_run_", irun, ".txt")),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE
        )
    }

    if (!is.null(bulk_counts_hgnc)) {
        tryCatch(
            {
                write.table(bulk_counts_hgnc,
                    file = file.path(output_dir, paste0("CIBERSORTx_bulks_hgnc_run_", irun, ".txt")),
                    quote = FALSE,
                    sep = "\t",
                    col.names = FALSE,
                    row.names = FALSE
                )
            },
            error = function(e) {
                print("hgnc not available")
            }
        )
    }
}

#' @export
write_harp_reference_input_cibersort <- function(reference) {
    gene_symbols <- rownames(reference)
    cell_names <- colnames(reference)
    reference_gene_names <- cbind("gene_name" = gene_symbols, reference)
    reference_matrix <- rbind(c("name", cell_names), reference_gene_names)
    write.table(reference_matrix,
        file = file.path(output_dir, paste0("CIBERSORTx_harpref_run_", irun, ".txt")),
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE
    )
    reference_matrix
}

#' @export
benchmark_cibersort <- function(output_dir) {
    cibersort_prp <- read.table(file.path(output_dir, "CIBERSORTx_Adjusted_run_1.txt"), sep = "\t", header = TRUE)
    proportions <- cibersort_prp %>%
        as_tibble() %>%
        dplyr::select(-P.value, -Correlation, -RMSE) %>%
        dplyr::rename(bulk_id = Mixture)
    return(proportions)
}

#' @export
process_lm22_output <- function(proportions, simulation_name) {
    if (simulation_name %in% c("roider_no_dist", "roider_dist")) {
        mapped <- c(
            "malignant_lymphoma",
            "malignant_lymphoma",
            "T_CD8",
            "T_CD4",
            "T_CD4",
            "T_CD4",
            "T_reg",
            "Tfh",
            "monocyte",
            "NK",
            "NK",
            "Plasmablasts"
        )
    } else {
        mapped <- c(
            "B",
            "B",
            "T_CD8",
            "T_CD4",
            "T_CD4",
            "T_CD4",
            "T_reg",
            "Tfh",
            "monocyte",
            "NK",
            "NK",
            "Plasmablasts"
        )
    }

    celltype_mapping <- tibble(
        celltype = c(
            "B.cells.naive",
            "B.cells.memory",
            "T.cells.CD8",
            "T.cells.CD4.naive",
            "T.cells.CD4.memory.resting",
            "T.cells.CD4.memory.activated",
            "T.cells.regulatory..Tregs.",
            "T.cells.follicular.helper",
            "Monocytes",
            "NK.cells.resting",
            "NK.cells.activated",
            "Plasma.cells"
        ),
        reduced_celltype = mapped
    )
    proportions <- proportions %>%
        right_join(celltype_mapping, by = "celltype") %>%
        dplyr::select(-celltype) %>%
        dplyr::rename(celltype = reduced_celltype) %>%
        group_by(celltype, bulk_id, run) %>%
        summarize(prp = sum(prp)) %>%
        relocate(bulk_id, prp) %>%
        add_column(algo = "cibersort_lm22")
}

#' @export
read_cibersort_proportions <- function(output_dir, cibersort_filename, algoname, irun) {
    cibersort_prp <- read.table(
        file.path(output_dir, paste0(cibersort_filename, "_", algoname, "_run_", irun, ".txt")),
        sep = "\t",
        header = TRUE
    )
    cibersort_proportions <- cibersort_prp %>%
        as_tibble() %>%
        dplyr::select(-P.value, -Correlation, -RMSE) %>%
        dplyr::rename(bulk_id = Mixture) %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
        add_column(run = irun, algo = paste0("cibersort_", algoname)) %>%
        mutate(bulk_id = as.character(bulk_id))
    return(cibersort_proportions)
}
#' @export
read_cibersort_proportions_result <- function(output_dir, cibersort_filename, algoname, irun) {
    cibersort_prp <- read.table(
        file.path(output_dir, cibersort_filename),
        sep = "\t",
        header = TRUE
    )
    cibersort_proportions <- cibersort_prp %>%
        as_tibble() %>%
        dplyr::select(-P.value, -Correlation, -RMSE) %>%
        dplyr::rename(bulk_id = Mixture) %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
        add_column(run = irun, algo = paste0("cibersort_", algoname)) %>%
        mutate(bulk_id = as.character(bulk_id))
    return(cibersort_proportions)
}

#' @export
add_cibersort_proportions <- function(output_dir, proportions_all, cibersort_filename, algonames, number_simulation_runs, simulation_name) {
    for (irun in 1:number_simulation_runs) {
        for (algo in algonames) {
            if (algo == "lm22") {
                cibersort_proportions <- process_lm22_output(
                    read_cibersort_proportions(output_dir, cibersort_filename, algo, irun), simulation_name
                )
            } else {
                (
                    cibersort_proportions <- read_cibersort_proportions(output_dir, cibersort_filename, algo, irun))
            }
            proportions_all <- rbind(
                proportions_all,
                cibersort_proportions
            )
        }
    }
    proportions_all <- proportions_all %>% arrange(run)
    return(proportions_all)
}

#' @export
read_cibersort_bulk_correlations_result <- function(output_dir, cibersort_filename, algoname, irun) {
    cibersort_prp <- read.table(
        file.path(output_dir, cibersort_filename),
        sep = "\t",
        header = TRUE
    )
    cibersort_bulk_correlations <- cibersort_prp %>%
        as_tibble() %>%
        dplyr::select(Correlation, Mixture) %>%
        dplyr::rename(bulk_id = Mixture) %>%
        pivot_longer(!bulk_id, values_to = "correlation") %>%
        add_column(run = irun, algo = paste0("cibersort_", algoname)) %>%
        mutate(bulk_id = as.character(bulk_id)) %>%
        select(-name)

    cibersort_bulk_correlations
}
