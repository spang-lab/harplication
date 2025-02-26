#' computes the proportions of a cell type within each bulk sample,
#' by comparing the celltype-wise total gene count to the bulk's total count
#' @export
compute_proportion <- function(bulk_counts, cse) {
    bulk_total <- compute_bulk_total(bulk_counts)
    proportion <- cse %>%
        # group_by(bulk_id) %>%
        # mutate(bulk_sum = sum(cse)) %>%
        group_by(bulk_id, celltype) %>%
        mutate(cse_total = sum(cse)) %>%
        ungroup() %>%
        dplyr::select(-cse) %>%
        dplyr::distinct() %>%
        inner_join(bulk_total, by = "bulk_id") %>%
        mutate(prp = cse_total / bulk_sum) %>%
        # TODO: positive select
        dplyr::select(-gene, -cse_total, -bulk_sum) %>%
        dplyr::distinct()
    return(proportion)
}

#' computes the cell type specific expression,
#' by weighting the reference according to the weights
#' resulting from the deconvolution
#' @export
compute_cse <- function(weights, reference) {
    cse <- weights %>%
        # celltype is our notion of cell-type
        inner_join(reference, by = "celltype") %>%
        pivot_longer(cols = -c(bulk_id, celltype, weight), names_to = "gene", values_to = "expression") %>%
        mutate(expression = expression * weight) %>%
        group_by(bulk_id, gene, celltype) %>%
        mutate(cse = sum(expression)) %>%
        ungroup() %>%
        # TODO: maybe positive select better
        dplyr::select(-expression, -weight) %>%
        dplyr::distinct()
    return(cse)
}

#' the total bulk count for a single bulk
#' is the sum over all gene counts within that bulk
#' @export
compute_bulk_total <- function(bulk_counts) {
    bulk_total <- bulk_counts %>%
        pivot_longer(cols = -bulk_id, names_to = "gene", values_to = "expression") %>%
        group_by(bulk_id) %>%
        mutate(bulk_sum = sum(expression)) %>%
        ungroup() %>%
        dplyr::select(bulk_id, bulk_sum) %>%
        dplyr::distinct()

    return(bulk_total)
}

#' computes the cell-type specific expression for the
#' simulated bulks by using the true cell weights stored
#' in the composition
#' This will yield our ground truth
#' @export
compute_bulk_cse <- function(composition, sclong) {
    scwide <- sclong %>% pivot_wider(names_from = gene, values_from = expression, values_fill = 0)
    celltype_comp <- composition %>%
        inner_join(scwide, by = "cell_id")
    excluded_names <- intersect(names(celltype_comp), c('bulk_id', 'cell_id', 'celltype', 'weight', 'Patient', 'testtrain', 'origin'))
    celltype_comp <- celltype_comp %>%
        pivot_longer(cols = -excluded_names, names_to = "gene", values_to = "expression") %>%
        # the cse equals each weighted sc profile cumulated over the respective cell type
        mutate(expression = expression * weight) %>%
        group_by(bulk_id, celltype, gene) %>%
        mutate(cse = sum(expression)) %>%
        ungroup()
    cse <- celltype_comp %>%
        dplyr::select(bulk_id, gene, cse, celltype) %>%
        dplyr::distinct()
    return(cse)
}

#' compute the actual proportions, by weighting the reference with c
#' and determine the proportions between the sum of cse counts
#' and the sum of bulk counts
#' @export
convert_weight_to_proportion <- function(weight, reference, bulk_counts) {
    ref_wide <- reference %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "expression") %>%
        pivot_wider(names_from = gene, values_from = expression)

    reshaped_weights <- tibble(celltype = rownames(weight)) %>%
        cbind(as_tibble(weight)) %>%
        pivot_longer(!celltype, names_to = "bulk_id", values_to = "weight") %>%
        relocate(bulk_id)

    reshaped_bulk_counts <- tibble(bulk_id = colnames(bulk_counts)) %>%
        cbind(as_tibble(t(bulk_counts))) %>% as_tibble()

    cse <- compute_cse(reshaped_weights, ref_wide)
    proportions <- compute_proportion(reshaped_bulk_counts, cse)

    cse <- cse %>%
        pivot_wider(names_from = "celltype", values_from = "cse")
    return(list(proportions = proportions, cse = cse))
}

#' distort the simulated bulks by multiplying actual counts in a fixed percentage of genes
#' with gene specific random bias
#' @export
distort_bulks <- function(bulks_sim, percentage_genes_to_distort, distortion_factor_mean, distortion_factor_std) {
    bulks <- bulks_sim$bulks
    genes <- names(bulks %>% dplyr::select(-bulk_id))
    ngenes <- length(genes)
    nbulks <- length(bulks %>% pull(bulk_id))
    genes_to_modify <- genes[sample(round(percentage_genes_to_distort * ngenes), replace = FALSE)]
    for (gene in genes_to_modify) {
        modifier <- function(x) {
            return(abs(rnorm(1, distortion_factor_mean, distortion_factor_std)) * x)
        }
        bulks <- bulks %>% mutate({{ gene }} := modifier(.data[[gene]]))
    }
    print(paste("Replaced", length(genes_to_modify), "genes with random but normally distributed counts"))
    return(list(bulks_sim = list(bulks = bulks, composition = bulks_sim$composition), modified_genes = genes_to_modify))
}

#' simulate cell death in FACS measurments. We artificially diminish the proportion
#' of each celltype in each bulk with a celltype specific death rate
#' @export
kill_cells <- function(proportions, celltypes, kill_cells_mean = 1, kill_cells_std = 0) {
    # the celltype specific death rate is obtained as normally distributed
    # weight centered around kill_cells_mean with std  kill_cells_std
    # so kill_cells_mean=0.8 would mean we keep 80 perfcent of the cells
    # but we limit the death to a lower bound, in order to keep at least 20 percent of the cells
    csdr_values <- pmax(0.2, pmin(abs(rnorm(length(celltypes), kill_cells_mean, kill_cells_std)), 1))
    csdr_keys <- celltypes
    cell_death_rates <- tibble(celltype = csdr_keys, rate  = csdr_values)
    csdr <- list()
    csdr[csdr_keys] <- csdr_values
    for (key in csdr_keys) {
        print(paste("Keeping only", round(csdr[[key]] * 100), "percent of", key, "cells due to cell specific death"))
        proportions <- proportions %>% mutate(prp = ifelse(celltype == key, csdr[[key]] * prp, prp))
    }
    return(list(proportions = proportions, cell_death_rates = cell_death_rates))
}

#' function to map gene ENSEMBL names to HGNC symbols
#' @export
map_feature_names <- function(
    vector.of.feature.names,
    input.type = "hgnc_symbol",
    output.type = "ensembl_gene_id",
    verbose = TRUE,
    undo.safety.check = FALSE,
    n.tries = Inf) {
    if (!is.numeric(n.tries)) {
        n.tries <- 10
    }

    err_msg <- paste0(
        "in 'map_feature_names': Your provided 'vector.of.feature.names'",
        " is invalid, empty or consists of some NAs; Please provide",
        " a vector of characters!"
    )
    # some safety checks
    try_err <- try(base::exists("vector.of.feature.names"), silent = TRUE)
    if (class(try_err) == "try-error") {
        stop(err_msg)
    } else {
        if (!all(is.character(vector.of.feature.names)) ||
            !is.vector(vector.of.feature.names) ||
            any(is.na(vector.of.feature.names))) {
            stop(err_msg)
        } else {
            library("biomaRt")

            if (undo.safety.check) {
                library(httr)
                set_config(config(ssl_verifypeer = FALSE))
            }

            # the function calls a mirror, and then downloads the mapping in real time.
            # due to server issues, this might fail.
            # Therefore, we call the function in a while loop.
            times <- 0

            while (times < n.tries) {
                times <- times + 1
                if (verbose) {
                    print(
                        paste0(
                            "trying to connect to ensembl at ",
                            strftime(Sys.time(), format = "%H:%M:%S")
                        )
                    )
                }
                mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")

                tmp <- try(
                    mapper <- getBM(
                        filters = input.type,
                        attributes = c(input.type, output.type),
                        values = vector.of.feature.names,
                        mart = mart,
                        useCache = FALSE
                    ),
                    silent = TRUE
                )

                if (class(tmp) == "try-error") {
                    warning(
                        paste0(
                            "At ", Sys.time(),
                            "In 'map_feature_names': calling getBM failed, with message: \n",
                            tmp[1]
                        ),
                        immediate. = TRUE
                    )
                } else {
                    break
                }
            }
            if (class(tmp) == "try-error") {
                mapper <- NULL
            }
            if (verbose) {
                if (times < n.tries) {
                    print(
                        paste0(
                            "In 'map_feature_names': mapping succeded at ", Sys.time(),
                            " after ", times, " tries."
                        )
                    )
                } else {
                    print(
                        paste0(
                            "In 'map_feature_names': mapping failed at ", Sys.time(),
                            " after ", times, " tries."
                        )
                    )
                }
            }

            if (undo.safety.check) {
                library(httr)
                reset_config()
            }
            return(mapper)
        }
    }
}

#' conversion from summarized experiment
#' @export
convert_summarized_experiment_sc_library <- function(sc_summarized_experiment) {
    sc_counts <- t(assays(sc_summarized_experiment)$counts)
    sc_library <- sc_counts %>%
        as_tibble() %>%
        add_column(
            cell_id = colData(sc_summarized_experiment)$cell_id,
            celltype = colData(sc_summarized_experiment)$celltype
        )
    return(sc_library)
}

#' conversion from summarized experiment
#' @export
convert_summarized_experiment_bulk_counts <- function(bulk_summarized_experiment) {
    bulk_counts <- t(assays(bulk_summarized_experiment)$counts)
    bulks <- bulk_counts %>%
        as_tibble() %>%
        add_column(
            bulk_id = colData(bulk_summarized_experiment)$bulk_id
        )
    return(bulks)
}

#' conversion from summarized experiment
#' @export
convert_summarized_experiment_bulk_pheno <- function(bulk_pheno_summarized_experiment) {
    bulk_pheno <- assays(bulk_pheno_summarized_experiment)$counts
    colnames(bulk_pheno) <- colData(bulk_pheno_summarized_experiment)$bulk_id
    return(bulk_pheno)
}

