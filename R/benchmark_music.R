#' give the bulk expression in the appropriate format
#' @export
compute_bulk_counts_music <- function(bulk_counts) {
    bulks_music <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bulks_music) <- bulk_counts %>% pull(bulk_id)
    bulks_music <- t(bulks_music)

    return(bulks_music)
}

#' @export
benchmark_music <- function(sc_summarized_experiment, bulk_test_summarized_experiment) {
    sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
    bulk_counts <- convert_summarized_experiment_bulk_counts(bulk_test_summarized_experiment)
    # music expects the sc library as a SingleCellExperiment
    # dataset, so prepare that here
    celltype <- sc_library %>% pull(celltype)
    cell_id <- sc_library %>% pull(cell_id)

    sc_counts <- sc_library %>%
        dplyr::select(-c(celltype, cell_id)) %>%
        as.matrix()
    gene_names <- colnames(sc_counts)
    rownames(sc_counts) <- cell_id
    sc_counts <- t(sc_counts)

    steen_sce <- SingleCellExperiment(list(counts = sc_counts),
        rowData = DataFrame(gene.name = gene_names),
        colData = DataFrame(sampleID = cell_id, cellType = celltype)
    )

    bulk_counts <- compute_bulk_counts_music(bulk_counts)

    music_output <- music_prop(
        bulk.mtx = bulk_counts, sc.sce = steen_sce,
        clusters = "cellType", samples = "sampleID"
    )

    # perform the music algorithm
    proportions <- music_output$Est.prop.weighted %>% as_tibble() %>%
        cbind(tibble(bulk_id = rownames(music_output$Est.prop.weighted))) %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp")

    return(proportions)
}