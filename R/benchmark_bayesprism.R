#' this subclustering refines cell types into cell states in order
#' to provide BayesPrism with more fine grained info
#' @export
compute_subclustering <- function(sc_summarized_experiment, cluster_resolution = 0.5) {
    sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
    clustering <- sc_library %>% pull(cell_id)
    names(clustering) <- clustering
    for (cl in levels(factor(sc_library %>% pull(celltype)))) {
        exprmat <- sc_library %>%
            filter(celltype == cl) %>%
            dplyr::select(-c(cell_id, celltype)) %>%
            as.matrix()
        rownames(exprmat) <- sc_library %>%
            filter(celltype == cl) %>%
            pull(cell_id)
        print(paste("clustering for cell-type", cl))
        print(dim(exprmat))
        # too few samples
        if (dim(exprmat)[1] < 2) {
            print(paste("BayesPrism subclustering: To few samples in sc library for cell-type", cl))
            print(paste("Thus, performing no subclustering."))
            cellids <- cl
            clustering[cellids] <- cl
            next
        }
        seurobj <- CreateSeuratObject(counts = t(exprmat), min.cells = 3, min.features = 5)
        print(paste("clustering done for cell-type", cl))
        # too small cluster size
        if (dim(seurobj)[2] < 100) {
            print(paste("BayesPrism subclustering: Cluster size to small for cell-type", cl))
            print(paste("Thus, performing no subclustering."))
            cellids <- colnames(seurobj)
            clustering[cellids] <- cl
            next
        }
        # Seurat v5 needs this preprocessing for PCA, due to somewhat new required layers
        # in SeuratObjects, see: https://satijalab.org/seurat/articles/seurat5_integration
        seurobj <- NormalizeData(seurobj)
        seurobj <- FindVariableFeatures(seurobj)
        seurobj <- ScaleData(seurobj)
        seurobj <- RunPCA(seurobj, features = rownames(seurobj))

        print("find neigbours")
        seurobj <- FindNeighbors(seurobj)
        print("find clusters")
        seurobj <- FindClusters(seurobj, resolution = cluster_resolution)

        clusters <- seurobj@meta.data %>%
            rownames_to_column() %>%
            as_tibble()
        cellids <- clusters %>% pull(rowname)
        clusternames <- paste0(cl, "_", clusters %>% pull(seurat_clusters))
        clustering[cellids] <- clusternames
    }
    return(clustering)
}

#' fit bulks using BayesPrism
#' @export
fit_bulks_bayesPrism <- function(
    sc_summarized_experiment, bulk_test_summarized_experiment, clustering,
    input.type = "count.matrix") {
    sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
    bulk_counts_test <- convert_summarized_experiment_bulk_counts(bulk_test_summarized_experiment)
    # deconvolute only test_samples of HARP
    bk.dat <- bulk_counts_test %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bk.dat) <- bulk_counts_test %>% pull(bulk_id)
    sc.dat <- sc_library %>%
        dplyr::select(-c(celltype, cell_id)) %>%
        as.matrix()
    rownames(sc.dat) <- sc_library %>% pull(cell_id)

    cell.type.labels <- sc_library %>% pull(celltype)
    names(cell.type.labels) <- sc_library %>% pull(cell_id)

    if (is.null(clustering)) {
        cell.state.labels <- cell.type.labels
    } else {
        cell.state.labels <- clustering[sc_library %>% pull(cell_id)]
    }


    myPrism <- new.prism(
        reference = sc.dat,
        mixture = bk.dat,
        input.type = input.type,
        cell.type.labels = cell.type.labels,
        cell.state.labels = cell.state.labels,
        key = NULL,
        outlier.cut = 0.01,
        outlier.fraction = 0.1
    )

    print("Run bayesPrism")
    bp.res <- run.prism(prism = myPrism, n.cores = 50)
    return(bp.res)
}

#' @export
get_bayesPrism_proportions <- function(bp.res) {
    # compute and reshape cell proportions
    cell_proportions <- get.fraction(
        bp = bp.res,
        which.theta = "final",
        state.or.type = "type"
    )
    # just reshape the proportions to have format
    # consistent with other algorithms
    cell_proportions <- cell_proportions %>%
        as_tibble(rownames = "bulk_id") %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp")
    return(cell_proportions)
}

#' this is a function calculates the bulk gene expressions predicted by bayesprism then calculates the pearson
#' correlation of them with real data
#' @param bp_res bayesprism object
#' @param true_bulk_expression matrix with dimention of genes  * samples
#' @param exclude_celltypes a vector of strings this is the names of the celltypes you want to remove form your prediction, for example
#' the unidentified cell types if you want to
#' @param algo string a name for the algorithm either "bp_subtypes" or "bp_harp"
#' @return a tibble
#' @export
calculate_bulk_correlation_bp <- function(bp_res, true_bulk_expression, exclude_celltypes = NULL, algo) {
    celltypes <- colnames(bp_res@posterior.initial.cellType@theta)

    if (!is.null(exclude_celltypes)) {
        celltypes <- celltypes[!celltypes %in% exclude_celltypes]
    }
    samples <- intersect(rownames(bp_res@posterior.initial.cellType@theta), colnames(true_bulk_expression))


    # compute cell-specific expression for each cell type
    cell_specific_expression <- lapply(celltypes, function(cell_type) {
        message(paste0("Computing specific expression for ", cell_type))
        get.exp(bp = bp_res, state.or.type = "type", cell.name = cell_type)
    })

    # name the list based on celltypes
    names(cell_specific_expression) <- celltypes

    # initialize bulk expression matrix by summing all cell-specific expression matrices
    bulk_expression <- Reduce("+", cell_specific_expression)

    # Set row and column names for the bulk expression matrix from the first cell type's expression data
    rownames(bulk_expression) <- rownames(cell_specific_expression[[1]])
    colnames(bulk_expression) <- colnames(cell_specific_expression[[1]])

    bulk_expression_bp <- t(bulk_expression)
    isec <- intersect(rownames(bulk_expression_bp), rownames(true_bulk_expression))


    correlations_bp <- calculate_bulk_correlation(
        predicted_bulk_expression = bulk_expression_bp[isec, samples],
        true_bulk_expression = true_bulk_expression[isec, samples]
    )

    correlation_bp <- correlations_bp %>%
        as_tibble(rownames = "bulk_id") %>%
        pivot_longer(!bulk_id, values_to = "correlation") %>%
        add_column(run = 1, algo = algo) %>%
        select(-name)
}
