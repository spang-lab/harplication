#' the profiles in the reference for harp are given as cell-type mean over sc_library
#' @export
compute_reference_harp <- function(sc_library) {
    reference <- tibble()
    clusters <- sc_library %>%
        pull(celltype) %>%
        unique()
    reference_tibble <- sc_library %>%
        pivot_longer(cols = -c(cell_id, celltype), names_to = "gene", values_to = "count") %>%
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
    rownames(reference) <- reference_tibble %>% pull(celltype)
    reference <- t(reference)
    # Normalize
    reference <- DTD::normalize_to_count(reference)
    return(reference)
}

#' give the bulk expression in the appropriate format
#' @export
compute_bulk_counts_harp <- function(bulk_counts) {
    bulks_harp <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bulks_harp) <- bulk_counts %>% pull(bulk_id)
    bulks_harp <- t(bulks_harp)
    return(bulks_harp)
}

#' give the ground truth proportions in the adequate format
#' this is basically the simulated FACS data where we give the option
#' of simulating cell death
#' cdeath controlls wether we introduce artificial cell death in the "FACS" data
#' @export
compute_bulk_pheno_harp <- function(true_prp, sc_library, kill_cells_mean = 1, kill_cells_std = 0) {
    # we only use the clusters available in the sc_library
    # as other clusters are not part of the deconvolution
    available_celltypes <- sc_library %>%
        pull(celltype) %>%
        unique()
    bulk_pheno_mat_true <- bulk_pheno_to_matrix(true_prp, available_celltypes)
    # in order to simulate cell death in FACS measurments we artificially kill
    # some cells by diminishing the proportions in the bulk composition
    death_prp <- kill_cells(true_prp, available_celltypes, kill_cells_mean, kill_cells_std)
    bulk_pheno_mat_cdeath <- bulk_pheno_to_matrix(death_prp$proportions, available_celltypes)
    return(list(true = bulk_pheno_mat_true, cdeath = bulk_pheno_mat_cdeath, cell_death_rates = death_prp$cell_death_rates))
}

#' this just converts the proportions tibble to a matrix
#' serving as valuable input for HARP
#' @export
bulk_pheno_to_matrix <- function(proportions, available_celltypes) {
    bulk_pheno <- proportions %>%
        filter(celltype %in% available_celltypes) %>%
        pivot_wider(names_from = bulk_id, values_from = prp) %>%
        replace(is.na(.), 0)
    bulk_pheno_mat <- bulk_pheno %>%
        dplyr::select(-celltype) %>%
        as.matrix()
    rownames(bulk_pheno_mat) <- bulk_pheno %>% pull(celltype)
    return(bulk_pheno_mat)
}

#' @export
preprocess_harp_input <- function(reference, bulk_counts, bulk_pheno, train_samples, test_samples) {
    print("Computing bulk counts for Harp")
    bulk_counts <- compute_bulk_counts_harp(bulk_counts)
    # filter data by predefined train and test samples
    training_data <- list(
        "mixtures" = bulk_counts[, train_samples],
        "quantities" = bulk_pheno[, train_samples]
    )
    test_data <- list(
        "mixtures" = bulk_counts[, test_samples],
        "quantities" = bulk_pheno[, test_samples]
    )
    return(list(
        reference = reference,
        train_samples = train_samples,
        test_samples = test_samples,
        training_data = training_data,
        test_data = test_data
    ))
}

#' @export
generate_input_harp <- function(prep_harp, train_samples, test_samples) {
    # filter data by predefined train and test samples
    training_data_true <- list(
        "mixtures" = prep_harp$bulk.counts[, train_samples],
        "quantities" = prep_harp$bulk.pheno[, train_samples]
    )
    test_data_true <- list(
        "mixtures" = prep_harp$bulk.counts[, test_samples],
        "quantities" = prep_harp$bulk.pheno[, test_samples]
    )
    training_data_cdeath <- list(
        "mixtures" = prep_harp$bulk.counts[, train_samples],
        "quantities" = prep_harp$bulk.pheno$cdeath[, train_samples]
    )
    test_data_cdeath <- list(
        "mixtures" = prep_harp$bulk.counts[, test_samples],
        "quantities" = prep_harp$bulk.pheno$cdeath[, test_samples]
    )
    return(list(
        reference = prep_harp$reference,
        train_samples = train_samples,
        test_samples = test_samples,
        training_data_true = training_data_true,
        test_data_true = test_data_true,
        training_data_cdeath = training_data_cdeath,
        test_data_cdeath = test_data_cdeath
    ))
}

convert_mixture_quantities <- function(bulk_counts_harp, bulk_pheno = NULL) {
    harp_data <- list(
        "mixtures" = bulk_counts_harp,
        "quantities" = bulk_pheno
    )
    return(harp_data)
}

#' @export
benchmark_harp <- function(sc_summarized_experiment,
                           bulk_train_summarized_experiment,
                           bulk_test_summarized_experiment,
                           bulk_pheno_train_summarized_experiment,
                           n_folds = 5,
                           lambda_seq,
                           ...) {
    sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
    bulk_counts_train <- convert_summarized_experiment_bulk_counts(bulk_train_summarized_experiment)
    bulk_counts_test <- convert_summarized_experiment_bulk_counts(bulk_test_summarized_experiment)
    bulk_pheno_train <- convert_summarized_experiment_bulk_pheno(bulk_pheno_train_summarized_experiment)
    print("Computing X_sc for Harp")
    reference <- compute_reference_harp(sc_library)
    print("Converting train/test counts")
    bulk_counts_train_harp <- compute_bulk_counts_harp(bulk_counts_train)
    bulk_counts_test_harp <- compute_bulk_counts_harp(bulk_counts_test)

    # normalize to count is essential for alpha pipeline to work
    bulk_counts_train_harp <- DTD::normalize_to_count(bulk_counts_train_harp)
    bulk_counts_test_harp <- DTD::normalize_to_count(bulk_counts_test_harp)


    training_data <- convert_mixture_quantities(bulk_counts_train_harp, bulk_pheno_train)


    celltype <- rownames(training_data$quantities)
    reference <- reference[, celltype]
    print(paste("Using", n_folds, "folds in cv"))
    output.harp.pipeline <- harp_pipeline(
        train_data = training_data,
        cell_reference_profile = reference,
        bulk_data = bulk_counts_test_harp,
        n_folds = n_folds,
        lambda_seq = lambda_seq,
        verbose = TRUE,
        ...
    )

    # bring to consistent format
    proportions_cse <- convert_weight_to_proportion(
        output.harp.pipeline$estimated_c$estimated_c_second,
        output.harp.pipeline$reference_profiles$estimated_reference_second,
        bulk_counts_test_harp
    )
    proportions <- proportions_cse$proportions
    cse <- proportions_cse$cse

    # for checking also uncorrected proportions (i.e., first run of alpha pipeline)
    proportions_cse_uncorrected <- convert_weight_to_proportion(
        output.harp.pipeline$estimated_c$estimated_c_first,
        output.harp.pipeline$reference_profiles$estimated_reference_first,
        bulk_counts_test_harp
    )
    proportions_uncorrected <- proportions_cse_uncorrected$proportions

    return(list(
        proportions = proportions,
        proportions_uncorrected = proportions_uncorrected,
        cse = cse,
        output_harp = output.harp.pipeline
    ))
}

#' @export
generate_input_dtd <- function(sc_library) {
    expr_data <- tibble()
    expr_data_tibble <- sc_library %>%
        pivot_longer(cols = -c(cell_id, celltype), names_to = "gene", values_to = "count") %>%
        dplyr::select(cell_id, gene, count)
    expr_data_tibble <- expr_data_tibble %>%
        pivot_wider(names_from = gene, values_from = count)
    expr_data <- expr_data_tibble %>%
        dplyr::select(-cell_id) %>%
        as.matrix()
    rownames(expr_data) <- expr_data_tibble %>% pull(cell_id)
    expr_data <- t(expr_data)

    clusters <- sc_library %>%
        dplyr::select(cell_id, celltype)


    indicator.list <- clusters %>% pull(celltype)
    names(indicator.list) <- clusters %>% pull(cell_id)

    include.in.X <- indicator.list %>% unique()
    include.in.X <- include.in.X[include.in.X != "extra"]
    percentage.of.all.cells <- 0.2
    set.seed(1)
    sample.X <- DTD::sample_random_X(
        included.in.X = include.in.X,
        pheno = indicator.list,
        expr.data = expr_data,
        percentage.of.all.cells = percentage.of.all.cells,
        normalize.to.count = TRUE
    )

    X.matrix <- sample.X$X.matrix
    samples.to.remove <- sample.X$samples.to.remove

    X.matrix <- X.matrix
    # remove the used samples
    train.expr <- expr_data[, -which(colnames(expr_data) %in% samples.to.remove)]
    indicator.train <- indicator.list[colnames(train.expr)]

    no_sc <- ncol(train.expr)
    if (no_sc < 300) {
        n.samples <- 50
    }
    if (1000 > no_sc && no_sc > 300) {
        n.samples <- 0.3 * ncol(train.expr)
    } else {
        n.samples <- nrow(X.matrix)
    }
    print(n.samples)
    n.per.mixture <- round(0.1 * nrow(sc_library))



    set.seed(1)
    training.data <- DTD::mix_samples(
        expr.data = train.expr,
        pheno = indicator.train,
        included.in.X = include.in.X,
        n.samples = n.samples,
        n.per.mixture = n.per.mixture,
        verbose = FALSE
    )

    return(list(reference = X.matrix, training_data = training.data))
}

#' @export
benchmark_dtd <- function(sc_summarized_experiment,
                          bulk_test_summarized_experiment) {
    sc_library <- convert_summarized_experiment_sc_library(sc_summarized_experiment)
    bulk_counts_test <- convert_summarized_experiment_bulk_counts(bulk_test_summarized_experiment)
    print("Constructing training bulks for dtd from sc library")
    input_dtd <- generate_input_dtd(sc_library)
    reference <- input_dtd$reference
    training_data <- input_dtd$training_data
    bulk_counts_test <- compute_bulk_counts_harp(bulk_counts_test)

    n.features <- nrow(reference)
    start.tweak <- rep(1, n.features)
    names(start.tweak) <- rownames(reference)
    print("DTD: Train DTD")
    DTD.DTD.x <- DTD::train_deconvolution_model(
        tweak = start.tweak,
        X.matrix = reference,
        train.data.list = training_data,
        test.data.list = NULL,
        estimate.c.type = "direct",
        use.implementation = "cxx",
        lambda.seq = "none",
        cv.verbose = TRUE,
        verbose = TRUE,
        NORM.FUN = "norm1"
    )

    bulk_counts_test <- DTD::normalize_to_count(bulk_counts_test)
    gene_order <- rownames(reference)
    bulk_counts_test <- bulk_counts_test[gene_order, ]

    # calculate c on test data
    print("DTD: Estimate c")
    esti.c.DTD.test <- estimate_c(
        X.matrix = reference,
        new.data = bulk_counts_test,
        DTD.model = DTD.DTD.x$best.model$Tweak,
        estimate.c.type = "direct"
    )
    esti.c.DTD.test[esti.c.DTD.test < 0] <- 0

    proportions_cse <- convert_weight_to_proportion(esti.c.DTD.test, reference, bulk_counts_test)
    proportions <- proportions_cse$proportions
    cse <- proportions_cse$cse

    return(list(proportions = proportions, cse = cse, model = DTD.DTD.x))
}

#' this is a function calculates the bulk gene expressions predicted by harp then calculates the pearson
#' correlation of them with real data
#' @param harp_output list
#' @param true_bulk_expression matrix with dimention of genes  * samples
#' @param exclude_celltypes a vector of strings this is the names of the celltypes you want to remove form your prediction, for example
#' the unidentified cell types if you want to
#' @return a tibble
#' @export

calculate_bulk_correlation_harp <- function(harp_output, true_bulk_expression, exclude_celltypes = NULL) {
    harp_ref <- harp_output$reference_profiles$estimated_reference_second
    harp_c <- harp_output$estimated_c$estimated_c_second
    celltype <- colnames(harp_ref)
    gene_names <- intersect(rownames(harp_ref), rownames(true_bulk_expression))
    samples <- intersect(colnames(harp_c), colnames(true_bulk_expression))
    true_bulk_expression <- true_bulk_expression[gene_names, samples]

    # predicted bulks from harp
    predicted_bulk_expression_harp <- harp_ref[gene_names, celltype] %*% harp_c[celltype, samples]

    correlations_harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp,
        true_bulk_expression = true_bulk_expression
    )

    correlations_bulk_harp <- correlations_harp %>%
        as_tibble(rownames = "bulk_id") %>%
        pivot_longer(!bulk_id, values_to = "correlation") %>%
        add_column(run = 1, algo = "harp") %>%
        select(-name)
}
