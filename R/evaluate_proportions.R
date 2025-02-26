# Here we evaluate the proportions output by the algorithms
# We use the metrics from the MuSic package (see their vignette)

#' this function transforms the proportions dataframe output
#' in our benchmarking to a matrix that can be input to the MuSic
#' evaluation function
#' @export
generate_proportion_matrices <- function(proportions) {
    algos <- proportions %>%
        pull(algo) %>%
        unique()
    proportions <- proportions %>% mutate(prp = ifelse(prp < 0, 0, prp))
    proportion_matrices <- list()
    for (alg in algos) {
        alg_proportions <- proportions %>%
            filter(algo == alg) %>%
            dplyr::select(bulk_id, celltype, prp) %>%
            pivot_wider(names_from = celltype, values_from = prp)
        bulk_ids <- alg_proportions %>% pull(bulk_id)
        alg_prp_mat <- alg_proportions %>%
            dplyr::select(-bulk_id) %>%
            as.matrix()
        rownames(alg_prp_mat) <- bulk_ids
        proportion_matrices <- append(proportion_matrices, list(alg_prp_mat))
    }

    names(proportion_matrices) <- algos
    return(proportion_matrices)
}

#' This function computes quality scores for each algorithm and returns
#' a proportions heatmap for each algorithm
#' @export
perform_evaluation <- function(proportion_matrices, output_dir, simulation) {
    true_proportions <- proportion_matrices$true
    estimated_proportions <- proportion_matrices[names(proportion_matrices) != "true"]

    eval_statistics <- Eval_multi(
        prop.real = true_proportions,
        prop.est = estimated_proportions,
        method.name = names(estimated_proportions)
    )
    print("Final evaluation statistics are: ")
    print(eval_statistics)

    algos <- names(proportion_matrices)

    proportion_matrices_sampled <- list()
    for (algo in algos) {
        prop_algo <- proportion_matrices[[algo]]
        proportion_matrices_sampled <- append(proportion_matrices_sampled, list(prop_algo))
    }
    names(proportion_matrices_sampled) <- algos

    true_proportions <- proportion_matrices_sampled$true
    estimated_proportions <- proportion_matrices_sampled[names(proportion_matrices) != "true"]

    plot <- Prop_comp_multi(
        prop.real = true_proportions,
        prop.est = estimated_proportions,
        eval = FALSE,
        method.name = names(estimated_proportions)
    )
    return(list(eval_statistics = eval_statistics, plot = plot))
}

#' @export
correlate_celltype <- function(proportions) {
    proportions <- proportions %>% filter(celltype != "extra")
    proportions_algo <- proportions %>%
        filter(algo != "true") %>%
        arrange(algo)
    proportions_true <- proportions %>% filter(algo == "true")
    algos <- proportions %>%
        pull(algo) %>%
        unique()

    correlations <- tibble()
    celltypes_true <- proportions_true %>%
        pull(celltype) %>%
        unique()
    for (alg in algos) {
        prop_alg <- proportions_algo %>%
            filter(algo == alg) %>%
            arrange(celltype)
        celltypes_alg <- prop_alg %>%
            pull(celltype) %>%
            unique()
        for (ct in intersect(celltypes_true, celltypes_alg)) {
            prop_true_ct <- proportions_true %>%
                filter(celltype == ct) %>%
                arrange(bulk_id) %>%
                pull(prp)
            if (abs(sd(prop_true_ct)) < 1e-6) {
                print(paste(
                    "The true proporions of celltype",
                    ct, "have no variance. So we do not compute metrics here."
                ))
                next
            }
            prop_alg_ct <- prop_alg %>%
                filter(celltype == ct) %>%
                arrange(bulk_id) %>%
                pull(prp)
            if (sd(prop_alg_ct) < 1e-6) {
                print(paste("Expression of algorithm", alg, "for celltype", ct, "is zero"))
            }
            correlations <- rbind(
                correlations,
                tibble(algo = alg, celltype = ct, correlation = cor(prop_true_ct, prop_alg_ct))
            )
        }
    }
    return(correlations)
}
#' @export
correlate_samples <- function(proportions) {
    proportions <- proportions %>% filter(celltype != "extra")
    proportions_algo <- proportions %>%
        filter(algo != "true") %>%
        arrange(algo)
    proportions_true <- proportions %>% filter(algo == "true")
    algos <- proportions %>%
        pull(algo) %>%
        unique()

    correlations <- tibble()
    # proportions <- proportions %>%
    #     mutate(prp = ifelse(prp < 0, 0, prp)) %>%
    #     mutate(prp = ifelse(prp > 1, 1, prp))
    celltypes_true <- proportions_true %>%
        pull(celltype) %>%
        unique()

    samples <- proportions_true %>%
        pull(bulk_id) %>%
        unique()
    for (alg in algos) {
        prop_alg <- proportions_algo %>%
            filter(algo == alg) %>%
            arrange(bulk_id)

        intersect_celltype <- intersect(celltypes_true, prop_alg$celltype)
        sample_alg <- prop_alg %>%
            pull(bulk_id) %>%
            unique()
        for (sample in sample_alg) {
            prop_true_ct <- proportions_true %>%
                filter(bulk_id == sample, celltype %in% intersect_celltype) %>%
                arrange(celltype) %>%
                pull(prp)
            if (abs(sd(prop_true_ct)) < 1e-6) {
                print(paste(
                    "The true proporions of celltype",
                    ct, "have no variance. So we do not compute metrics here."
                ))
                next
            }
            prop_alg_ct <- prop_alg %>%
                filter(bulk_id == sample, celltype %in% intersect_celltype) %>%
                arrange(celltype) %>%
                pull(prp)
            cor_sample <- cor(prop_true_ct, prop_alg_ct)
            correlations <- rbind(
                correlations,
                tibble(algo = alg, bulk_id = sample, correlation = cor_sample)
            )
        }
    }
    return(correlations)
}


#' @export
plot_correlation_heatmap <- function(correlations) {
    # generate heatmap plot
    cor_plot <- ggplot(data = correlations, aes(x = celltype, y = algo, fill = correlation)) +
        geom_tile() +
        geom_text(aes(label = round(correlation, 2)), vjust = 1) +
        ylab("algo") +
        xlab("celltype") +
        ggtitle("Cell proportion correlation wrt to ground truth for different algos") +
        scale_color_gradient2(
            low = "orange",
            mid = "#c7eb96",
            high = "#087108",
            midpoint = 0,
            limits = c(-1, 1),
            guide = "none"
        ) +
        scale_fill_gradient2(
            low = "orange",
            mid = "#c7eb96",
            high = "#087108",
            midpoint = 0,
            limits = c(-1, 1)
        ) +
        theme_classic()
    return(cor_plot)
}

#' @export
plot_eval_statistics_bar <- function(eval_statistics_all) {
    eval_statistics_summary <- eval_statistics_all %>%
        group_by(algo) %>%
        mutate(
            mean_RMSD = mean(RMSD),
            mean_mAD = mean(mAD),
            mean_R = mean(R),
            sd_RMSD = sd(RMSD),
            sd_mAD = sd(mAD),
            sd_R = sd(R)
        ) %>%
        ungroup() %>%
        dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
        dplyr::distinct() %>%
        arrange(desc(mean_R))

    eval_summary_mean <- eval_statistics_summary %>%
        dplyr::select(algo, starts_with("mean_")) %>%
        pivot_longer(-algo, names_to = "metric", values_to = "value") %>%
        mutate(metric = str_extract(metric, "[^_]+$"))

    eval_summary_sd <- eval_statistics_summary %>%
        dplyr::select(algo, starts_with("sd_")) %>%
        pivot_longer(-algo, names_to = "metric", values_to = "sd") %>%
        mutate(metric = str_extract(metric, "[^_]+$"))

    eval_statistics_plot <- inner_join(eval_summary_mean, eval_summary_sd, by = c("algo", "metric"))
    metrics <- c("R", "RMSD", "mAD")
    metric_plots <- list()
    for (met in metrics) {
        bar_summary <- ggplot(eval_statistics_plot %>% filter(metric == met), aes(x = algo, y = value, fill = algo)) +
            geom_bar(
                stat = "identity", color = "black",
                position = position_dodge()
            ) +
            geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                width = .2,
                position = position_dodge(.9)
            ) +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            ) +
            ggtitle(met)
        metric_plots[[met]] <- bar_summary
    }
    return(metric_plots)
}

#' @export
plot_correlation_bar <- function(correlations) {
    eval_statistics_celltype_summary <- correlations %>%
        group_by(algo, celltype) %>%
        mutate(
            mean_R = mean(correlation),
            sd_R = sd(correlation)
        ) %>%
        ungroup() %>%
        dplyr::select(algo, celltype, starts_with("mean_"), starts_with("sd_")) %>%
        dplyr::distinct()
    bar_summary <- ggplot(eval_statistics_celltype_summary, aes(x = celltype, y = mean_R, fill = algo)) +
        geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
        ) +
        geom_errorbar(aes(ymin = mean_R - sd_R, ymax = mean_R + sd_R),
            width = .2,
            position = position_dodge(.9)
        ) +
        ggtitle("Celltype-wise correlation")
    return(bar_summary)
}


#' plot_all_individual_eval_statistics_bars
#' @description thi is a function to plot all of the quality scores, R, RMSD,AAD and  cell type wise correlations. Required packages "patchwork", "ggplot2",
#' "tiydyverse","Stringer"
#' @param proportions data frame  with bulk_id, celltye,prp, algo in its column
#' @param data_names string the name of the data set that have been used for the comparison
#' @param file_type string either pdf or png
#' @param combined_width single float the width of the combined plot for music quality
#' @param combined_height single float the height of the combined plot for music quality
#' @param individual_width single float the width of single plots
#' @param individual_height single float the height of single plots
#' @export
plot_all_individual_eval_statistics_bars <- function(proportions,
                                                     output_dir,
                                                     data_name,
                                                     file_type = "pdf",
                                                     combined_width = 16, combined_height = 6,
                                                     individual_width = 8, individual_height = 8,
                                                     y_min = 0.04,
                                                     y_min_ct = 0.1,
                                                     y_min_r = 0.5,
                                                     angle_celltype_names = 0,
                                                     hjust_celltype_names = 0.5,
                                                     simulation = FALSE) {
    library(patchwork)
    # these are methods we are benchmarking for the paper
    mapping_model_name <- c(
        "true" = "true",
        "harp" = "Harp",
        "harp_true" = "Harp",
        "harp_cdeath" = "Harp_(Distorted)",
        "harp_cdeath_uncorrected" = "Harp_(Cell Death Uncorrected)",
        "music" = "MuSiC",
        "bp_subtypes" = "BayesPrism",
        "bp_harp" = "BayesPrism_(Harp)",
        "cibersort_lm22" = "CIBERSORT_(LM22)",
        "cibersort_sc" = "CIBERSORTx_(scRNA-seq)",
        "cibersort_rna" = "CIBERSORTx_(RNA-seq)",
        "cibersort_harp" = "CIBERSORTx_(Harp)",
        "cibersort_harp_lm22" = "CIBERSORTx_(Harp)"
    )

    # mapping the model names
    proportions <- proportions %>%
        mutate(algo = ifelse(algo %in% names(mapping_model_name), mapping_model_name[algo], algo))

    proportions$algo <- gsub("_", " ", proportions$algo)

    # calculation of quality matrices
    proportion_matrices <- generate_proportion_matrices(proportions = proportions)
    music_metrics <- perform_evaluation(
        proportion_matrices = proportion_matrices,
        output_dir = output_dir,
        simulation = data_name
    )
    correlations <- correlate_celltype(proportions = proportions)
    sample_correlations <- correlate_samples(proportions = proportions)
    # to have the same cell type names in the paper
    map_celltype_names <- function(X) {
        if (stringr::str_detect(X, "CD4")) {
            return("CD4 T cell")
        }
        if (stringr::str_detect(X, "CD8")) {
            if (simulation) {
                return("CD8 T")
            }
            return("CD8 T cell")
        }
        if (stringr::str_detect(X, "NKT_")) {
            return("NKT cells")
        }
        if (stringr::str_detect(X, "memory_B_cell")) {
            return(" Memory B cell ")
        }
        if (stringr::str_detect(X, "naive_B_cell")) {
            return(" Naive B cell ")
        }
        if (stringr::str_detect(X, "B_")) {
            return("B cell")
        }
        if (stringr::str_detect(X, "Mono|mono")) {
            return("Monocytes")
        }
        if (stringr::str_detect(X, "NK")) {
            return("NK cells")
        } else {
            return(gsub("_", " ", X))
        }
    }

    cell_type_map <- Vectorize(map_celltype_names) # mapping cell type names
    correlations$celltype <- cell_type_map(correlations$celltype)
    # Summarize statistics for the R,RMSD and AAD plot
    eval_statistics_summary <- music_metrics$eval_statistics %>%
        as_tibble(rownames = "algo") %>%
        group_by(algo) %>%
        mutate(
            mean_RMSD = mean(RMSD),
            mean_mAD = mean(mAD),
            mean_R = mean(R),
            sd_RMSD = sd(RMSD),
            sd_mAD = sd(mAD),
            sd_R = sd(R)
        ) %>%
        ungroup() %>%
        dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
        dplyr::distinct() %>%
        arrange(desc(mean_R))
    # print(eval_statistics_summary)

    # Mean and SD
    eval_summary_mean <- eval_statistics_summary %>%
        dplyr::select(algo, starts_with("mean_")) %>%
        pivot_longer(-algo, names_to = "metric", values_to = "value") %>%
        mutate(metric = str_extract(metric, "[^_]+$"))

    eval_summary_sd <- eval_statistics_summary %>%
        dplyr::select(algo, starts_with("sd_")) %>%
        pivot_longer(-algo, names_to = "metric", values_to = "sd") %>%
        mutate(metric = str_extract(metric, "[^_]+$"))

    # Combine means and SDs
    eval_statistics_plot <- inner_join(eval_summary_mean, eval_summary_sd, by = c("algo", "metric"))
    eval_statistics_plot$algo <- gsub("_", " ", eval_statistics_plot$algo)

    # data for cell type wise correlation plot
    eval_statistics_celltype_summary <- correlations %>%
        group_by(algo, celltype) %>%
        mutate(
            mean_R = mean(correlation),
            sd_R = sd(correlation)
        ) %>%
        ungroup() %>%
        dplyr::select(algo, celltype, starts_with("mean_"), starts_with("sd_")) %>%
        dplyr::distinct()
    # define the order of all algorithms that we have in the paper
    full_order_of_algos <- c(
        "Harp", "Harp (True Proportion)", "Harp (Distorted)", "Harp (Cell Death Uncorrected)",
        "CIBERSORT (LM22)", "CIBERSORTx (scRNA-seq)", "CIBERSORTx (RNA-seq)", "CIBERSORTx (Harp)",
        "BayesPrism", "BayesPrism (Harp)",
        "MuSiC", "DTD"
    )

    present_algos <- unique(eval_statistics_plot$algo)
    order_of_algos <- intersect(full_order_of_algos, present_algos)

    # handle empty order_of_algos
    if (length(order_of_algos) == 0) {
        stop("No matching algorithms found in the data.")
    }
    # "#BBBBBB"
    # define color palette and filter it to match present algorithms
    custom.col <- c(
        "#52854C", "#52854C", "#009E73", "#C3D7A4",
        "#FFDB6D", "#F4EDCA", "#C4961A", "#D16103",
        "#293352", "#4E84C4",
        "#CC6677", "#882255"
    )
    # to ensure evey algorithm gets the same color every time
    filtered_colors <- custom.col[match(order_of_algos, full_order_of_algos)]
    eval_statistics_celltype_summary <- eval_statistics_celltype_summary %>%
        filter(algo %in% order_of_algos)

    # order the algorithms for the plots
    eval_statistics_plot$algo <- factor(eval_statistics_plot$algo, levels = order_of_algos)

    # mapping of metric names
    mapping_quality_scores <- c("R" = "R", "RMSD" = "RMSD", "mAD" = "mAD")
    eval_statistics_plot$metric <- recode(eval_statistics_plot$metric, !!!mapping_quality_scores)

    # metrics
    metrics <- c("R", "RMSD", "mAD")
    metric_plots <- list()

    # loop over each metric and generate plots
    for (met in metrics) {
        bar_summary <- ggplot(
            eval_statistics_plot %>% filter(metric == met),
            aes(x = algo, y = value, fill = algo)
        ) +
            geom_bar(stat = "identity", position = position_dodge()) +
            geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                width = .2, position = position_dodge(.9)
            ) +
            theme(
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(color = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 24),
                axis.title.y = element_text(size = 24),
                plot.title = element_text(face = "bold", size = 24),
                legend.text = element_text(size = 24),
                legend.title = element_text(size = 24)
            ) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(fill = "Method", y = "") +
            ggtitle(met) +
            scale_fill_manual(values = filtered_colors)

        if (met == "R") {
            bar_summary <- bar_summary + coord_cartesian(ylim = c(y_min_r, 1))
        }
        if (met %in% c("RMSD", "mAD")) {
            bar_summary <- bar_summary + coord_cartesian(ylim = c(y_min, NA))
        }
        # Store plot with legend for separate use
        metric_plots[[paste0(met, "_with_legend")]] <- bar_summary

        # Remove legend for individual metric plots, except for R
        if (met != "RSMD") {
            bar_summary <- bar_summary + theme(legend.position = "none")
        }

        metric_plots[[met]] <- bar_summary
    }

    # Combine the three plots side by side with one shared legend

    combined_plot <- wrap_plots(metric_plots[["R"]], metric_plots[["RMSD"]], metric_plots[["mAD"]],
        guides = "collect"
    ) &
        theme(legend.position = "none")

    plot_filenames <- c("combined_quality_scores", "R_with_legend", "RMSD_with_legend", "mAD_with_legend")

    for (i in seq_along(plot_filenames)) {
        ggsave(file.path(output_dir, paste0("bar_plot_", plot_filenames[i], "_", data_name, ".", file_type)),
            plot = if (i == 1) combined_plot else metric_plots[[plot_filenames[i]]],
            width = if (i == 1) combined_width else individual_width,
            height = if (i == 1) combined_height else individual_height
        )
    }

    # order the algorithms for the plots
    eval_statistics_celltype_summary$algo <- gsub("_", " ", eval_statistics_celltype_summary$algo)
    eval_statistics_celltype_summary$algo <- factor(eval_statistics_celltype_summary$algo, levels = order_of_algos)

    bar_summary_celltype_wise <- ggplot(eval_statistics_celltype_summary, aes(x = celltype, y = mean_R, fill = algo)) +
        geom_bar(
            stat = "identity",
            position = position_dodge(), # Control the spacing within the group of bars
            width = 0.6
        ) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 24, angle = angle_celltype_names, hjust = hjust_celltype_names),
            axis.text.y = element_text(size = 24),
            axis.title.y = element_text(size = 24),
            plot.title = element_text(face = "bold", size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(fill = "Method", y = "") +
        ggtitle(expression(bold(R[c]) * "  ")) +
        coord_cartesian(ylim = c(y_min_ct, NA)) +
        scale_fill_manual(values = filtered_colors) +
        theme(legend.position = "right")
    ggsave(file.path(output_dir, paste0("bar_plot_celltypewise_correlation_", data_name, ".", file_type)),
        plot = bar_summary_celltype_wise, width = individual_width + 2, height = individual_height
    )

    # plot_boxplot_sample_wise_correlatio
    sample_correlations <- sample_correlations %>% filter(algo %in% order_of_algos)
    sample_correlations$algo <- gsub("_", " ", sample_correlations$algo)
    # order the algorithms for the plots
    sample_correlations$algo <- factor(sample_correlations$algo, levels = order_of_algos)

    # Plot a box plot of correlation by algorithm
    order_of_algos <- gsub(" ", "\n", order_of_algos)
    box_plot <- ggplot(sample_correlations, aes(x = algo, y = correlation, fill = algo)) +
        geom_boxplot() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 24), # , angle = 45,hjust = 1
            axis.text.y = element_text(size = 24),
            axis.title.y = element_text(size = 24),
            axis.line = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold", size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24)
        ) +
        labs(
            x = "",
            y = ""
        ) +
        ggtitle(expression(bold(R[s]) * "  ")) +
        scale_fill_manual(values = filtered_colors) +
        scale_x_discrete(labels = order_of_algos) +
        theme(legend.position = "none")
    # theme(axis.text.x = element_blank())
    ggsave(file.path(output_dir, paste0("box_plot_samplewise_correlation_with_names", data_name, ".", file_type)),
        plot = box_plot, bg = "white"
    )
    ggsave(file.path(output_dir, paste0("box_plot_samplewise_correlation_", data_name, ".", file_type)),
        plot = box_plot + theme(axis.text.x = element_blank()), bg = "white"
    )
    box_plot <- box_plot + theme(axis.text.x = element_blank())

    ret_list <- list(
        "combined_plot" = combined_plot,
        "individual_plots" = metric_plots,
        "celltype_wise" = bar_summary_celltype_wise,
        "sample_wise" = box_plot
    )
}

#' plot_umap_signature
#' @description this function first calculates umap then plots a umap of the signature references, reference from harp with optimal lambda
#' and lambda equal to zero. Required packages "uwot", "ggplot2","tiydyverse"
#' @export
plot_umap_signature <- function(signature_data,
                                data_name,
                                output_dir,
                                lambda_seq,
                                lambda_optimal,
                                file_type = "pdf",
                                n_neighbors = 10,
                                min_dist = 2,
                                spread = 2) {
    # check inputs
    if (!is.list(signature_data)) {
        stop("signature_data must be a list")
    }
    if (!dir.exists(output_dir)) {
        stop("output_dir does not exist.")
    }
    if (!file_type %in% c("png", "pdf", "jpeg", "tiff")) {
        stop("Invalid file type. Choose 'png', 'pdf', 'jpeg', or 'tiff'.")
    }

    # check if ordering of genes is consistent

    # convert signature_data to a data frame
    signature.df <- as.data.frame(signature_data)

    # UMAP Embedding
    set.seed(1)
    message("Starting UMAP embedding...")
    umap.embedded <- uwot::umap(
        X = t(signature.df),
        metric = "correlation",
        ret_model = TRUE,
        n_neighbors = n_neighbors,
        min_dist = min_dist,
        spread = spread
    )

    # prepare UMAP results for plotting
    umap.frame.monaco <- data.frame(
        "UMAP1" = umap.embedded$embedding[, 1],
        "UMAP2" = umap.embedded$embedding[, 2],
        "cell.type" = colnames(signature.df)
    )

    # clean up cell.type and group
    umap.frame.monaco <- umap.frame.monaco %>%
        dplyr::mutate(
            cell.type = gsub("^([^.]+)\\.|\\..*", "", cell.type) %>%
                gsub("_", " ", .),
            group = gsub("\\..*", "", rownames(.)) %>%
                gsub("_", " ", .)
        )
    lambda_optimal_index <- match(lambda_optimal, lambda_seq)
    # Create a new column called broad_group
    # Check if the value in group starts with "lambda" AND is not "lambda 1"
    # If true, assign "lambda"
    # If false, keep the original value from group
    umap.frame.monaco <- umap.frame.monaco %>%
        dplyr::mutate(broad_group = case_when(
            startsWith(group, "lambda") &
                group != "lambda 1" &
                group != paste("lambda", lambda_optimal_index) ~ "lambda",
            TRUE ~ group
        )) %>%
        dplyr::mutate(broad_group = ifelse(broad_group == "lambda 1", "naive", broad_group)) %>%
        dplyr::mutate(broad_group = ifelse(broad_group == paste("lambda", lambda_optimal_index), "lambda optimal", broad_group))

    # prepare seperate dataframes
    umap.frame.monaco <- subset(umap.frame.monaco, cell.type != "extra")
    umap.frame.monaco_sc <- subset(umap.frame.monaco, broad_group == "Single cell Signature")
    umap.frame.monaco_lambda_all <- subset(umap.frame.monaco, broad_group != "Single cell Signature")
    umap.frame.monaco_naive <- subset(umap.frame.monaco_lambda_all, broad_group == "naive")
    umap.frame.monaco_optimal <- subset(umap.frame.monaco_lambda_all, broad_group == "lambda optimal")
    umap.frame.monaco_lambda <- subset(umap.frame.monaco_lambda_all, broad_group != "naive" & broad_group != "lambda optimal")

    # add the actual lambda values to a new column
    # this is used to increase the size of points in the UMAP plot depending on lambda
    for (i in seq_along(lambda_seq)) {
        umap.frame.monaco_lambda$lambda_values[umap.frame.monaco_lambda$group == paste("lambda", i)] <- lambda_seq[i]
    }

    custom.col <- c(
        "#009E73", "#D16103", "#293352", "#4E84C4", "#FFDB6D",
        "#882255", "#CC6677"
    )
    # create UMAP plot (for sc, for optimal lambda, for naive and all other)
    umap_plot <- ggplot(data = umap.frame.monaco_sc, aes(x = UMAP1, y = UMAP2, col = cell.type)) +
        geom_point(size = 1, alpha = 0.1, shape = 16) +
        # stroke sets the thickness of the boundary of a point
        geom_point(
            data = umap.frame.monaco_naive,
            aes(x = UMAP1, y = UMAP2),
            shape = 1, size = 7, stroke = 1, alpha = 1
        ) +
        guides(shape = guide_legend(override.aes = list(size = 5))) +
        geom_point(
            data = umap.frame.monaco_optimal,
            aes(x = UMAP1, y = UMAP2),
            shape = 17, size = 7, stroke = 1, alpha = 1
        ) +
        guides(shape = guide_legend(override.aes = list(size = 5))) +
        geom_point(
            data = umap.frame.monaco_lambda,
            aes(x = UMAP1, y = UMAP2, size = lambda_values),
            shape = 22, stroke = 0.5, alpha = 1
        ) +
        # The squares get bigger with increasing lambda
        scale_size_continuous(
            range = c(6, 15)
        ) +
        guides(size = "none") +
        # Here we manually create the legend:
        # the circle is naive, the triangle lambda optimal, the squares the other lambdas (getting bigger with bigger lambda values)
        geom_point(
            data = umap.frame.monaco_lambda_all,
            aes(x = -Inf, y = -Inf, shape = broad_group),
        ) +
        scale_shape_manual(
            name = "Regularization value",
            values = c(
                "lambda optimal" = 17,
                "naive" = 1
                # "lambda" = 22
            ),
            labels = c(
                "lambda optimal" = expression("Optimal " * lambda * ""),
                "naive" = expression("Naive Solution (" * lambda * " = 0)")
                # "lambda" = "Lambda"
            )
        ) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold"),
            text = element_text(size = 20)
        ) +
        scale_color_manual(values = custom.col, name = "Cell type")

    # save plot
    output_path <- file.path(output_dir, paste0("umap_signatures", data_name, "n_neighbors_", n_neighbors, ".", file_type))
    ggsave(output_path, plot = umap_plot, device = cairo_pdf, width = 250, height = 250, units = "mm")

    message("UMAP plot saved to: ", output_path)

    # return the UMAP plot
    umap_plot
}

#' plot_lambda_trajectories
#' @description this plots the sequence of lambda regularization provided in HARP against their mean corellation
#' @param mean_cor a vector storing for each lambda in lambda_seq the mean cor over folds and cell-types
#' @param lambda_seq the sequence of lambdas in [0,1] used as regularization in HARP
#' @param output_dir the directory to store the plot
#' @param harp_type string either cellDeath or true
#' @param filter_lambda sequence of lambda values to be used
#' @param file_name output file name including file format
#' @param lim optional 2d vector for y-axis lim
#' @export
plot_lambda_trajectories <- function(mean_cor,
                                     lambda_seq,
                                     output_dir,
                                     harp_type,
                                     filter_lambda = NULL,
                                     file_name = "lambda_trajectories.png",
                                     lim = NULL) {
    data <- tibble(correlation = mean_cor, lambda = lambda_seq, run = "first") %>%
        mutate(run = factor(run, levels = c("first", "second"), labels = c("pre alpha", "post alpha")))
    # in order to avoid log(0) warning
    if (!is.null(filter_lambda)) {
        data <- data %>% filter(lambda %in% filter_lambda)
    }
    data$lambda[data$lambda == 0] <- 1e-1 * data %>%
        filter(lambda > 0) %>%
        pull(lambda) %>%
        min()
    data$correlation <- as.numeric(data$correlation)
    maximum <- data %>%
        slice_max(correlation, n = 1)

    traj_plot <- ggplot(data = data, aes(x = lambda, y = correlation)) +
        geom_line(linewidth = 1) +
        # scale_x_log10() +
        geom_point(data = maximum, aes(x = lambda, y = correlation)) +
        ggtitle(expression(bold(paste("mean ", R[c])))) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            text = element_text(size = 30)
        ) +
        labs(
            x = "Regularization",
            y = ""
        ) +
        scale_color_brewer(palette = "Set1", name = "") +
        xlim(0, 2.1) +
        if (!is.null(lim)) {
            ylim(lim[1], lim[2])
        }
    output_path <- file.path(output_dir, paste0(harp_type, "_", file_name))
    ggsave(output_path, plot = traj_plot, dpi = 300, device = cairo_pdf)
    message("Trajectory plot saved to: ", output_path)
    return(traj_plot)
}

#' plot_predicted_bulks_scatter
#' @description #' this is a function that plots the scatter plots of the predicted bulk expression  from HARP and real_data
#' @param predicted_bulk_expression  matrix composed of Harp reference and facs data Y = XC
#' @param true_bulk_expression matrix
#' @param predicted_bulk matrix
#' @param file.type string
#' @param data_name string
#' @param  reference_name string the type of the data that was used to reconstruct the bulks, scRNA-seq, RNA-seq,...
#' @param output_dir string
#' @export
plot_predicted_bulks_scatter <- function(
    predicted_bulk_expression,
    true_bulk_expression,
    predicted_bulk_from_experiment,
    file_type = "png",
    reference_name = "RNA-seq data",
    data_name,
    output_dir) {
    # correlation
    cor.predicted.bulks <- harp::estimated_c_correlation(
        true_c = t(true_bulk_expression),
        estimated_c = t(predicted_bulk_expression)
    )

    cor.predicted.bulks.experiments <- harp::estimated_c_correlation(
        true_c = t(true_bulk_expression),
        estimated_c = t(predicted_bulk_from_experiment)
    )
    mean.cor <- mean(cor.predicted.bulks)
    closest.to.mean.index <- which.min(abs(cor.predicted.bulks - mean.cor))
    cor.value.closest.mean <- cor.predicted.bulks[closest.to.mean.index]
    closest.sample.name <- names(cor.predicted.bulks)[closest.to.mean.index]

    # data prepration for predicted bulk from harp
    df_harp_plot <- data.frame(
        "Recounstrcuted_bulk" = log10(predicted_bulk_expression[, closest.sample.name] + 1),
        "Bulk_RNAseq_expression" = log10(true_bulk_expression[, closest.sample.name] + 1)
    )

    custom_theme <- theme_minimal() + theme(
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold")
    )
    plot_harp_predicted_bulk <- ggplot(df_harp_plot, aes(x = Recounstrcuted_bulk, y = Bulk_RNAseq_expression)) +
        geom_point(color = "#52854C") +
        custom_theme +
        labs(x = "Bulk RNA-seq expression", y = "Recounstructed Bulk (Harp)") +
        scale_y_continuous(expand = c(0.008, 0.008)) +
        scale_x_continuous(expand = c(0.008, 0.088))

    output_path <- file.path(output_dir, paste0("harp_scatter_plot", data_name, ".", file_type))

    ggsave(
        filename = output_path,
        plot = plot_harp_predicted_bulk,
        width = 12,
        height = 12
    )

    df_experiment_plot <- data.frame(
        "Recounstrcuted_bulk" = log10(predicted_bulk_from_experiment[, closest.sample.name] + 1),
        "Bulk_RNAseq_expression" = log10(true_bulk_expression[, closest.sample.name] + 1)
    )

    plot_predicted_bulks_experiment <- ggplot(
        data = df_experiment_plot, aes(x = Recounstrcuted_bulk, y = Bulk_RNAseq_expression)
    ) +
        geom_point(color = "#D16103") +
        custom_theme +
        labs(x = "Bulk RNA-seq expression", y = paste0("Recounstructed Bulk (", reference_name, "& Flow Cytometry)")) +
        scale_y_continuous(expand = c(0.008, 0.008)) +
        scale_x_continuous(expand = c(0.008, 0.008))

    output_path <- file.path(output_dir, paste0("experiment_scatter_plot", data_name, ".", file_type))
    ggsave(
        filename = output_path,
        plot = plot_predicted_bulks_experiment,
        width = 12,
        height = 12
    )
}

#' plot_hist_predicted_bulks
#' @description  this is a function that plots the histogram of the pearson correlation of predicted bulk expression from HARP and real_data with bulk rna-seq data
#' @param predicted_bulk_expression  matrix composed of Harp reference and facs data Y = XC
#' @param true_bulk_expression matrix
#' @param predicted_bulk matrix
#' @param reference_name string the type of the data that was used to reconstruct the bulks, e.g. scRNA-seq , RNA-seq, LM22, miroarray
#' @param file.type string
#' @param data_name string
#' @param output_dir string
#' @export
plot_box_predicted_bulks <- function(predicted_bulk_expression,
                                     true_bulk_expression,
                                     predicted_bulk_from_experiment,
                                     reference_name = "RNA-seq reference",
                                     file_type = "png",
                                     data_name,
                                     output_dir,
                                     width = 9,
                                     height = 8) {
    # Compute Pearson correlations
    cor.predicted.bulks <- harp::estimated_c_correlation(
        true_c = t(true_bulk_expression),
        estimated_c = t(predicted_bulk_expression)
    )

    cor.predicted.bulks.experiments <- harp::estimated_c_correlation(
        true_c = t(true_bulk_expression),
        estimated_c = t(predicted_bulk_from_experiment)
    )

    # Combine data for plotting
    df_harp <- data.frame(value = cor.predicted.bulks, group = "Harp reference")
    df_experiment <- data.frame(value = cor.predicted.bulks.experiments, group = reference_name)

    combined_df <- rbind(df_harp, df_experiment)

    # custom theme for plot
    custom_theme <- theme_minimal() + theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 23),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 36)
        # ,plot.margin = unit(c(0, 0, 0, 0), "cm")
    )

    plot_box <- ggplot(combined_df, aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 1, outlier.shape = NA) +
        custom_theme +
        theme(legend.position = "none") +
        labs(
            x = " ",
            y = "",
            fill = ""
        ) +
        ggtitle(expression(bold(rho))) +
        scale_fill_manual(values = c("#52854C", "#882255")) + #
        scale_x_discrete(labels = function(x) {
            sapply(x, function(label) {
                bquote(atop(.(label), "& Flow cytometry"))
            })
        })

    # Save plot
    output_path <- file.path(output_dir, paste0("box_predicted_bulks_", data_name, ".", file_type))
    ggsave(output_path, plot = plot_box, device = file_type, width = width, height = height)

    message("Histogram plot saved to: ", output_path)
    plot_box
}


#' function to map the proportions after the deconvolution for different data sets
#' TODO: will be updated for all data sets that we are using in the paper
#' @export
mapping_proportions_data_celltypes <- function(proportions, map_name) {
    # function for mapping cell types to the coarse cell types
    if (map_name == "zimmermann_monaco") {
        mapping_monaco_cells_to_coarse_celltypes <- function(X) {
            case_when(
                str_detect(X, "CD4|Th|TFH|Treg|CD8|VD|MAIT") ~ "T_cells",
                str_detect(X, "mono") ~ "Monocytes",
                str_detect(X, "B_") ~ "B_cells",
                str_detect(X, "NK") ~ "NK_cells",
                str_detect(X, "Plasmablasts") ~ "Plasmablasts",
                str_detect(X, "Progenitor|Neut|mDC|pDC|Basophil") ~ "extra",
                TRUE ~ "NA"
            )
        }
        proportions <- proportions %>%
            # Map cell types using the function
            mutate(coarse_celltype = mapping_monaco_cells_to_coarse_celltypes(celltype))
    }

    if (map_name == "zimmermann_lm22") {
        mapping_lm22_cells_to_coarse_celltypes <- function(X) {
            case_when(
                is.na(X) ~ NA_character_, # Return NA for NA input
                str_detect(X, "CD4|follicular|regulatory|CD8|T\\.cells\\.gamma\\.delta") ~ "T_cells",
                str_detect(X, "NK") ~ "NK_cells",
                str_detect(X, "B\\.") ~ "B_cells",
                str_detect(X, "Mono") ~ "Monocytes",
                # str_detect(X, "Plasma") ~ "Plasmablasts",
                TRUE ~ "extra" # Return "extra" for unmatched cases
            )
        }
        proportions <- proportions %>%
            # Map cell types using the function
            mutate(coarse_celltype = mapping_lm22_cells_to_coarse_celltypes(celltype))
    }
    if (map_name == "GSE65133_lm22") {
        mapping_lm22_cells_to_coarse_celltypes_GSE65133 <- function(X) {
            case_when(
                is.na(X) ~ NA_character_, # Return NA for NA input
                str_detect(X, "CD4|follicular|regulatory") ~ "T_cells_CD4",
                str_detect(X, "T\\.cells\\.gamma\\.delta") ~ "gamma_delta_T_cell",
                str_detect(X, "CD8") ~ "T_cells_CD8",
                str_detect(X, "NK") ~ "NK",
                str_detect(X, "B.cells.naive") ~ "naive_B_cell",
                str_detect(X, "B.cells.memory") ~ "memory_B_cell",
                str_detect(X, "Mono") ~ "Monocytes",
                TRUE ~ "extra" # Return "extra" for unmatched cases
            )
        }
        proportions <- proportions %>%
            # Map cell types using the function
            mutate(coarse_celltype = mapping_lm22_cells_to_coarse_celltypes_GSE65133(celltype))
    }
    proportions <- proportions %>%
        # Group by bulk_id and the newly mapped celltype
        group_by(bulk_id, coarse_celltype, algo, run) %>%
        # Sum the prp values for each group
        summarize(prp = sum(prp), .groups = "drop")

    # Select and rename columns in base R
    proportions <- proportions[, c("bulk_id", "coarse_celltype", "prp", "run", "algo")]
    names(proportions)[names(proportions) == "coarse_celltype"] <- "celltype"
    return(proportions)
}

#' @export
calculate_metrics <- function(proportions_all, output_dir, file_name = "sample_test_statistcs.rds") {
    print("Computing quality metrics for comparison")
    # Filter for test samples
    test_ids <- proportions_all %>%
        group_by(algo) %>%
        summarise(bulk_ids = list(unique(bulk_id))) %>%
        pull(bulk_ids) %>%
        Reduce(intersect, .)
    proportions_all <- proportions_all %>% filter(bulk_id %in% test_ids)

    proportion_matrices <- generate_proportion_matrices(proportions_all)
    ev <- perform_evaluation(proportion_matrices, output_dir, simulation_name)
    eval_statistics <- ev$eval_statistics %>%
        as_tibble(rownames = "algo")
    saveRDS(eval_statistics, file.path(output_dir, file_name))
    ggsave(file.path(output_dir, file_name), tableGrob(eval_statistics))
    return(eval_statistics)
}

#' calculate_bulk_correlation
#' This function computes the Pearson correlation coefficients between
#' corresponding columns of predicted and true bulk expression.
#'
#' @param predicted_bulk_expression  matrix with rows as genes and columns as samples; predicted from each algorithm,
#' their reference %*% their predicted cell counts
#' @param true_bulk_expression matrix with rows as genes and columns as samples,
#' from true bulk gene expression data, bulkrna-seq /microarray
#' @return a vector
#' @export
calculate_bulk_correlation <- function(predicted_bulk_expression, true_bulk_expression) {
    # ensure both matrices have the same column names
    if (!all(colnames(predicted_bulk_expression) %in% colnames(true_bulk_expression))) {
        stop("Both matrices must have the same column names")
    }

    correlations <- numeric(length(colnames(predicted_bulk_expression)))

    # Iterate through the column names
    for (sample in colnames(predicted_bulk_expression)) {
        # Calculate the correlation between the columns of both matrices
        sample_index <- which(colnames(predicted_bulk_expression) == sample)
        correlations[sample_index] <- cor(predicted_bulk_expression[, sample], true_bulk_expression[, sample])
    }
    names(correlations) <- colnames(predicted_bulk_expression)
    return(correlations)
}


#' box_plot_bulk_correlation
#' @description this is a function to plot correlation of bulks - required packages "patchwork", "ggplot2",
#' "tiydyverse","Stringer"
#' @param bulk_correlations_all data frame  with bulk_id, correlation, run, algo  in its column
#' @param data_names string the name of the data set that have been used for the comparison
#' @param file_type string either pdf or png
#' @param individual_width single float the width of single plots
#' @param individual_height single float the height of single plots
#' @export
box_plot_bulk_correlation <- function(
    bulk_correlations_all,
    data_name,
    file_type = "png",
    individual_width = 8,
    individual_height = 8,
    output_dir) {
    mapping_model_name <- c(
        "true" = "true",
        "harp" = "Harp",
        "harp_true" = "Harp_(True Proportion)",
        "harp_cdeath" = "Harp_(Cell Death)",
        "harp_cdeath_uncorrected" = "Harp_(Cell Death Uncorrected)",
        "music" = "MuSiC",
        "bp_subtypes" = "BayesPrism",
        "bp_harp" = "BayesPrism_(Harp)",
        "cibersort_lm22" = "CIBERSORT_(LM22)",
        "cibersort_sc" = "CIBERSORTx_(scRNA-seq)",
        "cibersort_rna" = "CIBERSORTx_(RNA-seq)",
        "cibersort_harp" = "CIBERSORTx_(Harp)", # only one of (cibersort_harp, cibersort_harp_lm22) must be in the main data
        "cibersort_harp_lm22" = "CIBERSORTx_(Harp)"
    )
    # mapping the model names

    bulk_correlations_all <- bulk_correlations_all %>%
        mutate(algo = ifelse(algo %in% names(mapping_model_name), mapping_model_name[algo], algo))

    bulk_correlations_all$algo <- gsub("_", " ", bulk_correlations_all$algo)
    full_order_of_algos <- c(
        "Harp", "Harp (True Proportion)", "Harp (Cell Death)", "Harp (Cell Death Uncorrected)",
        "CIBERSORT (LM22)", "CIBERSORTx (scRNA-seq)", "CIBERSORTx (RNA-seq)", "CIBERSORTx (Harp)",
        "BayesPrism", "BayesPrism (Harp)",
        "MuSiC", "DTD"
    )
    present_algos <- unique(bulk_correlations_all$algo)
    order_of_algos <- intersect(full_order_of_algos, present_algos)

    # handle empty order_of_algos
    if (length(order_of_algos) == 0) {
        stop("No matching algorithms found in the data.")
    }
    # "#BBBBBB"
    # define color palette and filter it to match present algorithms
    custom.col <- c(
        "#52854C", "#52854C", "#009E73", "#C3D7A4",
        "#FFDB6D", "#F4EDCA", "#C4961A", "#D16103",
        "#293352", "#4E84C4",
        "#CC6677", "#882255"
    )
    # to ensure evey alog get the same color every time
    filtered_colors <- custom.col[match(order_of_algos, full_order_of_algos)]
    bulk_correlations_all <- bulk_correlations_all %>%
        filter(algo %in% order_of_algos) %>%
        mutate(algo = factor(algo, levels = order_of_algos))
    # Create the box plot
    order_of_algos <- gsub(" ", "\n", order_of_algos)
    box_plot <- ggplot(bulk_correlations_all, aes(x = algo, y = correlation, fill = algo)) +
        geom_boxplot() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #  axis.text.x = element_text(size = 24), # , angle = 45,hjust = 1
            axis.text.y = element_text(size = 24),
            axis.title.y = element_text(size = 24),
            axis.line = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(face = "bold", size = 36),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24)
        ) +
        labs(
            x = "",
            y = "",
            fill = "Methods"
        ) +
        ggtitle(expression(bold(rho) * "  ")) +
        scale_fill_manual(values = filtered_colors) +
        scale_x_discrete(labels = order_of_algos) +
        theme(legend.position = "none")
    theme(axis.text.x = element_blank())
    ggsave(file.path(output_dir, paste0("box_plot_bulk_correlation_", data_name, ".", file_type)),
        plot = box_plot, bg = "white", width = individual_width, height = individual_height
    )
    box_plot
}


# plot predicted bulks from Harp reference * harp_c, x_sc*harp_c, harp_ref * experiment_c, x_sc * experiment_c
#' @export
box_plot_predicted_bulks_harp_experiment <- function(harp_model,
                                                     true_bulk_expression,
                                                     bulk_pheno_test,
                                                     ref_xsc,
                                                     reference_name = "RNA-seq reference",
                                                     file_type = "png",
                                                     data_name,
                                                     output_dir,
                                                     width = 8,
                                                     height = 8) {
    harp_ref <- harp_model$reference_profiles$estimated_reference_second
    harp_c <- harp_model$estimated_c$estimated_c_second
    celltype <- colnames(harp_ref)
    gene_names <- intersect(rownames(harp_ref), rownames(true_bulk_expression))
    samples <- intersect(colnames(harp_c), colnames(true_bulk_expression))
    true_bulk_expression <- true_bulk_expression[gene_names, samples]

    # predicted bulks from harp
    predicted_bulk_expression_harp <- harp_ref[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_harp_dtdc <- harp_ref[, celltype] %*% harp_c[celltype, samples]
    # cell counts from experiment
    predicted_bulk_expression_experiment <- ref_xsc[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_experiment_dtdc <- ref_xsc[gene_names, celltype] %*% harp_c[celltype, samples]

    cor.predicted.bulks <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.harp.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp_dtdc,
        true_bulk_expression = true_bulk_expression
    )
    cor.predicted.bulks.experiments <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.experiments.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment_dtdc,
        true_bulk_expression = true_bulk_expression
    )

    # Create data frame with separate columns for reference type and method
    df <- data.frame(
        value = c(
            cor.predicted.bulks,
            cor.predicted.bulks.experiments,
            cor.predicted.bulks.experiments.harp,
            cor.predicted.bulks.harp.harp
        ),
        reference = c(
            rep("Harp reference", length(cor.predicted.bulks)),
            rep(reference_name, length(cor.predicted.bulks.experiments)),
            rep(reference_name, length(cor.predicted.bulks.experiments.harp)),
            rep("Harp reference", length(cor.predicted.bulks.harp.harp))
        ),
        method = c(
            rep("Flow cytometry", length(cor.predicted.bulks)),
            rep("Flow cytometry", length(cor.predicted.bulks.experiments)),
            rep("Harp cell proportions", length(cor.predicted.bulks.experiments.harp)),
            rep("Harp cell proportions", length(cor.predicted.bulks.harp.harp))
        )
    )

    # Custom theme
    custom_theme <- theme_minimal() +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(color = "black", fill = NA),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_blank(), # Remove x-axis text
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 12),
            strip.background = element_blank(), # Remove strip background
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.text.y.right = element_text(angle = 270, margin = margin(l = 5)),
            # Remove borders from facet labels
            strip.background.x = element_blank(),
            strip.background.y = element_blank(),
            strip.text.y = element_text(margin = margin(r = 10))
        )

    # Create the faceted plot
    plot_box <- ggplot(df, aes(x = method, y = value, fill = reference)) +
        geom_boxplot(aes(x = 0), width = 0.5, alpha = 1, outlier.shape = NA) + # Set x = 0 and control width
        facet_grid(reference ~ method) +
        custom_theme +
        labs(
            title = expression(bold(rho)),
            y = "Score",
            x = ""
        ) +
        scale_fill_manual(values = c("#52854C", "#D16103")) +
        scale_x_continuous(limits = c(-1, 1)) +
        theme(
            legend.position = "none",
            strip.placement = "outside",
            strip.text.y = element_text(angle = 0),
            strip.text.x = element_text(angle = 0)
        )

    # Save plot
    output_path <- file.path(output_dir, paste0("all_box_predicted_bulks_harp_experiment_", data_name, ".", file_type))
    ggsave(output_path, plot = plot_box, device = file_type, width = width, height = height)

    message("Box plot saved to: ", output_path)
}

# Function to create correlation plots
#' @export
# create_correlation_scatter_plots <- function(
#     estimated_c, true_proportions,
#     output_dir,
#     data_name,
#     file_type = "png",
#     width = 9,
#     height = 9) {
#     # Verify matrices have same dimensions
#     if (!all(dim(estimated_c) == dim(true_proportions))) {
#         stop("Matrices must have the same dimensions")
#     }

#     # Create list to store plots
#     plot_list <- list()

#     # Custom color palette
#     custom.col <- c(
#         "#52854C",
#         "#FFDB6D",
#         "#D16103",
#         "#293352",
#         "#4E84C4",
#         "#882255",
#         "#CC6677"
#     )

#     # Ensure we have enough colors for all cell types
#     if (nrow(estimated_c) > length(custom.col)) {
#         warning("More cell types than available colors. Colors will be recycled.")
#         colors <- rep(custom.col, ceiling(nrow(estimated_c) / length(custom.col)))
#     } else {
#         colors <- custom.col[1:nrow(estimated_c)]
#     }

#     # Create plots for each cell type
#     for (i in 1:nrow(estimated_c)) {
#         # Extract data for current cell type
#         cell_type <- rownames(estimated_c)[i]

#         # Rename "extra" to "Unidentified"
#         if (cell_type == "extra") {
#             cell_type <- "Unidentified"
#         }

#         estimated_vals <- estimated_c[i, ]
#         true_vals <- true_proportions[i, ]

#         # Calculate Pearson correlation
#         correlation <- cor(estimated_vals, true_vals, method = "pearson")

#         # Create data frame for plotting
#         plot_data <- data.frame(
#             Estimated = estimated_vals,
#             True = true_vals
#         )

#         # Create individual plot
#         p <- ggplot(plot_data, aes(x = True, y = Estimated)) +
#             geom_point(color = colors[i], alpha = 1) +
#             geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#             labs(
#                 title = gsub("_", " ", cell_type),
#                 x = "True Proportions",
#                 y = "Estimated Proportions"
#             ) +
#             annotate(
#                 "text",
#                 x = max(true_vals) * 0.2,
#                 y = max(estimated_vals) * 0.9,
#                 label = bquote(R[c] == .(round(correlation, 3))),
#                 size = 4
#             ) +
#             theme_bw() +
#             theme(
#                 plot.title = element_text(hjust = 0.5),
#                 aspect.ratio = 1
#             ) +
#             coord_cartesian(
#                 xlim = c(0, max(true_vals)),
#                 ylim = c(0, max(estimated_vals))
#             )

#         plot_list[[i]] <- p
#     }

#     # Calculate number of rows and columns for grid arrangement
#     n_plots <- length(plot_list)
#     n_cols <- ceiling(sqrt(n_plots))
#     n_rows <- ceiling(n_plots / n_cols)

#     # Arrange plots in a grid
#     arranged_plots <- do.call(
#         grid.arrange,
#         c(plot_list,
#             ncol = n_cols,
#             nrow = n_rows
#         )
#     )

#     ggsave(file.path(output_dir, paste0("scatter_plots_harp_", data_name, ".", file_type)),
#         plot = arranged_plots, bg = "white", width = width, height = height
#     )
#     arranged_plots
# }
create_correlation_scatter_plots <- function(
    estimated_c, true_proportions,
    output_dir,
    data_name,
    file_type = "png",
    width = 9,
    height = 9,
    map_names = TRUE) {
    # Verify matrices have same dimensions
    if (!all(dim(estimated_c) == dim(true_proportions))) {
        stop("Matrices must have the same dimensions")
    }

    # Filter out rows with "extra" or "Unidentified" cells
    keep_rows <- !grepl("extra", rownames(estimated_c), ignore.case = TRUE)
    estimated_c <- estimated_c[keep_rows, ]
    true_proportions <- true_proportions[keep_rows, ]

    # Create list to store plots
    plot_list <- list()
    map_celltype_names <- function(X) {
        if (stringr::str_detect(X, "CD4")) {
            return("CD4 T cell")
        }
        if (stringr::str_detect(X, "CD8")) {
            return("CD8 T cell")
        }
        if (stringr::str_detect(X, "NKT")) {
            return("NKT cells")
        }
        if (stringr::str_detect(X, "memory_B_cell")) {
            return(" Memory B cell ")
        }
        if (stringr::str_detect(X, "naive_B_cell")) {
            return(" Naive B cell ")
        }
        if (stringr::str_detect(X, "B")) {
            return("B cell")
        }
        if (stringr::str_detect(X, "Mono|mono")) {
            return("Monocytes")
        }
        if (stringr::str_detect(X, "NK")) {
            return("NK cells")
        } else {
            return(gsub("_", " ", X))
        }
    }

    # Custom color palette
    custom.col <- c(
        "#52854C",
        "#FFDB6D",
        "#D16103",
        "#293352",
        "#4E84C4",
        "#882255",
        "#CC6677"
    )

    # Ensure we have enough colors for all cell types
    if (nrow(estimated_c) > length(custom.col)) {
        warning("More cell types than available colors. Colors will be recycled.")
        colors <- rep(custom.col, ceiling(nrow(estimated_c) / length(custom.col)))
    } else {
        colors <- custom.col[1:nrow(estimated_c)]
    }

    # Create plots for each cell type
    for (i in 1:nrow(estimated_c)) {
        # Extract data for current cell type
        cell_type <- rownames(estimated_c)[i]

        # Apply cell type mapping if enabled
        if (map_names) {
            cell_type <- map_celltype_names(cell_type)
        }

        estimated_vals <- estimated_c[i, ]
        true_vals <- true_proportions[i, ]

        # Calculate Pearson correlation
        correlation <- cor(estimated_vals, true_vals, method = "pearson")

        # Create data frame for plotting
        plot_data <- data.frame(
            Estimated = estimated_vals,
            True = true_vals
        )

        # Create individual plot
        p <- ggplot(plot_data, aes(x = True, y = Estimated)) +
            geom_point(color = colors[i], alpha = 1, size = 3) + # Increased point size
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            labs(
                title = cell_type,
                x = "True Proportions",
                y = "Estimated Proportions"
            ) +
            annotate(
                "text",
                x = max(true_vals) * 0.2,
                y = max(estimated_vals) * 0.9,
                label = bquote(R[c] == .(round(correlation, 3))),
                size = 5 # Increased correlation text size
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5, size = 20), # Increased title size
                axis.title = element_text(size = 18), # Increased axis titles size
                axis.text = element_text(size = 16), # Increased axis text size
                aspect.ratio = 1
            ) +
            coord_cartesian(
                xlim = c(0, max(true_vals)),
                ylim = c(0, max(estimated_vals))
            )

        plot_list[[i]] <- p
    }

    # Calculate number of rows and columns for grid arrangement
    n_plots <- length(plot_list)
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)

    # Arrange plots in a grid
    arranged_plots <- do.call(
        grid.arrange,
        c(plot_list,
            ncol = n_cols,
            nrow = n_rows
        )
    )

    ggsave(file.path(output_dir, paste0("scatter_plots_harp_", data_name, ".", file_type)),
        plot = arranged_plots, bg = "white", width = width, height = height
    )
    arranged_plots
}

#' @export
bulk_expression_correlations_data_frame <- function(harp_model,
                                                    true_bulk_expression,
                                                    bulk_pheno_test,
                                                    ref_xsc) {
    harp_ref <- harp_model$reference_profiles$estimated_reference_second
    harp_c <- harp_model$estimated_c$estimated_c_second
    celltype <- colnames(harp_ref)
    gene_names <- intersect(rownames(harp_ref), rownames(true_bulk_expression))
    samples <- intersect(colnames(harp_c), colnames(true_bulk_expression))
    true_bulk_expression <- true_bulk_expression[gene_names, samples]

    # predicted bulks from harp
    predicted_bulk_expression_harp <- harp_ref[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_harp_dtdc <- harp_ref[, celltype] %*% harp_c[celltype, samples]

    # cell counts from experiment
    predicted_bulk_expression_experiment <- ref_xsc[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_experiment_dtdc <- ref_xsc[gene_names, celltype] %*% harp_c[celltype, samples]

    # Calculate correlations
    cor.predicted.bulks <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.harp.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp_dtdc,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.experiments <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.experiments.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment_dtdc,
        true_bulk_expression = true_bulk_expression
    )

    # Create data frame with correlation statistics
    correlation_stats <- data.frame(
        Method = c("HARP_Experiment", "HARP-HARP", "Experiment-Experiment", "Experiment-HARP"),
        Mean = c(
            mean(cor.predicted.bulks),
            mean(cor.predicted.bulks.harp.harp),
            mean(cor.predicted.bulks.experiments),
            mean(cor.predicted.bulks.experiments.harp)
        ),
        SD = c(
            sd(cor.predicted.bulks),
            sd(cor.predicted.bulks.harp.harp),
            sd(cor.predicted.bulks.experiments),
            sd(cor.predicted.bulks.experiments.harp)
        ),
        N = c(
            length(cor.predicted.bulks),
            length(cor.predicted.bulks.harp.harp),
            length(cor.predicted.bulks.experiments),
            length(cor.predicted.bulks.experiments.harp)
        )
    )
    correlation_stats[, c("Mean", "SD")] <-
        round(correlation_stats[, c("Mean", "SD")], 4)

    return(correlation_stats)
}
#' @export
create_correlation_scatter_plots_new <- function(
    estimated_c, true_proportions,
    output_dir,
    data_name,
    file_type = "png",
    width = 9,
    height = 9,
    map_names = TRUE) {
    # Verify matrices have same dimensions
    if (!all(dim(estimated_c) == dim(true_proportions))) {
        stop("Matrices must have the same dimensions")
    }

    # Create list to store plots
    plot_list <- list()
    map_celltype_names <- function(X) {
        if (stringr::str_detect(X, "CD4")) {
            return("CD4 T cell")
        }
        if (stringr::str_detect(X, "CD8")) {
            return("CD8 T cell")
        }
        if (stringr::str_detect(X, "NKT")) {
            return("NKT cells")
        }
        if (stringr::str_detect(X, "memory_B_cell")) {
            return(" Memory B cell ")
        }
        if (stringr::str_detect(X, "naive_B_cell")) {
            return(" Naive B cell ")
        }
        if (stringr::str_detect(X, "B")) {
            return("B cell")
        }
        if (stringr::str_detect(X, "Mono|mono")) {
            return("Monocytes")
        }
        if (stringr::str_detect(X, "NK")) {
            return("NK cells")
        }
        if (stringr::str_detect(X, "extra")) {
            return("Unidentified")
        } else {
            return(gsub("_", " ", X))
        }
    }

    # Custom color palette
    custom.col <- c(
        "#52854C",
        "#FFDB6D",
        "#D16103",
        "#293352",
        "#4E84C4",
        "#882255",
        "#CC6677"
    )

    # Ensure we have enough colors for all cell types
    if (nrow(estimated_c) > length(custom.col)) {
        warning("More cell types than available colors. Colors will be recycled.")
        colors <- rep(custom.col, ceiling(nrow(estimated_c) / length(custom.col)))
    } else {
        colors <- custom.col[1:nrow(estimated_c)]
    }

    # Function to calculate RMSD
    calculate_rmsd <- function(true, estimated) {
        sqrt(mean((true - estimated)^2))
    }

    # Function to calculate MAD
    calculate_mad <- function(true, estimated) {
        mean(abs(true - estimated))
    }

    # Create plots for each cell type
    for (i in 1:nrow(estimated_c)) {
        # Extract data for current cell type
        cell_type <- rownames(estimated_c)[i]

        # Apply cell type mapping if enabled
        if (map_names) {
            cell_type <- map_celltype_names(cell_type)
        }

        estimated_vals <- estimated_c[i, ]
        true_vals <- true_proportions[i, ]

        # Calculate metrics
        correlation <- cor(estimated_vals, true_vals, method = "pearson")
        rmsd <- calculate_rmsd(true_vals, estimated_vals)
        mad <- calculate_mad(true_vals, estimated_vals)

        # Create data frame for plotting
        plot_data <- data.frame(
            Estimated = estimated_vals,
            True = true_vals
        )

        # Create individual plot
        p <- ggplot(plot_data, aes(x = True, y = Estimated)) +
            geom_point(color = colors[i], alpha = 1) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            labs(
                title = cell_type,
                x = "True Proportions",
                y = "Estimated Proportions"
            ) +
            annotate(
                "text",
                x = max(true_vals) * 0.2,
                y = max(estimated_vals) * 0.95,
                label = bquote(R[c] == .(round(correlation, 3))),
                size = 4
            ) +
            annotate(
                "text",
                x = max(true_vals) * 0.2,
                y = max(estimated_vals) * 0.85,
                label = bquote(RMSD == .(round(rmsd, 3))),
                size = 4
            ) +
            annotate(
                "text",
                x = max(true_vals) * 0.2,
                y = max(estimated_vals) * 0.75,
                label = bquote(MAD == .(round(mad, 3))),
                size = 4
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5),
                aspect.ratio = 1
            ) +
            coord_cartesian(
                xlim = c(0, max(true_vals)),
                ylim = c(0, max(estimated_vals))
            )

        plot_list[[i]] <- p
    }

    # Calculate number of rows and columns for grid arrangement
    n_plots <- length(plot_list)
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)

    # Arrange plots in a grid
    arranged_plots <- do.call(
        grid.arrange,
        c(plot_list,
            ncol = n_cols,
            nrow = n_rows
        )
    )

    ggsave(file.path(output_dir, paste0("scatter_plots_harp_", data_name, ".", file_type)),
        plot = arranged_plots, bg = "white", width = width, height = height
    )
    arranged_plots
}

bulk_expression_correlations_data_frame_first <- function(harp_model,
                                                          true_bulk_expression,
                                                          bulk_pheno_test,
                                                          ref_xsc) {
    harp_ref <- harp_model$reference_profiles$estimated_reference_first
    harp_c <- harp_model$estimated_c$estimated_c_first
    celltype <- colnames(harp_ref)
    gene_names <- intersect(rownames(harp_ref), rownames(true_bulk_expression))
    samples <- intersect(colnames(harp_c), colnames(true_bulk_expression))
    true_bulk_expression <- true_bulk_expression[gene_names, samples]

    # predicted bulks from harp
    predicted_bulk_expression_harp <- harp_ref[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_harp_dtdc <- harp_ref[, celltype] %*% harp_c[celltype, samples]

    # cell counts from experiment
    predicted_bulk_expression_experiment <- ref_xsc[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]
    predicted_bulk_expression_experiment_dtdc <- ref_xsc[gene_names, celltype] %*% harp_c[celltype, samples]

    # Calculate correlations
    cor.predicted.bulks <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.harp.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_harp_dtdc,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.experiments <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment,
        true_bulk_expression = true_bulk_expression
    )

    cor.predicted.bulks.experiments.harp <- calculate_bulk_correlation(
        predicted_bulk_expression = predicted_bulk_expression_experiment_dtdc,
        true_bulk_expression = true_bulk_expression
    )

    # Create data frame with correlation statistics
    correlation_stats <- data.frame(
        Method = c("HARP_Experiment", "HARP-HARP", "Experiment-Experiment", "Experiment-HARP"),
        Mean = c(
            mean(cor.predicted.bulks),
            mean(cor.predicted.bulks.harp.harp),
            mean(cor.predicted.bulks.experiments),
            mean(cor.predicted.bulks.experiments.harp)
        ),
        SD = c(
            sd(cor.predicted.bulks),
            sd(cor.predicted.bulks.harp.harp),
            sd(cor.predicted.bulks.experiments),
            sd(cor.predicted.bulks.experiments.harp)
        ),
        N = c(
            length(cor.predicted.bulks),
            length(cor.predicted.bulks.harp.harp),
            length(cor.predicted.bulks.experiments),
            length(cor.predicted.bulks.experiments.harp)
        )
    )
    correlation_stats[, c("Mean", "SD")] <-
        round(correlation_stats[, c("Mean", "SD")], 4)

    return(correlation_stats)
}


# this is to recounstruct bulk samples from real data and harp
#' @export
predict_bulk_expression <- function(harp_output, data) {
    # Extract required components from inputs
    harp_ref <- harp_output$reference_profiles$estimated_reference_second
    harp_c <- harp_output$estimated_c$estimated_c_second

    # Extract cell types and gene names
    celltype <- colnames(harp_ref)
    gene_names <- rownames(harp_ref)

    # Get test samples from data
    bulk_pheno_test <- data$bulk_pheno_test_true
    samples <- colnames(bulk_pheno_test)

    # Calculate true bulk expression
    true_bulk_expression <- normalize_to_count(as.matrix(data$bulk_counts_test[gene_names, samples]))

    # Calculate predicted bulk expression using HARP reference profiles
    predicted_bulk_expression_harp <- harp_ref %*% bulk_pheno_test[celltype, samples]

    # Compute experimental reference profiles
    ref_xsc <- compute_reference_harp(sc_library = data$sc_library)

    # Calculate predicted bulk expression using experimental reference profiles
    predicted_bulk_expression_experiment <- ref_xsc[gene_names, celltype] %*% bulk_pheno_test[celltype, samples]

    # Return the results as a list
    return(list(
        true_bulk_expression = true_bulk_expression,
        predicted_bulk_expression_harp = predicted_bulk_expression_harp,
        predicted_bulk_expression_experiment = predicted_bulk_expression_experiment
    ))
}
