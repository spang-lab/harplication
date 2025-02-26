#' prepares steen data for simulation
#' @export
preprocess_sc_bulk <- function(single_cell_file, sc_library_patients, bulk_patients, simulation_name) {
    # Use sc file to prepare for simulating bulks and consituting the reference matrix
    if (!file.exists(single_cell_file)) {
        stop(paste0(
            "Single cell file is not available in ",
            single_cell_file,
            ". Please download from https://doi.org/10.5281/zenodo.10568550."
        ))
    }
    sc <- readRDS(single_cell_file)
    if (grepl("roider", simulation_name, fixed = TRUE)) {
        if (!grepl("roider_steen", simulation_name, fixed = TRUE)) {
            # in the pure Roider case only use Roider
            sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Roider")
            sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
            scpheno <- as_tibble(sc$sc.pheno) %>%
                dplyr::select(celltype, colnames, origin, Patient, Sample) %>%
                mutate(Patient = ifelse(is.na(Patient), Sample, Patient)) %>%
                dplyr::select(celltype, colnames, origin, Patient) %>%
                dplyr::rename(cell_id = colnames)
        } else {
            # In the Roider Steen case map malignant
            # lymphoma to B cells in ordr to be compatible with steen
            scpheno <- as_tibble(sc$sc.pheno) %>%
                dplyr::select(celltype, colnames, origin, Patient, Sample) %>%
                mutate(Patient = ifelse(is.na(Patient), Sample, Patient)) %>%
                dplyr::select(celltype, colnames, origin, Patient) %>%
                dplyr::rename(cell_id = colnames) %>%
                mutate(celltype = ifelse(celltype == "malignant_lymphoma", "B", celltype))
        }
    } else {
        # pure Steen case
        sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Steen")
        sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
        scpheno <- as_tibble(sc$sc.pheno) %>%
            dplyr::select(celltype, colnames, Patient) %>%
            dplyr::rename(cell_id = colnames)
    }
    sccounts <- tibble(cell_id = colnames(sc$sc.counts)) %>% cbind(as_tibble(t(sc$sc.counts)))

    # split by patients into a dataset for bulk creation and a single cell library:
    patients <- scpheno %>%
        pull(Patient) %>%
        unique()
    if (!grepl("roider_steen", simulation_name, fixed = TRUE)) {
            scpheno <- scpheno %>%
                filter(Patient %in% c(sc_library_patients, bulk_patients)) %>%
                mutate(testtrain = ifelse(Patient %in% sc_library_patients, "train", "test")) %>%
                mutate(testtrain = as_factor(testtrain))
            npatients_sc <- length(sc_library_patients)
    } else {
        scpheno <- scpheno %>%
            mutate(testtrain = ifelse(origin == "Steen", "train", "test")) %>%
            mutate(testtrain = as_factor(testtrain))
        npatients_sc <- scpheno %>% filter(testtrain == "train") %>% pull(Patient) %>% unique() %>% length()
    }
    sclong <- sccounts %>%
        pivot_longer(!cell_id, names_to = "gene", values_to = "expression") %>%
        inner_join(scpheno)
    return(sclong)
}


#' @export
filter_sc_bulk_common_celltypes <- function(sclong) {
    pc <- get_patient_celltype(sclong)
    celltypes_sc_bulk <- pc$celltypes_sc_bulk
    celltypes_sc <- celltypes_sc_bulk %>% filter(testtrain == "train") %>% pull(celltype)
    # Here we determine the celltypes that are expressed in at least one patient, celltypes that are not expressed
    # in any of the patients are discarded from the simulation
    celltypes_bulks <- celltypes_sc_bulk %>% filter(testtrain == "test") %>% pull(celltype)
    # both in sc_library as well as in bulks keep the intersect of celltypes
    celltypes_isec <- intersect(celltypes_sc, celltypes_bulks)

    print("The downstream analyzed celltypes are: ")
    print(celltypes_isec)
    print("The single cell library will be composed of the patients:")
    print(pc$patient_celltypes %>% filter(testtrain == "train") %>% pull(Patient) %>% unique())
    print("The bulks will be simulated from the patients:")
    print(pc$patient_celltypes %>%
        filter(testtrain == "test") %>%
        pull(Patient) %>% unique())
    sclong <- sclong %>% filter(celltype %in% celltypes_isec)
    return(sclong)
}

#' @export
get_patient_celltype <- function(sclong) {
    patient_celltypes <- sclong %>%
        dplyr::select(Patient, celltype, testtrain) %>%
        dplyr::distinct()
    celltypes_sc_bulk <- patient_celltypes %>%
        dplyr::select(celltype, testtrain) %>%
        dplyr::distinct()
    return(list(patient_celltypes = patient_celltypes, celltypes_sc_bulk = celltypes_sc_bulk))
}

#' @export
filter_patient_common_celltypes_steen_roider <- function(sclong) {
    origins <- c("Roider", "Steen")
    celltypes_isec <- c()
    count <- 0
    for (og in origins) {
        celltypes_origin <- sclong %>%
            filter(origin == og) %>%
            pull(celltype) %>%
            unique()
        if (count > 0) {
            celltypes_isec <- intersect(celltypes_origin, celltypes_isec)
        } else {
            celltypes_isec <- celltypes_origin
        }
        count <- count + 1
    }
    sclong <- sclong %>% filter(celltype %in% celltypes_isec)
    return(sclong)
}