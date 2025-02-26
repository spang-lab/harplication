# DEPENDENCIES
library(harplication)
# set parameters
source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(extrafont)
loadfonts()


sc <- readRDS(single_cell_file)

# preprocessing for Steen
sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Steen")
sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
print(paste("Number of cells in Steen", ncol(sc$sc.counts)))
scpheno <- as_tibble(sc$sc.pheno) %>%
    dplyr::select(celltype, colnames, Patient) %>%
    dplyr::rename(cell_id = colnames)

sccounts <- tibble(cell_id = colnames(sc$sc.counts)) %>% cbind(as_tibble(t(sc$sc.counts)))
sclong <- sccounts %>%
    pivot_longer(!cell_id, names_to = "gene", values_to = "expression") %>%
    inner_join(scpheno)

# getting cell proportions
cell_distribution_patients <- patient_overview(sclong)
ggsave(
    file.path(
        output_dir,
        paste0("cell_distribution_patients_Steen.pdf")
    ),
    tableGrob(cell_distribution_patients),
    width = 10, height = 25
)

pie <- ggplot(cell_distribution_patients, aes(x = "", y = pct, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(~Patient, ncol = 4) +
    theme_void() +
    labs(
        fill = "Celltype"
    ) + theme(text = element_text(size = 20, family = "DejaVu Sans"))
ggsave(file.path(output_dir, "pie_steen.pdf"), width = 8, height = 5, device = cairo_pdf)

sc <- readRDS(single_cell_file)

# preprocessing for Roider
sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Roider")
sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
print(paste("Number of cells in Roider", ncol(sc$sc.counts)))
scpheno <- as_tibble(sc$sc.pheno) %>%
    dplyr::select(celltype, colnames, origin, Patient, Sample) %>%
    mutate(Patient = ifelse(is.na(Patient), Sample, Patient)) %>%
    dplyr::select(celltype, colnames, origin, Patient) %>%
    dplyr::rename(cell_id = colnames)

sccounts <- tibble(cell_id = colnames(sc$sc.counts)) %>% cbind(as_tibble(t(sc$sc.counts)))
sclong <- sccounts %>%
    pivot_longer(!cell_id, names_to = "gene", values_to = "expression") %>%
    inner_join(scpheno)

# getting cell proportions
cell_distribution_patients <- patient_overview(sclong)
ggsave(
    file.path(
        output_dir,
        paste0("cell_distribution_patients_Roider.pdf")
    ),
    tableGrob(cell_distribution_patients),
    width = 10, height = 25
)

pie <- ggplot(cell_distribution_patients, aes(x = "", y = pct, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(~Patient, ncol = 4) +
    theme_void() +
    labs(
        fill = "Celltype"
    ) + theme(text = element_text(size = 20, family = "DejaVu Sans"))
ggsave(file.path(output_dir, "pie_roider.pdf"), last_plot(), width = 8, height = 5, device = cairo_pdf)
