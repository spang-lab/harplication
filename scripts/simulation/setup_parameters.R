config <- config::get(file = system.file("R/inst/config", "config.yml", package = "harplication"), use_parent = FALSE)
# PARAMETERS
base_dir <- config$base_dir
simulation_name <- config$simulation_name
output_dir <- file.path(base_dir, simulation_name)
dir.create(file.path(output_dir))
single_cell_file <- config$single_cell_file
mapping_file <- config$mapping_file
amount_train_bulks <- config$amount_train_bulks
amount_test_bulks <- config$amount_test_bulks
subset_bulks <- config$subset_bulks
sc_library_patients <- config$sc_library_patients
bulk_patients <- config$bulk_patients
single_cell_bulk_split <- config$single_cell_bulk_split
bulks_simulation_train_split <- config$bulks_simulation_train_split
train_sc_patients <- config$train_sc_patients
test_sc_patients <- config$test_sc_patients
bulk_train_split <- config$bulk_train_split
kill_cells_mean <- config$kill_cells_mean
kill_cells_std <- config$kill_cells_std
lambda_seq <- config$lambda_seq
algorithms <- config$algorithms
percentage_genes_to_distort <- config$percentage_genes_to_distort
distortion_factor_mean <- config$distortion_factor_mean
distortion_factor_std <- config$distortion_factor_std
number_simulation_runs <- config$number_simulation_runs
proportions_filename <- config$proportions_filename
cibersort_filename <- config$cibersort_filename
cibersort_algorithms <- config$cibersort_algorithms
commit_resluts <- config$commit_results