# A helper script to combine proportions tibbles from several algorithms

library(harplication)
library(tidyverse)

config <- config::get(file = system.file("R/inst/config", "config.yml", package = "harplication"), use_parent = FALSE)
if (config$standalone) {
    source(system.file("scripts/simulation", "setup_parameters.R", package = "harplication"))
}

# note CIBERSORT proportions are added in evaluate_proportions.R
algorithms <- c("harp_true",
                "harp_cdeath",
                "harp_noise",
                "bayesPrism",
                "bayesPrism_X1",
                "bayesPrism_harp",
                "music")
proportions_all <- tibble()
for (irun in 1:number_simulation_runs) {
    proportions_true <- readRDS(file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))
    proportions_all <- rbind(proportions_all, proportions_true)
}
proportions_all <- proportions_all %>% mutate(bulk_id = as.character(bulk_id))
for (alg in algorithms) {
    for (irun in 1:number_simulation_runs) {
        proportions_alg <- readRDS(file.path(output_dir, paste0("proportions_", alg, "_run_", irun , ".rds")))
        proportions_all <- rbind(proportions_all, proportions_alg)
    }
    print(paste("Appended", alg))
}
saveRDS(proportions_all, "/home/phuettl/zahra/harplication/data/generated/roider_steen_dist/proportions_all.rds")
