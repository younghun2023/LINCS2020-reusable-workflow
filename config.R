# Portable configuration. Run scripts from this package folder.
project_dir <- getwd()
data_dir <- file.path(project_dir,"data","LINCS2020")
run_root <- file.path(project_dir,"results")
# Required CSV: symbol,direction, with direction Up or Down.
# Optional gene_id column supplies user-verified Entrez IDs for legacy symbols.
query_file <- file.path(project_dir,"inputs","gene_signature.csv")
# Optional CSV: Drug_name,Network_Z,Network_adjusted_P.
# These must be independently computed for YOUR disease/gene set.
network_file <- ""
download_data <- FALSE
download_only <- FALSE
install_missing_packages <- FALSE
require_qc_pass <- TRUE
require_is_hiq <- FALSE
chunk_size <- 250L
minimum_signatures <- 5L
rank_percentile_cutoff <- 5
network_z_cutoff <- -0.15
