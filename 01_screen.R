# Generic bidirectional human-gene LINCS2020 screen. See README.md.
source("config.R")
dir.create(data_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(run_root,recursive=TRUE,showWarnings=FALSE)
# DEPENDENCIES
# Packages are checked before analysis. Installation is disabled by default.
# Run 00_install_packages.R explicitly if required. Package versions are saved
# in session information; they are not pinned by this package.
required_packages <- c("cmapR","dplyr","readxl","openxlsx","ggplot2","ggrepel")
missing <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing) && install_missing_packages) {
  cran <- setdiff(missing,"cmapR")
  if(length(cran)) install.packages(cran,repos="https://cloud.r-project.org")
  if("cmapR" %in% missing) {
    if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
    BiocManager::install("cmapR",ask=FALSE,update=FALSE)
  }
}
missing <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing packages: ",paste(missing,collapse=", "), ". Run 00_install_packages.R if needed.")
library(dplyr)
library(ggplot2)
if(!download_only) for(f in c(query_file)) {
  if(!file.exists(f)) stop("Check input path before downloading: ",f)
}

# 1B. Download release-matched data; retain partial downloads ------------------
base_url <- "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/"
remote_names <- c("level5/level5_beta_trt_cp_n720216x12328.gctx",
                  "siginfo_beta.txt","geneinfo_beta.txt","compoundinfo_beta.txt",
                  "cellinfo_beta.txt","README.txt")
curl_bin <- Sys.which(if(.Platform$OS.type=="windows") "curl.exe" else "curl")
if(download_data && !nzchar(curl_bin)) stop("curl command is needed for resumable downloads.")
# REFERENCE FILES
# With download_data=FALSE, nonempty cached files are required and used locally.
# With download_data=TRUE, remote size/ETag metadata are checked, completed files
# are reused, and interrupted .part downloads are resumed. A changed remote
# fingerprint triggers a stop to avoid silently mixing reference versions.
# ETags can identify a release object but multipart ETags are not MD5 hashes.
download_one <- function(remote) {
  dest <- file.path(data_dir,basename(remote))
  url <- paste0(base_url,remote)
  if(!download_data) {
    if(!file.exists(dest) || file.info(dest)$size==0) stop("Missing cached file: ",dest)
    return(data.frame(URL=url,File=basename(dest),Bytes=file.info(dest)$size,ETag=NA_character_))
  }
  hdr <- tempfile()
  status <- system2(curl_bin,c("--fail","--silent","--show-error","--location","--head",
    "--retry","3","--connect-timeout","30","--max-time","120",
    "--output",shQuote(hdr),shQuote(url)))
  if(status!=0) stop("Cannot verify download URL: ",url)
  h <- readLines(hdr,warn=FALSE)
  get_header <- function(key) {
    z <- grep(paste0("^",key,":"),h,value=TRUE,ignore.case=TRUE)
    if(!length(z)) return(NA_character_)
    trimws(sub("^[^:]+:","",tail(z,1)))
  }
  size <- as.numeric(get_header("Content-Length")); etag <- get_header("ETag")
  if(!is.finite(size) || size<=0) stop("Missing remote file length: ",url)
  fingerprint <- list(url=url,bytes=size,etag=etag)
  stamp <- paste0(dest,".remote.rds")
  if(file.exists(stamp) && !identical(readRDS(stamp),fingerprint))
    stop("Remote release file changed. Use a new data_dir to preserve the previous release: ",dest)
  if(file.exists(dest)) {
    if(file.info(dest)$size!=size) stop("Existing file size mismatch; preserve it and use a new data_dir: ",dest)
  } else {
    part <- paste0(dest,".part")
    if(file.exists(part) && !file.exists(stamp)) stop("Unverified partial file: ",part)
    saveRDS(fingerprint,stamp)
    if(!file.exists(part) || file.info(part)$size!=size) {
      message("Downloading/resuming ",basename(dest)," (",round(size/1e9,3)," GB)")
      status <- system2(curl_bin,c("--fail","--location","--retry","5","--retry-delay","10",
        "--connect-timeout","30","--continue-at","-","--output",shQuote(part),shQuote(url)))
      if(status!=0) stop("Download interrupted; rerun to resume: ",part)
    }
    if(file.info(part)$size!=size) stop("Incomplete download: ",part)
    if(!file.rename(part,dest)) stop("Could not finalize download: ",part)
  }
  saveRDS(fingerprint,stamp)
  data.frame(URL=url,File=basename(dest),Bytes=size,ETag=etag)
}
manifest <- bind_rows(lapply(remote_names,download_one))
write.csv(manifest,file.path(data_dir,"download_manifest.csv"),row.names=FALSE)
# Length checks and ETags establish release identity; multipart ETags are not MD5 checksums.
if(download_only) stop("Downloads complete. Set download_only <- FALSE and rerun for analysis.")
matrix_file <- file.path(data_dir,basename(remote_names[1]))
sig_file <- file.path(data_dir,"siginfo_beta.txt")
gene_file <- file.path(data_dir,"geneinfo_beta.txt")
compound_file <- file.path(data_dir,"compoundinfo_beta.txt")
# Input/settings hash separates corrected reruns and makes checkpoints reusable.
# RUN IDENTIFICATION AND RESUMPTION
# A hash of inputs/settings defines a separate output folder for each run.
# Matching inputs/settings permit reuse of completed score chunks. Moving files
# or changing their timestamps may produce a new run folder. Do not rename an
# old checkpoint directory to force its reuse with a different gene query.
# Input hashes support provenance; the large matrix is recorded by file metadata.
# Generated manifests and RDS files may contain the executing user's paths;
# inspect these separately before sharing outputs.
run_spec <- list(version="LINCS2020_generic_v1.0",inputs=tools::md5sum(c(query_file,sig_file,gene_file,compound_file)),
                 matrix_bytes=file.info(matrix_file)$size,matrix_mtime=file.info(matrix_file)$mtime,
                 qc=require_qc_pass,hiq=require_is_hiq,chunk=chunk_size,min_n=minimum_signatures,
                 rank_cutoff=rank_percentile_cutoff,z_cutoff=network_z_cutoff,network_input=if(nzchar(network_file)) tools::md5sum(network_file) else NULL)
spec_file <- tempfile(); saveRDS(run_spec,spec_file)
run_id <- substr(unname(tools::md5sum(spec_file)),1,12)
output_dir <- file.path(run_root,paste0("run_",run_id))
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
saveRDS(run_spec,file.path(output_dir,"run_specification.rds"))
writeLines(normalizePath(output_dir,winslash="/"),file.path(run_root,"LATEST_RUN.txt"))
write.csv(manifest,file.path(output_dir,"download_manifest.csv"),row.names=FALSE)
message("Run outputs: ",output_dir)

write_output <- function(x, name) write.csv(x, file.path(output_dir, name), row.names=FALSE, na="NA", fileEncoding="UTF-8")
require_columns <- function(x, cols, label) {
  absent <- setdiff(cols, names(x))
  if (length(absent)) stop(label, " missing: ", paste(absent, collapse=", "))
}
clean_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}
numeric_checked <- function(x, label) {
  v <- suppressWarnings(as.numeric(x))
  missing <- is.na(x) | trimws(as.character(x)) %in% c("", "NA", "NaN")
  if (any(!missing & !is.finite(v))) stop("Invalid numeric values: ", label)
  v[missing] <- NA_real_
  v
}
collapse_text <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse="; ")

# 2. Read the user's gene signature -----------------------------------------
# No predefined module, gene count, disease, tissue, or significance criterion.
# Define your signature upstream and document the biological contrast and filters.
# Up means higher expression in the condition to reverse; Down means lower.
query <- read.csv(query_file,colClasses="character",check.names=FALSE)
require_columns(query,c("symbol","direction"),"Query")
query$symbol <- trimws(query$symbol)
if(anyNA(query$symbol) || any(!nzchar(query$symbol)) || anyDuplicated(query$symbol))
  stop("Query symbols must be nonempty and unique.")
d <- tolower(trimws(query$direction))
if(anyNA(d) || any(!d %in% c("up","down"))) stop("Direction must be Up or Down.")
query$direction <- ifelse(d=="up","Up","Down")
if(!all(c("Up","Down") %in% query$direction)) stop("This implementation requires both directions.")
if(!"gene_id" %in% names(query)) query$gene_id <- NA_character_
query$gene_id <- trimws(query$gene_id)
query$gene_id[query$gene_id %in% c("","NA")] <- NA_character_
names(query)[names(query)=="gene_id"] <- "supplied_gene_id"
write_output(query,"query_before_mapping.csv")
# 3. Map to actual matrix IDs; retain mapping losses for review ---------------
# QUERY MAPPING
# Only query IDs present in the matrix enter scoring. Document mapping losses
# rather than replacing missing values or inventing expression measurements.
# Measured landmark, best-inferred and inferred features are all eligible here.
# The complete matrix is still the ranked enrichment background; it is not
# reduced to the query genes. Review the exported feature-space annotations.
genes <- read.delim(gene_file, colClasses="character", check.names=FALSE)
require_columns(genes, c("gene_id", "gene_symbol", "feature_space"), "Gene annotations")
map <- genes %>% transmute(symbol=trimws(gene_symbol), gene_id=trimws(gene_id),
                           feature_space=feature_space) %>% distinct()
map <- map %>% filter(!is.na(gene_id), nzchar(gene_id), !is.na(symbol), nzchar(symbol))
ambiguous <- map %>% filter(symbol %in% query$symbol[is.na(query$supplied_gene_id)]) %>% count(symbol) %>% filter(n>1)
if (nrow(ambiguous)) {
  write_output(ambiguous, "ambiguous_gene_mappings.csv")
  stop("Ambiguous gene mappings exported; resolve before screening.")
}
row_ids <- as.character(cmapR::read_gctx_ids(matrix_file, dim="row"))
column_ids <- as.character(cmapR::read_gctx_ids(matrix_file, dim="column"))
if (anyDuplicated(row_ids) || anyDuplicated(column_ids)) stop("Duplicate matrix identifiers.")
# User-supplied IDs take precedence; record both input symbol and mapped symbol.
# No study-specific alias mappings are embedded in this reusable workflow.
exact_query <- query %>% filter(is.na(supplied_gene_id))
provided_query <- query %>% filter(!is.na(supplied_gene_id))
mapped <- bind_rows(
  left_join(exact_query,map,by="symbol") %>%
    mutate(mapping_method=ifelse(is.na(gene_id),"Unmapped","Exact symbol")),
  provided_query %>% mutate(gene_id=supplied_gene_id,
                            mapping_method="User-supplied Entrez ID")
)
annotation_match <- match(mapped$gene_id,genes$gene_id)
mapped$LINCS_annotation_symbol <- genes$gene_symbol[annotation_match]
mapped$feature_space <- genes$feature_space[annotation_match]
mapped <- mapped %>%
  mutate(mapping_status=case_when(is.na(gene_id) ~ "Symbol not mapped",
                                 !gene_id %in% row_ids ~ "Absent from matrix",
                                 TRUE ~ "Included"))
write_output(mapped, "query_mapping_audit.csv")
included <- mapped %>% filter(mapping_status=="Included")
conflicts <- included %>% group_by(gene_id) %>% summarise(n=n_distinct(direction), .groups="drop") %>% filter(n>1)
if (nrow(conflicts)) stop("A matrix gene maps to opposite query directions.")
up_ids <- unique(included$gene_id[included$direction=="Up"])
down_ids <- unique(included$gene_id[included$direction=="Down"])
if (!length(up_ids) || !length(down_ids)) stop("Both query directions must map to the matrix.")
writeLines(up_ids, file.path(output_dir, "query_up_Entrez.txt"))
writeLines(down_ids, file.path(output_dir, "query_down_Entrez.txt"))
write_output(mapped %>% filter(mapping_status!="Included"), "query_excluded_genes.csv")
message("Matrix-supported query: ", length(up_ids), " up; ", length(down_ids), " down.")
if (any(mapped$mapping_status!="Included")) warning("Some genes were excluded. Review query_mapping_audit.csv before finalizing the analysis.")

# 4. Select small-molecule signatures and record coverage --------------------
meta <- read.delim(sig_file, colClasses="character", check.names=FALSE)
require_columns(meta,c("cmap_name","cell_iname","qc_pass","is_hiq"),"LINCS2020 metadata")
meta$pert_iname <- meta$cmap_name
meta$cell_id <- meta$cell_iname
# QUALITY-FLAG POLICY
# This reproduces the chosen exact-1 filter: 1/true passes; 0, 2 and 3 do not.
# Recognizing codes 2 and 3 prevents parser failure, but does not assign them a
# biological meaning. Their definitions must be checked against source metadata.
# is_hiq is optional and is not required by the default configuration.
# is_ncs_sig is not used as a P value or as evidence that our query is significant.
parse_qc <- function(x) {
  x <- tolower(trimws(as.character(x)))
  # Preserve the configured exact-1 rule. Codes 2 and 3 are observed
  # in this release; do not infer that nonzero codes all mean pass.
  allowed <- c("1","0","2","3","true","false","","na","-666")
  unexpected <- unique(x[!is.na(x) & !x %in% allowed])
  if(length(unexpected)) stop("Unexpected QC values: ",paste(unexpected,collapse=", "))
  !is.na(x) & x %in% c("1","true")
}
write_output(as.data.frame(table(meta$qc_pass,useNA="ifany")),"qc_pass_value_counts.csv")
write_output(as.data.frame(table(meta$is_hiq,useNA="ifany")),"is_hiq_value_counts.csv")
meta$passes_configured_qc <- (!require_qc_pass | parse_qc(meta$qc_pass)) &
                              (!require_is_hiq | parse_qc(meta$is_hiq))
require_columns(meta, c("sig_id", "pert_id", "pert_iname", "pert_type", "cell_id"), "Signature metadata")
if (anyDuplicated(meta$sig_id)) stop("Duplicate signature metadata IDs.")
meta$in_matrix <- meta$sig_id %in% column_ids
meta$valid_identity <- !is.na(meta$pert_id) & nzchar(meta$pert_id) &
                       !is.na(meta$pert_iname) & nzchar(meta$pert_iname)
selected <- meta %>% filter(pert_type=="trt_cp", in_matrix, valid_identity, passes_configured_qc)
if (!nrow(selected)) stop("No trt_cp signatures with matrix coverage.")
write_output(meta %>% filter(pert_type=="trt_cp", !in_matrix | !valid_identity | !passes_configured_qc), "excluded_signature_metadata.csv")

# 5. Custom unweighted enrichment (full matrix gene universe) ----------------
# This is the unweighted running-sum statistic used in the original script.
# No weighted-CMap implementation or reference-query calibration is implied.
# SCORE DEFINITION
# Within each reference signature, genes are ordered from highest to lowest.
# Query hits increase the running sum; non-hits decrease it. The signed largest
# absolute deviation is the unweighted enrichment score. Expression ties retain
# matrix row order. A constant or nonfinite profile is not scored.
# Negative CS indicates opposition to the up/down query. This is a local
# unweighted procedure, not an implementation of official CMap calibration.
unweighted_es <- function(ranked_ids, set_ids) {
  hits <- ranked_ids %in% set_ids
  nh <- sum(hits); n <- length(hits)
  if (nh==0L || nh==n) return(NA_real_)
  running <- cumsum(hits/nh) - cumsum((!hits)/(n-nh))
  running[which.max(abs(running))]
}
score_one <- function(values, ids) {
  if (any(!is.finite(values)) || max(values)==min(values)) return(c(NA_real_,NA_real_,NA_real_))
  ranked <- ids[order(values, decreasing=TRUE, method="radix")]
  eu <- unweighted_es(ranked, up_ids)
  ed <- unweighted_es(ranked, down_ids)
  cs <- if (!is.finite(eu) || !is.finite(ed)) NA_real_ else if (sign(eu)!=sign(ed)) (eu-ed)/2 else 0
  c(eu, ed, cs)
}

# 6. Score in chunks --------------------------------------------------------
# COMPUTATION AND CHECKPOINTS
# Each chunk reads all gene rows for a limited number of signatures.
# Reducing chunk_size can reduce memory use; it changes the run specification.
# Valid completed chunks are reused on rerun. Partial temporary checkpoints do
# not replace completed RDS files. Do not interrupt R while files are being moved.
# Processing the full resource can take many hours; no fixed runtime is assumed.
chunks <- split(selected$sig_id, ceiling(seq_len(nrow(selected))/chunk_size))
scores <- vector("list", length(chunks))
checkpoint_dir <- file.path(output_dir,"checkpoints")
dir.create(checkpoint_dir,showWarnings=FALSE)
for (i in seq_along(chunks)) {
  checkpoint <- file.path(checkpoint_dir,sprintf("chunk_%05d.rds",i))
  if(file.exists(checkpoint)) {
    cached <- readRDS(checkpoint)
    if(!setequal(cached$sig_id,chunks[[i]]) || anyDuplicated(cached$sig_id)) stop("Invalid checkpoint: ",checkpoint)
    scores[[i]] <- cached
    next
  }
  message("Scoring chunk ", i, "/", length(chunks), " (", length(chunks[[i]]), " signatures)")
  gct <- cmapR::parse_gctx(matrix_file, cid=as.character(chunks[[i]]), matrix_only=TRUE)
  mat <- gct@mat
  if (is.null(rownames(mat)) || is.null(colnames(mat))) stop("Matrix returned without gene/signature IDs.")
  if (!setequal(rownames(mat), row_ids) || !setequal(colnames(mat), chunks[[i]])) stop("Unexpected matrix subset identifiers.")
  values <- vapply(seq_len(ncol(mat)), function(j) score_one(mat[,j], rownames(mat)), numeric(3))
  scores[[i]] <- data.frame(sig_id=colnames(mat), ES_up=values[1,], ES_down=values[2,], Local_CS=values[3,])
  tmp_checkpoint <- paste0(checkpoint,".tmp")
  saveRDS(scores[[i]],tmp_checkpoint)
  if(!file.rename(tmp_checkpoint,checkpoint)) stop("Could not save checkpoint.")
  rm(mat, gct, values)
  invisible(gc(FALSE))
}
results <- left_join(selected, bind_rows(scores), by="sig_id")
write_output(results %>% filter(!is.finite(Local_CS)), "unscorable_signatures.csv")
results <- results %>% filter(is.finite(Local_CS))
if (!nrow(results)) stop("No scorable signatures.")
denominator <- mean(abs(results$Local_CS))
if (!is.finite(denominator) || denominator==0) stop("All custom connectivity scores are zero.")
results$Local_scaled_CS <- results$Local_CS/denominator
# Global mean-absolute rescaling preserves order; this is NOT CMap NCS.
write_output(results, "LINCS_signature_custom_scores.csv")

# 7. Compound aggregation ---------------------------------------------------
# Equal weighting of available signatures, preserving the original approach.
# This gives unevenly sampled conditions uneven influence; inspect cell-level
# results and compare aggregation rules before publication.
# COMPOUND SUMMARY
# Summarize by perturbagen ID, not drug name, to preserve distinct LINCS records.
# Equal weighting across signatures means heavily sampled cell/dose/time
# conditions have more influence. Cell-level summaries are exported to help
# examine heterogeneity. A pooled mean is not a cell-specific efficacy estimate.
# Minimum coverage is five signatures by default; no minimum number of cells is
# imposed. This eligibility rule is distinct from the quality filter above.
compound <- results %>% group_by(pert_id) %>% summarise(
  Drug_name=first(pert_iname), All_names=collapse_text(pert_iname),
  Mean_local_CS=mean(Local_CS), Median_local_CS=median(Local_CS),
  Mean_local_scaled_CS=mean(Local_scaled_CS),
  n_signatures=n(), n_cells=n_distinct(cell_id),
  Fraction_negative=mean(Local_CS<0), .groups="drop")
compound$Eligible <- compound$n_signatures>=minimum_signatures
compound$Rank_percentile <- NA_real_
k <- which(compound$Eligible)
if (!length(k)) stop("No compound meets minimum_signatures.")
compound$Rank_percentile[k] <- 100*rank(compound$Mean_local_CS[k],ties.method="average")/length(k)
compound$Local_screen_selected <- compound$Eligible & compound$Mean_local_CS<0 &
                                  !is.na(compound$Rank_percentile) & compound$Rank_percentile<=rank_percentile_cutoff
cell_summary <- results %>% group_by(pert_id, cell_id) %>% summarise(
  Mean_local_CS=mean(Local_CS), n_signatures=n(), .groups="drop")
write_output(cell_summary, "compound_cell_summary.csv")
info <- read.delim(compound_file, colClasses="character",check.names=FALSE)
require_columns(info,c("pert_id","moa","target"),"Compound annotation")
info <- info %>% group_by(pert_id) %>% summarise(MOA=collapse_text(moa), Targets=collapse_text(target), .groups="drop")
compound <- left_join(compound,info,by="pert_id") %>% arrange(Mean_local_CS)
write_output(compound,"LINCS_compound_custom_scores.csv")


# 8. Save run summaries, then generate the report ----------------------------
counts <- data.frame(Metric=c("Input query genes","Mapped up genes","Mapped down genes",
  "Matrix signatures","Selected signatures","Scorable signatures","Scored perturbagens",
  "Eligible perturbagens"),Count=c(nrow(query),length(up_ids),length(down_ids),
  length(column_ids),nrow(selected),nrow(results),nrow(compound),sum(compound$Eligible)))
write_output(counts,"screening_counts.csv")
capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
saveRDS(list(minimum_signatures=minimum_signatures,rank_percentile_cutoff=rank_percentile_cutoff,
 network_z_cutoff=network_z_cutoff,require_qc_pass=require_qc_pass,require_is_hiq=require_is_hiq),
 file.path(output_dir,"report_settings.rds"))
source("02_report.R")
