# Regenerate reports from saved scores. No downloads or matrix rescoring.
# For standalone use, set completed_run to a run folder before sourcing.
source("config.R")
library(dplyr)
library(ggplot2)
if(exists("completed_run") && nzchar(completed_run)) {
  output_dir <- completed_run
} else {
  pointer <- file.path(run_root,"LATEST_RUN.txt")
  if(!file.exists(pointer)) stop("Run 01_screen.R first or set completed_run.")
  output_dir <- readLines(pointer,warn=FALSE)[1]
}
s <- readRDS(file.path(output_dir,"report_settings.rds"))
report_dir <- file.path(output_dir,"report")
dir.create(report_dir,showWarnings=FALSE)
read_result <- function(n) read.csv(file.path(output_dir,n),check.names=FALSE)
write_result <- function(x,n) write.csv(x,file.path(report_dir,n),row.names=FALSE,na="NA")
compound <- read_result("LINCS_compound_custom_scores.csv")
# This linear rank rescaling is not official CMap Tau or a P value.
compound$Tau_like_score <- 2*compound$Rank_percentile-100
write_result(compound,"Compound_results.csv")
openxlsx::write.xlsx(list(Compound_results=compound,
 Query_mapping=read_result("query_mapping_audit.csv"),
 Counts=read_result("screening_counts.csv"),
 Settings=data.frame(Parameter=names(s),Value=as.character(unlist(s)))),
 file.path(report_dir,"Transcriptomic_screen.xlsx"),overwrite=TRUE)

# Optional network integration. Network scores must concern the user's own module.
# No assumed target set, interaction network, null model or P-value correction.
if(nzchar(network_file)) {
  if(!file.exists(network_file)) stop("Network CSV not found: ",network_file)
  network <- read.csv(network_file,check.names=FALSE)
  required <- c("Drug_name","Network_Z","Network_adjusted_P")
  if(!all(required %in% names(network))) stop("Network CSV needs: ",paste(required,collapse=", "))
  for(nm in required[-1]) {
    raw <- trimws(as.character(network[[nm]]))
    value <- suppressWarnings(as.numeric(raw))
    if(any(!is.na(raw) & !raw %in% c("","NA") & !is.finite(value))) stop("Invalid network values: ",nm)
    network[[nm]] <- value
  }
  if(any(network$Network_adjusted_P<0 | network$Network_adjusted_P>1,na.rm=TRUE)) stop("Adjusted P outside [0,1].")
  clean <- function(x) {
    y <- tolower(trimws(as.character(x)))
    y <- gsub("[[:space:]]+"," ",y)
    y[is.na(y) | !nzchar(y)] <- NA_character_
    y
  }
  compound$key <- clean(compound$Drug_name)
  network$key <- clean(network$Drug_name)
  ambiguous <- unique(c(compound$key[duplicated(compound$key)],network$key[duplicated(network$key)]))
  write_result(compound %>% filter(is.na(key) | key %in% ambiguous),"Ambiguous_compounds.csv")
  write_result(network %>% filter(is.na(key) | key %in% ambiguous),"Ambiguous_network.csv")
  # Prefix all matching keys to avoid collisions with unmatched IDs.
  compound$join_key <- ifelse(is.na(compound$key) | compound$key %in% ambiguous,
    paste0("C_",seq_len(nrow(compound))),paste0("MATCH_",compound$key))
  network$join_key <- ifelse(is.na(network$key) | network$key %in% ambiguous,
    paste0("N_",seq_len(nrow(network))),paste0("MATCH_",network$key))
  network$Network_name <- network$Drug_name
  combined <- full_join(compound,network %>% select(join_key,Network_name,Network_Z,Network_adjusted_P),by="join_key") %>%
    mutate(Drug_name=coalesce(Drug_name,Network_name),
      Joint_exploratory=ifelse(is.finite(Tau_like_score) & is.finite(Network_Z),
                               Local_screen_selected & Network_Z<=s$network_z_cutoff,NA),
      Evidence_category=case_when(
        Joint_exploratory %in% TRUE ~ "Joint exploratory",
        Network_Z<0 & Network_adjusted_P<.05 ~ "Significant proximity",
        Local_screen_selected %in% TRUE ~ "Transcriptomic selection",
        Network_Z<=s$network_z_cutoff ~ "Exploratory proximity",
        TRUE ~ "Other"))
  write_result(combined,"Integrated_results.csv")
  write_result(network,"Network_results.csv")
  plot_data <- combined %>% filter(is.finite(Tau_like_score),is.finite(Network_Z))
  if(nrow(plot_data)) {
    labels <- plot_data %>% filter(Joint_exploratory %in% TRUE |
                                    (Tau_like_score<= -50 & Network_Z<=s$network_z_cutoff))
    p <- ggplot(plot_data,aes(Network_Z,Tau_like_score))+
      geom_rect(data=data.frame(xmin=-Inf,xmax=s$network_z_cutoff,ymin=-100,
         ymax=2*s$rank_percentile_cutoff-100,Region="Joint selection region"),
         aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Region),
         inherit.aes=FALSE,alpha=.25,color=NA,
         show.legend=c(fill=TRUE,colour=FALSE,shape=FALSE))+
      geom_point(aes(color=Evidence_category,shape=Evidence_category),alpha=.75,size=2)+
      geom_vline(xintercept=s$network_z_cutoff,linetype="dashed")+
      geom_hline(yintercept=2*s$rank_percentile_cutoff-100,linetype="dashed")+
      scale_fill_manual(values=c("Joint selection region"="#F0E442"))+
      scale_color_manual(values=c("Joint exploratory"="black","Significant proximity"="#CC79A7",
        "Transcriptomic selection"="#D55E00","Exploratory proximity"="#0072B2","Other"="grey65"))+
      scale_shape_manual(values=c("Joint exploratory"=18,"Significant proximity"=16,
        "Transcriptomic selection"=16,"Exploratory proximity"=16,"Other"=16))+
      ggrepel::geom_text_repel(data=labels,aes(label=Drug_name),seed=123,max.overlaps=Inf,size=3,show.legend=FALSE)+
      labs(x="Network proximity Z-score",y="Tau-like connectivity score",color=NULL,shape=NULL,fill=NULL)+
      theme_classic(base_size=12)+theme(legend.position="bottom",legend.text=element_text(size=12))
    temp <- tempfile(fileext=".png")
    ggsave(temp,p,width=14,height=9,dpi=600,bg="white")
    if(!file.copy(temp,file.path(report_dir,"Integrated_screen.png"),overwrite=TRUE)) stop("PNG copy failed.")
    write_result(plot_data,"Plotted_data.csv")
    write_result(labels,"Labeled_candidates.csv")
  }
}
message("Report saved: ",normalizePath(report_dir))
