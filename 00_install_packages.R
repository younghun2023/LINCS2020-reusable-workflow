# Optional explicit dependency installation. Does not download reference data.
cran <- c("dplyr","readxl","openxlsx","ggplot2","ggrepel","BiocManager")
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
if(!requireNamespace("cmapR",quietly=TRUE)) BiocManager::install("cmapR",ask=FALSE,update=FALSE)
