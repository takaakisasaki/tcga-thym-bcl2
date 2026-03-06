cat("=== Installing required R packages ===\n\n")
cran_pkgs <- c("ggplot2","ggpubr","survival","survminer","dplyr","tidyr",
               "tibble","readr","RColorBrewer","circlize","gridExtra")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing CRAN:", pkg, "\n")
    install.packages(pkg, repos="https://cran.r-project.org", quiet=TRUE)
  } else cat("OK:", pkg, "\n")
}
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cran.r-project.org", quiet=TRUE)
bioc_pkgs <- c("TCGAbiolinks","SummarizedExperiment","ComplexHeatmap")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("Installing Bioc:", pkg, "\n")
    BiocManager::install(pkg, update=FALSE, ask=FALSE, quiet=TRUE)
  } else cat("OK:", pkg, "\n")
}
cat("\nVerifying...\n")
all_ok <- TRUE
for (pkg in c(cran_pkgs, bioc_pkgs)) {
  ok <- requireNamespace(pkg, quietly=TRUE)
  if (!ok) { cat("MISSING:", pkg, "\n"); all_ok <- FALSE }
}
if (all_ok) cat("All packages OK.\n") else cat("Some packages failed.\n")
