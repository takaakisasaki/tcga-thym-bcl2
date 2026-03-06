suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(dplyr)
})
OUTPUT_DIR <- "data"
dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat("=== TCGA-THYM Data Download ===\n\n")

# 1. RNA-seq
cat("[1/3] RNA-seq expression...\n")
query_exp <- GDCquery(project="TCGA-THYM",
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts")
GDCdownload(query_exp, method="api", files.per.chunk=10)
data_exp <- GDCprepare(query_exp, save=TRUE,
  save.filename=file.path(OUTPUT_DIR,"THYM_RNAseq.rda"))
tpm_matrix <- assay(data_exp, "tpm_unstrand")
gene_info <- as.data.frame(rowData(data_exp))
saveRDS(list(tpm=tpm_matrix, gene_info=gene_info),
  file.path(OUTPUT_DIR,"THYM_tpm_matrix.rds"))
cat("  Samples:", ncol(data_exp), "\n\n")

# 2. Copy Number
cat("[2/3] Copy number...\n")
query_cnv <- tryCatch({
  GDCquery(project="TCGA-THYM",
    data.category="Copy Number Variation",
    data.type="Gene Level Copy Number")
}, error=function(e) {
  cat("  Gene Level CN not available, trying Copy Number Segment...\n")
  GDCquery(project="TCGA-THYM",
    data.category="Copy Number Variation",
    data.type="Copy Number Segment")
})
GDCdownload(query_cnv, method="api", files.per.chunk=10)
data_cnv <- GDCprepare(query_cnv)
saveRDS(data_cnv, file.path(OUTPUT_DIR,"THYM_cnv_data.rds"))
cat("  CNV data saved.\n\n")

# 3. Clinical
cat("[3/3] Clinical data...\n")
clinical <- GDCquery_clinic(project="TCGA-THYM", type="clinical")
saveRDS(clinical, file.path(OUTPUT_DIR,"THYM_clinical.rds"))
cat("  Patients:", nrow(clinical), "\n\n")

# 4. Subtypes
cat("[Subtypes] Trying TCGAquery_subtype...\n")
subtypes <- tryCatch({
  PanCancerAtlas_subtypes() %>% filter(cancer.type == "THYM")
}, error=function(e) {
  tryCatch(TCGAquery_subtype(tumor="THYM"), error=function(e2) NULL)
})
if (!is.null(subtypes) && nrow(subtypes) > 0) {
  saveRDS(subtypes, file.path(OUTPUT_DIR,"THYM_subtypes.rds"))
  cat("  Subtypes:", nrow(subtypes), "samples\n")
} else {
  cat("  Subtypes not available via API. Will use clinical data.\n")
}

# 5. Extract BCL2 family
cat("\n[BCL2 family extraction...]\n")
bcl2_genes <- c("BCL2","MCL1","BCL2L1","BCL2L2","BCL2L11","BAX","BAK1",
  "BBC3","PMAIP1","BID","BAD","TP53","CDKN2A","DNMT3A","DNMT1","TET2")
gene_idx <- which(gene_info$gene_name %in% bcl2_genes)
if (length(gene_idx) > 0) {
  bcl2_tpm <- tpm_matrix[gene_idx, ]
  rownames(bcl2_tpm) <- gene_info$gene_name[gene_idx]
  bcl2_tpm <- bcl2_tpm[!duplicated(rownames(bcl2_tpm)), ]
  saveRDS(bcl2_tpm, file.path(OUTPUT_DIR,"THYM_bcl2_family_tpm.rds"))
  cat("  Found:", paste(rownames(bcl2_tpm), collapse=", "), "\n")
}

cat("\n=== Download complete ===\n")
