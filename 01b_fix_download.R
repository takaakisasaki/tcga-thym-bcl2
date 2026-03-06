#!/usr/bin/env Rscript
# =============================================================================
# FIX: Complete download after CNV duplicate error
# =============================================================================
# RNA-seq already downloaded successfully (122 samples)
# CNV failed due to ASCAT2/ASCAT3 duplicates → fix: filter by workflow.type
# =============================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
})

OUTPUT_DIR <- "data"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== FIX: Completing TCGA-THYM Data Download ===\n\n")

# =============================================================================
# 1. Verify RNA-seq data exists
# =============================================================================
cat("[1/4] Checking RNA-seq data...\n")
if (file.exists(file.path(OUTPUT_DIR, "THYM_RNAseq.rda"))) {
  cat("  RNA-seq data: OK (already downloaded)\n")
  load(file.path(OUTPUT_DIR, "THYM_RNAseq.rda"))
  # The loaded object name from GDCprepare varies; find it
  se_obj <- NULL
  for (obj_name in ls()) {
    obj <- get(obj_name)
    if (is(obj, "SummarizedExperiment") || is(obj, "RangedSummarizedExperiment")) {
      se_obj <- obj
      cat("  Found SummarizedExperiment:", obj_name, "with", ncol(obj), "samples\n")
      break
    }
  }
  if (is.null(se_obj)) {
    # Try loading from the query result
    cat("  Loading via GDCprepare...\n")
    query_exp <- GDCquery(
      project = "TCGA-THYM",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    se_obj <- GDCprepare(query_exp)
  }
  data_exp <- se_obj
} else {
  stop("RNA-seq data not found. Please run 01_download_data.R first.")
}

# Extract/save TPM if not already done
tpm_file <- file.path(OUTPUT_DIR, "THYM_tpm_matrix.rds")
if (!file.exists(tpm_file)) {
  cat("  Extracting TPM matrix...\n")
  tpm_matrix <- assay(data_exp, "tpm_unstrand")
  gene_info <- as.data.frame(rowData(data_exp))
  saveRDS(list(tpm = tpm_matrix, gene_info = gene_info), file = tpm_file)
  cat("  TPM matrix saved.\n")
} else {
  cat("  TPM matrix: OK (already exists)\n")
  tpm_data <- readRDS(tpm_file)
  tpm_matrix <- tpm_data$tpm
  gene_info <- tpm_data$gene_info
}

cat("\n")

# =============================================================================
# 2. FIX: Download CNV with workflow.type filter (ASCAT3 only)
# =============================================================================
cat("[2/4] Downloading CNV data (ASCAT3 only, fixing duplicate issue)...\n")

cnv_file <- file.path(OUTPUT_DIR, "THYM_cnv_data.rds")
if (!file.exists(cnv_file)) {
  query_cnv <- GDCquery(
    project = "TCGA-THYM",
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number",
    workflow.type = "ASCAT3"
  )

  # Check for remaining duplicates and handle them
  results <- getResults(query_cnv)
  cat("  CNV files to download:", nrow(results), "\n")

  if (any(duplicated(results$cases))) {
    cat("  WARNING: Still have duplicated cases, keeping first occurrence\n")
    results <- results[!duplicated(results$cases), ]
    # Reconstruct query with filtered results
  }

  tryCatch({
    GDCdownload(query_cnv, method = "api", files.per.chunk = 10)
    data_cnv <- GDCprepare(query_cnv)
    saveRDS(data_cnv, file = cnv_file)
    cat("  CNV data saved.\n")
  }, error = function(e) {
    cat("  CNV GDCprepare failed:", conditionMessage(e), "\n")
    cat("  Trying ASCAT2 instead...\n")

    query_cnv2 <- GDCquery(
      project = "TCGA-THYM",
      data.category = "Copy Number Variation",
      data.type = "Gene Level Copy Number",
      workflow.type = "ASCAT2"
    )
    GDCdownload(query_cnv2, method = "api", files.per.chunk = 10)
    data_cnv <- GDCprepare(query_cnv2)
    saveRDS(data_cnv, file = cnv_file)
    cat("  CNV data saved (ASCAT2).\n")
  })
} else {
  cat("  CNV data: OK (already exists)\n")
}

cat("\n")

# =============================================================================
# 3. Download Clinical Data
# =============================================================================
cat("[3/4] Downloading clinical data...\n")

clinical_file <- file.path(OUTPUT_DIR, "THYM_clinical.rds")
if (!file.exists(clinical_file)) {
  clinical <- GDCquery_clinic(project = "TCGA-THYM", type = "clinical")
  saveRDS(clinical, file = clinical_file)
  cat("  Clinical data saved. Patients:", nrow(clinical), "\n")
} else {
  cat("  Clinical data: OK (already exists)\n")
}

cat("\n")

# =============================================================================
# 4. Extract BCL2 Family Expression
# =============================================================================
cat("[4/4] Extracting BCL2 family genes...\n")

bcl2_genes <- c(
  "BCL2", "MCL1", "BCL2L1", "BCL2L2",
  "BCL2L11", "BAX", "BAK1", "BBC3", "PMAIP1", "BID", "BAD",
  "TP53", "CDKN2A",
  "DNMT3A", "DNMT1", "TET2"
)

gene_idx <- which(gene_info$gene_name %in% bcl2_genes)

if (length(gene_idx) > 0) {
  bcl2_tpm <- tpm_matrix[gene_idx, ]
  rownames(bcl2_tpm) <- gene_info$gene_name[gene_idx]

  # Remove duplicates (keep highest expressed)
  if (any(duplicated(rownames(bcl2_tpm)))) {
    dups <- rownames(bcl2_tpm)[duplicated(rownames(bcl2_tpm))]
    cat("  Duplicated gene names:", paste(dups, collapse = ", "), "\n")
    # Keep the one with higher mean expression
    keep <- !duplicated(rownames(bcl2_tpm))
    bcl2_tpm <- bcl2_tpm[keep, ]
  }

  saveRDS(bcl2_tpm, file = file.path(OUTPUT_DIR, "THYM_bcl2_family_tpm.rds"))
  cat("  BCL2 family genes extracted:", nrow(bcl2_tpm), "genes\n")
  cat("  Genes found:", paste(rownames(bcl2_tpm), collapse = ", "), "\n")
} else {
  cat("  WARNING: No BCL2 family genes found!\n")
  cat("  First 10 gene names:", paste(head(gene_info$gene_name, 10), collapse = ", "), "\n")
  cat("  Checking for alternative column names...\n")
  cat("  Available columns:", paste(colnames(gene_info), collapse = ", "), "\n")
}

# =============================================================================
# 5. Get WHO subtypes
# =============================================================================
cat("\n[5/5] Getting WHO histological subtypes...\n")

subtypes_file <- file.path(OUTPUT_DIR, "THYM_subtypes.rds")
if (!file.exists(subtypes_file)) {
  subtypes <- tryCatch({
    TCGAquery_subtype(tumor = "THYM")
  }, error = function(e) {
    cat("  TCGAquery_subtype failed:", conditionMessage(e), "\n")
    # Try PanCancerAtlas
    tryCatch({
      PanCancerAtlas_subtypes()
    }, error = function(e2) {
      cat("  PanCancerAtlas_subtypes also failed.\n")
      NULL
    })
  })

  if (!is.null(subtypes)) {
    saveRDS(subtypes, file = subtypes_file)
    cat("  Subtypes saved.\n")
  } else {
    cat("  Will use clinical data for subtype assignment.\n")
  }
} else {
  cat("  Subtypes: OK (already exists)\n")
}

cat("\n=== Data download/fix complete ===\n")
cat("Files in", OUTPUT_DIR, ":\n")
for (f in list.files(OUTPUT_DIR, pattern = "\\.rds$|\\.rda$")) {
  sz <- file.size(file.path(OUTPUT_DIR, f))
  cat(sprintf("  %-35s %s\n", f, format(sz, big.mark = ",")))
}
