#!/usr/bin/env Rscript
# =============================================================================
# Alternative: cBioPortal API Approach
# =============================================================================
# Use this if TCGAbiolinks has issues or to get WHO subtype mapping
# cBioPortal has curated TCGA-THYM data with WHO histological classification
# =============================================================================

suppressPackageStartupMessages({
  library(cBioPortalData)
  library(dplyr)
  library(readr)
})

DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE)

cat("=== cBioPortal Data Retrieval for TCGA-THYM ===\n\n")

# --- Connect to cBioPortal ---
cbio <- cBioPortal()

# --- Study: thym_tcga ---
study_id <- "thym_tcga"

cat("[1] Downloading clinical data with WHO subtypes...\n")
clin <- clinicalData(cbio, studyId = study_id)
cat("  Patients:", nrow(clin), "\n")

# WHO subtype column varies: check available columns
cat("  Available clinical columns:\n")
cat("  ", paste(colnames(clin), collapse = ", "), "\n\n")

# Look for histological type column
hist_cols <- grep("histol|subtype|WHO|type", colnames(clin), ignore.case = TRUE, value = TRUE)
if (length(hist_cols) > 0) {
  cat("  Histology-related columns:", paste(hist_cols, collapse = ", "), "\n")
  for (hc in hist_cols) {
    cat("  ", hc, ":\n")
    print(table(clin[[hc]], useNA = "ifany"))
  }
}

saveRDS(clin, file.path(DATA_DIR, "THYM_cbio_clinical.rds"))

cat("\n[2] Downloading molecular profiles...\n")
profiles <- molecularProfiles(cbio, studyId = study_id)
cat("  Available profiles:\n")
print(profiles[, c("molecularProfileId", "name", "molecularAlterationType")])

# --- Get BCL2 family expression ---
cat("\n[3] Downloading BCL2 family expression...\n")

bcl2_genes <- c("BCL2", "MCL1", "BCL2L1", "BCL2L2", "BCL2L11",
                "BAX", "BAK1", "BBC3", "PMAIP1", "BID", "BAD",
                "TP53", "CDKN2A", "DNMT3A", "DNMT1", "TET2")

# Get mRNA expression (z-scores or RSEM)
mrna_profile <- grep("mrna|rna_seq", profiles$molecularProfileId,
                     ignore.case = TRUE, value = TRUE)
cat("  mRNA profiles:", paste(mrna_profile, collapse = ", "), "\n")

if (length(mrna_profile) > 0) {
  # Get all patient IDs
  samples <- allSamples(cbio, studyId = study_id)

  for (prof in mrna_profile) {
    cat("  Fetching from:", prof, "...\n")
    expr_data <- tryCatch({
      getDataByGenes(
        api = cbio,
        studyId = study_id,
        molecularProfileIds = prof,
        genes = bcl2_genes,
        by = "hugoGeneSymbol"
      )
    }, error = function(e) {
      cat("    Error:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(expr_data) && length(expr_data) > 0) {
      saveRDS(expr_data, file.path(DATA_DIR,
              paste0("THYM_cbio_expr_", gsub("[^a-zA-Z0-9]", "_", prof), ".rds")))
      cat("    Saved.\n")
    }
  }
}

# --- Get Copy Number data ---
cat("\n[4] Downloading BCL2 copy number data...\n")
cn_profile <- grep("gistic|cna|copy", profiles$molecularProfileId,
                   ignore.case = TRUE, value = TRUE)
cat("  CN profiles:", paste(cn_profile, collapse = ", "), "\n")

if (length(cn_profile) > 0) {
  samples <- allSamples(cbio, studyId = study_id)

  for (prof in cn_profile) {
    cat("  Fetching from:", prof, "...\n")
    cn_data <- tryCatch({
      getDataByGenes(
        api = cbio,
        studyId = study_id,
        molecularProfileIds = prof,
        genes = bcl2_genes,
        by = "hugoGeneSymbol"
      )
    }, error = function(e) {
      cat("    Error:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(cn_data) && length(cn_data) > 0) {
      saveRDS(cn_data, file.path(DATA_DIR,
              paste0("THYM_cbio_cn_", gsub("[^a-zA-Z0-9]", "_", prof), ".rds")))
      cat("    Saved.\n")
    }
  }
}

# --- Get Mutation data ---
cat("\n[5] Downloading BCL2 mutation data...\n")
mut_profile <- grep("mutation", profiles$molecularProfileId,
                    ignore.case = TRUE, value = TRUE)
if (length(mut_profile) > 0) {
  for (prof in mut_profile) {
    cat("  Fetching from:", prof, "...\n")
    mut_data <- tryCatch({
      getDataByGenes(
        api = cbio,
        studyId = study_id,
        molecularProfileIds = prof,
        genes = bcl2_genes,
        by = "hugoGeneSymbol"
      )
    }, error = function(e) {
      cat("    Error:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(mut_data) && length(mut_data) > 0) {
      saveRDS(mut_data, file.path(DATA_DIR,
              paste0("THYM_cbio_mut_", gsub("[^a-zA-Z0-9]", "_", prof), ".rds")))
      cat("    Saved.\n")
    }
  }
}

cat("\n=== cBioPortal data retrieval complete ===\n")
cat("Files saved in:", normalizePath(DATA_DIR), "\n")
