#!/usr/bin/env Rscript
# =============================================================================
# TCGA-THYM BCL2 Analysis - V2 (fixed for actual colData structure)
# =============================================================================
# Generates 4 publication-quality figures for thymic carcinoma case report
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(SummarizedExperiment)
})

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

cat("=== TCGA-THYM BCL2 Analysis V2 ===\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================
cat("[1] Loading data...\n")

bcl2_tpm <- readRDS("data/THYM_bcl2_family_tpm.rds")
clinical  <- readRDS("data/THYM_clinical.rds")
cnv_se    <- readRDS("data/THYM_cnv_data.rds")
tpm_data  <- readRDS("data/THYM_tpm_matrix.rds")

cat("  BCL2 TPM:", nrow(bcl2_tpm), "genes x", ncol(bcl2_tpm), "samples\n")
cat("  Clinical:", nrow(clinical), "samples\n")
cat("  CNV SE:", ncol(cnv_se), "samples\n")

# =============================================================================
# 2. Assign WHO subtypes from primary_diagnosis
# =============================================================================
cat("\n[2] Assigning WHO subtypes...\n")

# clinical is colData from SE; barcode = column names of TPM matrix
meta <- data.frame(
  barcode = clinical$barcode,
  patient = clinical$patient,
  primary_diagnosis = clinical$primary_diagnosis,
  vital_status = clinical$vital_status,
  days_to_death = as.numeric(clinical$days_to_death),
  days_to_last_follow_up = as.numeric(clinical$days_to_last_follow_up),
  age_at_diagnosis = as.numeric(clinical$age_at_diagnosis),
  gender = clinical$gender,
  masaoka_stage = if ("masaoka_stage" %in% colnames(clinical)) clinical$masaoka_stage else NA,
  stringsAsFactors = FALSE
)

# Parse WHO subtype from primary_diagnosis
parse_who <- function(dx) {
  if (is.na(dx) || dx == "" || dx == "Not Reported") return("Unknown")
  dx <- tolower(dx)
  if (grepl("carcinoma|thymic carcinoma|squamous cell", dx)) return("TC")
  if (grepl("type a,|type a[^b]|\\btype a\\b", dx)) return("A")
  if (grepl("type ab", dx)) return("AB")
  if (grepl("type b1", dx)) return("B1")
  if (grepl("type b2", dx)) return("B2")
  if (grepl("type b3", dx)) return("B3")
  if (grepl("thymoma", dx)) return("Unknown")
  return("Unknown")
}

meta$who_subtype <- sapply(meta$primary_diagnosis, parse_who)
meta$who_subtype <- factor(meta$who_subtype,
                           levels = c("A", "AB", "B1", "B2", "B3", "TC"))

cat("  Subtype counts:\n")
print(table(meta$who_subtype, useNA = "ifany"))

# =============================================================================
# 3. Merge BCL2 expression with metadata
# =============================================================================
cat("\n[3] Merging expression data...\n")

# BCL2 TPM colnames should match clinical$barcode
bcl2_df <- as.data.frame(t(bcl2_tpm))
bcl2_df$barcode <- rownames(bcl2_df)

merged <- inner_join(bcl2_df, meta, by = "barcode")
cat("  Merged samples:", nrow(merged), "\n")

# Filter to known subtypes
merged_known <- merged %>% filter(!is.na(who_subtype))
cat("  With known subtype:", nrow(merged_known), "\n")

# =============================================================================
# Figure A: BCL2 mRNA expression by WHO subtype
# =============================================================================
cat("\n[Fig A] BCL2 expression by WHO subtype...\n")

# Kruskal-Wallis test
kw_test <- kruskal.test(BCL2 ~ who_subtype, data = merged_known)
cat("  Kruskal-Wallis p =", format(kw_test$p.value, digits = 3), "\n")

# Pairwise Wilcoxon: TC vs each
tc_data <- merged_known %>% filter(who_subtype == "TC")
pw_results <- data.frame()
for (sub in c("A", "AB", "B1", "B2", "B3")) {
  sub_data <- merged_known %>% filter(who_subtype == sub)
  if (nrow(sub_data) > 1 && nrow(tc_data) > 1) {
    wt <- wilcox.test(tc_data$BCL2, sub_data$BCL2)
    pw_results <- rbind(pw_results, data.frame(
      comparison = paste("TC vs", sub),
      p.value = wt$p.value,
      TC_median = median(tc_data$BCL2),
      other_median = median(sub_data$BCL2)
    ))
  }
}
pw_results$p.adj <- p.adjust(pw_results$p.value, method = "BH")
cat("  Pairwise comparisons:\n")
print(pw_results)

# Color palette
subtype_colors <- c("A" = "#4DAF4A", "AB" = "#377EB8", "B1" = "#FF7F00",
                     "B2" = "#E41A1C", "B3" = "#984EA3", "TC" = "#A65628")

# Create comparison list for stat_compare_means
comparisons_list <- lapply(c("A", "AB", "B1", "B2", "B3"),
                           function(x) c("TC", x))

p_a <- ggplot(merged_known, aes(x = who_subtype, y = log2(BCL2 + 1), fill = who_subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  scale_fill_manual(values = subtype_colors) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test", label = "p.signif",
                     size = 3, tip.length = 0.01,
                     step.increase = 0.08) +
  labs(
    title = "BCL2 mRNA Expression in TCGA-THYM by WHO Subtype",
    subtitle = sprintf("Kruskal-Wallis p = %s  (n = %d)", 
                       format(kw_test$p.value, digits = 3), nrow(merged_known)),
    x = "WHO Histological Subtype",
    y = expression(log[2](TPM + 1)),
    fill = "Subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "gray40"))

ggsave("figures/FigA_BCL2_expression_by_subtype.pdf", p_a, width = 7, height = 6)
ggsave("figures/FigA_BCL2_expression_by_subtype.png", p_a, width = 7, height = 6, dpi = 300)
cat("  Saved FigA\n")

# =============================================================================
# Figure B: BCL2 Copy Number by subtype
# =============================================================================
cat("\n[Fig B] BCL2 copy number by subtype...\n")

# Extract BCL2 CN from the CNV SummarizedExperiment
cnv_gene_info <- as.data.frame(rowData(cnv_se))
cnv_coldata <- as.data.frame(colData(cnv_se))

# Find BCL2 in CNV data
bcl2_row <- which(cnv_gene_info$gene_name == "BCL2")
if (length(bcl2_row) == 0) {
  # Try alternative column
  for (col in c("gene_name", "gene_id", "Gene.Symbol", "Hugo_Symbol")) {
    if (col %in% colnames(cnv_gene_info)) {
      bcl2_row <- which(cnv_gene_info[[col]] == "BCL2")
      if (length(bcl2_row) > 0) break
    }
  }
}

if (length(bcl2_row) > 0) {
  bcl2_cn <- assay(cnv_se, "copy_number")[bcl2_row[1], ]
  
  cnv_df <- data.frame(
    barcode = colnames(cnv_se),
    patient = cnv_coldata$patient,
    bcl2_cn = as.numeric(bcl2_cn),
    stringsAsFactors = FALSE
  )
  
  # Categorize CN: normal diploid = 2
  cnv_df$cn_cat <- cut(cnv_df$bcl2_cn,
                       breaks = c(-Inf, 0.5, 1.5, 2.5, 3.5, Inf),
                       labels = c("Deep Del", "Shallow Del", "Diploid", "Gain", "Amp"))
  
  # Match patients to WHO subtype
  patient_subtype <- meta %>%
    select(patient, who_subtype) %>%
    distinct(patient, .keep_all = TRUE)
  
  cnv_df <- left_join(cnv_df, patient_subtype, by = "patient")
  cnv_known <- cnv_df %>% filter(!is.na(who_subtype))
  
  cat("  CNV samples with subtype:", nrow(cnv_known), "\n")
  
  # Stacked bar plot
  cn_colors <- c("Deep Del" = "#0571B0", "Shallow Del" = "#92C5DE",
                 "Diploid" = "#F7F7F7", "Gain" = "#F4A582", "Amp" = "#CA0020")
  
  cn_summary <- cnv_known %>%
    count(who_subtype, cn_cat) %>%
    group_by(who_subtype) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  
  p_b <- ggplot(cn_summary, aes(x = who_subtype, y = pct, fill = cn_cat)) +
    geom_bar(stat = "identity", width = 0.7, color = "gray30", linewidth = 0.3) +
    scale_fill_manual(values = cn_colors, name = "Copy Number") +
    labs(
      title = "BCL2 Copy Number Alterations by WHO Subtype",
      subtitle = sprintf("TCGA-THYM (n = %d)", nrow(cnv_known)),
      x = "WHO Histological Subtype",
      y = "Percentage of Samples (%)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, color = "gray40"))
  
  ggsave("figures/FigB_BCL2_copynumber_by_subtype.pdf", p_b, width = 7, height = 5)
  ggsave("figures/FigB_BCL2_copynumber_by_subtype.png", p_b, width = 7, height = 5, dpi = 300)
  cat("  Saved FigB\n")
  
  # TC-specific stats
  tc_cnv <- cnv_known %>% filter(who_subtype == "TC")
  cat("  TC CNV breakdown:\n")
  print(table(tc_cnv$cn_cat))
  
} else {
  cat("  WARNING: BCL2 not found in CNV data.\n")
  cat("  Available gene columns:", paste(colnames(cnv_gene_info)[1:10], collapse = ", "), "\n")
  cat("  First genes:", paste(head(cnv_gene_info[[1]], 5), collapse = ", "), "\n")
  # Create placeholder
  p_b <- NULL
}

# =============================================================================
# Figure C: Kaplan-Meier survival by BCL2 expression
# =============================================================================
cat("\n[Fig C] Survival analysis...\n")

# Prepare survival data
surv_df <- merged_known %>%
  mutate(
    os_time = ifelse(!is.na(days_to_death) & days_to_death > 0,
                     days_to_death,
                     days_to_last_follow_up),
    os_event = ifelse(vital_status == "Dead", 1, 0),
    bcl2_group = ifelse(BCL2 >= median(BCL2, na.rm = TRUE), "BCL2-High", "BCL2-Low")
  ) %>%
  filter(!is.na(os_time) & os_time > 0)

cat("  Survival samples:", nrow(surv_df), "\n")
cat("  Events:", sum(surv_df$os_event), "\n")
cat("  BCL2-High:", sum(surv_df$bcl2_group == "BCL2-High"),
    " BCL2-Low:", sum(surv_df$bcl2_group == "BCL2-Low"), "\n")

if (nrow(surv_df) > 10 && sum(surv_df$os_event) >= 3) {
  # Kaplan-Meier
  fit <- survfit(Surv(os_time / 365.25, os_event) ~ bcl2_group, data = surv_df)
  
  # Log-rank test
  logrank <- survdiff(Surv(os_time / 365.25, os_event) ~ bcl2_group, data = surv_df)
  logrank_p <- 1 - pchisq(logrank$chisq, 1)
  cat("  Log-rank p =", format(logrank_p, digits = 3), "\n")
  
  # Cox regression (univariate)
  cox_uni <- coxph(Surv(os_time, os_event) ~ BCL2, data = surv_df)
  cox_summary <- summary(cox_uni)
  cat("  Cox HR (continuous BCL2):", round(cox_summary$conf.int[1, 1], 3),
      " (", round(cox_summary$conf.int[1, 3], 3), "-",
      round(cox_summary$conf.int[1, 4], 3), ")\n")
  
  # Cox multivariate (adjusted for subtype)
  cox_multi <- tryCatch({
    coxph(Surv(os_time, os_event) ~ BCL2 + who_subtype, data = surv_df)
  }, error = function(e) NULL)
  
  # KM plot
  p_c <- ggsurvplot(
    fit,
    data = surv_df,
    pval = TRUE,
    pval.coord = c(0.1, 0.1),
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = c("#E41A1C", "#377EB8"),
    legend.labs = c("BCL2-High", "BCL2-Low"),
    legend.title = "",
    xlab = "Time (years)",
    ylab = "Overall Survival Probability",
    title = "Overall Survival by BCL2 Expression (TCGA-THYM)",
    subtitle = sprintf("Median split | n = %d, events = %d",
                       nrow(surv_df), sum(surv_df$os_event)),
    ggtheme = theme_classic(base_size = 12),
    font.title = c(13, "bold"),
    font.subtitle = c(10, "plain", "gray40")
  )
  
  pdf("figures/FigC_BCL2_survival_KM.pdf", width = 7, height = 6)
  print(p_c)
  dev.off()
  
  png("figures/FigC_BCL2_survival_KM.png", width = 7, height = 6, units = "in", res = 300)
  print(p_c)
  dev.off()
  
  cat("  Saved FigC\n")
} else {
  cat("  WARNING: Too few events for survival analysis.\n")
  p_c <- NULL
}

# =============================================================================
# Figure D: BCL2 family heatmap
# =============================================================================
cat("\n[Fig D] BCL2 family heatmap...\n")

# Prepare heatmap data (Z-score of log2 TPM)
hm_data <- log2(bcl2_tpm + 1)
hm_zscore <- t(scale(t(hm_data)))

# Order genes by functional category
gene_categories <- data.frame(
  gene = c("BCL2", "MCL1", "BCL2L1", "BCL2L2",
           "BCL2L11", "BAX", "BAK1", "BBC3", "PMAIP1", "BID", "BAD",
           "TP53", "CDKN2A",
           "DNMT3A", "DNMT1", "TET2"),
  category = c(rep("Anti-apoptotic", 4),
               rep("Pro-apoptotic", 7),
               rep("Tumor Suppressor", 2),
               rep("Epigenetic", 3)),
  stringsAsFactors = FALSE
)

# Match available genes
avail_genes <- intersect(gene_categories$gene, rownames(hm_zscore))
gene_categories <- gene_categories %>% filter(gene %in% avail_genes)
hm_zscore <- hm_zscore[gene_categories$gene, ]

cat("  Heatmap genes:", nrow(hm_zscore), "\n")

# Sample annotation (WHO subtype)
sample_meta <- meta %>% filter(barcode %in% colnames(hm_zscore))
# Reorder to match heatmap columns
sample_meta <- sample_meta[match(colnames(hm_zscore), sample_meta$barcode), ]

# Top annotation
ha_top <- HeatmapAnnotation(
  `WHO Subtype` = sample_meta$who_subtype,
  col = list(`WHO Subtype` = subtype_colors),
  na_col = "gray90",
  annotation_name_side = "left",
  show_legend = TRUE
)

# Row annotation (gene category)
cat_colors <- c("Anti-apoptotic" = "#E41A1C", "Pro-apoptotic" = "#377EB8",
                "Tumor Suppressor" = "#4DAF4A", "Epigenetic" = "#FF7F00")

ha_row <- rowAnnotation(
  Category = gene_categories$category,
  col = list(Category = cat_colors),
  show_legend = TRUE
)

# Cap extreme z-scores for visualization
hm_zscore[hm_zscore > 3] <- 3
hm_zscore[hm_zscore < -3] <- -3

# Column order: group by subtype
col_order <- order(sample_meta$who_subtype, na.last = TRUE)

# Create heatmap
col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

ht <- Heatmap(
  hm_zscore[, col_order],
  name = "Z-score",
  col = col_fun,
  top_annotation = HeatmapAnnotation(
    `WHO Subtype` = sample_meta$who_subtype[col_order],
    col = list(`WHO Subtype` = subtype_colors),
    na_col = "gray90",
    annotation_name_side = "left"
  ),
  right_annotation = ha_row,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = "BCL-2 Family Expression Profile in TCGA-THYM",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  row_split = factor(gene_categories$category,
                     levels = c("Anti-apoptotic", "Pro-apoptotic",
                                "Tumor Suppressor", "Epigenetic")),
  row_gap = unit(2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Z-score\n(log2 TPM)",
    legend_width = unit(4, "cm"),
    direction = "horizontal"
  )
)

pdf("figures/FigD_BCL2_family_heatmap.pdf", width = 10, height = 6)
draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE,
     padding = unit(c(5, 5, 5, 15), "mm"))
dev.off()

png("figures/FigD_BCL2_family_heatmap.png", width = 10, height = 6, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE,
     padding = unit(c(5, 5, 5, 15), "mm"))
dev.off()

cat("  Saved FigD\n")

# =============================================================================
# Summary Tables
# =============================================================================
cat("\n[Tables] Saving summary statistics...\n")

# Table 1: BCL2 expression summary by subtype
expr_summary <- merged_known %>%
  group_by(who_subtype) %>%
  summarise(
    n = n(),
    BCL2_median = round(median(BCL2, na.rm = TRUE), 2),
    BCL2_mean = round(mean(BCL2, na.rm = TRUE), 2),
    BCL2_sd = round(sd(BCL2, na.rm = TRUE), 2),
    BCL2_min = round(min(BCL2, na.rm = TRUE), 2),
    BCL2_max = round(max(BCL2, na.rm = TRUE), 2),
    MCL1_median = round(median(MCL1, na.rm = TRUE), 2),
    BCL2L1_median = round(median(BCL2L1, na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(expr_summary, "tables/BCL2_expression_by_subtype.csv", row.names = FALSE)

# Table 2: Pairwise comparisons
write.csv(pw_results, "tables/BCL2_pairwise_wilcoxon_TC_vs_each.csv", row.names = FALSE)

# Table 3: BCL2 family expression (all genes, all subtypes)
family_summary <- merged_known %>%
  select(who_subtype, all_of(gene_categories$gene)) %>%
  pivot_longer(-who_subtype, names_to = "gene", values_to = "TPM") %>%
  group_by(who_subtype, gene) %>%
  summarise(
    median_TPM = round(median(TPM, na.rm = TRUE), 2),
    mean_TPM = round(mean(TPM, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = who_subtype, values_from = c(median_TPM, mean_TPM))
write.csv(family_summary, "tables/BCL2_family_expression_all_subtypes.csv", row.names = FALSE)

# Table 4: Survival statistics
if (exists("cox_uni")) {
  surv_stats <- data.frame(
    test = c("Log-rank", "Cox univariate (BCL2)", "Cox multivariate (BCL2 + subtype)"),
    p_value = c(
      format(logrank_p, digits = 4),
      format(cox_summary$coefficients[1, "Pr(>|z|)"], digits = 4),
      if (!is.null(cox_multi)) format(summary(cox_multi)$coefficients["BCL2", "Pr(>|z|)"], digits = 4) else "N/A"
    ),
    HR = c(
      "N/A",
      round(cox_summary$conf.int[1, 1], 3),
      if (!is.null(cox_multi)) round(summary(cox_multi)$conf.int["BCL2", 1], 3) else "N/A"
    ),
    HR_95CI = c(
      "N/A",
      sprintf("%.3f-%.3f", cox_summary$conf.int[1, 3], cox_summary$conf.int[1, 4]),
      if (!is.null(cox_multi)) sprintf("%.3f-%.3f",
        summary(cox_multi)$conf.int["BCL2", 3],
        summary(cox_multi)$conf.int["BCL2", 4]) else "N/A"
    )
  )
  write.csv(surv_stats, "tables/BCL2_survival_statistics.csv", row.names = FALSE)
}

cat("  Tables saved.\n")

# =============================================================================
# Final Summary
# =============================================================================
cat("\n")
cat("=============================================\n")
cat(" ANALYSIS COMPLETE\n")
cat("=============================================\n")
cat("\nFigures:\n")
for (f in list.files("figures/")) cat("  ", f, "\n")
cat("\nTables:\n")
for (f in list.files("tables/")) cat("  ", f, "\n")
cat("\n")

# Key findings for the paper
cat("=== KEY FINDINGS FOR DISCUSSION ===\n")
tc_expr <- merged_known %>% filter(who_subtype == "TC")
cat(sprintf("TC samples: n = %d\n", nrow(tc_expr)))
cat(sprintf("BCL2 median TPM in TC: %.2f\n", median(tc_expr$BCL2)))
cat(sprintf("BCL2 median TPM overall: %.2f\n", median(merged_known$BCL2)))
cat(sprintf("MCL1 median TPM in TC: %.2f\n", median(tc_expr$MCL1)))
cat(sprintf("BCL2L1 (BCL-xL) median TPM in TC: %.2f\n", median(tc_expr$BCL2L1)))
cat(sprintf("Kruskal-Wallis p = %s\n", format(kw_test$p.value, digits = 4)))
if (exists("logrank_p")) cat(sprintf("Log-rank survival p = %s\n", format(logrank_p, digits = 4)))
