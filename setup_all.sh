#!/bin/bash
# =============================================================================
# SHIROKANE Setup: このスクリプトをrc001で実行すると全ファイルが生成されます
# Usage: bash setup_all.sh
# =============================================================================
set -e

WORK_DIR="${HOME}/tcga_thym_bcl2"
mkdir -p "${WORK_DIR}"/{logs,data,figures,tables}
cd "${WORK_DIR}"

echo "Creating pipeline files in ${WORK_DIR}..."

# =============================================================================
# check_env.sh
# =============================================================================
cat > check_env.sh << 'SCRIPT_END'
#!/bin/bash
echo "=== SHIROKANE Environment Check ==="
module use /usr/local/package/modulefiles 2>/dev/null
echo ""
echo "[1] Available R versions:"
module avail 2>&1 | grep -i "R/" | head -10
echo ""
echo "[2] Loading R..."
module load R 2>/dev/null && echo "  Loaded: $(which R)" || echo "  Failed. Try: module load R/4.x.x"
echo ""
echo "[3] R version:"
R --version 2>/dev/null | head -1 || echo "  R not found"
echo ""
echo "[4] Checking R packages..."
Rscript -e '
pkgs <- c("TCGAbiolinks","ggplot2","survival","ComplexHeatmap","survminer","ggpubr","circlize","dplyr")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(sprintf("  %-20s %s\n", p, ifelse(ok,"OK","MISSING")))
}' 2>/dev/null || echo "  Could not run Rscript"
echo ""
echo "[5] GDC API connectivity:"
curl -s --connect-timeout 5 https://api.gdc.cancer.gov/status | head -c 100 && echo -e "\n  GDC API: OK" || echo "  GDC API: BLOCKED"
echo ""
echo "[6] qsub:"
which qsub 2>/dev/null && echo "  OK" || echo "  NOT FOUND"
echo ""
echo "=== Next Steps ==="
echo "  1. Rscript 00_install_packages.R   # if MISSING above"
echo "  2. qsub run_shirokane.sh           # submit job"
echo "  Or interactive: qlogin -l s_vmem=32G,mem_req=32G"
SCRIPT_END

# =============================================================================
# 00_install_packages.R
# =============================================================================
cat > 00_install_packages.R << 'SCRIPT_END'
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
SCRIPT_END

# =============================================================================
# 01_download_data.R
# =============================================================================
cat > 01_download_data.R << 'SCRIPT_END'
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
SCRIPT_END

# =============================================================================
# 02_analysis_visualization.R
# =============================================================================
cat > 02_analysis_visualization.R << 'SCRIPT_END'
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(survival); library(survminer)
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
  library(dplyr); library(tidyr); library(tibble)
})

DATA_DIR <- "data"; FIG_DIR <- "figures"; TABLE_DIR <- "tables"
dir.create(FIG_DIR, showWarnings=FALSE); dir.create(TABLE_DIR, showWarnings=FALSE)

WHO_COLORS <- c(A="#4DAF4A", AB="#377EB8", B1="#FF7F00", B2="#E41A1C",
                B3="#984EA3", TC="#000000", C="#000000")

theme_pub <- theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold",size=14),
        axis.title=element_text(face="bold"), plot.margin=margin(10,15,10,10))

cat("=== TCGA-THYM BCL2 Analysis ===\n\n")

# --- Load data ----------------------------------------------------------------
bcl2_tpm <- readRDS(file.path(DATA_DIR,"THYM_bcl2_family_tpm.rds"))
clinical <- readRDS(file.path(DATA_DIR,"THYM_clinical.rds"))
subtypes <- tryCatch(readRDS(file.path(DATA_DIR,"THYM_subtypes.rds")), error=function(e) NULL)
cnv_data <- tryCatch(readRDS(file.path(DATA_DIR,"THYM_cnv_data.rds")), error=function(e) NULL)

# --- Build metadata -----------------------------------------------------------
meta <- data.frame(sample_id=colnames(bcl2_tpm),
                   patient_id=substr(colnames(bcl2_tpm),1,12), stringsAsFactors=FALSE)
if ("submitter_id" %in% colnames(clinical))
  meta <- meta %>% left_join(clinical, by=c("patient_id"="submitter_id"))

# WHO subtype assignment
assign_who <- function(meta, subtypes) {
  if (!is.null(subtypes)) {
    sub_cols <- intersect(colnames(subtypes),
      c("patient","pan.samplesID","Histology","histology","Subtype","subtype",
        "WHO.type","WHO_Classification","histological_type"))
    if (length(sub_cols) >= 1) {
      id_col <- sub_cols[grep("patient|sample", sub_cols, ignore.case=TRUE)[1]]
      type_col <- setdiff(sub_cols, id_col)[1]
      if (!is.na(id_col) && !is.na(type_col)) {
        mapping <- subtypes[, c(id_col, type_col)]
        colnames(mapping) <- c("match_id", "who_raw")
        mapping$match_id_short <- substr(mapping$match_id, 1, 12)
        meta <- meta %>% left_join(mapping, by=c("patient_id"="match_id_short"))
        if ("who_raw" %in% colnames(meta) && !all(is.na(meta$who_raw))) {
          meta$who_subtype <- meta$who_raw
          cat("  WHO subtype from subtypes data.\n")
          return(meta)
        }
      }
    }
  }
  # Fallback: infer from primary_diagnosis
  if ("primary_diagnosis" %in% colnames(meta)) {
    meta$who_subtype <- dplyr::case_when(
      grepl("type A,|type A$|thymoma, type A\\b", meta$primary_diagnosis, ignore.case=TRUE) ~ "A",
      grepl("type AB", meta$primary_diagnosis, ignore.case=TRUE) ~ "AB",
      grepl("type B1", meta$primary_diagnosis, ignore.case=TRUE) ~ "B1",
      grepl("type B2", meta$primary_diagnosis, ignore.case=TRUE) ~ "B2",
      grepl("type B3", meta$primary_diagnosis, ignore.case=TRUE) ~ "B3",
      grepl("carcinoma|NOS|squamous", meta$primary_diagnosis, ignore.case=TRUE) ~ "TC",
      TRUE ~ "Unknown"
    )
    cat("  WHO subtype inferred from primary_diagnosis.\n")
  } else {
    meta$who_subtype <- "Unknown"
    cat("  WARNING: No subtype info.\n")
  }
  meta
}
meta <- assign_who(meta, subtypes)
cat("  Subtype counts:\n"); print(table(meta$who_subtype, useNA="ifany")); cat("\n")

# =============================================================================
# FIGURE A: BCL2 expression by WHO subtype
# =============================================================================
cat("[Fig A] BCL2 by subtype...\n")
bcl2_expr <- data.frame(sample_id=colnames(bcl2_tpm),
                        BCL2_TPM=as.numeric(bcl2_tpm["BCL2",])) %>%
  left_join(meta %>% select(sample_id, who_subtype), by="sample_id") %>%
  filter(!is.na(who_subtype), who_subtype != "Unknown") %>%
  mutate(log2_TPM=log2(BCL2_TPM+1),
         who_subtype=factor(who_subtype, levels=c("A","AB","B1","B2","B3","TC")))

kw <- kruskal.test(log2_TPM ~ who_subtype, data=bcl2_expr)
cat("  Kruskal-Wallis p:", format(kw$p.value, digits=3), "\n")

comparisons_tc <- lapply(setdiff(levels(bcl2_expr$who_subtype), "TC"),
                         function(x) c(x, "TC"))

fig_a <- ggplot(bcl2_expr, aes(x=who_subtype, y=log2_TPM, fill=who_subtype)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.6) +
  geom_jitter(width=0.15, size=1.5, alpha=0.6) +
  scale_fill_manual(values=WHO_COLORS, guide="none") +
  stat_compare_means(comparisons=comparisons_tc, method="wilcox.test",
                     label="p.signif", step.increase=0.08, tip.length=0.01) +
  stat_compare_means(label.y=max(bcl2_expr$log2_TPM,na.rm=TRUE)*1.4,
                     method="kruskal.test", size=3.5) +
  labs(title="BCL2 mRNA Expression in Thymic Epithelial Tumors (TCGA-THYM)",
       x="WHO Histological Subtype", y=expression(BCL2~mRNA~(log[2]~TPM+1)),
       caption="Data: TCGA-THYM (Radovich et al., Cancer Cell 2018)") +
  theme_pub

ggsave(file.path(FIG_DIR,"FigA_BCL2_by_subtype.pdf"), fig_a, width=8, height=6, dpi=300)
ggsave(file.path(FIG_DIR,"FigA_BCL2_by_subtype.png"), fig_a, width=8, height=6, dpi=300)

bcl2_summary <- bcl2_expr %>% group_by(who_subtype) %>%
  summarise(n=n(), mean_TPM=mean(BCL2_TPM), median_TPM=median(BCL2_TPM),
            sd_TPM=sd(BCL2_TPM), .groups="drop")
write.csv(bcl2_summary, file.path(TABLE_DIR,"BCL2_expression_summary.csv"), row.names=FALSE)
cat("  Fig A saved.\n\n")

# =============================================================================
# FIGURE B: BCL2 Copy Number
# =============================================================================
cat("[Fig B] BCL2 copy number...\n")
fig_b_done <- FALSE

if (!is.null(cnv_data)) {
  bcl2_cn <- tryCatch({
    if (is(cnv_data, "SummarizedExperiment")) {
      cn_mat <- assay(cnv_data)
      cn_genes <- rowData(cnv_data)$gene_name
      idx <- which(cn_genes == "BCL2")
      if (length(idx) > 0) {
        data.frame(sample_id=colnames(cn_mat), BCL2_CN=as.numeric(cn_mat[idx[1],]))
      } else NULL
    } else if (is.data.frame(cnv_data) || is(cnv_data, "data.frame")) {
      # Segmented data - need to check BCL2 locus (chr18:63,123,346-63,320,128 hg38)
      bcl2_seg <- cnv_data %>%
        filter(Chromosome == "chr18" | Chromosome == "18") %>%
        filter(Start <= 63320128 & End >= 63123346)
      if (nrow(bcl2_seg) > 0) {
        bcl2_seg %>% select(Sample, Segment_Mean) %>%
          rename(sample_id=Sample, BCL2_CN_seg=Segment_Mean)
      } else NULL
    } else NULL
  }, error=function(e) { cat("  CN extraction error:", conditionMessage(e), "\n"); NULL })

  if (!is.null(bcl2_cn) && nrow(bcl2_cn) > 0) {
    bcl2_cn$patient_id <- substr(bcl2_cn$sample_id, 1, 12)
    bcl2_cn <- bcl2_cn %>%
      left_join(meta %>% select(patient_id, who_subtype) %>% distinct(), by="patient_id") %>%
      filter(!is.na(who_subtype), who_subtype != "Unknown")

    if ("BCL2_CN" %in% colnames(bcl2_cn)) {
      bcl2_cn$CN_status <- factor(dplyr::case_when(
        bcl2_cn$BCL2_CN <= -2 ~ "Deep Deletion",
        bcl2_cn$BCL2_CN == -1 ~ "Shallow Deletion",
        bcl2_cn$BCL2_CN == 0  ~ "Diploid",
        bcl2_cn$BCL2_CN == 1  ~ "Gain",
        bcl2_cn$BCL2_CN >= 2  ~ "Amplification"),
        levels=c("Deep Deletion","Shallow Deletion","Diploid","Gain","Amplification"))
    } else if ("BCL2_CN_seg" %in% colnames(bcl2_cn)) {
      bcl2_cn$CN_status <- factor(dplyr::case_when(
        bcl2_cn$BCL2_CN_seg < -0.5 ~ "Loss",
        bcl2_cn$BCL2_CN_seg > 0.3  ~ "Gain/Amp",
        TRUE ~ "Neutral"),
        levels=c("Loss","Neutral","Gain/Amp"))
    }

    cn_freq <- bcl2_cn %>% group_by(who_subtype, CN_status) %>%
      summarise(count=n(), .groups="drop") %>%
      group_by(who_subtype) %>% mutate(freq=count/sum(count)*100)

    cn_colors <- c("Deep Deletion"="#2166AC", "Shallow Deletion"="#92C5DE",
                   "Diploid"="#F7F7F7", "Gain"="#F4A582", "Amplification"="#B2182B",
                   "Loss"="#2166AC", "Neutral"="#F7F7F7", "Gain/Amp"="#B2182B")

    fig_b <- ggplot(cn_freq, aes(x=who_subtype, y=freq, fill=CN_status)) +
      geom_bar(stat="identity", color="grey30", linewidth=0.3) +
      scale_fill_manual(values=cn_colors, name="BCL2 CN") +
      labs(title="BCL2 Copy Number Alterations (TCGA-THYM)",
           x="WHO Subtype", y="Frequency (%)") + theme_pub

    ggsave(file.path(FIG_DIR,"FigB_BCL2_copynumber.pdf"), fig_b, width=8, height=6, dpi=300)
    ggsave(file.path(FIG_DIR,"FigB_BCL2_copynumber.png"), fig_b, width=8, height=6, dpi=300)

    cn_summary <- bcl2_cn %>% group_by(who_subtype, CN_status) %>%
      summarise(n=n(), .groups="drop")
    write.csv(cn_summary, file.path(TABLE_DIR,"BCL2_CN_summary.csv"), row.names=FALSE)
    cat("  Fig B saved.\n\n")
    fig_b_done <- TRUE
  }
}
if (!fig_b_done) cat("  CNV data unavailable or BCL2 not found. Skipped.\n\n")

# =============================================================================
# FIGURE C: Kaplan-Meier Survival
# =============================================================================
cat("[Fig C] Survival analysis...\n")

surv_data <- data.frame(sample_id=colnames(bcl2_tpm),
                        patient_id=substr(colnames(bcl2_tpm),1,12),
                        BCL2_TPM=as.numeric(bcl2_tpm["BCL2",])) %>%
  left_join(clinical, by=c("patient_id"="submitter_id")) %>%
  left_join(meta %>% select(patient_id, who_subtype) %>% distinct(), by="patient_id") %>%
  mutate(
    OS_status = ifelse(vital_status=="Dead", 1, ifelse(vital_status=="Alive", 0, NA)),
    OS_time = dplyr::coalesce(as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    OS_months = OS_time / 30.44
  )

bcl2_med <- median(surv_data$BCL2_TPM, na.rm=TRUE)
surv_data$BCL2_group <- factor(ifelse(surv_data$BCL2_TPM >= bcl2_med, "High", "Low"),
                               levels=c("Low","High"))
cat("  BCL2 median TPM:", round(bcl2_med, 2), "\n")

surv_c <- surv_data %>% filter(!is.na(OS_time), OS_time > 0, !is.na(OS_status), !is.na(BCL2_group))

if (nrow(surv_c) >= 10) {
  fit <- survfit(Surv(OS_months, OS_status) ~ BCL2_group, data=surv_c)
  lr <- survdiff(Surv(OS_months, OS_status) ~ BCL2_group, data=surv_c)
  lr_p <- 1 - pchisq(lr$chisq, 1)
  cat("  Log-rank p:", format(lr_p, digits=4), "\n")

  cox <- coxph(Surv(OS_months, OS_status) ~ BCL2_group, data=surv_c)
  cs <- summary(cox)
  cat("  HR:", round(cs$conf.int[1,1],2), "(95%CI:",
      round(cs$conf.int[1,3],2), "-", round(cs$conf.int[1,4],2), ")\n")

  fig_c <- ggsurvplot(fit, data=surv_c, pval=TRUE, pval.method=TRUE,
    risk.table=TRUE, risk.table.col="strata", risk.table.height=0.25,
    palette=c("#2166AC","#B2182B"),
    legend.labs=c(paste0("BCL2-Low (n=",sum(surv_c$BCL2_group=="Low"),")"),
                  paste0("BCL2-High (n=",sum(surv_c$BCL2_group=="High"),")")),
    legend.title="BCL2 Expression", xlab="Time (months)",
    title="Overall Survival by BCL2 Expression (TCGA-THYM)",
    subtitle=paste0("Median split: TPM=",round(bcl2_med,1)),
    ggtheme=theme_pub, risk.table.y.text=FALSE, surv.median.line="hv")

  pdf(file.path(FIG_DIR,"FigC_BCL2_survival.pdf"), width=8, height=7)
  print(fig_c); dev.off()
  png(file.path(FIG_DIR,"FigC_BCL2_survival.png"), width=8, height=7, units="in", res=300)
  print(fig_c); dev.off()
  cat("  Fig C saved.\n\n")
} else {
  cat("  Insufficient survival data (n<10). Skipped.\n\n")
}

# =============================================================================
# FIGURE D: BCL2 Family Heatmap
# =============================================================================
cat("[Fig D] BCL2 family heatmap...\n")

hm_genes <- c("BCL2","MCL1","BCL2L1","BCL2L2","BCL2L11","BAX","BAK1",
              "BBC3","PMAIP1","BID","BAD","TP53","CDKN2A","DNMT3A","DNMT1","TET2")
avail <- intersect(hm_genes, rownames(bcl2_tpm))
cat("  Genes:", length(avail), "/", length(hm_genes), "\n")

if (length(avail) >= 4) {
  hm_log <- log2(bcl2_tpm[avail,] + 1)
  hm_z <- t(scale(t(hm_log)))
  hm_z[hm_z > 3] <- 3; hm_z[hm_z < -3] <- -3

  col_meta <- meta %>% filter(sample_id %in% colnames(hm_z)) %>%
    arrange(match(sample_id, colnames(hm_z)))
  hm_z <- hm_z[, col_meta$sample_id]
  who_v <- col_meta$who_subtype; who_v[is.na(who_v)] <- "Unknown"

  ord <- order(factor(who_v, levels=c("A","AB","B1","B2","B3","TC","Unknown")))
  hm_z <- hm_z[, ord]; who_v <- who_v[ord]

  ha <- HeatmapAnnotation(
    WHO = who_v,
    BCL2 = anno_barplot(log2(as.numeric(bcl2_tpm["BCL2",col_meta$sample_id[ord]])+1),
                        height=unit(2,"cm"),
                        gp=gpar(fill=ifelse(who_v=="TC","black","grey70"))),
    col=list(WHO=WHO_COLORS), annotation_name_side="left")

  func <- dplyr::case_when(
    avail %in% c("BCL2","MCL1","BCL2L1","BCL2L2") ~ "Anti-apoptotic",
    avail %in% c("BCL2L11","BAX","BAK1","BBC3","PMAIP1","BID","BAD") ~ "Pro-apoptotic",
    avail %in% c("TP53","CDKN2A") ~ "Tumor Suppressor",
    avail %in% c("DNMT3A","DNMT1","TET2") ~ "Epigenetic", TRUE ~ "Other")
  ha_r <- rowAnnotation(Function=func,
    col=list(Function=c("Anti-apoptotic"="#B2182B","Pro-apoptotic"="#2166AC",
                        "Tumor Suppressor"="#4DAF4A","Epigenetic"="#FF7F00","Other"="grey70")))

  labels <- c(BCL2="BCL-2",MCL1="MCL-1",BCL2L1="BCL-xL",BCL2L2="BCL-W",
    BCL2L11="BIM",BAX="BAX",BAK1="BAK",BBC3="PUMA",PMAIP1="NOXA",BID="BID",
    BAD="BAD",TP53="TP53",CDKN2A="CDKN2A",DNMT3A="DNMT3A",DNMT1="DNMT1",TET2="TET2")
  rl <- labels[avail]; rl[is.na(rl)] <- avail[is.na(rl)]

  hm <- Heatmap(hm_z, name="Z-score",
    col=colorRamp2(c(-3,0,3), c("#2166AC","white","#B2182B")),
    top_annotation=ha, right_annotation=ha_r, row_labels=rl,
    show_column_names=FALSE, cluster_columns=FALSE, cluster_rows=TRUE,
    row_names_gp=gpar(fontsize=10, fontface="bold"),
    column_title="BCL-2 Family Expression in Thymic Epithelial Tumors (TCGA-THYM)",
    column_title_gp=gpar(fontsize=14, fontface="bold"))

  pdf(file.path(FIG_DIR,"FigD_BCL2_family_heatmap.pdf"), width=14, height=8)
  draw(hm, padding=unit(c(10,10,10,10),"mm")); dev.off()
  png(file.path(FIG_DIR,"FigD_BCL2_family_heatmap.png"), width=14, height=8, units="in", res=300)
  draw(hm, padding=unit(c(10,10,10,10),"mm")); dev.off()
  cat("  Fig D saved.\n\n")
}

# --- Save results summary ---
results <- list(
  bcl2_median = bcl2_med,
  kruskal_p = kw$p.value,
  logrank_p = if(exists("lr_p")) lr_p else NA,
  cox_hr = if(exists("cs")) cs$conf.int[1,] else NA,
  n_samples = ncol(bcl2_tpm),
  n_TC = sum(meta$who_subtype == "TC", na.rm=TRUE),
  subtype_counts = table(meta$who_subtype)
)
saveRDS(results, file.path(DATA_DIR,"analysis_results.rds"))

cat("\n=== COMPLETE ===\n")
cat("Figures:", paste(list.files(FIG_DIR), collapse=", "), "\n")
cat("Tables:", paste(list.files(TABLE_DIR), collapse=", "), "\n")
SCRIPT_END

# =============================================================================
# run_shirokane.sh (UGE)
# =============================================================================
cat > run_shirokane.sh << 'SCRIPT_END'
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N tcga_bcl2
#$ -o logs/tcga_bcl2_$JOB_ID.out
#$ -e logs/tcga_bcl2_$JOB_ID.err
#$ -l s_vmem=32G,mem_req=32G
#$ -pe def_slot 4
#$ -l d_rt=04:00:00,s_rt=04:00:00

set -euo pipefail
cd "${HOME}/tcga_thym_bcl2"
mkdir -p logs data figures tables

echo "=== TCGA-THYM BCL2 Pipeline ==="
echo "Date: $(date) | Node: $(hostname) | Job: ${JOB_ID:-local}"

module use /usr/local/package/modulefiles
module load R 2>/dev/null || module load R/4.3.2 2>/dev/null || true
echo "R: $(which R) ($(R --version 2>/dev/null | head -1))"

echo "[0] Packages..."
Rscript 00_install_packages.R 2>&1 | tail -5
echo "[1] Download..."
Rscript 01_download_data.R 2>&1
echo "[2] Analysis..."
Rscript 02_analysis_visualization.R 2>&1

echo ""
echo "=== Done: $(date) ==="
ls -la figures/ tables/
SCRIPT_END

chmod +x check_env.sh run_shirokane.sh

echo ""
echo "=== Setup complete ==="
echo "Files created in ${WORK_DIR}:"
ls -la *.R *.sh
echo ""
echo "Next steps:"
echo "  1. bash check_env.sh"
echo "  2. Rscript 00_install_packages.R  (if packages missing)"
echo "  3. qsub run_shirokane.sh"
