# TCGA-THYM BCL2 Family Analysis Pipeline

## Purpose

Supplementary analysis for a thymic carcinoma case report (CR achieved with venetoclax + 5-AZA).
This pipeline analyzes BCL2 family expression, copy number, and survival data from the TCGA-THYM cohort.

## Generated Figures

| Figure | Description | Filename |
|--------|-------------|----------|
| **Fig A** | BCL2 mRNA expression box plot by WHO histological subtype | `FigA_BCL2_expression_by_subtype.pdf` |
| **Fig B** | BCL2 copy number alteration frequency by subtype | `FigB_BCL2_copynumber_by_subtype.pdf` |
| **Fig C** | Kaplan-Meier survival curve by BCL2 high/low expression | `FigC_BCL2_survival_KM.pdf` |
| **Fig D** | Heatmap of 16 BCL2 family genes | `FigD_BCL2_family_heatmap.pdf` |

## Directory Structure

```
tcga_thym_bcl2/
├── 00_install_packages.R           # Package installation
├── 01_download_data.R              # Data download via TCGAbiolinks
├── 01b_fix_download.R              # Download fix script
├── 02_analysis_v2.R                # Main analysis (v2)
├── 02_analysis_visualization.R     # Analysis & visualization
├── 02b_cbioportal_alternative.R    # cBioPortal alternative (WHO subtype)
├── run_shirokane.sh                # Job submission script
├── setup_all.sh                    # Full setup script
├── README.md
├── data/                           # Downloaded data (auto-generated)
├── figures/                        # Output figures (auto-generated)
├── tables/                         # Summary statistics CSV (auto-generated)
└── logs/                           # Job logs (auto-generated)
```

## Usage

### Running on SHIROKANE

```bash
cd ~/tcga_thym_bcl2
bash setup_all.sh        # Install packages + download data
bash run_shirokane.sh    # Submit analysis job
```

### Local Execution

```bash
cd tcga_thym_bcl2
Rscript 00_install_packages.R
Rscript 01_download_data.R
Rscript 02_analysis_visualization.R
```

## Target Genes

### BCL-2 Family (Apoptosis-related)
| Gene | Protein | Function | Relevance |
|------|---------|----------|-----------|
| `BCL2` | BCL-2 | Anti-apoptotic | **Venetoclax target** |
| `MCL1` | MCL-1 | Anti-apoptotic | Major survival factor in TETs |
| `BCL2L1` | BCL-xL | Anti-apoptotic | Major survival factor in TETs |
| `BCL2L2` | BCL-W | Anti-apoptotic | - |
| `BCL2L11` | BIM | Pro-apoptotic | Escape upon MCL-1 inhibition |
| `BAX` | BAX | Effector | Apoptosis execution |
| `BAK1` | BAK | Effector | Apoptosis execution |
| `BBC3` | PUMA | Pro-apoptotic | p53 downstream |
| `PMAIP1` | NOXA | Pro-apoptotic | MCL-1 antagonist |
| `BID` | BID | Pro-apoptotic | - |
| `BAD` | BAD | Pro-apoptotic | BCL-2/BCL-xL antagonist |

### Other Key Genes
| Gene | Relevance |
|------|-----------|
| `TP53` | Apoptosis / cell cycle |
| `CDKN2A` | Frequently deleted in thymic carcinoma (p16) |
| `DNMT3A` | **5-azacytidine target** (DNA methyltransferase) |
| `DNMT1` | **5-azacytidine target** (DNA methyltransferase) |
| `TET2` | Epigenetic regulation |

## TCGA-THYM Dataset

- **Publication**: Radovich et al., "The Integrated Genomic Landscape of Thymic Epithelial Tumors", Cancer Cell 2018
- **Samples**: 117 TETs (Type A: ~11, AB: ~39, B1: ~10, B2: ~18, B3: ~17, TC: ~11)
- **Data types**: RNA-seq, SNP array (CN), somatic mutations, DNA methylation
- **GDC Project ID**: TCGA-THYM

## References

1. Radovich M, et al. The Integrated Genomic Landscape of Thymic Epithelial Tumors. Cancer Cell. 2018;33(2):244-258.e10.
2. Petrini I, et al. Copy number aberrations of BCL2 and CDKN2A/B identified by array-CGH in thymic epithelial tumors. Cell Death Dis. 2012;3(7):e351.
3. Mueller D, et al. Functional apoptosis profiling identifies MCL-1 and BCL-xL as prognostic markers and therapeutic targets in advanced thymomas and thymic carcinomas. BMC Med. 2021;19(1):300.

## Requirements

- R >= 4.2.0
- Bioconductor >= 3.17
- Memory: 16 GB or more recommended
- Disk: ~5 GB (for TCGA data download)
- Internet connection (for data download)
