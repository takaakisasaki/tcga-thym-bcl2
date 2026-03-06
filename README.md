# TCGA-THYM BCL2 Family Analysis Pipeline

## Purpose

胸腺がん症例報告（Venetoclax + 5-AZA によるCR達成）の補足データとして、
TCGA-THYMデータベースからBCL2ファミリーの発現・コピー数・予後データを解析する。

## 生成されるFigure

| Figure | 内容 | ファイル名 |
|--------|------|------------|
| **Fig A** | BCL2 mRNA発現のWHO組織型別ボックスプロット | `FigA_BCL2_expression_by_subtype.pdf` |
| **Fig B** | BCL2コピー数変異の組織型別頻度 | `FigB_BCL2_copynumber_by_subtype.pdf` |
| **Fig C** | BCL2高/低発現のKaplan-Meier生存曲線 | `FigC_BCL2_survival_KM.pdf` |
| **Fig D** | BCL2ファミリー16遺伝子のヒートマップ | `FigD_BCL2_family_heatmap.pdf` |

## ディレクトリ構造

```
tcga_thym_bcl2/
├── 00_install_packages.R           # パッケージインストール
├── 01_download_data.R              # TCGAbiolinksでデータダウンロード
├── 02_analysis_visualization.R     # メイン解析・可視化
├── 02b_cbioportal_alternative.R    # cBioPortal代替（WHO subtype取得）
├── run_shirokane.sh                # SLURM投入スクリプト
├── README.md                       # 本ファイル
├── data/                           # ダウンロードデータ（自動生成）
├── figures/                        # 出力Figure（自動生成）
├── tables/                         # 統計サマリーCSV（自動生成）
└── logs/                           # SLURMログ（自動生成）
```

## 使用方法

### SHIROKANE上での実行

```bash
# 1. ファイルをSHIROKANEに転送
scp -r tcga_thym_bcl2/ username@shirokane:/home/username/

# 2. SHIROKANEにログイン
ssh username@shirokane

# 3. ジョブ投入
cd ~/tcga_thym_bcl2
sbatch run_shirokane.sh

# 4. ジョブ確認
squeue -u $USER
cat logs/tcga_bcl2_*.out
```

### ローカル実行（テスト用）

```bash
cd tcga_thym_bcl2
Rscript 00_install_packages.R
Rscript 01_download_data.R
Rscript 02_analysis_visualization.R
```

## 解析対象遺伝子

### BCL-2ファミリー（アポトーシス関連）
| 遺伝子 | タンパク質 | 機能 | 本症例との関連 |
|---------|-----------|------|---------------|
| `BCL2` | BCL-2 | 抗アポトーシス | **Venetoclaxの標的** |
| `MCL1` | MCL-1 | 抗アポトーシス | TETの主要survival factor |
| `BCL2L1` | BCL-xL | 抗アポトーシス | TETの主要survival factor |
| `BCL2L2` | BCL-W | 抗アポトーシス | - |
| `BCL2L11` | BIM | アポトーシス促進 | MCL-1阻害時のescape |
| `BAX` | BAX | エフェクター | アポトーシス実行 |
| `BAK1` | BAK | エフェクター | アポトーシス実行 |
| `BBC3` | PUMA | アポトーシス促進 | p53下流 |
| `PMAIP1` | NOXA | アポトーシス促進 | MCL-1と拮抗 |
| `BID` | BID | アポトーシス促進 | - |
| `BAD` | BAD | アポトーシス促進 | BCL-2/BCL-xLと拮抗 |

### その他重要遺伝子
| 遺伝子 | 関連 |
|---------|------|
| `TP53` | アポトーシス/cell cycle |
| `CDKN2A` | 胸腺がんで欠失頻出（p16） |
| `DNMT3A` | **5-azacytidineの標的**（DNA methyltransferase） |
| `DNMT1` | **5-azacytidineの標的**（DNA methyltransferase） |
| `TET2` | エピジェネティック制御 |

## TCGA-THYMデータセット概要

- **論文**: Radovich et al., "The Integrated Genomic Landscape of Thymic Epithelial Tumors", Cancer Cell 2018
- **サンプル数**: 117 TETs
  - Type A: ~11, AB: ~39, B1: ~10, B2: ~18, B3: ~17, TC: ~11
- **データ種類**: RNA-seq, SNP array (CN), 体細胞変異, DNAメチレーション
- **GDC Project ID**: TCGA-THYM

## 論文での使用方法

### Discussionへの追記例

> To further validate BCL-2 expression in thymic carcinoma, we analyzed
> the TCGA-THYM cohort (n=117). BCL2 mRNA expression was significantly
> elevated in thymic carcinoma compared to type A and AB thymomas
> (Kruskal-Wallis test, p = X.XX; Supplementary Figure X). BCL2 copy
> number gain was observed in X% (X/11) of thymic carcinomas.
> Kaplan-Meier analysis suggested that higher BCL2 expression was
> associated with [shorter/comparable] overall survival in the TCGA
> cohort (log-rank p = X.XX; Supplementary Figure X).
>
> Interestingly, analysis of the BCL-2 family expression profile
> revealed that MCL-1 and BCL-xL were also highly expressed in thymic
> carcinoma (Supplementary Figure X), consistent with recent BH3
> profiling studies (Müller et al., BMC Med 2021). The synergistic
> effect of venetoclax and 5-azacytidine observed in our case may
> reflect epigenetic priming by 5-azacytidine that shifts the apoptotic
> balance, enhancing the cytotoxic effect of BCL-2 inhibition even in
> tumors with concurrent MCL-1/BCL-xL expression.

## 追加引用すべき文献

1. Radovich M, et al. The Integrated Genomic Landscape of Thymic
   Epithelial Tumors. Cancer Cell. 2018;33(2):244-258.e10.
   doi:10.1016/j.ccell.2018.01.003

2. Petrini I, et al. Copy number aberrations of BCL2 and CDKN2A/B
   identified by array-CGH in thymic epithelial tumors. Cell Death Dis.
   2012;3(7):e351. doi:10.1038/cddis.2012.92

3. Müller D, et al. Functional apoptosis profiling identifies MCL-1
   and BCL-xL as prognostic markers and therapeutic targets in advanced
   thymomas and thymic carcinomas. BMC Med. 2021;19(1):300.
   doi:10.1186/s12916-021-02158-3

## トラブルシューティング

### TCGAbiolinksのダウンロードが失敗する場合
- SHIROKANEのcompute nodeからインターネットにアクセスできるか確認
- login nodeで`01_download_data.R`を実行（qloginで対話セッション使用）
- またはcBioPortal代替スクリプト(`02b_cbioportal_alternative.R`)を使用

### WHO subtypeが取得できない場合
1. cBioPortalの[TCGA-THYM study](https://www.cbioportal.org/study/summary?id=thym_tcga)にアクセス
2. "Clinical Data" タブからCSVをダウンロード
3. `data/` ディレクトリに `THYM_subtypes_manual.csv` として保存
4. 02_analysis_visualization.R の subtype assignment 部分を修正

### ComplexHeatmapのフォント問題
```bash
# SHIROKANEでフォントが見つからない場合
export FONTCONFIG_PATH=/etc/fonts
```

## 環境要件

- R >= 4.2.0
- Bioconductor >= 3.17
- メモリ: 16GB以上推奨
- ディスク: ~5GB（TCGAデータダウンロード用）
- インターネット接続（データダウンロード時）
