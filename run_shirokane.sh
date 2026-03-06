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
