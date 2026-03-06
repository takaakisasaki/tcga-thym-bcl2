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
