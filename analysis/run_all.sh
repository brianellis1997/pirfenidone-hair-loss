#!/bin/zsh
# Rebuilds every result from scratch: data, then each analysis track in order. Logs go to results/logs/.
set -e
cd "$(dirname "$0")"
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null && conda activate 394_HW1 2>/dev/null || true
./00_download_data.sh
mkdir -p ../results/logs
for s in 01_expression_signatures 02_drug_signature_calibration 03_genetics_druggability 04_fibrosis_program 05_tgfb_fibroblasts_pirfenidone 06_faers 07_novelty_literature 08_scalp_cell_localization 09_llm_evidence_agent 10_genetic_natural_experiment 11_cxcl12_link 12_beard_vs_scalp 13_community_scan 14_dose_feasibility 15_report_figures; do
  echo "== $s"; python -u $s.py > ../results/logs/$s.log 2>&1 || { echo "FAILED: $s"; tail -20 ../results/logs/$s.log; exit 1; }; tail -3 ../results/logs/$s.log
done
