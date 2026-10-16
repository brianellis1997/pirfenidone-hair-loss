#!/bin/zsh
# Fetches every public dataset this project uses into data/. Safe to re-run; skips files that exist.
set -e
cd "$(dirname "$0")/.."
mkdir -p data/geo data/drugs data/ref
G=https://ftp.ncbi.nlm.nih.gov/geo/series
P=https://ftp.ncbi.nlm.nih.gov/geo/platforms
fetch() { [ -s "$2" ] || curl -s -L -f -o "$2" "$1"; }
cd data/geo
fetch $G/GSE36nnn/GSE36169/matrix/GSE36169_series_matrix.txt.gz GSE36169_series_matrix.txt.gz
fetch $G/GSE90nnn/GSE90594/matrix/GSE90594_series_matrix.txt.gz GSE90594_series_matrix.txt.gz
fetch $G/GSE78nnn/GSE78722/matrix/GSE78722_series_matrix.txt.gz GSE78722_series_matrix.txt.gz
fetch $G/GSE66nnn/GSE66664/matrix/GSE66664_series_matrix.txt.gz GSE66664_series_matrix.txt.gz
fetch $G/GSE267nnn/GSE267664/suppl/GSE267664_DHT-vs-control-all.gene.xls.gz GSE267664_DHT-vs-control-all.gene.xls.gz
fetch $G/GSE184nnn/GSE184501/suppl/GSE184501_2D_cuffnorm_fpkm_quantilenormalised.txt.gz GSE184501_2D_cuffnorm_fpkm_quantilenormalised.txt.gz
fetch $G/GSE178nnn/GSE178518/suppl/GSE178518_ipf_voom_normalized.txt.gz GSE178518_ipf_voom_normalized.txt.gz
fetch $G/GSE226nnn/GSE226249/suppl/GSE226249_counts.tsv.gz GSE226249_counts.tsv.gz
fetch $G/GSE189nnn/GSE189684/suppl/GSE189684_RAW.tar GSE189684_RAW.tar && mkdir -p GSE189684 && tar -xf GSE189684_RAW.tar -C GSE189684
fetch $G/GSE178nnn/GSE178374/suppl/GSE178374_RAW.tar GSE178374_RAW.tar && mkdir -p GSE178374 && tar -xf GSE178374_RAW.tar -C GSE178374
fetch $P/GPLnnn/GPL96/annot/GPL96.annot.gz GPL96.annot.gz
fetch $P/GPL10nnn/GPL10558/annot/GPL10558.annot.gz GPL10558.annot.gz
# GPL17077 family.soft is 3.7 GB (every sample ever run on the chip). The view=data table is 27 MB and has GENE_SYMBOL.
fetch "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17077&targ=self&form=text&view=data" GPL17077_table.txt
cd ../drugs
for lib in LINCS_L1000_Chem_Pert_Consensus_Sigs Old_CMAP_up Old_CMAP_down RummaGEO_DrugPerturbations_2025; do
  fetch "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=$lib" $lib.gmt
done
fetch https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt repurposing_drugs_20200324.txt
# SIDER needs https; the http URLs return a redirect page, not the file.
fetch https://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz meddra_all_se.tsv.gz
fetch https://sideeffects.embl.de/media/download/drug_names.tsv drug_names.tsv
cd ../ref
fetch https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt hgnc_complete_set.txt
mkdir -p ../scalp_atlas && cd ../scalp_atlas
# Ober-Reynolds 2023 scalp atlas: the authors' shared Drive folder (github.com/GreenleafLab/scScalpChromatin).
D="https://drive.usercontent.google.com/download?export=download&confirm=t&id="
fetch "${D}1Mk2WBU18BTGi1KPfJzhme5_AQDzVS6s8" pbFineClust_log2CP10kRNA.tsv
fetch "${D}1tLX1JEL07-r0OQeGY2gtFkBkAy5P1LDV" Supplemental_Tables.xlsx
mkdir -p ../gwas && cd ../gwas
# Yap et al. 2018 male-pattern baldness GWAS (GCST007020), 582 MB. Track J reads gene-activity data remotely (eQTL Catalogue tabix).
W=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007020
fetch $W/mpb_bolt_lmm_aut_x.tab.zip mpb_bolt_lmm_aut_x.tab.zip
fetch $W/MPB_GWAS_summary_statistics_README.pdf MPB_GWAS_summary_statistics_README.pdf
echo "all data present"
