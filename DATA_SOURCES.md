# Data sources

Everything is public, and `analysis/00_download_data.sh` fetches all of it into `data/`. Accessed Oct 3 to 4, 2026.

## Gene expression (NCBI GEO)
| Accession | What it is | n | Used for |
|---|---|---|---|
| GSE36169 | Bald vs haired scalp, the same men (Garza et al., 2012, *Sci Transl Med*), Affymetrix HG-U133A (GPL96) | 5 pairs | In-vivo bald signature |
| GSE90594 | AGA vs healthy scalp biopsies (Michel et al., 2017, *Br J Dermatol*), Agilent G3 (GPL17077) | 14 vs 14 | In-vivo bald signature; composition check |
| GSE78722 | Frontal, vertex and occipital scalp, 6 men with early AGA, Agilent G3 | 18 | Replication attempt (did not replicate) |
| GSE267664 | DHT organ culture of AGA frontal follicles (authors' DESeq2 table) | 3 vs 3 | Tried; no signal |
| GSE178374 | DHT on primary dermal papilla cells from balding scalp, 2D and 3D (RNA-seq) | 12 | DHT-induced brakes |
| GSE66664 | Immortalized balding vs non-balding DPCs, DHT 1 or 10 nM, 0 to 48 h, Illumina HT-12 (GPL10558) | 140 | DHT time course |
| GSE184501 | Balding vs non-balding DPCs on 10 nM DHT, 2D (FPKM) | 6 | Intrinsic differences |
| GSE189684 | Beard, frontal, occipital and armpit follicles from one AGA patient set (TPM) | 4 | Paradox contrast (not yet used) |
| GSE178518 | Normal human lung fibroblasts ± TGF-beta, with pirfenidone (20 ug/mL), nintedanib or saracatinib (voom) | 30 | TGF-beta drive; drug effects |
| GSE226249 | Lung fibroblasts from 3 donors, serum ± pirfenidone 1 mM (counts) | 27 | Pirfenidone reversal |
| GSE212447/8 | Scalp single-cell RNA and ATAC atlas (Ober-Reynolds et al., 2023, *Nat Genet*): 6 healthy + 4 alopecia areata for RNA, mostly women | 10 + 13 | Cell localization (Track H). We use the authors' processed files, below, not the raw matrices |

The atlas's processed files come from the authors' shared Google Drive folder, linked from
github.com/GreenleafLab/scScalpChromatin:
- `pbFineClust_log2CP10kRNA.tsv`: pseudobulk expression, log2(CP10k + 1), for 43 annotated fine clusters (8 MB).
- `Supplemental_Tables.xlsx`: Tables S1 to S12. We use S4 (cluster names) and S12 (gkmSVM-prioritized hair-loss risk
  variants, with linked genes and motifs) (36 MB).

Platform maps:
- GPL96 and GPL10558: GEO `annot` files.
- GPL17077: the `view=data` table, not the 3.7 GB family file.

## Genetics and targets
- **Open Targets Platform GraphQL API**, disease MONDO_0005339 (androgenetic alopecia). Gives associated targets,
  datatype scores (genetic association pulls in GWAS such as UK Biobank), tractability, and drugs/clinical candidates.
- **HGNC complete set:** Ensembl to gene symbol.

## Drugs
- **LINCS L1000 Chem Pert Consensus Signatures**, via Enrichr: 5,425 drugs.
- **Original Connectivity Map** (Old_CMAP up/down), via Enrichr: 1,308 drugs.
- **RummaGEO Drug Perturbations 2025**, via Enrichr.
- **Broad Drug Repurposing Hub** (2020-03-24): clinical phase, mechanism, targets.
- **SIDER 4.1** (https only): side effects, used for the calibration labels.

## Genetics of baldness and gene activity (Track J)
- **Baldness GWAS:** Yap et al., 2018, *Nat Commun*, 205,327 UK Biobank men. Balding pattern score, BOLT-LMM, autosomes + X.
  GWAS Catalog GCST007020, `mpb_bolt_lmm_aut_x.tab.zip` (582 MB). Positions are GRCh37; we match on rsID.
- **Gene-activity variants (eQTLs):** GTEx v8 through the eQTL Catalogue's tabix files, read remotely one gene region at a
  time and cached in `data/eqtl/`:
  - QTD000216: cultured skin fibroblasts, n = 483.
  - QTD000311: suprapubic skin (not sun-exposed), n = 517.
  - QTD000316: lower-leg skin (sun-exposed), n = 602.
  The REST API was retired; the file paths are listed in eQTL-Catalogue-resources `tabix/tabix_ftp_paths.tsv`.
- **Gene coordinates:** Ensembl REST (GRCh38).

## LLM evidence agent (Track I)
- **PubMed E-utilities:** abstracts for 17 claims (5 hypothesis links, plus 2 claims each for 6 drugs). Frozen on Oct 4 in
  `data/llm_agent/pubmed_frozen.json`.
- **OpenRouter:** Llama 3.3 70B and DeepSeek V3.1 as independent reviewers, Gemini 2.5 Flash as the tiebreaker, and
  GPT-4o-mini as a shadow with no vote, for comparison. Every reply
  is cached in `data/llm_agent/llm_calls.jsonl`, so reruns are exact and free. A fresh run needs `OPENROUTER_API_KEY`.

## Real-world and literature
- **openFDA drug adverse events (FAERS):** reporting odds ratios and the active comparator.
- **PubMed E-utilities** and **ClinicalTrials.gov API v2:** novelty counts.
- **Enrichr API:** MSigDB Hallmark 2020 and Ligand Perturbations from GEO, for enrichment.
