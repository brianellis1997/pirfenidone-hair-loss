# Report figure captions (draft, Oct 10)

Every figure is built by `analysis/15_report_figures.py` from `results/`. Each one comes as a PDF (vector, for the ACM
report) and a PNG (300 dpi, for slides).

**Figure 1. Workflow.** Four questions, the public data behind each one, and what came out. The negative results are
reported alongside.

**Figure 2. Signature screening fails its own calibration.** Drugs were ranked by how well they reverse the bald-scalp
or DHT gene signature. Known hair-growth drugs do not rank above known hair-loss drugs: AUROC 0.51 to 0.66, and two
signatures are inverted. This is why the project changed engines.

**Figure 3. The scarring program survives a follicle-loss adjustment.** The data are 14 bald vs 14 healthy scalp
biopsies (GSE90594). A gray dot is a gene's change before adjusting for how much follicle tissue each biopsy holds; a
colored dot is the change after. Filled means still significant. 15 of 21 scarring genes stay significant and about 90%
of the effect remains. Immune genes hold up less well (6 of 11, about 75%).

**Figure 4. Pirfenidone in human fibroblasts.**
- (a) At 1 mM, pirfenidone lowers 48% of the 23 TGF-beta-driven bald-scalp genes. At 0.11 mM, about the pill's blood
  level, it lowers 4%. Nintedanib and saracatinib lower 9%.
- (b) Pirfenidone at 1 mM cuts CXCL12 3.8-fold. CXCL12 is the signal whose blockade regrew hair in a 2025 mouse study.

**Figure 5. Where the program lives in human scalp** (Ober-Reynolds 2023 single-cell atlas, 41 annotated cell
clusters). The TGF-beta-driven genes, the subset pirfenidone lowers, and the scarring genes are concentrated in the
follicle's own connective tissue (sheath and papilla) and nearby dermal fibroblasts. Genes with strong hair-loss genetics
peak in the papilla. Genes lost in bald scalp sit in follicle epithelium, as a sanity check. * marks q < 0.05 against
expression-matched random gene sets.

**Figure 6. FDA adverse event reports (20.7M reports).**
- (a) The method picks up known hair-growth drugs (minoxidil, diazoxide, prostaglandins) and known hair-loss drugs
  (chemotherapy, leflunomide).
- (b) Lung-fibrosis patients on pirfenidone report hair loss about 70% less often than patients on nintedanib, the
  other drug for the same disease. This holds in men, women and ages 60 to 85. These are voluntary reports, not a trial.

**Figure 7. Genetic natural experiment.** The question is whether the DNA variants that change a gene's activity in skin
or fibroblasts (GTEx) are the same ones that change baldness risk (UK Biobank, 205k men). For the strongest known
baldness genes they often are (10 of 23 testable pairs). For the TGF-beta and scarring genes, none are (0 of 8). We
report this as a negative result.

**Figure 8. Blind audits of the LLM evidence agent.**
- (a) After two rounds of fixes, the pipeline labels 89% of unseen papers correctly. GPT-4o-mini
  alone manages 64%.
- (b) When the two reviewer models agree they are right 95% of the time; when they disagree, 50%. Disagreement is the
  cue for human review.
- The reference labels come from Claude Opus 5.5 reading each abstract without seeing any model's label. They are not
  human expert labels.

**Figure 9. Dose.**
- Cells needed about 1 mM, but the pill already slows lung fibrosis at 0.035 to 0.068 mM, so cell culture overstates the
  bar.
- A skin gel averages a similar level across the whole dermis.
- Tissue right at an application site reaches tens of mM (rabbit eye).
- The proposed lab test therefore includes a dose-response.
