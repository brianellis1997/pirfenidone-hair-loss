"""Builds HRP_Code_Walkthrough.ipynb: a guided tour of the 15 analysis scripts, with each track's core code, saved
results and figure, plus three live recomputations. Execute it afterwards with nbconvert (see the bottom of this file)."""
import pathlib

import nbformat as nbf

ROOT = pathlib.Path(__file__).parent
OUT = ROOT / "HRP_Code_Walkthrough.ipynb"
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip()))


md(r"""
# High-risk Project code: testing pirfenidone for pattern hair loss with an audited LLM evidence agent
**Brian Ellis · A I 380H AI in Healthcare · Fall 2026**

### AI disclosure
I wrote and tested the code together with Claude (Opus 5.5), and it also helped me come up with ideas and draft this
notebook. The question and the calls on where to take it were mine, and I went through the results myself before
anything went in. The evidence agent itself runs on Llama 3.3 70B, DeepSeek V3.1, Gemini 2.5 Flash and GPT-4o-mini
through OpenRouter.

### What this notebook does
The analysis is 15 numbered scripts in `analysis/`, one per research track, plus shared loaders in `src/`. This notebook
walks through them in order. For each track it states the question, prints the core code straight from the script, and
shows the saved results and the report figure. Three results are recomputed live from the intermediate tables (the FDA
comparator ratio, the evidence agent's grades and the dose arithmetic), to show the saved numbers come from the code.

It doesn't re-query the live APIs. PubMed counts drift as papers are added, so a live re-query would no longer match the
report. To rebuild everything from public downloads, run `analysis/run_all.sh` (about 10 minutes, needs internet). The
evidence agent replays every LLM reply from its cache, so a rebuild needs no API key.

| Path | What's there |
|---|---|
| `analysis/00_download_data.sh` | Downloads all public data into `data/` |
| `analysis/01_*.py` to `15_*.py` | One script per track (A to N), then the report figures |
| `src/signatures.py`, `src/drug_scoring.py` | GEO loaders, differential expression, drug scoring |
| `data/llm_agent/` | The frozen PubMed set and every cached LLM reply |
| `results/` | Tables, plus `NN_summary.json` with each track's key numbers |
| `figs/report/` | The nine report figures |

Environment: conda env `394_HW1` (Python 3.11, numpy 1.26, pandas, scipy, statsmodels, scikit-learn, pysam, matplotlib).
""")

code(r"""
import ast, hashlib, json, pathlib, platform, runpy, sys
import numpy as np, pandas as pd, scipy
from IPython.display import Image, display

ROOT = pathlib.Path.cwd()
RESULTS, FIGS = ROOT / "results", ROOT / "figs" / "report"
sys.path.insert(0, str(ROOT / "analysis"))
pd.set_option("display.width", 170)
pd.set_option("display.max_colwidth", 70)


def show_code(path, names):
    src = (ROOT / path).read_text()
    for node in ast.parse(src).body:
        name = getattr(node, "name", None)
        if name is None and isinstance(node, ast.Assign):
            name = getattr(node.targets[0], "id", None)
        if name in names:
            print(ast.get_source_segment(src, node), end="\n\n")


def show_from(path, marker):
    lines = (ROOT / path).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(marker))
    print("\n".join(lines[start:]))


def summary(n):
    return json.load(open(RESULTS / f"{n}_summary.json"))


def fig(name, width=900):
    display(Image(filename=str(FIGS / f"{name}.png"), width=width))


print("python", platform.python_version(), "| numpy", np.__version__, "| pandas", pd.__version__, "| scipy", scipy.__version__)
""")

md(r"""
## A. Is there a reproducible bald-scalp signature?
Two independent studies compare bald with haired scalp: GSE36169 (5 men, paired biopsies) and GSE90594 (14 AGA vs 14
healthy). They are combined into one consensus signature (Stouffer z), keeping only genes that move the same way in both.
""")
code(r"""
show_code("analysis/common.py", ["consensus_bald_signature"])
s = summary("01")
a = s["agreement"]["GSE36169 vs GSE90594"]
print(f"Agreement between the two in-vivo studies: rho {a['rho_sig_either']:.2f}; {a['same_direction']} of {a['both_sig']} genes significant in both move the same way")
print(f"Consensus signature: {s['consensus_genes']:,} genes")
""")

md(r"""
## B. Does drug-signature screening pass its own calibration?
Every drug in LINCS L1000 and the original Connectivity Map is scored on how strongly it reverses the bald-scalp or
DHT-response signature. The check: do drugs known to grow hair (from SIDER side effects) rank above drugs known to cause
hair loss? It fails, so signature screening is reported as a negative result.
""")
code(r"""
show_code("src/drug_scoring.py", ["reversal_z", "calibration"])
cal = pd.read_csv(RESULTS / "02_calibration.csv").dropna(subset=["auroc"])
display(cal[["library", "signature", "drugs", "n_growth", "n_loss", "auroc", "p_perm"]].round(3))
fig("fig2_calibration", 480)
""")

md(r"""
## C. Which genes carry the genetic risk, and can any drug act on them?
Open Targets lists every gene with evidence for androgenetic alopecia; we count how many strong genetic hits have any drug.
""")
code(r"""
s = summary("03")
print(f"Associated targets: {s['associated_targets']:,}; genetic score >= 0.5: {s['genetic_ge_0.5']}; with any drug or candidate: {s['with_any_drug']}; with an approved drug: {s['with_approved_drug']}")
for t in s["approved_drug_targets"]:
    print(f"  {t['symbol']:7} {t['approved_drugs'][:90]}")
""")

md(r"""
## D. Is there a scarring program in bald scalp, and is it just lost follicles?
Bald-up genes are tested against the MSigDB Hallmark sets. Then, in GSE90594, each scarring and immune gene is refit with
a follicle-content score as a covariate, so a gene can't look raised just because the biopsy holds more dermis.
""")
code(r"""
show_from("analysis/04_fibrosis_program.py", "FOLLICLE = ")
s = summary("04")
print("\nTop Hallmark sets among bald-up genes:")
for r in s["hallmark_top"]["bald_up"][:3]:
    print(f"  {r['term']:38} adjusted p {r['adj_p']:.1e} ({r['n_genes']} genes)")
for k, v in s["composition_GSE90594"].items():
    print(f"{k}: {v['sig_adjusted']} of {v['genes']} genes stay significant after adjustment, median effect kept {v['median_effect_kept']:.0%}")
fig("fig3_scarring_program")
""")

md(r"""
## E. Is part of the program TGF-beta driven, and does pirfenidone switch it off?
TGF-beta-induced genes in lung fibroblasts (GSE178518) are intersected with bald-up genes. Then the drugs' effects on those
genes are measured: pirfenidone at 0.11 mM, nintedanib and saracatinib (GSE178518), and pirfenidone at 1 mM (GSE226249).
""")
code(r"""
print((ROOT / "analysis/05_tgfb_fibroblasts_pirfenidone.py").read_text())
s = summary("05")
e = s["tgfb_enrichment"]
print(f"TGF-beta-driven bald genes: {e['overlap']} of {e['bald_up_tested']} (odds ratio {e['odds_ratio']:.1f}, p {e['p']:.1e})")
for k, v in s["drugs"].items():
    print(f"  {k:28} lowers {v['driven_lowered_sig']:.0%} of them (median log2 change {v['driven_median_change']:+.2f})")
fig("fig4_pirfenidone_fibroblasts")
""")

md(r"""
## F. Real-world reports: FDA adverse events
Reporting odds ratios from openFDA, first on drugs with known hair effects (calibration), then pirfenidone against
nintedanib, the other drug for the same lung disease. **Live check:** the comparator ratio is recomputed below from the saved
report counts with the same formula the script uses.
""")
code(r"""
show_code("analysis/06_faers.py", ["ror"])
comp = pd.read_csv(RESULTS / "06_faers_pirfenidone_vs_nintedanib.csv")
r = (comp.pirf_alopecia / (comp.pirf_reports - comp.pirf_alopecia)) / (comp.nint_alopecia / (comp.nint_reports - comp.nint_alopecia))
se = np.sqrt(1 / comp.pirf_alopecia + 1 / (comp.pirf_reports - comp.pirf_alopecia) + 1 / comp.nint_alopecia + 1 / (comp.nint_reports - comp.nint_alopecia))
comp["recomputed"] = r
comp["lo95"], comp["hi95"] = np.exp(np.log(r) - 1.96 * se), np.exp(np.log(r) + 1.96 * se)
display(comp[["stratum", "pirf_alopecia", "pirf_reports", "nint_alopecia", "nint_reports", "ratio", "recomputed", "lo95", "hi95"]].round(3))
print("recomputed ratios match the saved ones:", bool(np.allclose(comp.ratio, comp.recomputed)))
fig("fig6_faers")
""")

md(r"""
## G. Has anyone tested this already?
PubMed and ClinicalTrials.gov counts per candidate drug and per mechanism (saved on the day they were queried).
""")
code(r"""
s = summary("07")
display(pd.DataFrame(s["drugs"])[["group", "drug", "pattern_hair_loss_papers", "hair_papers", "all_papers"]])
print(s["mechanism_papers"])
""")

md(r"""
## H. Which scalp cells carry the program?
Using the authors' pseudobulk expression for 41 annotated clusters of the Ober-Reynolds 2023 scalp atlas, each gene set
is scored per cell type against 10,000 expression-matched random gene sets.
""")
code(r"""
show_code("analysis/08_scalp_cell_localization.py", ["module_scores"])
s = summary("08")
for name in ["TGF-beta-driven bald genes", "Pirfenidone-lowered (1 mM)", "Hair-loss genetics, all strong (>=0.5)"]:
    top = ", ".join(f"{c['cluster']} {c['score']:+.2f} (q {c['q']:.3f})" for c in s["sets"][name]["top_clusters"][:4])
    print(f"{name}: {top}")
fig("fig5_scalp_localization")
""")

md(r"""
## I. The LLM evidence agent
Two reviewer models label every claim-abstract pair independently (stance, evidence type, verbatim quote); a tiebreaker
decides disagreements; grades come from rules fixed in code. **Live check:** the grades are recomputed below from the saved
labels with the script's own `grade` function.
""")
code(r"""
show_code("analysis/09_llm_evidence_agent.py", ["REVIEWERS", "ADJUDICATOR", "SHADOW", "grade"])
src = (ROOT / "analysis/09_llm_evidence_agent.py").read_text()
node = next(n for n in ast.parse(src).body if getattr(n, "name", None) == "grade")
ns = {}
exec(ast.get_source_segment(src, node), ns)
labels = pd.read_csv(RESULTS / "09_evidence_labels.csv")
saved = pd.read_csv(RESULTS / "09_claim_grades.csv", keep_default_na=False).set_index("claim_id")
live = pd.DataFrame({cid: ns["grade"](rows) for cid, rows in labels.groupby("claim_id")}).T
check = saved[["grade", "supporting_primary", "supporting_human", "contradicting"]].join(live[["grade"]], rsuffix="_live")
display(check)
print("live grades match the saved ones:", bool((check.grade == check.grade_live).all()))
""")
code(r"""
a = summary("09")["audit"]
acc = pd.DataFrame(a["accuracy"]).T[["population_weighted", "agree_stratum", "disagree_stratum", "clear_items_only"]]
display(acc.round(3))
print("held-out pairs:", a["n"], "| reference labels:", a["auditor"])
fig("fig8_agent_audit")
""")

md(r"""
## J. A genetic natural experiment
Colocalization between the UK Biobank baldness GWAS and gene-activity variants (GTEx v8 skin and fibroblasts): does one
causal variant drive both? Reference genes with strong baldness genetics show the method can find a signal.
""")
code(r"""
show_code("analysis/10_genetic_natural_experiment.py", ["labf", "coloc"])
p = summary("10")["panels"]
for k in ("hypothesis", "reference"):
    print(f"{k}: {p[k]['shared_variant_H4']} of {p[k]['testable']} testable gene-tissue pairs share a causal variant")
print("Fisher p, hypothesis vs reference:", round(p["fisher_hypothesis_vs_reference"]["p"], 3))
fig("fig7_genetics", 420)
""")

md(r"""
## K and L. The CXCL12 link, and beard vs scalp
""")
code(r"""
c = summary("11")
b, f = c["bald_scalp"]["CXCL12"], c["fibroblast_drugs"]["CXCL12"]
print(f"CXCL12 in bald scalp: consensus z {b['z']:.2f} ({b['z_rank_pct']:.1f}th percentile)")
print(f"Pirfenidone 1 mM: log2 {f['pirf_1mM_lfc']:+.2f} (p {f['pirf_1mM_p']:.4f}); 0.11 mM: {f['pirf_0p1mM_lfc']:+.3f}; nintedanib: {f['nint_lfc']:+.2f}")
print(f"Share of atlas expression in dermal clusters: {c['atlas_share_in_dermal_clusters']['CXCL12']:.0%}")
l = pd.DataFrame(summary("12")["sets"])
display(l[l.gene_set.isin(["TGF-beta-driven bald genes", "Scarring genes up in bald scalp"])][["gene_set", "contrast", "mean_log2_diff", "perm_p_two_sided", "share_of_genes_up"]].round(3))
""")

md(r"""
## M and N. The community scan, and dose (**live**: the dose arithmetic reruns here)
""")
code(r"""
m = summary("13")
print(f"YouTube scan ({m['scanned_on']}): {m['videos']} videos, {m['videos_on_pirfenidone']} on pirfenidone; earliest: {m['earliest_pirfenidone_video']}")
print()
_ = runpy.run_path(str(ROOT / "analysis/14_dose_feasibility.py"), run_name="__main__")
fig("fig9_dose", 520)
""")

md(r"""
## Reproducibility
Fingerprints of every summary file the report quotes. On Oct 4, full rebuilds with `analysis/run_all.sh` reproduced the
summaries for tracks 01 to 12 byte for byte, apart from PubMed counts that drift as new papers are indexed.
""")
code(r"""
rows = [(p.name, hashlib.sha256(p.read_bytes()).hexdigest()[:12]) for p in sorted(RESULTS.glob("*_summary.json"))]
display(pd.DataFrame(rows, columns=["file", "sha256 (first 12)"]))
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, OUT)
print("wrote", OUT, len(cells), "cells")
# Execute: jupyter nbconvert --to notebook --execute --inplace HRP_Code_Walkthrough.ipynb
