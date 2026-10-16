"""Track H: which scalp cell types carry the TGF-beta / scarring program, the hair-loss genetics and the androgen receptor?
Ober-Reynolds 2023 scalp atlas (GSE212447): the authors' pseudobulk expression, log2(CP10k + 1), for 43 annotated cell
clusters from 10 scalp samples (6 healthy, 4 alopecia areata), plus their table of fine-mapped hair-loss risk variants."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from statsmodels.stats.multitest import multipletests

from common import RESULTS, ROOT, consensus_bald_signature, save_summary

ATLAS = ROOT / "data" / "scalp_atlas"
RNA_COMPARTMENTS = {
    "Dermal fibroblasts": ["D.Fib_1", "D.Fib_2", "D.Fib_3", "D.Fib_4", "D.Fib_5"],
    "Follicle dermis": ["D.Sheath", "D.Papilla"],
    "Muscle and pericytes": ["Muscle", "Pericytes"],
    "Vessels": ["Vas.Endo_1", "Vas.Endo_2", "Vas.Endo_3", "Vas.Endo_4", "Lymph.Endo_1"],
    "Epidermis": ["Basal.Kc_1", "Spinous.Kc_1", "Spinous.Kc_2", "Spinous.Kc_3", "Cyc.Kc_1"],
    "Follicle epithelium": ["Infundibulum", "Inf.Segment", "Isthmus"],
    "Glands": ["Sebaceous", "Eccrine"],
    "Melanocytes": ["Melanocytes_1", "Melanocytes_2"],
    "Immune": ["Plasma", "Mast", "cDC2", "M2.macs_1", "M2.macs_2", "CLEC9a.DC", "M1.macs", "TCR.macs", "TREM2.macs",
               "CD4.Tc_1", "CD4.Tc_2", "CD8_Tc_1", "Treg", "NK", "Cyc.Tc"],
}
ATAC_COMPARTMENT_PREFIX = {"D.Fib": "Dermal fibroblasts", "D.Sheath": "Follicle dermis", "D.Papilla": "Follicle dermis",
                           "Muscle": "Muscle and pericytes", "Pericytes": "Muscle and pericytes", "Vas.Endo": "Vessels",
                           "Lymph.Endo": "Vessels", "Unknown_2": "Vessels", "Basal.Kc": "Epidermis", "Spinous.Kc": "Epidermis",
                           "Infundibulum": "Follicle epithelium", "Inf.Segment": "Follicle epithelium", "Isthmus": "Follicle epithelium",
                           "Matrix": "Follicle epithelium", "Sebaceous": "Glands", "Eccrine": "Glands", "Unknown_1": "Epidermis",
                           "Melanocytes": "Melanocytes"}
KEY_GENES = ["AR", "SRD5A1", "SRD5A2", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD3", "TWIST1", "TWIST2", "PDGFA", "PDGFRA",
             "PDGFRB", "DKK1", "SFRP2", "IL6", "COL1A1", "COL3A1", "LOX", "SPARC", "TIMP1", "TGFBI", "ACTA2", "CORIN", "SOX2"]


def module_scores(z, genes, bins, n_perm=10000, seed=0):
    genes = [g for g in genes if g in z.index]
    obs = z.loc[genes].mean()
    rng = np.random.default_rng(seed)
    pools = {b: np.flatnonzero(bins.values == b) for b in np.unique(bins.values)}
    zv = z.values.astype(np.float32)
    exceed = np.zeros(z.shape[1])
    for chunk in np.array_split(np.arange(n_perm), 10):
        draws = np.stack([rng.choice(pools[bins[g]], len(chunk)) for g in genes], axis=1)
        exceed += (zv[draws].mean(axis=1) >= obs.values).sum(axis=0)
    p = (exceed + 1) / (n_perm + 1)
    return pd.DataFrame({"score": obs, "p": p, "q": multipletests(p, method="fdr_bh")[1]}), genes


pb = pd.read_csv(ATLAS / "pbFineClust_log2CP10kRNA.tsv", sep="\t", index_col=0).drop(columns=["Unknown", "Plasma.contam"])
compartment = pd.Series({c: k for k, v in RNA_COMPARTMENTS.items() for c in v})
assert set(compartment.index) == set(pb.columns), set(pb.columns) ^ set(compartment.index)
pb = pb[compartment.index]
expressed = pb[pb.max(axis=1) >= 1]
z = expressed.sub(expressed.mean(axis=1), axis=0).div(expressed.std(axis=1), axis=0)
bins = pd.Series(pd.qcut(expressed.mean(axis=1), 25, labels=False), index=expressed.index)
print(f"atlas: {pb.shape[1]} clusters, {len(expressed)} expressed genes (max log2CP10k >= 1)")

t5 = json.load(open(RESULTS / "05_summary.json"))
driven = t5["tgfb_enrichment"]["driven_genes"]
fx = pd.read_csv(RESULTS / "05_fibroblast_tgfb_drug_effects.csv", index_col=0)
pirf_lowered = [g for g in driven if g in fx.index and fx.loc[g, "pirf_1mM_lfc"] < 0 and fx.loc[g, "pirf_1mM_p"] < .05]
comp = pd.read_csv(RESULTS / "04_composition_check_GSE90594.csv")
fib_up = comp[(comp.group == "fibrosis") & (comp.bald_effect > 0) & (comp.p_raw < .05)].gene.tolist()
imm_up = comp[(comp.group == "immune") & (comp.bald_effect > 0) & (comp.p_raw < .05)].gene.tolist()
gfib = pd.read_csv(RESULTS / "03_fibrosis_gene_genetics.csv").dropna()
genetics_tgfb = gfib[gfib.genetic >= .3].symbol.tolist()
gen_all = pd.read_csv(RESULTS / "03_genetic_targets_druggability.csv").symbol.tolist()
inv = consensus_bald_signature()
bald_up = inv[inv.z > 0].nlargest(200, "z").index.tolist()
bald_down = inv[inv.z < 0].nsmallest(200, "z").index.tolist()
sets = {"TGF-beta-driven bald genes": driven, "Pirfenidone-lowered (1 mM)": pirf_lowered, "Scarring genes up in bald scalp": fib_up,
        "Immune genes up in bald scalp": imm_up, "Hair-loss genetics, TGF-beta/scarring": genetics_tgfb,
        "Hair-loss genetics, all strong (>=0.5)": gen_all, "Bald-up top 200": bald_up, "Bald-down top 200": bald_down}

table = pd.DataFrame({"compartment": compartment})
summary = {"clusters": int(pb.shape[1]), "expressed_genes": int(len(expressed)), "sets": {}}
for name, genes in sets.items():
    res, used = module_scores(z, genes, bins)
    table[f"{name} | score"] = res.score
    table[f"{name} | q"] = res.q
    peak = expressed.loc[used].idxmax(axis=1).map(compartment).value_counts()
    by_comp = res.score.groupby(compartment).mean().sort_values(ascending=False)
    top = res.sort_values("score", ascending=False).head(5)
    summary["sets"][name] = {"genes_in_atlas": f"{len(used)}/{len(genes)}", "top_clusters": [
        {"cluster": c, "score": round(r.score, 2), "q": float(f"{r.q:.2g}")} for c, r in top.iterrows()],
        "clusters_q_lt_0.05": res.index[res.q < .05].tolist(), "compartment_mean_score": by_comp.round(2).to_dict(),
        "peak_cluster_compartment_counts": peak.to_dict()}
    print(f"\n{name} ({len(used)}/{len(genes)} genes in atlas)")
    print("  top clusters: " + ", ".join(f"{c} {r.score:+.2f} (q {r.q:.1g})" for c, r in top.iterrows()))
    print("  genes peak in: " + ", ".join(f"{k} {v}" for k, v in peak.items()))
table.round(4).to_csv(RESULTS / "08_cluster_scores.csv")

summary["robustness"] = {}
for label, cols, cutoff in (("no immune clusters", compartment.index[compartment != "Immune"], 1.0),
                            ("looser cutoff (max >= 0.5)", compartment.index, 0.5)):
    e2 = pb.loc[pb[cols].max(axis=1) >= cutoff, cols]
    z2 = e2.sub(e2.mean(axis=1), axis=0).div(e2.std(axis=1), axis=0)
    b2 = pd.Series(pd.qcut(e2.mean(axis=1), 25, labels=False), index=e2.index)
    summary["robustness"][label] = {}
    for name in ("TGF-beta-driven bald genes", "Pirfenidone-lowered (1 mM)", "Hair-loss genetics, all strong (>=0.5)"):
        r2, _ = module_scores(z2, sets[name], b2)
        top = r2.sort_values("score", ascending=False).head(4)
        summary["robustness"][label][name] = [{"cluster": c, "score": round(r.score, 2), "q": float(f"{r.q:.2g}")} for c, r in top.iterrows()]
        print(f"[{label}] {name}: " + ", ".join(f"{c} {r.score:+.2f} (q {r.q:.1g})" for c, r in top.iterrows()))

key = pb.loc[[g for g in KEY_GENES if g in pb.index]].T
key.insert(0, "compartment", compartment)
key.round(3).to_csv(RESULTS / "08_key_genes_by_cluster.csv")
print("\nkey genes, top 3 clusters each:")
for g in key.columns[1:]:
    print(f"  {g:7} " + ", ".join(f"{c} {v:.1f}" for c, v in key[g].nlargest(3).items()))
ar_rho, ar_p = stats.spearmanr(key["AR"], table["TGF-beta-driven bald genes | score"])
summary["AR_top_clusters"] = key["AR"].nlargest(5).round(2).to_dict()
summary["AR_vs_tgfb_driven_module"] = {"spearman_rho": ar_rho, "p": ar_p}
summary["key_gene_peak_cluster"] = {g: key[g].idxmax() for g in key.columns[1:]}
print(f"\nAR vs TGF-beta-driven module across clusters: rho {ar_rho:+.2f} (p {ar_p:.2g})")

x = pd.ExcelFile(ATLAS / "Supplemental_Tables.xlsx")
t4 = x.parse("TableS4", header=None)
start = t4.index[t4[0] == "scATAC High-Resolution Clusters"][0] + 2
atac = {}
for _, r in t4.iloc[start:].iterrows():
    if pd.isna(r[0]):
        break
    atac[r[0]] = r[1]
v = x.parse("TableS12", header=1)
v = v[v.trait == "AGA"].copy()
v["cluster"] = v.top_model_cluster.map(atac)
v["compartment"] = v.cluster.map(lambda c: next((k for p, k in ATAC_COMPARTMENT_PREFIX.items() if str(c).startswith(p)), "Immune"))
tgfb_terms = set(driven) | set(fib_up) | set(genetics_tgfb) | {"SMAD3", "SMAD4", "TGFB2", "TGFBR1"}
v["tgfb_link"] = v.apply(lambda r: "SMAD" in str(r.top_motifs) or bool(set(str(r.linkedGenes).split(";")) & tgfb_terms), axis=1)
v.to_csv(RESULTS / "08_aga_variants.csv", index=False)
summary["aga_variants"] = {"n": int(len(v)), "by_compartment": v.compartment.value_counts().to_dict(),
                           "tgfb_or_scarring_linked": v[v.tgfb_link][["snp", "cluster", "fm_probs", "top_motifs", "linkedGenes"]].to_dict(orient="records")}
print(f"\nAGA high-effect variants (n={len(v)}) by compartment of the best-predicting chromatin model:", v.compartment.value_counts().to_dict())
print(v[v.tgfb_link][["snp", "cluster", "fm_probs", "top_motifs", "linkedGenes"]].to_string(index=False))
save_summary("08", summary)

shown = {"TGF-beta-driven bald genes": "TGF-beta-driven bald-scalp genes (22)", "Pirfenidone-lowered (1 mM)": "...of those, lowered by pirfenidone (10)",
         "Scarring genes up in bald scalp": "Scarring genes up in bald scalp (16)", "Hair-loss genetics, all strong (>=0.5)": "Hair-loss genetics, strong (66)",
         "Bald-down top 200": "Genes lost in bald scalp (111)"}
cols = [c for c in compartment.index if compartment[c] != "Immune"]
mat = pd.DataFrame({lab: table[f"{k} | score"] for k, lab in shown.items()}).T
qs = pd.DataFrame({lab: table[f"{k} | q"] for k, lab in shown.items()}).T
imm = compartment.index[compartment == "Immune"]
mat = pd.concat([mat[cols], mat[imm].mean(axis=1).rename("Immune (15, mean)")], axis=1)
qs = pd.concat([qs[cols], qs[imm].min(axis=1).rename("Immune (15, mean)") * np.nan], axis=1)
groups = [compartment[c] for c in cols] + ["Immune"]
cmap = LinearSegmentedColormap.from_list("div", ["#3B6FB6", "#F2F2EF", "#C4572A"])
fig, ax = plt.subplots(figsize=(13, 4.6))
lim = float(np.nanmax(np.abs(mat.values)))
im = ax.imshow(mat.values, cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if qs.iloc[i, j] < .05 and mat.iloc[i, j] > 0:
            ax.text(j, i, "*", ha="center", va="center", fontsize=11, color="#1f2a33")
ax.set_xticks(range(mat.shape[1]))
ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=8)
ax.set_yticks(range(mat.shape[0]))
ax.set_yticklabels(mat.index, fontsize=9)
edges = [j for j in range(1, len(groups)) if groups[j] != groups[j - 1]]
for j in edges:
    ax.axvline(j - .5, color="white", lw=2.5)
starts = [0] + edges
for s0, s1 in zip(starts, edges + [len(groups)]):
    short = {"Follicle dermis": "Follicle\ndermis", "Muscle and pericytes": "Muscle,\npericytes", "Follicle epithelium": "Follicle\nepithelium",
             "Melanocytes": "Melano-\ncytes", "Immune": "Imm."}
    ax.text((s0 + s1 - 1) / 2, -0.6, short.get(groups[s0], groups[s0]), ha="center", va="bottom", fontsize=8, color="#44525c")
for spine in ax.spines.values():
    spine.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("mean z across cell types", fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.set_title("Where each gene set is expressed in human scalp (Ober-Reynolds 2023 atlas). * = enriched, q < 0.05", fontsize=10, loc="left", pad=30)
fig.tight_layout()
(ROOT / "figs").mkdir(exist_ok=True)
fig.savefig(ROOT / "figs" / "08_scalp_localization.png", dpi=200)
print("saved figs/08_scalp_localization.png")
