"""Track L: the beard-versus-scalp paradox. DHT grows beard hair and shrinks frontal scalp hair. GSE189684 has surgically
removed follicles from one man with AGA: frontal (balding), occipital (spared), beard and armpit. If the scarring arm
is part of why frontal follicles shrink, the scarring and TGF-beta genes should run higher in frontal follicles than in
occipital and beard follicles from the same person.
This is one person, one sample per site. The permutation p-values say whether a gene set moves together more than random
gene sets do in that one comparison. They are not evidence about people in general."""
import json

import numpy as np
import pandas as pd

from common import RESULTS, S, consensus_bald_signature, save_summary

fs = S.follicle_sites()
expressed = fs[fs.max_log_tpm >= 1]
bins = pd.Series(pd.qcut(expressed.max_log_tpm, 20, labels=False, duplicates="drop"), index=expressed.index)
t5 = json.load(open(RESULTS / "05_summary.json"))
fx = pd.read_csv(RESULTS / "05_fibroblast_tgfb_drug_effects.csv", index_col=0)
driven = t5["tgfb_enrichment"]["driven_genes"]
comp = pd.read_csv(RESULTS / "04_composition_check_GSE90594.csv")
inv = consensus_bald_signature()
sets = {
    "TGF-beta-driven bald genes": driven,
    "Pirfenidone-lowered (1 mM)": [g for g in driven if g in fx.index and fx.loc[g, "pirf_1mM_lfc"] < 0 and fx.loc[g, "pirf_1mM_p"] < .05],
    "Scarring genes up in bald scalp": comp[(comp.group == "fibrosis") & (comp.bald_effect > 0) & (comp.p_raw < .05)].gene.tolist(),
    "Immune genes up in bald scalp": comp[(comp.group == "immune") & (comp.bald_effect > 0) & (comp.p_raw < .05)].gene.tolist(),
    "Bald-up top 200": inv[inv.z > 0].nlargest(200, "z").index.tolist(),
    "Bald-down top 200": inv[inv.z < 0].nsmallest(200, "z").index.tolist(),
}
CONTRASTS = {"frontal_vs_occipital": "frontal minus occipital", "beard_vs_frontal": "beard minus frontal", "armpit_vs_frontal": "armpit minus frontal"}
rng = np.random.default_rng(0)
pools = {b: np.flatnonzero(bins.values == b) for b in np.unique(bins.values)}
rows = []
for name, genes in sets.items():
    used = [g for g in genes if g in expressed.index]
    draws = np.stack([rng.choice(pools[bins[g]], 10000) for g in used], axis=1)
    for col, label in CONTRASTS.items():
        vals = expressed[col].values
        obs = float(expressed.loc[used, col].mean())
        null = vals[draws].mean(axis=1)
        p = float((np.abs(null - null.mean()) >= abs(obs - null.mean())).mean())
        rows.append({"gene_set": name, "genes": len(used), "contrast": label, "mean_log2_diff": obs, "null_mean": float(null.mean()),
                     "perm_p_two_sided": max(p, 1e-4), "share_of_genes_up": float((expressed.loc[used, col] > 0).mean())})
res = pd.DataFrame(rows)
res.to_csv(RESULTS / "12_beard_vs_scalp.csv", index=False)
pd.set_option("display.width", 200)
print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
single = ["CXCL12", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD3", "TWIST1", "TWIST2", "AR", "SRD5A2", "COL1A1", "COL3A1", "LOX", "SPARC", "TIMP1", "DKK1", "PTGDS"]
genes = fs.reindex(single)[list(CONTRASTS) + ["max_log_tpm"]].round(2)
print("\nsingle genes (log2 TPM differences):\n", genes.to_string())
save_summary("12", {"design": "GSE189684, one man with AGA, one surgical follicle sample per site", "sets": rows,
                    "single_genes": genes.replace({np.nan: None}).to_dict(orient="index")})
