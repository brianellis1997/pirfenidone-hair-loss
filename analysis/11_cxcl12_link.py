"""Track K: the CXCL12 link. The evidence agent (Track I) surfaced a 2025 study (PMID 40724819, testosterone-induced mouse
model): dermal fibroblasts make CXCL12, CXCL12 switches on TGF-beta and scarring, and blocking CXCL12 reversed the remodeling
and regrew hair. Three checks in human data we already have:
1. Is CXCL12 (and its receptors ACKR3, CXCR4) up in bald human scalp? (GSE36169, GSE90594, consensus)
2. Does pirfenidone lower CXCL12 in human fibroblasts, and at which dose? (GSE226249 1 mM, GSE178518 about 0.1 mM)
3. Which scalp cells make CXCL12 and carry its receptors? (Ober-Reynolds atlas)"""
import numpy as np
import pandas as pd

from common import RESULTS, ROOT, consensus_bald_signature, save_summary

GENES = ["CXCL12", "ACKR3", "CXCR4"]
inv = consensus_bald_signature()
bald = inv.reindex(GENES)[["lfc_g", "p_g", "lfc_m", "p_m", "z"]]
bald["z_rank_pct"] = [float((inv.z < inv.loc[g, "z"]).mean() * 100) if g in inv.index else np.nan for g in GENES]
print("bald scalp (g = GSE36169 paired, m = GSE90594):\n", bald.round(4).to_string())

fx = pd.read_csv(RESULTS / "05_fibroblast_tgfb_drug_effects.csv", index_col=0)
drug = fx.reindex(GENES)[["tgfb_lfc", "tgfb_p", "pirf_1mM_lfc", "pirf_1mM_p", "pirf_0p1mM_lfc", "pirf_0p1mM_p", "nint_lfc", "nint_p"]]
print("\nfibroblasts (lfc = log2 change):\n", drug.round(4).to_string())

pb = pd.read_csv(ROOT / "data" / "scalp_atlas" / "pbFineClust_log2CP10kRNA.tsv", sep="\t", index_col=0).drop(columns=["Unknown", "Plasma.contam"])
atlas = {g: pb.loc[g].nlargest(5).round(2).to_dict() for g in GENES if g in pb.index}
for g, top in atlas.items():
    print(f"\natlas, {g} top clusters: " + ", ".join(f"{c} {v}" for c, v in top.items()))
fib_cols = [c for c in pb.columns if c.startswith("D.")]
share = {g: float((2 ** pb.loc[g, fib_cols] - 1).sum() / (2 ** pb.loc[g] - 1).sum()) for g in GENES if g in pb.index}
print("\nshare of summed pseudobulk expression in dermal fibroblast, sheath and papilla clusters (7 of 41):", {k: round(v, 2) for k, v in share.items()})

out = bald.join(drug)
out.to_csv(RESULTS / "11_cxcl12_link.csv")
save_summary("11", {"source_paper": "PMID 40724819 (2025, Int J Mol Sci): CXCL12 drives reversible fibroimmune remodeling in AGA (mouse)",
                    "bald_scalp": bald.replace({np.nan: None}).to_dict(orient="index"), "fibroblast_drugs": drug.replace({np.nan: None}).to_dict(orient="index"),
                    "atlas_top_clusters": atlas, "atlas_share_in_dermal_clusters": share})
