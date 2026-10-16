"""Track A: gene-expression signatures of pattern hair loss and of DHT, and how well the studies agree."""
import json

import numpy as np
import pandas as pd
from scipy import stats

from common import RESULTS, S, consensus_bald_signature, save_summary

b1 = S.bald_vs_haired().rename(columns={"lfc_bald": "lfc", "p_bald": "p"})
b2 = S.frontal_vertex_vs_occiput()
b3 = S.alopecia_vs_healthy()
oc = S.dht_organ_culture()
dp = S.dht_on_dermal_papilla()
tc = S.dht_time_course()
bb = S.balding_vs_nonbalding_under_dht()

signal = {
    "GSE36169 bald vs haired, 5 paired men": (len(b1), int((b1.p < .05).sum())),
    "GSE90594 AGA vs healthy, 14 vs 14": (len(b3), int((b3.p < .05).sum())),
    "GSE78722 front+crown vs back, 6 men": (len(b2), int((b2.p < .05).sum())),
    "GSE267664 DHT organ culture, 3 vs 3": (len(oc), int((oc.p_dht < .05).sum())),
    "GSE178374 DHT on balding DP cells 2D": (len(dp), int((dp.p_2d < .05).sum())),
    "GSE178374 DHT on balding DP cells 3D": (len(dp), int((dp.p_3d < .05).sum())),
    "GSE66664 DHT time course, balding cells": (len(tc), int((tc.p_bab < .05).sum())),
    "GSE184501 balding vs non-balding DP cells on DHT": (len(bb), int((bb.p < .05).sum())),
}
print("signal per study (genes tested, genes p<.05, share):")
for k, (n, s) in signal.items():
    print(f"   {k:48} {n:6} {s:6} {s / n:5.0%}")


def agree(x, y, cx, cy, px, py):
    m = x.join(y, lsuffix="_x", rsuffix="_y", how="inner")
    s = m[(m[px + "_x"] < .05) | (m[py + "_y"] < .05)]
    both = m[(m[px + "_x"] < .05) & (m[py + "_y"] < .05)]
    return {"shared": len(m), "rho_sig_either": float(stats.spearmanr(s[cx + "_x"], s[cy + "_y"])[0]),
            "both_sig": len(both), "same_direction": int((np.sign(both[cx + "_x"]) == np.sign(both[cy + "_y"])).sum())}


agreement = {
    "GSE36169 vs GSE90594": agree(b1, b3, "lfc", "lfc", "p", "p"),
    "GSE36169 vs GSE78722": agree(b1, b2, "lfc", "lfc", "p", "p"),
    "GSE90594 vs GSE78722": agree(b3, b2, "lfc", "lfc", "p", "p"),
}
x = dp.join(tc, how="inner", lsuffix="_dp", rsuffix="_tc")
s2 = x[(x.p_2d < .05) | (x.p_bab < .05)]
s3 = x[(x.p_3d < .05) | (x.p_bab < .05)]
agreement["GSE178374 2D vs GSE66664 balding"] = {"rho": float(stats.spearmanr(s2.lfc_2d, s2.resp_bab)[0])}
agreement["GSE178374 3D vs GSE66664 balding"] = {"rho": float(stats.spearmanr(s3.lfc_3d, s3.resp_bab)[0])}
print("\nagreement:", json.dumps(agreement, indent=1))

inv = consensus_bald_signature()
y = inv.join(x, how="inner")
sel = y[np.abs(y.z) > 2.5]
invivo_vs_dht = {c: float(stats.spearmanr(sel.z, sel[c])[0]) for c in ("lfc_2d", "lfc_3d", "resp_bab", "resp_diff")}
print("\nin-vivo bald signature vs DHT responses (rho):", invivo_vs_dht)

markers = ["PTGDS", "LEF1", "HOXC13", "MSX2", "DKK1", "SFRP2", "TGFB1", "IL6", "IGF1", "SRD5A2", "AR", "WNT10A", "FGF5"]
rows = []
for g in markers:
    rows.append({"gene": g,
                 **{f"{n}_lfc": (d.loc[g, "lfc"] if g in d.index else np.nan) for n, d in (("GSE36169", b1), ("GSE78722", b2), ("GSE90594", b3))},
                 **{f"{n}_p": (d.loc[g, "p"] if g in d.index else np.nan) for n, d in (("GSE36169", b1), ("GSE78722", b2), ("GSE90594", b3))},
                 "dht_dp_2d": dp.loc[g, "lfc_2d"] if g in dp.index else np.nan, "dht_dp_3d": dp.loc[g, "lfc_3d"] if g in dp.index else np.nan})
pd.DataFrame(rows).to_csv(RESULTS / "01_marker_genes.csv", index=False)
inv.to_csv(RESULTS / "01_consensus_bald_signature.csv")
dp.to_csv(RESULTS / "01_dht_on_dermal_papilla.csv")

universe = set().union(*[u | d for u, d in __import__("drug_scoring").lincs_consensus().values()])
both = dp[(dp.p_2d < .05) & (dp.p_3d < .05) & (np.sign(dp.lfc_2d) == np.sign(dp.lfc_3d))].copy()
both["eff"] = (both.lfc_2d + both.lfc_3d) / 2
sigs = {}
for n in (100, 200):
    iu = inv[inv.index.isin(universe)]
    sigs[f"invivo_{n}"] = [sorted(iu[iu.z > 0].nlargest(n, "z").index), sorted(iu[iu.z < 0].nsmallest(n, "z").index)]
    bu = both[both.index.isin(universe)]
    sigs[f"dht_{n}"] = [sorted(bu[bu.eff > 0].nlargest(n, "eff").index), sorted(bu[bu.eff < 0].nsmallest(n, "eff").index)]
json.dump(sigs, open(RESULTS / "disease_signatures.json", "w"))
save_summary("01", {"signal": signal, "agreement": agreement, "invivo_vs_dht_rho": invivo_vs_dht,
                    "consensus_genes": len(inv), "dht_consistent_genes": len(both)})
