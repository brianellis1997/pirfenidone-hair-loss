"""Track E: is part of the bald-scalp program TGF-beta-driven, and does pirfenidone switch it off in human fibroblasts?
GSE178518: lung fibroblasts +/- TGF-beta with pirfenidone (20 ug/mL, about 0.1 mM), nintedanib or saracatinib.
GSE226249: lung fibroblasts from 3 donors, serum-activated, +/- pirfenidone 1 mM."""
import numpy as np
import pandas as pd
from scipy import stats

from common import RESULTS, S, consensus_bald_signature, save_summary

inv = consensus_bald_signature()
fib, sizes = S.tgfb_and_antifibrotics()
sr = S.pirfenidone_serum_fibroblasts()
tg = fib["tgfb"]
up200 = set(inv[inv.z > 0].nlargest(200, "z").index)
universe = set(inv.index) & set(tg.index)
tg_up = set(tg[(tg.p < .05) & (tg.lfc > .5)].index) & universe
a = len(up200 & tg_up); b = len((up200 & universe) - tg_up); c = len(tg_up - up200); d = len(universe - up200 - tg_up)
odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
driven = sorted(up200 & tg_up)
print(f"bald-up genes that TGF-beta induces: {a}/{a + b} vs background {len(tg_up) / len(universe):.1%}; OR {odds:.2f}, p {p:.1e}")
print("TGF-beta-driven bald genes:", ", ".join(driven))
sel = inv[np.abs(inv.z) > 2.5]
res = {"tgfb_enrichment": {"overlap": a, "bald_up_tested": a + b, "background_rate": len(tg_up) / len(universe), "odds_ratio": odds, "p": p,
                           "driven_genes": driven, "rho_all_bald_genes": float(stats.spearmanr(sel.join(tg).dropna().z, sel.join(tg).dropna().lfc)[0])},
       "group_sizes": sizes, "drugs": {}}
for drug in ("pirfenidone", "nintedanib", "saracatinib"):
    dd = fib[drug]
    s2 = sel.join(dd, how="inner")
    rho, pr = stats.spearmanr(s2.z, s2.lfc)
    dv = dd.loc[[g for g in driven if g in dd.index]]
    res["drugs"][f"{drug}_GSE178518"] = {"rho": rho, "p": pr, "driven_lowered_sig": float(((dv.lfc < 0) & (dv.p < .05)).mean()), "driven_median_change": float(dv.lfc.median())}
pf = sr["pirfenidone"]
s3 = sel.join(pf, how="inner")
dv = pf.loc[[g for g in driven if g in pf.index]]
rho, pr = stats.spearmanr(s3.z, s3.lfc)
res["drugs"]["pirfenidone_1mM_GSE226249"] = {"rho": rho, "p": pr, "driven_lowered_sig": float(((dv.lfc < 0) & (dv.p < .05)).mean()), "driven_median_change": float(dv.lfc.median())}
for k, v in res["drugs"].items():
    print(f"{k:28} bald-vs-drug rho {v['rho']:+.3f} (p {v['p']:.1e}) | TGF-beta-driven bald genes lowered {v['driven_lowered_sig']:.0%}, median {v['driven_median_change']:+.2f}")
out = pd.DataFrame({"bald_z": inv.z}).join(tg.add_prefix("tgfb_"), how="outer").join(fib["pirfenidone"].add_prefix("pirf_0p1mM_"))
out = out.join(fib["nintedanib"].add_prefix("nint_")).join(pf.add_prefix("pirf_1mM_"))
out.to_csv(RESULTS / "05_fibroblast_tgfb_drug_effects.csv")
save_summary("05", res)
