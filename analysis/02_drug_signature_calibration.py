"""Track B: score every drug in three public drug-response libraries for reversing the signatures, and test whether
the scores separate drugs known to grow hair from drugs known to cause hair loss (SIDER side effects)."""
import json
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from common import RESULTS, save_summary
import drug_scoring as D

sigs = {k: (set(v[0]), set(v[1])) for k, v in json.load(open(RESULTS / "disease_signatures.json")).items()}
labels = D.sider_labels()
print(f"SIDER: {sum(labels.values())} hair-growth drugs, {len(labels) - sum(labels.values())} hair-loss drugs")
known = ["finasteride", "dutasteride", "minoxidil", "latanoprost", "bimatoprost", "spironolactone", "cyproterone",
         "ketoconazole", "tofacitinib", "ruxolitinib", "baricitinib", "valproicacid", "ciclosporin", "diazoxide", "bicalutamide", "flutamide"]
cal_rows, all_scores = [], {}
for lib in ("lincs", "cmap"):
    for name, (u, d) in sigs.items():
        sc = D.score_library(u, d, lib)
        cal = D.calibration(sc, labels, n_perm=2000)
        pct = {k: round(100 * (1 - sc.index[sc.key == D.norm(k)][0] / len(sc)), 1) for k in known if (sc.key == D.norm(k)).any()}
        cal_rows.append({"library": lib, "signature": name, "drugs": len(sc), "auroc": cal["auroc"], "p_perm": cal["p_perm"],
                         "n_growth": cal["n_pos"], "n_loss": cal["n_neg"], "top10": "; ".join(sc.drug.head(10)), **{f"pct_{k}": v for k, v in pct.items()}})
        all_scores[f"{lib}|{name}"] = sc.set_index("drug").z
        print(f"[{lib} | {name}] AUROC {cal['auroc']:.3f} (p={cal['p_perm']:.3f}) | top: {', '.join(sc.drug.head(6))}")
rg = D.read_gmt(D.DRUGS / "RummaGEO_DrugPerturbations_2025.gmt")
pairs = defaultdict(dict)
for term, genes in rg.items():
    m = re.match(r"^(.*) (HUMAN|MOUSE) (UP|DN)$", term)
    if m and m.group(2) == "HUMAN":
        pairs[m.group(1)][m.group(3)] = genes
contrasts = {k: (v["UP"], v["DN"]) for k, v in pairs.items() if "UP" in v and "DN" in v and len(v["UP"]) >= 10 and len(v["DN"]) >= 10}
uni = set().union(*[u | d for u, d in contrasts.values()])
names = {k: re.sub(r"^\S+ \S+ V \S+ ", "", k) for k in contrasts}
for name in ("invivo_200", "dht_200"):
    u, d = sigs[name]
    df = pd.DataFrame([(names[k], D.reversal_z(u, d, cu, cd, uni)) for k, (cu, cd) in contrasts.items()], columns=["drug", "z"])
    df = df.groupby("drug").z.mean().reset_index()
    df["key"] = df.drug.map(D.norm)
    cal = D.calibration(df, labels, n_perm=500)
    cal_rows.append({"library": "rummageo", "signature": name, "drugs": df.drug.nunique(), "auroc": cal["auroc"] if cal else np.nan,
                     "p_perm": cal["p_perm"] if cal else np.nan, "top10": "; ".join(df.sort_values("z", ascending=False).drug.head(10))})
    print(f"[rummageo | {name}] {df.drug.nunique()} perturbagens, calibration " + (f"AUROC {cal['auroc']:.3f}" if cal else "not possible (too few labeled drugs)"))
pd.DataFrame(cal_rows).to_csv(RESULTS / "02_calibration.csv", index=False)
pd.concat(all_scores, axis=1).to_csv(RESULTS / "02_drug_scores.csv")
save_summary("02", {"calibration": cal_rows})
