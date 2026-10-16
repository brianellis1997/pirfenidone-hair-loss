"""Track D: the scarring (fibrosis) program in bald scalp: pathway enrichment, and whether it survives adjusting for
how much follicle tissue each biopsy contains."""
import json
import time
import urllib.request
import uuid

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from common import RESULTS, S, consensus_bald_signature, save_summary


def enrichr(genes, library, desc="list"):
    boundary = uuid.uuid4().hex
    body = (f"--{boundary}\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) +
            f"\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\n{desc}\r\n--{boundary}--\r\n").encode()
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body,
                                 headers={"Content-Type": f"multipart/form-data; boundary={boundary}", "User-Agent": "Mozilla/5.0"})
    uid = json.loads(urllib.request.urlopen(req, timeout=120).read())["userListId"]
    time.sleep(0.7)
    r = json.loads(urllib.request.urlopen(urllib.request.Request(
        f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={library}", headers={"User-Agent": "Mozilla/5.0"}), timeout=120).read())
    return pd.DataFrame([{"term": x[1], "p": x[2], "adj_p": x[6], "n_genes": len(x[5]), "genes": ",".join(x[5])} for x in r[library]])


inv = consensus_bald_signature()
up = inv[inv.z > 0].nlargest(200, "z").index.tolist()
dn = inv[inv.z < 0].nsmallest(200, "z").index.tolist()
gen = pd.read_csv(RESULTS / "03_genetic_targets_druggability.csv")
enr = {}
for label, genes in (("bald_up", up), ("bald_down", dn), ("aga_genetics", gen.symbol.tolist())):
    e = enrichr(genes, "MSigDB_Hallmark_2020", label)
    e.to_csv(RESULTS / f"04_hallmark_{label}.csv", index=False)
    enr[label] = e.head(6)[["term", "adj_p", "n_genes"]].to_dict(orient="records")
    print(f"\n{label}:\n", e.head(6)[["term", "adj_p", "n_genes"]].to_string(index=False))
    time.sleep(1)

FOLLICLE = ["KRT25", "KRT27", "KRT71", "KRT72", "KRT73", "TCHH", "KRT31", "KRT32", "KRT33A", "KRT33B", "KRT35", "KRT85",
            "KRT86", "LHX2", "MSX2", "HOXC13", "KRT75", "DSG4"]
FIB = ["COL1A1", "COL1A2", "COL3A1", "SPARC", "FBLN1", "FBLN2", "FBLN5", "TIMP1", "PDGFRB", "DCN", "LOX", "TGFBI", "FN1",
       "ACTA2", "VIM", "S100A4", "LUM", "POSTN", "TWIST1", "TGFB1", "SERPINE1", "CCN2"]
IMM = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DMA", "C1S", "C1R", "CFH", "VCAM1", "CD3E", "PTPRC", "CCL19"]
expr, meta = S.series_matrix_with_meta(S.GEO / "GSE90594_series_matrix.txt.gz")
expr = S.collapse(expr, S.probe_map("GPL17077"))
bald = (meta["disease state"] == "alopecia").astype(int).values
fol = [g for g in FOLLICLE if g in expr.index]
fscore = expr.loc[fol].sub(expr.loc[fol].mean(axis=1), axis=0).div(expr.loc[fol].std(axis=1), axis=0).mean().values
rows = []
for g in [x for x in FIB + IMM if x in expr.index]:
    y = expr.loc[g].values
    m0 = sm.OLS(y, sm.add_constant(bald)).fit()
    m1 = sm.OLS(y, sm.add_constant(np.column_stack([bald, fscore]))).fit()
    rows.append({"gene": g, "group": "fibrosis" if g in FIB else "immune", "bald_effect": m0.params[1], "p_raw": m0.pvalues[1],
                 "bald_effect_adj": m1.params[1], "p_adj_for_follicles": m1.pvalues[1], "follicle_coef": m1.params[2]})
comp = pd.DataFrame(rows)
comp.to_csv(RESULTS / "04_composition_check_GSE90594.csv", index=False)
comp_summary = {}
for grp in ("fibrosis", "immune"):
    x = comp[comp.group == grp]
    comp_summary[grp] = {"genes": len(x), "sig_raw": int((x.p_raw < .05).sum()), "sig_adjusted": int((x.p_adj_for_follicles < .05).sum()),
                         "median_effect_kept": float(np.nanmedian(x.bald_effect_adj / x.bald_effect.replace(0, np.nan)))}
print("\ncomposition check GSE90594:", comp_summary)
print(f"follicle score: hair loss {fscore[bald == 1].mean():+.2f} vs healthy {fscore[bald == 0].mean():+.2f}")
save_summary("04", {"hallmark_top": enr, "composition_GSE90594": comp_summary,
                    "follicle_score_means": {"hair_loss": float(fscore[bald == 1].mean()), "healthy": float(fscore[bald == 0].mean())},
                    "note_GSE36169": "all 5 men lost follicle content almost equally (-1.5 to -2.1 SD), so adjustment is collinear and uninformative there"})
