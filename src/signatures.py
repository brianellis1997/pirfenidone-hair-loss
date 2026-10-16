"""
Scalp DHT damage signature: the DHT organ-culture study (GSE267664), bald vs haired scalp (GSE36169), and the
beard/scalp/armpit follicle contrast (GSE189684).
"""
import gzip
import io
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent.parent
GEO = HERE / "data" / "geo"


def dht_organ_culture():
    df = pd.read_csv(GEO / "GSE267664_DHT-vs-control-all.gene.xls.gz", sep="\t", low_memory=False)
    df = df.rename(columns={"gene_id": "gene", "log2FoldChange": "lfc_dht", "p-value": "p_dht", "q-value": "q_dht"})
    df = df[(df.coding_type == "protein_coding") & (df.BaseMean_DESeq2 >= 10)]
    return df[["gene", "lfc_dht", "p_dht", "q_dht", "BaseMean_DESeq2"]].drop_duplicates("gene").set_index("gene")


def read_series_matrix(path):
    lines = gzip.open(path, "rt").read().splitlines()
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end")
    titles = next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]
    table = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    table.columns = [t.strip('"') for t in titles]
    return table


def bald_vs_haired():
    expr = read_series_matrix(GEO / "GSE36169_series_matrix.txt.gz")
    if expr.max().max() > 100:
        expr = np.log2(expr.clip(lower=1))
    annot = read_annot(GEO / "GPL96.annot.gz")[["ID", "Gene symbol"]]
    annot = annot.dropna()
    annot = annot[~annot["Gene symbol"].str.contains("///")]
    expr = expr.join(annot.set_index("ID")["Gene symbol"], how="inner")
    expr["mean"] = expr.drop(columns="Gene symbol").mean(axis=1)
    expr = expr.sort_values("mean", ascending=False).drop_duplicates("Gene symbol").set_index("Gene symbol").drop(columns="mean")
    bald = expr[[c for c in expr.columns if "Bald" in c]].values
    haired = expr[[c for c in expr.columns if "Haired" in c]].values
    diff = bald - haired
    t, p = stats.ttest_rel(bald, haired, axis=1)
    out = pd.DataFrame({"lfc_bald": diff.mean(axis=1), "p_bald": p, "mean_expr": expr.mean(axis=1)}, index=expr.index)
    return out[out.mean_expr >= np.percentile(out.mean_expr, 25)]


def follicle_sites():
    names = {"AB": "beard", "AF": "frontal", "AO": "occipital", "AX": "armpit"}
    tables = []
    for f in sorted((GEO / "GSE189684").glob("*.txt.gz")):
        site = names[f.name.split("_")[1][:2]]
        t = pd.read_csv(f, sep="\t", index_col=0)
        tables.append(t.iloc[:, 0].rename(site))
    tpm = pd.concat(tables, axis=1).fillna(0)
    lg = np.log2(tpm + 1)
    return pd.DataFrame({"frontal_vs_occipital": lg.frontal - lg.occipital, "beard_vs_frontal": lg.beard - lg.frontal,
                         "armpit_vs_frontal": lg.armpit - lg.frontal, "max_log_tpm": lg.max(axis=1)})


def series_matrix_with_meta(path):
    lines = gzip.open(path, "rt", errors="ignore").read().splitlines()
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end")
    meta = {}
    for l in lines[:start]:
        if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1") or l.startswith("!Sample_source_name_ch1"):
            vals = [v.strip('"').strip() for v in l.split("\t")[1:]]
            key = l.split("\t")[0]
            if key == "!Sample_characteristics_ch1" and vals and ":" in vals[0]:
                key = vals[0].split(":")[0].strip()
                vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            meta[key] = vals
    table = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    meta = pd.DataFrame(meta, index=table.columns)
    return table, meta


def read_annot(path):
    lines = gzip.open(path, "rt", errors="ignore").read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    j = next(k for k, l in enumerate(lines) if l.startswith("!platform_table_end"))
    return pd.read_csv(io.StringIO("\n".join(lines[i:j])), sep="\t", low_memory=False)


def probe_map(platform):
    if platform == "GPL17077":
        t = open(GEO / "GPL17077_table.txt").read().splitlines()
        i = next(k for k, l in enumerate(t) if l.startswith("!platform_table_begin")) + 1
        j = next(k for k, l in enumerate(t) if l.startswith("!platform_table_end"))
        df = pd.read_csv(io.StringIO("\n".join(t[i:j])), sep="\t", low_memory=False)
        return df.set_index("ID")["GENE_SYMBOL"].dropna()
    if platform == "GPL10558":
        df = read_annot(GEO / "GPL10558.annot.gz").dropna(subset=["Gene symbol"])
        return df[~df["Gene symbol"].str.contains("///")].set_index("ID")["Gene symbol"]
    raise ValueError(platform)


def collapse(expr, mapping):
    expr = expr.join(mapping.rename("symbol"), how="inner")
    expr["_m"] = expr.drop(columns="symbol").mean(axis=1)
    expr = expr.sort_values("_m", ascending=False).drop_duplicates("symbol").set_index("symbol").drop(columns="_m")
    if expr.max().max() > 100:
        expr = np.log2(expr.clip(lower=1))
    return expr


def frontal_vertex_vs_occiput():
    expr, meta = series_matrix_with_meta(GEO / "GSE78722_series_matrix.txt.gz")
    expr = collapse(expr, probe_map("GPL17077"))
    subj = meta["!Sample_title"].str.split().str[0]
    region = meta["region"]
    affected = np.vstack([expr[[c for c in expr.columns if subj[c] == s and region[c] in ("frontal", "vertex")]].mean(axis=1) for s in sorted(subj.unique())]).T
    occiput = np.vstack([expr[[c for c in expr.columns if subj[c] == s and region[c] == "occiput"]].mean(axis=1) for s in sorted(subj.unique())]).T
    t, p = stats.ttest_rel(affected, occiput, axis=1)
    out = pd.DataFrame({"lfc": (affected - occiput).mean(axis=1), "p": p, "mean_expr": expr.mean(axis=1)}, index=expr.index)
    return out[out.mean_expr >= np.percentile(out.mean_expr, 25)]


def alopecia_vs_healthy():
    expr, meta = series_matrix_with_meta(GEO / "GSE90594_series_matrix.txt.gz")
    expr = collapse(expr, probe_map("GPL17077"))
    state = meta["disease state"]
    a = expr[[c for c in expr.columns if state[c] == "alopecia"]]
    h = expr[[c for c in expr.columns if state[c] == "healthy control"]]
    t, p = stats.ttest_ind(a, h, axis=1, equal_var=False)
    out = pd.DataFrame({"lfc": a.mean(axis=1) - h.mean(axis=1), "p": p, "mean_expr": expr.mean(axis=1)}, index=expr.index)
    return out[out.mean_expr >= np.percentile(out.mean_expr, 25)]


def ensembl_to_symbol():
    h = pd.read_csv(HERE / "data" / "ref" / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], low_memory=False)
    return h.dropna().drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]


def dht_on_dermal_papilla():
    cols = {}
    for f in sorted((GEO / "GSE178374").glob("*.txt.gz")):
        name = f.name.split("_", 1)[1].replace(".txt.gz", "")
        cols[name] = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
    expr = pd.DataFrame(cols)
    sym = ensembl_to_symbol()
    expr = expr.join(sym.rename("symbol"), how="inner").groupby("symbol").max()
    lg = np.log2(expr + 1)
    rows = {}
    for fmt in ("2D", "3D"):
        d = lg[[c for c in lg.columns if c.startswith("DHT" + fmt)]]
        e = lg[[c for c in lg.columns if c.startswith("empty" + fmt)]]
        t, p = stats.ttest_ind(d, e, axis=1)
        rows[fmt] = pd.DataFrame({"lfc": d.mean(axis=1) - e.mean(axis=1), "p": p})
    out = rows["2D"].join(rows["3D"], lsuffix="_2d", rsuffix="_3d")
    out["mean_expr"] = lg.mean(axis=1)
    return out[out.mean_expr >= 1]


LATE = ("12h", "16h", "20h", "24h", "36h", "48h")


def dht_time_course():
    """GSE66664: balding (BAB) vs non-balding (BAN) dermal papilla cells, DHT 1 or 10 nM, 0h to 48h, triplicates."""
    expr, meta = series_matrix_with_meta(GEO / "GSE66664_series_matrix.txt.gz")
    expr = collapse(expr, probe_map("GPL10558"))
    t = meta["time (treatment duration)"]
    line = meta["cell line"]
    out = {}
    late_vals = {}
    for ln in ("BAB", "BAN"):
        base = expr[[c for c in expr.columns if line[c] == ln and t[c] == "0h"]]
        late = expr[[c for c in expr.columns if line[c] == ln and t[c] in LATE]]
        stat, pv = stats.ttest_ind(late, base, axis=1, equal_var=False)
        out[f"resp_{ln.lower()}"] = late.mean(axis=1) - base.mean(axis=1)
        out[f"p_{ln.lower()}"] = pd.Series(pv, index=expr.index)
        late_vals[ln] = late.sub(base.mean(axis=1), axis=0)
    stat, pv = stats.ttest_ind(late_vals["BAB"], late_vals["BAN"], axis=1, equal_var=False)
    out["resp_diff"] = out["resp_bab"] - out["resp_ban"]
    out["p_diff"] = pd.Series(pv, index=expr.index)
    out["mean_expr"] = expr.mean(axis=1)
    df = pd.DataFrame(out)
    return df[df.mean_expr >= np.percentile(df.mean_expr, 30)]


def balding_vs_nonbalding_under_dht():
    """GSE184501: immortalised balding vs non-balding dermal papilla cells, both on 10 nM DHT for 24 h, triplicates."""
    df = pd.read_csv(GEO / "GSE184501_2D_cuffnorm_fpkm_quantilenormalised.txt.gz", sep="\t")
    g = df.groupby("gene_short_name")[[c for c in df.columns if c.startswith(("BAN", "BAB"))]].sum()
    lg = np.log2(g + 1)
    bab, ban = lg[[c for c in lg.columns if c.startswith("BAB")]], lg[[c for c in lg.columns if c.startswith("BAN")]]
    stat, pv = stats.ttest_ind(bab, ban, axis=1, equal_var=False)
    out = pd.DataFrame({"lfc": bab.mean(axis=1) - ban.mean(axis=1), "p": pv, "mean_expr": lg.mean(axis=1)}, index=lg.index)
    return out[out.mean_expr >= 1]


def _gsm_blocks(acc):
    import re as _re
    import time as _time
    import urllib.request as _ur
    _time.sleep(0.6)
    txt = _ur.urlopen(_ur.Request(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief",
                                  headers={"User-Agent": "Mozilla/5.0"}), timeout=90).read().decode("utf-8", "ignore")
    rows = []
    for block in txt.split("^SAMPLE")[1:]:
        title = _re.search(r"^!Sample_title = (.*)$", block, _re.M).group(1).strip()
        ch = dict(c.split(":", 1) for c in _re.findall(r"^!Sample_characteristics_ch1 = (.*)$", block, _re.M) if ":" in c)
        rows.append({"title": title, **{k.strip(): v.strip() for k, v in ch.items()}})
    return pd.DataFrame(rows).set_index("title")


def _welch(a, b):
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    return pd.DataFrame({"lfc": a.mean(axis=1) - b.mean(axis=1), "p": p}, index=a.index)


def tgfb_and_antifibrotics():
    """GSE178518: normal human lung fibroblasts +/- TGF-beta, with pirfenidone, nintedanib, saracatinib or vehicle (24 h)."""
    meta = _gsm_blocks("GSE178518")
    expr = pd.read_csv(GEO / "GSE178518_ipf_voom_normalized.txt.gz", sep="\t", index_col=0)
    expr = expr.join(ensembl_to_symbol().rename("symbol"), how="inner").groupby("symbol").mean()
    expr = expr[expr.mean(axis=1) > 2]
    stim, inh = meta["stimulation"], meta["inhibitor"]
    cols = lambda s, i: [c for c in expr.columns if c in meta.index and stim[c] == s and inh[c] == i]
    out = {"tgfb": _welch(expr[cols("TGFb", "Vehicle A")], expr[cols("Vehicle B", "Vehicle A")])}
    for drug, label in (("Pirfenidone", "pirfenidone"), ("Nintedanib", "nintedanib"), ("Saracatanib", "saracatinib")):
        out[label] = _welch(expr[cols("TGFb", drug)], expr[cols("TGFb", "Vehicle A")])
    return out, {k: len(v) for k, v in {"tgfb_veh": cols("TGFb", "Vehicle A"), "veh_veh": cols("Vehicle B", "Vehicle A"),
                                        "tgfb_pirf": cols("TGFb", "Pirfenidone"), "tgfb_nint": cols("TGFb", "Nintedanib")}.items()}


def pirfenidone_serum_fibroblasts():
    """GSE226249: lung fibroblasts from 3 donors, serum-starved, 20% FBS, or FBS + 1 mM pirfenidone, triplicates."""
    counts = pd.read_csv(GEO / "GSE226249_counts.tsv.gz", sep="\t", index_col=0)
    counts = counts.join(ensembl_to_symbol().rename("symbol"), how="inner").groupby("symbol").sum()
    cpm = np.log2(counts / counts.sum() * 1e6 + 1)
    cpm = cpm[cpm.mean(axis=1) > 1]
    donors = sorted({c.split("_")[0] for c in cpm.columns})
    mean = lambda d, cond: cpm[[c for c in cpm.columns if c.startswith(f"{d}_{cond}_")]].mean(axis=1)
    pf = np.vstack([mean(d, "PFD") - mean(d, "FBS") for d in donors]).T
    act = np.vstack([mean(d, "FBS") - mean(d, "SF") for d in donors]).T
    res = {}
    for name, diff in (("pirfenidone", pf), ("serum_activation", act)):
        t, p = stats.ttest_1samp(diff, 0, axis=1)
        res[name] = pd.DataFrame({"lfc": diff.mean(axis=1), "p": p}, index=cpm.index)
    return res
