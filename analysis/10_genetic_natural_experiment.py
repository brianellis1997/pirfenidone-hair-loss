"""Track J: a genetic natural experiment. People inherit DNA variants that turn a gene's activity up or down, more or less at
random. If variants that raise a TGF-beta / scarring gene's activity in skin or fibroblasts are the same variants that raise
baldness risk, that supports the gene as a cause, and lowering it as a treatment.
Two tests per gene and tissue:
- colocalization (coloc ABF, Giambartolomei 2014): do the gene-activity signal and the baldness signal share one causal
  variant (H4) or come from different variants (H3)?
- Mendelian randomization (Wald ratio at the strongest activity variant): baldness change per unit of genetically
  higher gene activity.
Data:
- baldness GWAS: Yap et al. 2018, 205,327 UK Biobank men, BOLT-LMM, GWAS Catalog GCST007020;
- gene activity (eQTL): GTEx v8 via the eQTL Catalogue: cultured skin fibroblasts, suprapubic skin and lower-leg skin.
AR was meant as a positive control, but it has no usable activity variants in these tissues. Instead, a reference panel of
the 36 genes with the strongest baldness genetics (Open Targets score >= 0.7, Track C) goes through the same test, to show
how often a real baldness gene passes it."""
import concurrent.futures as cf
import gzip
import io
import json
import threading
import time
import urllib.error
import urllib.request
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import logsumexp

from common import RESULTS, ROOT, save_summary

GWAS_ZIP = ROOT / "data" / "gwas" / "mpb_bolt_lmm_aut_x.tab.zip"
GWAS_SUBSET = ROOT / "data" / "gwas" / "mpb_cis_subset.tsv.gz"
EQTL = ROOT / "data" / "eqtl"
EQTL.mkdir(parents=True, exist_ok=True)
GWAS_N = 205_327
DATASETS = {"QTD000216": "fibroblasts (cultured, skin)", "QTD000311": "skin, not sun-exposed", "QTD000316": "skin, sun-exposed"}
FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{d}/{d}.all.tsv.gz"
EQTL_COLS = ["molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type",
             "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"]
GENES = {"TGFB1": "TGF-beta ligand", "TGFB2": "TGF-beta ligand", "TGFBR1": "TGF-beta receptor", "TGFBR2": "TGF-beta receptor",
         "SMAD3": "TGF-beta signal carrier", "TWIST1": "scarring master switch", "HDAC9": "neighbor of TWIST1, approved-drug target",
         "TWIST2": "scarring master switch", "PDGFA": "fibroblast growth signal", "AR": "positive control (hormone receptor)"}
P1, P2, P12 = 1e-4, 1e-4, 1e-5


def gene_coords(symbols):
    path = EQTL / "genes_grch38.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for g in symbols:
        if g in out:
            continue
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?content-type=application/json"
        for attempt in range(6):
            try:
                d = json.loads(urllib.request.urlopen(urllib.request.Request(url), timeout=60).read())
                out[g] = {"id": d["id"], "chr": d["seq_region_name"], "start": d["start"], "end": d["end"]}
                break
            except urllib.error.HTTPError as e:
                if e.code in (400, 404):
                    out[g] = None
                    break
                time.sleep(2 + 2 * attempt)
        time.sleep(0.4)
    path.write_text(json.dumps(out, indent=1))
    return {g: out[g] for g in symbols if out.get(g)}


_idx_lock = threading.Lock()


def eqtl_region(dataset, gene, c):
    path = EQTL / f"{dataset}_{gene}.tsv.gz"
    if not path.exists():
        import pysam
        idx = EQTL / "tbi_cache" / f"{dataset}.all.tsv.gz.tbi"
        with _idx_lock:
            if not idx.exists():
                idx.parent.mkdir(exist_ok=True)
                urllib.request.urlretrieve(FTP.format(d=dataset) + ".tbi", idx)
        tb = pysam.TabixFile(FTP.format(d=dataset), index=str(idx))
        rows = [r.split("\t") for r in tb.fetch(c["chr"], max(0, c["start"] - 1_000_000), c["end"] + 1_000_000)]
        df = pd.DataFrame([r for r in rows if r[16] == c["id"]], columns=EQTL_COLS)
        df.to_csv(path, sep="\t", index=False)
    df = pd.read_csv(path, sep="\t")
    df = df[(df.molecular_trait_id == c["id"]) & df.rsid.astype(str).str.startswith("rs")]
    df = df.astype({"maf": float, "pvalue": float, "beta": float, "se": float, "an": float})
    return df.sort_values("pvalue").drop_duplicates("rsid")


def gwas_subset(rsids):
    saved = GWAS_SUBSET.parent / "mpb_cis_subset_rsids.json"
    if GWAS_SUBSET.exists() and saved.exists() and set(json.loads(saved.read_text())) >= rsids:
        return pd.read_csv(GWAS_SUBSET, sep="\t")
    keep = []
    with zipfile.ZipFile(GWAS_ZIP) as z:
        name = [n for n in z.namelist() if not n.endswith("/")][0]
        with z.open(name) as f:
            for chunk in pd.read_csv(io.TextIOWrapper(f), sep=r"\s+", chunksize=1_000_000):
                keep.append(chunk[chunk.SNP.isin(rsids)])
    sub = pd.concat(keep)
    sub.to_csv(GWAS_SUBSET, sep="\t", index=False)
    (GWAS_SUBSET.parent / "mpb_cis_subset_rsids.json").write_text(json.dumps(sorted(rsids)))
    return sub


def labf(beta, se, prior_sd):
    z2 = (beta / se) ** 2
    r = prior_sd ** 2 / (prior_sd ** 2 + se ** 2)
    return 0.5 * np.log(1 - r) + 0.5 * r * z2


def coloc(l1, l2):
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    h3_inner = s1 + s2 + np.log1p(-np.exp(min(s12 - (s1 + s2), -1e-12)))
    lh = np.array([0.0, np.log(P1) + s1, np.log(P2) + s2, np.log(P1) + np.log(P2) + h3_inner, np.log(P12) + s12])
    return np.exp(lh - logsumexp(lh))


gen = pd.read_csv(RESULTS / "03_genetic_targets_druggability.csv")
REFERENCE = [g for g in gen[gen.genetic >= 0.7].sort_values("genetic", ascending=False).symbol if g not in GENES]
coords = gene_coords(list(GENES) + REFERENCE)
jobs = [(d, g) for d in DATASETS for g in coords]
with cf.ThreadPoolExecutor(4) as ex:
    eq = dict(zip(jobs, ex.map(lambda j: eqtl_region(j[0], j[1], coords[j[1]]), jobs)))
var = pd.read_csv(RESULTS / "08_aga_variants.csv")
var = var[var.tgfb_link]
DIRECTION_CHECK = "rs6152:66765627:G:A"
rsids = set().union(*[set(v.rsid) for v in eq.values()]) | set(var.snp) | {DIRECTION_CHECK}
gw = gwas_subset(rsids).drop_duplicates("SNP").set_index("SNP")
print(f"eQTL variants: {len(rsids):,}; found in the baldness GWAS: {gw.index.isin(list(rsids)).sum():,}")
maf = np.minimum(gw.A1FREQ, 1 - gw.A1FREQ)
sdy_gwas = float(np.sqrt(np.median(gw.SE ** 2 * 2 * GWAS_N * maf * (1 - maf))))
print(f"baldness score SD estimated from the GWAS standard errors: {sdy_gwas:.3f}")
dc = gw.loc[DIRECTION_CHECK]
g_beta = dc.BETA if dc.ALLELE1 == "G" else -dc.BETA
print(f"direction check, AR rs6152: the known risk allele G has beta {g_beta:+.3f}, so a higher score means more balding")
assert g_beta > 0

rows = []
for (d, g), e in eq.items():
    m = e.set_index("rsid").join(gw, how="inner")
    m = m[np.isfinite(m.SE) & (m.SE > 0) & np.isfinite(m.se) & (m.se > 0)]
    if len(m) < 50:
        rows.append({"gene": g, "panel": "hypothesis" if g in GENES else "reference", "role": GENES.get(g, "strong baldness genetics"),
                     "tissue": DATASETS[d], "dataset": d, "n_snps": len(m)})
        continue
    n_eqtl = float(np.median(m.an)) / 2
    sdy_e = float(np.sqrt(np.median(m.se ** 2 * 2 * n_eqtl * m.maf * (1 - m.maf))))
    pp = coloc(labf(m.beta.values, m.se.values, 0.15 * sdy_e), labf(m.BETA.values, m.SE.values, 0.15 * sdy_gwas))
    lead = m.sort_values("pvalue").iloc[0]
    lead_rs = m.sort_values("pvalue").index[0]
    same = (lead.alt == lead.ALLELE1) and (lead.ref == lead.ALLELE0)
    flip = (lead.alt == lead.ALLELE0) and (lead.ref == lead.ALLELE1)
    palin = {lead.ref, lead.alt} in ({"A", "T"}, {"C", "G"}) and lead.maf > 0.4
    mr = {}
    if (same or flip) and not palin and lead.pvalue < 1e-5:
        bg = lead.BETA * (1 if same else -1)
        b = bg / lead.beta
        s = lead.SE / abs(lead.beta)
        mr = {"mr_variant": lead_rs, "mr_beta_per_sd_expression": b / sdy_gwas, "mr_se": s / sdy_gwas, "mr_p": 2 * stats.norm.sf(abs(b / s)),
              "mr_direction": "higher activity, more balding" if b > 0 else "higher activity, less balding"}
    gw_lead = m.sort_values("P_BOLT_LMM").iloc[0]
    verdict = ("shared variant (H4)" if pp[4] > 0.8 else "different variants (H3)" if pp[3] > 0.8 else
               "no gene-activity signal" if pp[1] + pp[3] + pp[4] < 0.5 else "no baldness signal here" if pp[2] + pp[3] + pp[4] < 0.5 else "unclear")
    rows.append({"gene": g, "panel": "hypothesis" if g in GENES else "reference", "role": GENES.get(g, "strong baldness genetics"),
                 "tissue": DATASETS[d], "dataset": d, "n_snps": len(m), "eqtl_lead": lead_rs,
                 "eqtl_lead_p": lead.pvalue, "eqtl_lead_beta": lead.beta, "gwas_lead_in_window": m.sort_values("P_BOLT_LMM").index[0],
                 "gwas_min_p_in_window": gw_lead.P_BOLT_LMM, "gwas_p_at_eqtl_lead": lead.P_BOLT_LMM,
                 **{f"PP.H{i}": round(float(v), 4) for i, v in enumerate(pp)}, "verdict": verdict, **mr})
res = pd.DataFrame(rows)
res["testable"] = res.eqtl_lead_p < 1e-5
res.to_csv(RESULTS / "10_coloc_mr.csv", index=False)
panel = {}
for name, grp in res.groupby("panel"):
    t = grp[grp.testable]
    panel[name] = {"gene_tissue_pairs": int(len(grp)), "testable": int(len(t)), "shared_variant_H4": int((t["PP.H4"] > 0.8).sum()),
                   "different_variants_H3": int((t["PP.H3"] > 0.8).sum()),
                   "genes_with_H4": sorted(t[t["PP.H4"] > 0.8].gene.unique().tolist())}
ph, pr = panel.get("hypothesis", {}), panel.get("reference", {})
if ph.get("testable") and pr.get("testable"):
    table = [[ph["shared_variant_H4"], ph["testable"] - ph["shared_variant_H4"]], [pr["shared_variant_H4"], pr["testable"] - pr["shared_variant_H4"]]]
    panel["fisher_hypothesis_vs_reference"] = {"odds_ratio": float(stats.fisher_exact(table)[0]), "p": float(stats.fisher_exact(table)[1]),
                                               "caveat": "reference genes were chosen by genetic scores that partly use colocalization evidence"}
print("\npanels (testable = activity variant p < 1e-5):", json.dumps(panel, indent=1))
show = ["gene", "panel", "tissue", "n_snps", "eqtl_lead_p", "gwas_min_p_in_window", "PP.H3", "PP.H4", "verdict", "mr_beta_per_sd_expression", "mr_p", "mr_direction"]
pd.set_option("display.width", 250)
print(res[[c for c in show if c in res.columns]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
vrows = []
for r in var.itertuples():
    ref, alt = r.linked_refalt.split(":")
    active = ref if r.active_allele == "ref" else alt
    row = {"snp": r.snp, "atlas_cluster": r.cluster, "linked_genes": r.linkedGenes, "motif": r.top_motifs, "fine_map_prob": r.fm_probs,
           "more_open_allele": active}
    if r.snp in gw.index:
        x = gw.loc[r.snp]
        b = x.BETA if active == x.ALLELE1 else -x.BETA if active == x.ALLELE0 else np.nan
        row.update({"gwas_beta_more_open_allele_sd": b / sdy_gwas, "gwas_p": x.P_BOLT_LMM,
                    "reading": "more open DNA, more balding" if b > 0 else "more open DNA, less balding" if b < 0 else "alleles do not match"})
    vrows.append(row)
vd = pd.DataFrame(vrows)
vd.to_csv(RESULTS / "10_atlas_variant_directions.csv", index=False)
print("\natlas variants tied to TGF-beta or the scarring genes, read against the baldness GWAS:")
print(vd.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
save_summary("10", {"gwas": "Yap 2018, GCST007020, 205,327 UK Biobank men", "baldness_score_sd": sdy_gwas, "priors": {"p1": P1, "p2": P2, "p12": P12},
                    "direction_check": {"variant": DIRECTION_CHECK, "risk_allele": "G", "beta_risk_allele": float(g_beta)},
                    "panels": panel, "results": res.replace({np.nan: None}).to_dict(orient="records"),
                    "atlas_variant_directions": vd.replace({np.nan: None}).to_dict(orient="records")})

t = res[res.testable].copy()
t["label"] = t.gene + "  (" + t.tissue.str.replace("fibroblasts (cultured, skin)", "fibroblasts").str.replace("skin, not sun-exposed", "skin").str.replace("skin, sun-exposed", "sun skin") + ")"
t = t.sort_values(["panel", "PP.H4"], ascending=[True, True])
fig, ax = plt.subplots(figsize=(7.5, 0.26 * len(t) + 1.2))
colors = t.panel.map({"hypothesis": "#C4572A", "reference": "#3B6FB6"})
ax.hlines(range(len(t)), 0, t["PP.H4"], color="#d5dde2", lw=1)
ax.scatter(t["PP.H4"], range(len(t)), c=colors, s=40, zorder=3, edgecolor="white", linewidth=1.5)
ax.axvline(0.8, color="#44525c", lw=1, ls="--")
ax.text(0.81, len(t) - 0.5, "shared variant\n(PP.H4 > 0.8)", fontsize=7, color="#44525c", va="top")
ax.set_yticks(range(len(t)))
ax.set_yticklabels(t.label, fontsize=7)
ax.set_xlim(0, 1.05)
ax.set_xlabel("probability that gene activity and baldness share one causal variant (coloc PP.H4)", fontsize=8)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#C4572A", label="TGF-beta / scarring genes"),
                   Line2D([], [], marker="o", ls="", color="#3B6FB6", label="reference: strongest baldness genes")],
          fontsize=7, loc="lower right", frameon=False)
ax.set_title("Genetic natural experiment: testable gene-tissue pairs (activity variant p < 1e-5)", fontsize=9, loc="left")
fig.tight_layout()
(ROOT / "figs").mkdir(exist_ok=True)
fig.savefig(ROOT / "figs" / "10_coloc_panels.png", dpi=200)
print("saved figs/10_coloc_panels.png")
