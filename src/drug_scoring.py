"""
Scores every drug in public drug-response libraries on how strongly it reverses a disease gene signature, and
calibrates the scores against drugs known to grow hair or cause hair loss (SIDER side effects).
"""
import gzip
import pathlib
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent.parent
DRUGS = HERE / "data" / "drugs"

SALTS = r"\b(hydrochloride|dihydrochloride|hcl|sodium|potassium|calcium|maleate|mesylate|mesilate|citrate|tartrate|" \
        r"sulfate|sulphate|acetate|phosphate|fumarate|succinate|bromide|chloride|hydrobromide|besylate|tosylate|" \
        r"nitrate|lactate|monohydrate|dihydrate|hemihydrate|anhydrous|free base|salt)\b"


def norm(name):
    n = name.lower().strip()
    n = re.sub(r"\(.*?\)", " ", n)
    n = re.sub(SALTS, " ", n)
    n = re.sub(r"[^a-z0-9]+", "", n)
    return n


def read_gmt(path):
    sets = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = {g.split(",")[0] for g in parts[2:] if g}
    return sets


def lincs_consensus():
    sets = read_gmt(DRUGS / "LINCS_L1000_Chem_Pert_Consensus_Sigs.gmt")
    drugs = defaultdict(dict)
    for term, genes in sets.items():
        name, direction = term.rsplit(" ", 1)
        drugs[name][direction.lower()] = genes
    return {d: (v["up"], v["down"]) for d, v in drugs.items() if "up" in v and "down" in v}


def old_cmap():
    up, down = read_gmt(DRUGS / "Old_CMAP_up.gmt"), read_gmt(DRUGS / "Old_CMAP_down.gmt")
    inst = defaultdict(list)
    for term in up:
        if term in down:
            inst[term.rsplit("-", 1)[0]].append((up[term], down[term]))
    return inst


def reversal_z(sig_up, sig_dn, drug_up, drug_dn, universe):
    sig_up, sig_dn = sig_up & universe, sig_dn & universe
    drug_up, drug_dn = drug_up & universe, drug_dn & universe
    n = len(universe)
    pairs_rev = [(sig_up, drug_dn), (sig_dn, drug_up)]
    pairs_mim = [(sig_up, drug_up), (sig_dn, drug_dn)]
    obs = sum(len(a & b) for a, b in pairs_rev) - sum(len(a & b) for a, b in pairs_mim)
    exp_each = [len(a) * len(b) / n for a, b in pairs_rev + pairs_mim]
    var = sum(e * (1 - e / max(len(a), 1)) for e, (a, b) in zip(exp_each, pairs_rev + pairs_mim))
    exp = sum(exp_each[:2]) - sum(exp_each[2:])
    return (obs - exp) / np.sqrt(max(var, 1e-9))


def score_library(sig_up, sig_dn, library="lincs"):
    if library == "lincs":
        lib = lincs_consensus()
        universe = set().union(*[u | d for u, d in lib.values()])
        rows = [(drug, reversal_z(sig_up, sig_dn, u, d, universe)) for drug, (u, d) in lib.items()]
    else:
        lib = old_cmap()
        universe = set().union(*[u | d for insts in lib.values() for u, d in insts])
        rows = [(drug, float(np.mean([reversal_z(sig_up, sig_dn, u, d, universe) for u, d in insts]))) for drug, insts in lib.items()]
    df = pd.DataFrame(rows, columns=["drug", "z"])
    df["key"] = df.drug.map(norm)
    return df.sort_values("z", ascending=False).drop_duplicates("key").reset_index(drop=True)


def sider_labels():
    names = pd.read_csv(DRUGS / "drug_names.tsv", sep="\t", header=None, names=["cid", "name"])
    se = pd.read_csv(DRUGS / "meddra_all_se.tsv.gz", sep="\t", header=None,
                     names=["cid", "cid_stereo", "umls", "type", "umls_meddra", "effect"])
    se = se[se.type == "PT"]
    se["effect"] = se.effect.str.lower()
    grow = {"hypertrichosis", "growth of eyelashes", "eyelash thickening"}
    loss = {"alopecia", "androgenetic alopecia", "diffuse alopecia"}
    g = set(se[se.effect.isin(grow)].cid)
    l = set(se[se.effect.isin(loss)].cid)
    names["key"] = names.name.map(norm)
    lab = {}
    for cid, key in zip(names.cid, names.key):
        if cid in g and cid not in l:
            lab[key] = 1
        elif cid in l and cid not in g:
            lab[key] = 0
    return lab


def calibration(scores, labels, n_perm=10000, seed=0):
    df = scores[scores.key.isin(labels)].copy()
    df["y"] = df.key.map(labels)
    if df.y.nunique() < 2:
        return None
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(df.y, df.z)
    rng = np.random.default_rng(seed)
    perm = np.array([roc_auc_score(rng.permutation(df.y.values), df.z.values) for _ in range(n_perm)])
    return {"auroc": auc, "p_perm": float((np.sum(perm >= auc) + 1) / (n_perm + 1)), "n_pos": int(df.y.sum()),
            "n_neg": int((df.y == 0).sum()), "positives": df[df.y == 1].sort_values("z", ascending=False)[["drug", "z"]]}
