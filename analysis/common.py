"""Shared helpers for the numbered analysis scripts: paths, the consensus bald-scalp signature, and result saving."""
import json
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

import signatures as S  # noqa: E402


def consensus_bald_signature():
    """Stouffer-combined bald-scalp signature from GSE36169 (5 paired men) and GSE90594 (14 AGA vs 14 healthy),
    keeping genes that move the same direction in both."""
    b1 = S.bald_vs_haired().rename(columns={"lfc_bald": "lfc", "p_bald": "p"})
    b3 = S.alopecia_vs_healthy()
    inv = b1.join(b3, lsuffix="_g", rsuffix="_m", how="inner")
    inv = inv[np.sign(inv.lfc_g) == np.sign(inv.lfc_m)].copy()
    inv["z"] = (np.sign(inv.lfc_g) * stats.norm.isf(inv.p_g / 2) + np.sign(inv.lfc_m) * stats.norm.isf(inv.p_m / 2)) / np.sqrt(2)
    return inv


def save_summary(name, summary):
    path = RESULTS / f"{name}_summary.json"
    json.dump(summary, open(path, "w"), indent=1, default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    print(f"saved {path.relative_to(ROOT)}")
