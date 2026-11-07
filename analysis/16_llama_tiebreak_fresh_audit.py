"""Track I, version 4: does letting Llama 3.3 break ties beat the Gemini tiebreaker, on papers nobody has looked at?

The v3 held-out audit found Llama right on 80% of the reviewer disagreements and the Gemini tiebreaker right on 50%.
Switching the rule on that same audit would be tuning on the test set, so this script draws a fresh audit instead.
It excludes every PMID seen during development (data/llm_agent/dev_seen_pmids.json) and every PMID in the v3 held-out
audit. No model is called: version 4 only changes the tie rule, so it is scored from the stored v3 labels.

First run writes results/16_fresh_audit_blind.csv (claim, title, abstract, no model labels). Reference labels go in
results/16_fresh_audit_labels.csv (audit_id, audit_stance, audit_type, audit_confidence, note), made by a Claude
(Opus 5.5) subagent that only ever saw the blind file. Second run scores every labeler against them.

AI disclosure: I wrote and tested this script together with Claude (Opus 5.5). The fresh audit's reference labels
were made by a Claude subagent reading only the claim and the abstract, so they are Claude's labels, not a human
expert's.
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from common import RESULTS, ROOT

AGENT = ROOT / "data" / "llm_agent"
LABELS = RESULTS / "09_evidence_labels.csv"
BLIND = RESULTS / "16_fresh_audit_blind.csv"
REF = RESULTS / "16_fresh_audit_labels.csv"

df = pd.read_csv(LABELS, dtype={"pmid": str})
df.index.name = "audit_id"
frozen = json.loads((AGENT / "pubmed_frozen.json").read_text())
articles = {}
for v in frozen.values():
    if isinstance(v, dict) and isinstance(v.get("articles"), dict):
        articles.update(v["articles"])

dev_seen = set(map(str, json.loads((AGENT / "dev_seen_pmids.json").read_text())["pmids"]))
v3_audit = set(pd.read_csv(RESULTS / "09_audit_comparison_heldout.csv", dtype={"pmid": str}).pmid)
seen = dev_seen | v3_audit
fresh = df[~df.pmid.isin(seen) & df.has_abstract.astype(str).eq("True")]
dis = fresh.index[fresh.resolution != "agree"].to_numpy()
agr = fresh.index[fresh.resolution == "agree"].to_numpy()
rng = np.random.default_rng(16)
pick = np.concatenate([rng.choice(dis, min(20, len(dis)), replace=False), rng.choice(agr, 20, replace=False)])
print(f"pairs {len(df)} | PMIDs excluded {len(seen)} (dev {len(dev_seen)}, v3 audit {len(v3_audit)}) | "
      f"fresh pairs {len(fresh)}: {len(dis)} disagreements, {len(agr)} agreements | sampled {len(pick)}")

if not BLIND.exists():
    rows = []
    for k in sorted(pick):
        r = df.loc[k]
        art = articles[r.pmid]
        rows.append({"audit_id": int(k), "claim": r.claim, "pmid": r.pmid, "year": art.get("year", ""), "title": art["title"],
                     "abstract": art["abstract"]})
    pd.DataFrame(rows).to_csv(BLIND, index=False)
    print(f"wrote {BLIND.name}; fill {REF.name} blind, then rerun")

if REF.exists():
    au = pd.read_csv(REF).merge(df.reset_index(), on="audit_id")
    au["audit_stance"] = au.audit_stance.str.strip().str.lower()
    au["v4"] = np.where(au.resolution == "agree", au.final, au.A_stance)
    w = (df.resolution == "agree").mean()
    strata = au.resolution.eq("agree")
    acc = {}
    for name, col in (("v3 pipeline (Gemini breaks ties)", "final"), ("v4 pipeline (Llama breaks ties)", "v4"),
                      ("Llama 3.3 alone", "A_stance"), ("DeepSeek V3.1 alone", "B_stance"), ("GPT-4o-mini (shadow)", "S_stance")):
        hit = au[col] == au.audit_stance
        acc[name] = {"agree_stratum": float(hit[strata].mean()), "disagree_stratum": float(hit[~strata].mean()),
                     "population_weighted": float(w * hit[strata].mean() + (1 - w) * hit[~strata].mean()),
                     "clear_items_only": float(hit[au.audit_confidence.eq("clear")].mean())}
    seen_c = au.C_stance.notna()
    d = au[~strata]
    v3_hit, v4_hit = d.final == d.audit_stance, d.v4 == d.audit_stance
    only_v4, only_v3 = int((v4_hit & ~v3_hit).sum()), int((v3_hit & ~v4_hit).sum())
    old = pd.read_csv(RESULTS / "09_audit_comparison_heldout.csv")
    old_d = old[old.resolution != "agree"]
    mcnemar = {"v4_right_v3_wrong": only_v4, "v3_right_v4_wrong": only_v3,
               "exact_p": float(binomtest(only_v4, only_v4 + only_v3, 0.5).pvalue) if only_v4 + only_v3 else None,
               "disagreements_right": {"v3": int(v3_hit.sum()), "v4": int(v4_hit.sum()), "n": int(len(d))},
               "pooled_with_v3_audit": {"llama_right": int((old_d.A_stance == old_d.audit_stance).sum() + v4_hit.sum()),
                                        "tiebreak_pipeline_right": int((old_d.final == old_d.audit_stance).sum() + v3_hit.sum()),
                                        "n": int(len(old_d) + len(d))}}
    out = {"question": "Does letting Llama break ties beat the Gemini tiebreaker on a fresh audit?",
           "auditor": "Claude Opus 5.5 subagent, reading only the claim and the abstract, never the model labels",
           "n": int(len(au)), "n_agree": int(strata.sum()), "n_disagree": int((~strata).sum()),
           "n_clear": int(au.audit_confidence.eq("clear").sum()), "full_run_agreement_rate": float(w), "accuracy": acc,
           "gemini_on_disagreements_it_saw": float((au.C_stance[seen_c] == au.audit_stance[seen_c]).mean()) if seen_c.any() else None,
           "excluded_pmids": len(seen), "disagreements_mcnemar": mcnemar}
    keep = ["audit_id", "claim_id", "pmid", "audit_stance", "audit_confidence", "resolution", "final", "v4", "A_stance", "B_stance",
            "C_stance", "S_stance", "note"]
    au[keep].to_csv(RESULTS / "16_fresh_audit_comparison.csv", index=False)
    (RESULTS / "16_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
