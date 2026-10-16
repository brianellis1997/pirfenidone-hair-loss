"""Track I: an LLM evidence agent. For each link in the hypothesis, and for each candidate drug, it pulls PubMed abstracts
and has two independent LLM reviewers from different model families label every abstract against the claim: relevance,
stance, evidence type, and a supporting sentence copied from the abstract. A third model breaks ties. Every quote is
checked against the abstract text, and each claim is graded with fixed rules written before the run.
Reliability is measured three ways: reviewer agreement (Cohen's kappa), the share of quotes found word for word, and
accuracy against a blind audit sample (results/09_audit_labels.csv, filled in by hand from results/09_audit_blind.csv).
PubMed results and LLM replies are frozen in data/llm_agent/ on the first run, so reruns are exact and free.
This is version 2. Version 1 (archived in results/09_v1/) read only the top 25 papers per claim, let "scar*" pull in lung
papers, and often scored "treatment X blocked the effect" as contradicting the effect. The PMIDs inspected while finding
those problems are excluded from the blind audit. On the 8 inspected papers, Llama 3.3 70B got 7 right, DeepSeek v3.1 6,
Gemini 2.5 Flash 5, GPT-4.1-mini 4 and GPT-4o-mini 2, so v2 uses Llama and DeepSeek as reviewers and Gemini as the
tiebreaker. GPT-4o-mini still labels every pair as a shadow reviewer with no vote, for comparison.
Version 3 (this file) adds a scope rule. The v2 blind audit (results/09_v2/) found the reviewers counting other organs (eye,
kidney, heart) as evidence for skin claims, and oral side effects as evidence against topical use. v3 is scored on a fresh
held-out audit of 40 pairs that excludes every PMID seen so far, and on the 60 v2 audit pairs as a development check."""
import concurrent.futures as cf
import hashlib
import json
import os
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from common import RESULTS, ROOT, save_summary

AGENT = ROOT / "data" / "llm_agent"
AGENT.mkdir(parents=True, exist_ok=True)
PUBMED_FROZEN = AGENT / "pubmed_frozen.json"
LLM_CACHE = AGENT / "llm_calls.jsonl"
AUDIT_BLIND = RESULTS / "09_audit_blind.csv"
AUDIT_LABELS = RESULTS / "09_audit_labels.csv"
REVIEWERS = ["meta-llama/llama-3.3-70b-instruct", "deepseek/deepseek-chat-v3.1"]
ADJUDICATOR = "google/gemini-2.5-flash"
SHADOW = "openai/gpt-4o-mini"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TGFB = '"TGF-beta"[tiab] OR "TGF-β"[tiab] OR "transforming growth factor"[tiab] OR TGFB1[tiab]'
SKIN = ('skin[tiab] OR dermal[tiab] OR cutaneous[tiab] OR keloid*[tiab] OR "hypertrophic scar*"[tiab] OR scleroderma[tiab] OR '
        '"systemic sclerosis"[tiab] OR "wound healing"[tiab]')
CHAIN = [
    ("C1 scarring in bald scalp", "Bald scalp in androgenetic alopecia shows fibrosis or scarring around hair follicles (perifollicular fibrosis, "
     "collagen deposition).", '("androgenetic alopecia"[tiab] OR "male pattern baldness"[tiab] OR "pattern hair loss"[tiab]) AND (fibrosis[tiab] '
     'OR fibrotic[tiab] OR perifollicular[tiab] OR collagen[tiab] OR microinflammation[tiab])', 60),
    ("C2 DHT raises TGF-beta", "Androgens such as DHT or testosterone increase TGF-beta production by dermal papilla cells or hair follicles.",
     f'(androgen*[tiab] OR dihydrotestosterone[tiab] OR testosterone[tiab]) AND ({TGFB}) AND ("dermal papilla"[tiab] OR "hair follicle"[tiab] '
     'OR alopecia[tiab])', 60),
    ("C3 TGF-beta stops hair growth", "TGF-beta signaling pushes hair follicles into regression (catagen) or inhibits hair growth.",
     f'({TGFB}) AND (catagen[tiab] OR "hair growth"[tiab] OR "hair cycle"[tiab] OR "hair follicle regression"[tiab])', 60),
    ("C4 pirfenidone is antifibrotic in skin", "Pirfenidone reduces TGF-beta-driven collagen production, fibroblast activation or fibrosis in "
     "skin, skin scars or dermal fibroblasts.", f'pirfenidone[tiab] AND ({SKIN})', 60),
    ("C5 topical pirfenidone is tolerated", "Pirfenidone applied topically to skin is safe and well tolerated.",
     'pirfenidone[tiab] AND (topical*[tiab] OR gel[tiab] OR cream[tiab] OR ointment[tiab])', 60),
]
CANDIDATES = ["pirfenidone", "tranilast", "telmisartan", "losartan", "nintedanib", "galunisertib"]
CAND_CLAIMS = [
    ("antifibrotic in skin", "{d} reduces TGF-beta signaling, collagen production or fibrosis in skin, skin scars or dermal fibroblasts.",
     f'{{d}}[tiab] AND ({TGFB} OR fibrosis[tiab] OR fibrotic[tiab] OR collagen[tiab]) AND ({SKIN})', 40),
    ("good for hair", "{d} promotes hair growth or protects against hair loss. (Evidence that it causes hair loss contradicts the claim.)",
     '{d}[tiab] AND (hair[tiab] OR alopecia[tiab] OR hypertrichosis[tiab] OR "hair follicle"[tiab])', 40),
]
STANCES = {"supports", "contradicts", "mixed", "none"}
TYPES = {"human_patients_or_tissue", "human_ex_vivo", "human_cells", "animal", "review", "other"}
SYSTEM = ("You are a careful biomedical evidence reviewer. You judge one PubMed abstract against one claim, using only what the "
          "abstract says. Reply with a single JSON object.")
PROMPT = """CLAIM: {claim}

ABSTRACT (PMID {pmid}, {year}, {journal}):
{title}
{abstract}

How to judge the stance:
- "supports": the abstract shows the effect in the claim happening, or states it as a finding or conclusion. This includes studies
  where a treatment blocks or reverses that effect, as long as they show the effect exists. Example: "drug X blocked the DHT-induced
  rise in TGF-beta" supports "DHT raises TGF-beta".
- "contradicts": the abstract reports that the effect is absent, or goes the opposite way. Example: TGF-beta promoting hair growth
  contradicts "TGF-beta inhibits hair growth".
- "mixed": the abstract reports both.
- "none": the abstract does not test or state the relationship in the claim, even if it uses the same words.
- Scope: the evidence has to match the claim's tissue and route. A study in another organ (lung, kidney, heart, eye) does not
  bear on a claim about skin, and effects of the drug taken by mouth do not bear on a claim about topical use. Use "none" for
  evidence outside the claim's scope. Animal and cell studies count when the tissue and route match; record their evidence type.

Return JSON with these keys:
"relevant": true if the abstract reports evidence (new data, or a review's stated conclusion) that bears directly on the claim, else false.
"stance": "supports", "contradicts", "mixed", or "none" (use "none" when not relevant).
"evidence_type": one of "human_patients_or_tissue" (people or biopsies), "human_ex_vivo" (organ-cultured human tissue), "human_cells" (cultured human cells), "animal", "review", "other".
"quote": one sentence copied exactly, word for word, from the abstract that best justifies the stance ("" if not relevant).
"reason": at most 20 words."""


def eutils(endpoint, **params):
    url = EU + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(6):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "hrp-evidence-agent"}), timeout=90) as r:
                data = r.read()
            time.sleep(0.45)
            return data
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503):
                time.sleep(2 + 2 * attempt)
                continue
            raise
    raise RuntimeError(f"PubMed failed: {url}")


def pubmed(query, n):
    res = json.loads(eutils("esearch.fcgi", db="pubmed", term=query, retmax=n, sort="relevance", retmode="json"))["esearchresult"]
    ids, arts = res["idlist"], {}
    if ids:
        root = ET.fromstring(eutils("efetch.fcgi", db="pubmed", id=",".join(ids), rettype="abstract", retmode="xml"))
        for a in root.findall(".//PubmedArticle"):
            t = a.find(".//ArticleTitle")
            date = a.find(".//JournalIssue/PubDate")
            arts[a.findtext(".//PMID")] = {
                "title": "".join(t.itertext()).strip() if t is not None else "",
                "abstract": " ".join("".join(x.itertext()).strip() for x in a.findall(".//Abstract/AbstractText")),
                "year": (date.findtext("Year") or (date.findtext("MedlineDate") or "")[:4]) if date is not None else "",
                "journal": a.findtext(".//Journal/ISOAbbreviation") or "",
                "types": [p.text for p in a.findall(".//PublicationType")]}
    return {"query": query, "count": int(res["count"]), "pmids": [i for i in ids if i in arts], "articles": arts}


def openrouter_key():
    if os.environ.get("OPENROUTER_API_KEY"):
        return os.environ["OPENROUTER_API_KEY"]
    raise RuntimeError("OPENROUTER_API_KEY is not set")


_cache, _lock = {}, threading.Lock()
if LLM_CACHE.exists():
    for line in LLM_CACHE.read_text().splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        _cache[rec["key"]] = rec


def parse_label(text):
    m = re.search(r"\{.*\}", text or "", re.S)
    if not m:
        return None
    try:
        d = json.loads(m.group(0))
    except json.JSONDecodeError:
        return None
    stance = str(d.get("stance", "")).strip().lower()
    etype = str(d.get("evidence_type", "")).strip().lower()
    rel = d.get("relevant")
    rel = rel if isinstance(rel, bool) else str(rel).strip().lower() == "true"
    if stance not in STANCES:
        return None
    return {"relevant": bool(rel) and stance != "none", "stance": stance if rel else "none",
            "evidence_type": etype if etype in TYPES else "other", "quote": str(d.get("quote", "")).strip(), "reason": str(d.get("reason", ""))[:200]}


def llm_label(model, claim, art, pmid):
    msgs = [{"role": "system", "content": SYSTEM}, {"role": "user", "content": PROMPT.format(
        claim=claim, pmid=pmid, year=art["year"], journal=art["journal"], title=art["title"], abstract=art["abstract"] or "(no abstract)")}]
    for seed in range(3):
        req = {"model": model, "messages": msgs, "temperature": 0, "seed": seed, "max_tokens": 1500 if model == ADJUDICATOR else 400,
               "response_format": {"type": "json_object"}, "provider": {"require_parameters": True}}
        key = hashlib.sha256(json.dumps(req, sort_keys=True).encode()).hexdigest()
        with _lock:
            rec = _cache.get(key)
        if rec is None:
            body = json.dumps(req).encode()
            for attempt in range(6):
                try:
                    r = urllib.request.Request("https://openrouter.ai/api/v1/chat/completions", data=body, headers={
                        "Content-Type": "application/json", "Authorization": "Bearer " + openrouter_key()})
                    out = json.loads(urllib.request.urlopen(r, timeout=120).read())
                    break
                except urllib.error.HTTPError as e:
                    if attempt == 0:
                        print(f"OpenRouter {model}: HTTP {e.code} {e.read()[:300]}")
                    time.sleep(3 + 3 * attempt)
                except (urllib.error.URLError, TimeoutError):
                    time.sleep(3 + 3 * attempt)
            else:
                raise RuntimeError(f"OpenRouter failed for {model}")
            rec = {"key": key, "model": model, "text": out["choices"][0]["message"].get("content") or "",
                   "usage": out.get("usage", {})}
            with _lock:
                with LLM_CACHE.open("a") as f:
                    f.write(json.dumps(rec) + "\n")
                _cache[key] = rec
        lab = parse_label(rec["text"])
        if lab:
            return lab
    return {"relevant": False, "stance": "invalid", "evidence_type": "other", "quote": "", "reason": "unparseable"}


def norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9%.,;:()/+\- ]", " ", s.lower())).strip(" .")


def grounded(quote, art):
    q = norm(quote)
    return bool(q) and q in norm(art["title"] + " " + art["abstract"])


def kappa(a, b):
    a, b = pd.Series(a), pd.Series(b)
    cats = sorted(set(a) | set(b))
    po = (a.values == b.values).mean()
    pe = sum((a == c).mean() * (b == c).mean() for c in cats)
    return (po - pe) / (1 - pe) if pe < 1 else 1.0


def grade(rows):
    prim = rows[rows.final.isin(["supports", "contradicts", "mixed"]) & ~rows.final_type.isin(["review", "other"])]
    n_sup, n_con, n_mix = (prim.final == "supports").sum(), (prim.final == "contradicts").sum(), (prim.final == "mixed").sum()
    n_hum = ((prim.final == "supports") & prim.final_type.isin(["human_patients_or_tissue", "human_ex_vivo"])).sum()
    n_rev = ((rows.final == "supports") & (rows.final_type == "review")).sum()
    tot = n_sup + n_con + n_mix
    con_share = (n_con + 0.5 * n_mix) / tot if tot else 0.0
    if tot >= 3 and con_share > 1 / 3:
        g = "Contested"
    elif n_sup >= 3 and n_hum >= 1 and con_share <= 0.2:
        g = "Strong"
    elif n_sup >= 2 and con_share <= 1 / 3:
        g = "Moderate"
    elif n_sup >= 1 or n_rev >= 1:
        g = "Weak"
    else:
        g = "None"
    return {"grade": g, "supporting_primary": int(n_sup), "supporting_human": int(n_hum), "contradicting": int(n_con), "mixed": int(n_mix),
            "supporting_reviews": int(n_rev), "relevant": int(rows.final.isin(["supports", "contradicts", "mixed"]).sum()), "screened": int(len(rows))}


claims = [(cid, text, q, n) for cid, text, q, n in CHAIN]
claims += [(f"{d} | {cid}", text.format(d=d.capitalize()), q.replace("{d}", d), n) for d in CANDIDATES for cid, text, q, n in CAND_CLAIMS]
frozen = json.loads(PUBMED_FROZEN.read_text()) if PUBMED_FROZEN.exists() else {}
changed = False
for cid, text, q, n in claims:
    if cid not in frozen:
        frozen[cid] = pubmed(q, n)
        changed = True
        print(f"PubMed {cid}: {frozen[cid]['count']} hits, kept {len(frozen[cid]['pmids'])}")
if changed:
    frozen.setdefault("_frozen_on", time.strftime("%Y-%m-%d"))
    PUBMED_FROZEN.write_text(json.dumps(frozen))
print(f"PubMed set frozen on {frozen['_frozen_on']}")

pairs = [(cid, text, pmid, frozen[cid]["articles"][pmid]) for cid, text, q, n in claims for pmid in frozen[cid]["pmids"]]
jobs = {(i, m): (m, p[1], p[3], p[2]) for i, p in enumerate(pairs) for m in REVIEWERS + [SHADOW]}
with cf.ThreadPoolExecutor(16) as ex:
    labels = dict(zip(jobs, ex.map(lambda j: llm_label(*j), jobs.values())))
rows = []
for i, (cid, text, pmid, art) in enumerate(pairs):
    a, b = labels[(i, REVIEWERS[0])], labels[(i, REVIEWERS[1])]
    rows.append({"claim_id": cid, "claim": text, "pmid": pmid, "year": art["year"], "title": art["title"], "has_abstract": bool(art["abstract"]),
                 "A_stance": a["stance"], "A_type": a["evidence_type"], "A_quote": a["quote"], "A_grounded": grounded(a["quote"], art) if a["quote"] else None,
                 "B_stance": b["stance"], "B_type": b["evidence_type"], "B_quote": b["quote"], "B_grounded": grounded(b["quote"], art) if b["quote"] else None,
                 "S_stance": labels[(i, SHADOW)]["stance"]})
df = pd.DataFrame(rows)
need = df.index[(df.A_stance != df.B_stance) | (df.A_stance.isin(["supports", "contradicts", "mixed"]) & (df.A_type != df.B_type))].tolist()
with cf.ThreadPoolExecutor(16) as ex:
    adj = dict(zip(need, ex.map(lambda k: llm_label(ADJUDICATOR, pairs[k][1], pairs[k][3], pairs[k][2]), need)))
df["C_stance"] = [adj[k]["stance"] if k in adj else None for k in df.index]
df["C_type"] = [adj[k]["evidence_type"] if k in adj else None for k in df.index]
df["C_quote"] = [adj[k]["quote"] if k in adj else None for k in df.index]
df["C_grounded"] = [grounded(adj[k]["quote"], pairs[k][3]) if k in adj and adj[k]["quote"] else None for k in df.index]


def resolve(r):
    if r.A_stance == r.B_stance:
        return r.A_stance, "agree", r.A_type if r.A_type == r.B_type else r.C_type, r.A_quote if r.A_grounded else r.B_quote
    if r.C_stance in (r.A_stance, r.B_stance):
        w = "A" if r.C_stance == r.A_stance else "B"
        return r.C_stance, "majority", r.C_type, r[f"{w}_quote"] if r[f"{w}_grounded"] else r.C_quote
    return "unresolved", "three-way split", r.C_type, ""


df[["final", "resolution", "final_type", "final_quote"]] = df.apply(lambda r: pd.Series(resolve(r)), axis=1)
df.to_csv(RESULTS / "09_evidence_labels.csv", index=False)

quotes = pd.concat([df[f"{m}_grounded"].dropna().astype(float) for m in "ABC"])
rel = {"pairs_screened": int(len(df)), "abstracts": int(df.pmid.nunique()), "missing_abstracts": int((~df.has_abstract).sum()),
       "reviewer_agreement": float((df.A_stance == df.B_stance).mean()), "cohen_kappa_stance": float(kappa(df.A_stance, df.B_stance)),
       "cohen_kappa_relevant": float(kappa(df.A_stance != "none", df.B_stance != "none")),
       "adjudicated": int(len(need)), "three_way_split": int((df.resolution == "three-way split").sum()),
       "quotes_found_verbatim": float(quotes.mean()), "quotes_checked": int(len(quotes)),
       "quotes_found_by_model": {REVIEWERS[0]: float(df.A_grounded.dropna().astype(float).mean()), REVIEWERS[1]: float(df.B_grounded.dropna().astype(float).mean()),
                                 ADJUDICATOR: float(df.C_grounded.dropna().astype(float).mean()) if df.C_grounded.notna().any() else None},
       "invalid_replies": int((df[["A_stance", "B_stance"]] == "invalid").sum().sum()),
       "shadow_agreement_with_final": float((df.S_stance == df.final).mean()), "shadow_kappa_vs_final": float(kappa(df.S_stance, df.final))}
print("\nagent reliability:", json.dumps(rel, indent=1))

grades = pd.DataFrame([{"claim_id": cid, "claim": text, "pubmed_hits": frozen[cid]["count"], **grade(df[df.claim_id == cid])}
                       for cid, text, q, n in claims])
best = (df[df.final.isin(["supports", "contradicts", "mixed"]) & (df.final_quote.str.len() > 0)]
        .assign(rank=lambda x: x.final_type.map({"human_patients_or_tissue": 0, "human_ex_vivo": 1, "human_cells": 2, "animal": 3, "review": 4, "other": 5}))
        .sort_values(["claim_id", "rank", "year"], ascending=[True, True, False]))
grades["key_citations"] = grades.claim_id.map(lambda c: "; ".join(
    f"PMID {r.pmid} ({r.year}, {r.final_type}, {r.final}): {r.final_quote[:160]}" for r in best[best.claim_id == c].head(3).itertuples()))
grades.to_csv(RESULTS / "09_claim_grades.csv", index=False)
print("\n" + grades[["claim_id", "grade", "supporting_primary", "supporting_human", "contradicting", "mixed", "supporting_reviews", "relevant",
                     "screened", "pubmed_hits"]].to_string(index=False))


def evaluate(labels_path, pool, tag):
    au = pd.read_csv(labels_path).merge(df.reset_index().rename(columns={"index": "audit_id"}), on="audit_id")
    au["audit_stance"] = au.audit_stance.str.strip().str.lower()
    stratum = au.resolution.eq("agree")
    acc = {}
    for name, col in (("pipeline", "final"), (REVIEWERS[0], "A_stance"), (REVIEWERS[1], "B_stance"), (SHADOW + " (shadow, no vote)", "S_stance")):
        hit = au[col] == au.audit_stance
        a = {"overall": float(hit.mean())}
        if pool is not None and stratum.any() and (~stratum).any():
            w = (pool.resolution == "agree").mean()
            a.update({"agree_stratum": float(hit[stratum].mean()), "disagree_stratum": float(hit[~stratum].mean())})
            a["population_weighted"] = w * a["agree_stratum"] + (1 - w) * a["disagree_stratum"]
        if "audit_confidence" in au:
            a["clear_items_only"] = float(hit[au.audit_confidence.eq("clear")].mean())
        acc[name] = a
    au["pipeline_error"] = np.where(au.final == au.audit_stance, "", au.final + " (audit: " + au.audit_stance + ")")
    keep = ["audit_id", "claim_id", "pmid", "audit_stance", "final", "A_stance", "B_stance", "C_stance", "S_stance", "resolution", "pipeline_error"]
    au[keep + [c for c in ("audit_confidence", "note") if c in au]].to_csv(RESULTS / f"09_audit_comparison_{tag}.csv", index=False)
    seen_c = au.C_stance.notna()
    return {"n": int(len(au)), "n_clear": int(au.audit_confidence.eq("clear").sum()) if "audit_confidence" in au else None,
            "pipeline_errors": int((au.pipeline_error != "").sum()), "accuracy": acc,
            "adjudicator_accuracy_on_items_it_saw": float((au.C_stance[seen_c] == au.audit_stance[seen_c]).mean()) if seen_c.any() else None}


rng = np.random.default_rng(1)
seen_path = AGENT / "dev_seen_pmids.json"
seen = set(json.loads(seen_path.read_text())["pmids"]) if seen_path.exists() else set()
eligible = df[~df.pmid.astype(str).isin(seen)]
dis = eligible.index[eligible.resolution != "agree"].to_numpy()
agr = eligible.index[eligible.resolution == "agree"].to_numpy()
pick = np.concatenate([rng.choice(dis, min(20, len(dis)), replace=False), rng.choice(agr, 20, replace=False)])
blind = pd.DataFrame([{"audit_id": int(k), "claim": pairs[k][1], "pmid": pairs[k][2], "title": pairs[k][3]["title"], "abstract": pairs[k][3]["abstract"]}
                      for k in sorted(pick)])
if not AUDIT_BLIND.exists():
    blind.to_csv(AUDIT_BLIND, index=False)
summary = {"models": {"reviewers": REVIEWERS, "adjudicator": ADJUDICATOR, "shadow": SHADOW}, "pubmed_frozen_on": frozen["_frozen_on"], "reliability": rel,
           "grades": grades.set_index("claim_id")[["grade", "supporting_primary", "supporting_human", "contradicting", "mixed",
                                                   "supporting_reviews", "relevant", "screened", "pubmed_hits"]].to_dict(orient="index")}
dev_labels = RESULTS / "09_v2" / "09_audit_labels.csv"
if dev_labels.exists():
    summary["dev_check_on_v2_audit_pairs"] = evaluate(dev_labels, None, "dev_v2_pairs")
    print("\ndevelopment check on the 60 v2 audit pairs (seen while designing v3, so optimistic):",
          json.dumps(summary["dev_check_on_v2_audit_pairs"]["accuracy"], indent=1))
if AUDIT_LABELS.exists():
    summary["audit"] = {**evaluate(AUDIT_LABELS, eligible, "heldout"), "auditor": "Claude Opus 5.5, reading each abstract blind to the model labels"}
    print("\nheld-out audit accuracy:", json.dumps(summary["audit"], indent=1))
else:
    print(f"\nblind audit sample written to {AUDIT_BLIND.relative_to(ROOT)}; fill {AUDIT_LABELS.name} (audit_id, audit_stance, audit_type, audit_confidence, note)")
save_summary("09", summary)
