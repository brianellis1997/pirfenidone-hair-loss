"""Track C: which genes have genetic evidence for androgenetic alopecia (Open Targets), and how many can existing drugs act on."""
import json
import time
import urllib.request

import pandas as pd

from common import RESULTS, save_summary

AGA = "MONDO_0005339"


def gql(query, variables=None):
    for attempt in range(4):
        try:
            req = urllib.request.Request("https://api.platform.opentargets.org/api/v4/graphql",
                                         data=json.dumps({"query": query, "variables": variables or {}}).encode(),
                                         headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                out = json.loads(r.read())
            if "errors" in out:
                raise RuntimeError(out["errors"][0]["message"][:300])
            return out["data"]
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 + 2 * attempt)


rows = []
for idx in range(3):
    d = gql("""query($i: Int!) { disease(efoId: "MONDO_0005339") { associatedTargets(page: {index: $i, size: 500}) { count rows { score
        target { id approvedSymbol approvedName } datatypeScores { id score } } } } }""", {"i": idx})["disease"]["associatedTargets"]
    for r in d["rows"]:
        ds = {x["id"]: x["score"] for x in r["datatypeScores"]}
        rows.append({"ensembl": r["target"]["id"], "symbol": r["target"]["approvedSymbol"], "name": r["target"]["approvedName"],
                     "overall": r["score"], "genetic": ds.get("genetic_association", 0), "literature": ds.get("literature", 0),
                     "known_drug": ds.get("known_drug", 0)})
    if len(rows) >= d["count"]:
        break
df = pd.DataFrame(rows)
df.to_csv(RESULTS / "03_aga_all_associated_targets.csv", index=False)
gen = df[df.genetic >= 0.5].sort_values("genetic", ascending=False).reset_index(drop=True)
out = []
for r in gen.itertuples():
    t = gql("""query($id: String!) { target(ensemblId: $id) { tractability { label modality value }
        drugAndClinicalCandidates { count rows { maxClinicalStage drug { name } } } } }""", {"id": r.ensembl})["target"]
    tr = t.get("tractability") or []
    dc = t.get("drugAndClinicalCandidates") or {}
    drugs = sorted({(x["drug"]["name"], x.get("maxClinicalStage")) for x in (dc.get("rows") or [])})
    out.append({"symbol": r.symbol, "genetic": round(r.genetic, 3), "name": r.name, "n_drugs": dc.get("count", 0) or 0,
                "approved_drugs": ";".join(sorted({n for n, s in drugs if s in ("APPROVAL", "PHASE_4")})),
                "sm_tractable": ";".join(sorted({x["label"] for x in tr if x["modality"] == "SM" and x["value"]})),
                "ab_tractable": ";".join(sorted({x["label"] for x in tr if x["modality"] == "AB" and x["value"]}))})
    time.sleep(0.15)
res = pd.DataFrame(out)
res.to_csv(RESULTS / "03_genetic_targets_druggability.csv", index=False)
fib_genes = ["TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TWIST1", "TWIST2", "PDGFA", "PDGFRB", "AGTR1", "ACE", "SMAD3", "COL1A1", "LOX", "SPARC"]
fib = df.sort_values("genetic", ascending=False).drop_duplicates("symbol").set_index("symbol").reindex(fib_genes)[["genetic", "literature", "overall"]]
fib.to_csv(RESULTS / "03_fibrosis_gene_genetics.csv")
summary = {"associated_targets": len(df), "genetic_ge_0.5": len(gen), "genetic_ge_0.7": int((df.genetic >= .7).sum()),
           "with_any_drug": int((res.n_drugs > 0).sum()), "with_approved_drug": int((res.approved_drugs != "").sum()),
           "approved_drug_targets": res[res.approved_drugs != ""][["symbol", "approved_drugs"]].to_dict(orient="records"),
           "fibrosis_gene_genetics": fib.fillna(0).round(2).to_dict(orient="index")}
print(json.dumps({k: v for k, v in summary.items() if k != "fibrosis_gene_genetics"}, indent=1))
print(fib.round(2))
save_summary("03", summary)
