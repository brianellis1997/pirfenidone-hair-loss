"""Track G: how much has already been published or trialed, per candidate drug and per mechanism (PubMed, ClinicalTrials.gov)."""
import json
import time
import urllib.parse
import urllib.request

import pandas as pd

from common import RESULTS, save_summary


def get(url):
    for attempt in range(5):
        try:
            time.sleep(0.45)
            return json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=90).read())
        except Exception:
            time.sleep(2 + 2 * attempt)


def pubmed(term):
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": term, "retmode": "json", "retmax": 0}))
    return int(r["esearchresult"]["count"])


HAIR = '(alopecia OR "hair loss" OR "hair growth" OR "hair follicle")'
AGA = '("androgenetic alopecia" OR "male pattern" OR "pattern hair loss")'
cands = {"TGF-beta / scarring": ["pirfenidone", "telmisartan", "losartan", "tranilast", "nintedanib", "galunisertib"],
         "IL-6": ["tocilizumab", "tofacitinib", "ruxolitinib", "baricitinib"],
         "PGD2": ["setipiprant", "fevipiprant", "ramatroban"],
         "Wnt brakes": ["lithium", "valproic acid", "cyclosporine", "WAY-316606", "romosozumab"],
         "reference": ["minoxidil", "finasteride", "dutasteride", "latanoprost"]}
rows = []
for group, drugs in cands.items():
    for d in drugs:
        rows.append({"group": group, "drug": d, "hair_papers": pubmed(f'"{d}"[tiab] AND {HAIR}[tiab]'),
                     "pattern_hair_loss_papers": pubmed(f'"{d}"[tiab] AND {AGA}[tiab]'), "all_papers": pubmed(f'"{d}"[tiab]')})
        print(rows[-1])
pd.DataFrame(rows).to_csv(RESULTS / "07_drug_novelty.csv", index=False)
mech = {"fibrosis + pattern hair loss": f'(fibrosis OR fibrotic OR "perifollicular fibrosis")[tiab] AND {AGA}[tiab]',
        "TGF-beta + pattern hair loss": f'("TGF-beta" OR "TGF-beta1" OR "transforming growth factor")[tiab] AND {AGA}[tiab]',
        "TGF-beta + catagen": '("TGF-beta" OR "transforming growth factor")[tiab] AND catagen[tiab]',
        "pirfenidone + skin/scar/topical": 'pirfenidone[tiab] AND (scar OR keloid OR skin OR topical)[tiab]'}
mech_counts = {k: pubmed(v) for k, v in mech.items()}
trials = {}
for q in ("pirfenidone alopecia", "pirfenidone hair"):
    r = get("https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode({"query.term": q, "countTotal": "true", "pageSize": 5}))
    trials[q] = {"count": r.get("totalCount"), "titles": [s["protocolSection"]["identificationModule"].get("briefTitle", "") for s in r.get("studies", [])]}
print(mech_counts, trials)
save_summary("07", {"drugs": rows, "mechanism_papers": mech_counts, "trials": trials})
