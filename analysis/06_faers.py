"""Track F: real-world evidence from the FDA adverse event reporting system (openFDA). Reporting odds ratios for hair loss
and for excess hair growth, calibrated on known drugs, plus an active comparator: pirfenidone vs nintedanib (both used
for lung fibrosis, so similar patients)."""
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

from common import RESULTS, save_summary

B = "https://api.fda.gov/drug/event.json?"
LOSS = 'patient.reaction.reactionmeddrapt:"alopecia"'
GROW = ('(patient.reaction.reactionmeddrapt:"hypertrichosis" OR patient.reaction.reactionmeddrapt:"hair growth abnormal" OR '
        'patient.reaction.reactionmeddrapt:"growth of eyelashes" OR patient.reaction.reactionmeddrapt:"hair growth increased")')


def total(search):
    for attempt in range(5):
        try:
            time.sleep(1.6)
            with urllib.request.urlopen(urllib.request.Request(B + urllib.parse.urlencode({"search": search, "limit": 1}),
                                                               headers={"User-Agent": "Mozilla/5.0"}), timeout=90) as r:
                return json.loads(r.read())["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 0
            time.sleep(3 + 3 * attempt)
        except Exception:
            time.sleep(3 + 3 * attempt)
    return 0


def ror(a, n_drug, n_event, n_all):
    b, c = n_drug - a, n_event - a
    d = n_all - n_drug - c
    if min(a, b, c, d) <= 0:
        return float("nan"), float("nan"), float("nan")
    r = (a / b) / (c / d)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return r, math.exp(math.log(r) - 1.96 * se), math.exp(math.log(r) + 1.96 * se)


ALL = total("receivedate:[20040101 TO 20261231]")
n_loss, n_grow = total(LOSS), total(GROW)
groups = {"growth_control": ["minoxidil", "cyclosporine", "diazoxide", "latanoprost", "bimatoprost"],
          "loss_control": ["paclitaxel", "doxorubicin", "docetaxel", "leflunomide"],
          "candidate": ["pirfenidone", "nintedanib", "telmisartan", "candesartan", "losartan", "valsartan", "irbesartan", "olmesartan"],
          "hormonal": ["finasteride", "dutasteride", "spironolactone"]}
rows = []
for group, names in groups.items():
    for d in names:
        q = f'patient.drug.openfda.generic_name:"{d.upper()}"'
        n_d = total(q)
        a_l, a_g = total(f"{q} AND {LOSS}"), total(f"{q} AND {GROW}")
        rl, ll, ul = ror(a_l, n_d, n_loss, ALL)
        rg, lg, ug = ror(a_g, n_d, n_grow, ALL)
        rows.append({"group": group, "drug": d, "reports": n_d, "alopecia": a_l, "ror_alopecia": rl, "ror_alopecia_lo": ll, "ror_alopecia_hi": ul,
                     "hair_growth": a_g, "ror_growth": rg, "ror_growth_lo": lg, "ror_growth_hi": ug})
        print(f"{group:15} {d:14} reports {n_d:7} | alopecia ROR {rl:6.2f} [{ll:.2f}, {ul:.2f}] | hair growth ROR {rg:6.2f} [{lg:.2f}, {ug:.2f}]")
pd.DataFrame(rows).to_csv(RESULTS / "06_faers_disproportionality.csv", index=False)
comp = []
for label, extra in [("all", ""), ("men", "patient.patientsex:1"), ("women", "patient.patientsex:2"),
                     ("age_60_85", "patient.patientonsetage:[60 TO 85] AND patient.patientonsetageunit:801")]:
    vals = {}
    for d in ("pirfenidone", "nintedanib"):
        q = f'patient.drug.openfda.generic_name:"{d.upper()}"' + (f" AND {extra}" if extra else "")
        vals[d] = (total(q), total(f"{q} AND {LOSS}"))
    (n1, a1), (n2, a2) = vals["pirfenidone"], vals["nintedanib"]
    r = (a1 / (n1 - a1)) / (a2 / (n2 - a2))
    se = math.sqrt(1 / a1 + 1 / (n1 - a1) + 1 / a2 + 1 / (n2 - a2))
    comp.append({"stratum": label, "pirf_reports": n1, "pirf_alopecia": a1, "nint_reports": n2, "nint_alopecia": a2,
                 "ratio": r, "lo": math.exp(math.log(r) - 1.96 * se), "hi": math.exp(math.log(r) + 1.96 * se)})
    print(f"pirfenidone vs nintedanib, {label:9}: {a1}/{n1} vs {a2}/{n2} -> ratio {r:.2f} [{comp[-1]['lo']:.2f}, {comp[-1]['hi']:.2f}]")
pd.DataFrame(comp).to_csv(RESULTS / "06_faers_pirfenidone_vs_nintedanib.csv", index=False)
save_summary("06", {"total_reports": ALL, "alopecia_reports": n_loss, "hair_growth_reports": n_grow, "drugs": rows, "active_comparator": comp})
