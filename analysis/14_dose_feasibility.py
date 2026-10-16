"""Track N: dose feasibility. In cells, pirfenidone reversed the bald-scalp program at 1 mM (GSE226249) but not at about
0.1 mM (GSE178518). Can a pill or a skin gel get scalp tissue anywhere near that? Every input below is a number quoted
from a PubMed abstract (PMID given), or a stated textbook assumption with a range. This is a back-of-envelope check,
not a pharmacokinetic model."""
import numpy as np

from common import save_summary

MW = 185.22
MG_PER_L_PER_MM = MW

oral_cmax_mg_l = {"fed, capsules": 6.56, "fed, tablet": 7.64, "fasted, capsules": 12.5, "fasted, tablet": 12.6}
oral_source = ("PMID 28808905 (Pan 2017, Table 2 of the PMC full text): single 801 mg dose, the clinical dose taken 3 times a "
               "day; geometric mean Cmax 6,560 and 7,640 ng/mL fed, 12,500 and 12,600 ng/mL fasted (capsules, tablet). "
               "PMID 17906160 (Shi 2007, 200 to 600 mg doses) reports 9.2 mg/L fed and 13.0 mg/L fasted, the same range")
cell_doses_mm = {"no effect (GSE178518, 20 ug/mL)": 20 / MG_PER_L_PER_MM, "reversal (GSE226249)": 1.0}
skin_flux_ug_cm2_h = {"liposomal hydrogel, pig ear skin (PMID 26114208)": 10.9,
                      "8% nanoemulsion gel, pig skin, 404.9 ug/cm2 over 24 h (PMID 33006045)": 404.9 / 24,
                      "plain pirfenidone, about 5x below the liposomal gel (PMID 26114208)": 10.9 / 5}
eye_tissue_cmax_mg_g = {"conjunctiva": 9.62, "cornea": 9.64, "sclera": 2.13}
eye_source = "PMID 21866212: one 50 uL drop of 0.5% pirfenidone on the rabbit eye (250 ug); tissue Cmax in mg/g"
mechanism_note = ("PMID 36585256 (Ma 2023, the paper behind GSE226249): pirfenidone inhibits MRTF activation in primary human lung "
                  "fibroblasts with IC50 50 to 150 uM and > 90% maximal inhibition; patients' peak concentration is < 100 uM")
clinical_topical = ["PMID 21659848: 8% gel 3x/day for 6 months improved pediatric hypertrophic burn scars, no serious adverse effects",
                    "PMID 25533576: pirfenidone gel in localized scleroderma, phase II",
                    "PMID 30115530: topical pirfenidone sped healing of skin graft donor sites"]
skin_blood_flow_ml_h_cm2 = (0.8, 2.5)
blood_flow_note = ("assumption: resting skin gets about 250 mL/min across about 1.8 m2, so about 0.8 mL/h per cm2; the scalp is "
                   "better perfused, so the range goes up to about 3x that")

oral_mm = {k: v / MG_PER_L_PER_MM for k, v in oral_cmax_mg_l.items()}
gap = 1.0 / np.array(list(oral_mm.values()))
print("Oral pill, peak blood level:", {k: f"{v:.3f} mM" for k, v in oral_mm.items()}, f"| that is {gap.min():.0f} to {gap.max():.0f}x below the 1 mM cell dose")
print("Cell doses:", {k: f"{v:.2f} mM" for k, v in cell_doses_mm.items()})
eye_mm = {k: v * 1000 / MG_PER_L_PER_MM for k, v in eye_tissue_cmax_mg_g.items()}
print("Topical on the rabbit eye, peak tissue level at the site:", {k: f"{v:.0f} mM" for k, v in eye_mm.items()})
rows = []
for name, flux in skin_flux_ug_cm2_h.items():
    lo, hi = flux / skin_blood_flow_ml_h_cm2[1], flux / skin_blood_flow_ml_h_cm2[0]
    rows.append({"formulation": name, "flux_ug_cm2_h": round(flux, 1), "dermis_average_mM_low": round(lo / MG_PER_L_PER_MM, 3),
                 "dermis_average_mM_high": round(hi / MG_PER_L_PER_MM, 3)})
    print(f"Skin gel, {name}: flux {flux:.1f} ug/cm2/h -> whole-dermis average if cleared by blood flow "
          f"{lo / MG_PER_L_PER_MM:.3f} to {hi / MG_PER_L_PER_MM:.3f} mM")
verdict = ("Neither route clearly reaches the 1 mM cell-culture level as a whole-dermis average (about 0.01 to 0.1 mM by "
           "flux / blood flow, the same order as the pill). But the cell-culture level is probably not the bar: the pill "
           "slows lung fibrosis in patients at 0.035 to 0.068 mM (the 801 mg dose, fed to fasted), 15 to 28x below it. A gel concentrates drug at the "
           "surface and in follicle openings, where the upper-follicle scarring sits (site tissue reached tens of mM on "
           "the rabbit eye), and 8% gel has clinical antifibrotic effects in human skin. The real test is a dose-response "
           "in human scalp follicles (about 0.03 to 3 mM), which the proposed organ-culture experiment should include.")
print("\n" + verdict)
save_summary("14", {"oral_cmax_mM": oral_mm, "oral_source": oral_source, "cell_doses_mM": cell_doses_mm,
                    "oral_gap_to_1mM_x": [float(gap.min()), float(gap.max())], "topical_eye_site_mM": eye_mm, "eye_source": eye_source,
                    "skin_gel_dermis_average_estimates": rows, "skin_blood_flow_ml_h_cm2": skin_blood_flow_ml_h_cm2,
                    "blood_flow_assumption": blood_flow_note, "clinical_topical_evidence": clinical_topical,
                    "antifibrotic_mechanism_at_clinical_levels": mechanism_note, "verdict": verdict})
