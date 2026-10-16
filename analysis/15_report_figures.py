"""Report figures, built only from results/ so they always match the numbers in the log. Vector PDF for the ACM report
and PNG (300 dpi) for slides, in figs/report/. Widths fit ACM two-column (3.35 in single, 7.0 in double)."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import textwrap

from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FuncFormatter, LogLocator

from common import RESULTS, ROOT

OUT = ROOT / "figs" / "report"
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GRAY, LIGHT, INK, MUTED = "#3B6FB6", "#C4572A", "#8A949C", "#D5DDE2", "#1F2A33", "#5B6770"
SINGLE, DOUBLE = 3.35, 7.0
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "pdf.fonttype": 42,
                     "axes.titlelocation": "left", "axes.titleweight": "bold"})


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def vgrid(ax):
    ax.grid(axis="x", color=LIGHT, linewidth=0.5)
    ax.set_axisbelow(True)


def plain_log(ax):
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def num(v, fmt="+.2f"):
    return format(0.0 if abs(v) < 0.005 else v, fmt)


def fig1_workflow():
    fig, ax = plt.subplots(figsize=(DOUBLE, 3.3))
    ax.set_xlim(-0.6, 101)
    ax.set_ylim(0, 50)
    ax.axis("off")
    stages = [
        ("1  What drives it,\n    besides DHT?", ["Bald vs haired scalp, 2 studies (GEO)", "DHT on follicle cells (GEO)", "Open Targets genetics"],
         "Scarring (TGF-beta) program, survives the follicle-loss adjustment"),
        ("2  Does a drug\n    reverse it?", ["Fibroblasts with TGF-beta and pirfenidone at 0.11 and 1 mM (GEO)", "CXCL12 check (2025 mouse study)"],
         "Pirfenidone reverses it at 1 mM and cuts CXCL12 3.8x"),
        ("3  Right cells,\n    right people?", ["Scalp single-cell atlas", "Beard vs scalp, one man", "FDA reports, 20.7M (FAERS)"],
         "Sheath, papilla and fibroblasts; 70% fewer hair-loss reports than nintedanib"),
        ("4  Does it hold up?", ["Baldness GWAS + eQTL", "LLM evidence agent, blind-audited", "Community scan; dose check"],
         "Mechanism Strong in the literature; genetics weak; dose uncertain; untested"),
    ]
    x0, w, gap = 0.4, 22.2, 3.6
    for i, (title, inputs, result) in enumerate(stages):
        x = x0 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 1.5), w, 46, boxstyle="round,pad=0.3,rounding_size=1.2", fc="white", ec=LIGHT, lw=0.8))
        ax.text(x + 1.2, 46, title, fontsize=6.7, fontweight="bold", va="top", color=INK, linespacing=1.15)
        y = 37.5
        for line in inputs:
            wrapped = textwrap.fill(line, 22, subsequent_indent="  ", break_on_hyphens=False)
            ax.text(x + 1.2, y, "\u2022 " + wrapped, fontsize=6.1, va="top", color=MUTED, linespacing=1.2)
            y -= 3.6 + 3.3 * wrapped.count("\n")
        ax.add_patch(FancyBboxPatch((x + 0.9, 2.8), w - 1.8, 12.4, boxstyle="round,pad=0.2,rounding_size=0.8", fc="#F3F5F7", ec="none"))
        ax.text(x + 1.6, 14.2, textwrap.fill(result, 22, break_on_hyphens=False), fontsize=6.1, va="top", color=INK, linespacing=1.2)
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + w + 0.8, 25), (x + w + gap - 0.8, 25), arrowstyle="-|>", mutation_scale=10,
                                         color=GRAY, lw=1.4, shrinkA=0, shrinkB=0))
    ax.text(0.5, -2.5, "Negative results reported alongside: drug-signature screening failed calibration; 92% of the non-hormonal "
            "genetic targets have no drug;\nblood-pressure drugs (ARBs) dropped; the genetic natural experiment did not confirm "
            "the pathway.", fontsize=6.1, color=MUTED, va="top", linespacing=1.25)
    save(fig, "fig1_workflow")


def fig2_calibration():
    cal = pd.read_csv(RESULTS / "02_calibration.csv").dropna(subset=["auroc"])
    names = {"lincs": "LINCS L1000", "cmap": "Connectivity Map"}
    sigs = {"invivo_100": "bald scalp, top 100", "invivo_200": "bald scalp, top 200", "dht_100": "DHT response, top 100", "dht_200": "DHT response, top 200"}
    cal = cal[cal.library.isin(names)].copy()
    cal["label"] = cal.library.map(names) + "\n" + cal.signature.map(sigs)
    cal = cal.iloc[::-1]
    fig, ax = plt.subplots(figsize=(SINGLE, 2.9))
    colors = [ORANGE if a < 0.5 else (BLUE if p < 0.05 else GRAY) for a, p in zip(cal.auroc, cal.p_perm)]
    ax.hlines(range(len(cal)), 0.5, cal.auroc, color=LIGHT, lw=1.2)
    ax.scatter(cal.auroc, range(len(cal)), color=colors, s=22, zorder=3, edgecolor="white", linewidth=1)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls="--")
    for y, (a, p) in enumerate(zip(cal.auroc, cal.p_perm)):
        ax.text(1.03, y, f"{a:.2f}" + ("  inverted" if a < 0.5 else f"  p {p:.2f}"), transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=6.3, color=ORANGE if a < 0.5 else INK)
    ax.set_yticks(range(len(cal)))
    ax.set_yticklabels(cal.label, fontsize=6.4, linespacing=1.1)
    ax.set_xlim(0.25, 0.75)
    ax.set_xticks([0.3, 0.5, 0.7])
    ax.set_xlabel("AUROC, hair-growth vs hair-loss drugs\n(0.5 = chance)")
    ax.set_title("Signature screening fails its own calibration", x=-0.66, ha="left")
    fig.subplots_adjust(left=0.36, right=0.80, top=0.91, bottom=0.17)
    vgrid(ax)
    save(fig, "fig2_calibration")


def fig3_scarring():
    comp = pd.read_csv(RESULTS / "04_composition_check_GSE90594.csv")
    summ = json.load(open(RESULTS / "04_summary.json"))["composition_GSE90594"]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.4), gridspec_kw={"width_ratios": [21, 13]})
    for ax, group, col, label in ((axes[0], "fibrosis", ORANGE, "Scarring genes"), (axes[1], "immune", BLUE, "Immune genes")):
        g = comp[comp.group == group].sort_values("bald_effect").reset_index(drop=True)
        n = summ[group]
        for i, r in g.iterrows():
            ax.plot([r.bald_effect, r.bald_effect_adj], [i, i], color=LIGHT, lw=1.6, zorder=1)
            ax.scatter(r.bald_effect, i, s=14, color=GRAY, zorder=2)
            sig = r.p_adj_for_follicles < 0.05
            ax.scatter(r.bald_effect_adj, i, s=26, facecolor=col if sig else "white", edgecolor=col, linewidth=1.2, zorder=3)
        ax.axvline(0, color=MUTED, lw=0.6)
        ax.set_yticks(range(len(g)))
        ax.set_yticklabels(g.gene, fontsize=7)
        ax.set_ylim(-0.8, len(g) - 0.2)
        ax.set_xlim(-0.15, 1.05)
        ax.set_xlabel("bald vs healthy scalp, log2 fold change")
        ax.set_title(f"{label}: {n['sig_adjusted']} of {n['genes']} stay significant,\nabout {n['median_effect_kept'] * 100:.0f}% of the effect kept",
                     fontsize=8)
        vgrid(ax)
    handles = [Line2D([], [], marker="o", ls="", color=GRAY, markersize=4.5, label="before adjusting for follicle content"),
               Line2D([], [], marker="o", ls="", color=INK, markersize=5.5, label="after adjusting, still p < 0.05"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=INK, markersize=5.5, label="after adjusting, not significant")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=3, frameon=False)
    fig.text(0.5, -0.05, "GSE90594: 14 bald vs 14 healthy scalp biopsies. Adjustment: an 18-gene follicle-content score as a covariate.",
             ha="center", fontsize=6.5, color=MUTED)
    fig.subplots_adjust(wspace=0.35)
    save(fig, "fig3_scarring_program")


def fig4_pirfenidone():
    s5 = json.load(open(RESULTS / "05_summary.json"))["drugs"]
    cx = pd.read_csv(RESULTS / "11_cxcl12_link.csv", index_col=0).loc["CXCL12"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 2.0), gridspec_kw={"width_ratios": [1.15, 1]})
    rows = [("Pirfenidone 0.11 mM\n(just above pill level)", s5["pirfenidone_GSE178518"]["driven_lowered_sig"], GRAY),
            ("Nintedanib 1 uM", s5["nintedanib_GSE178518"]["driven_lowered_sig"], GRAY),
            ("Saracatinib 0.3 uM", s5["saracatinib_GSE178518"]["driven_lowered_sig"], GRAY),
            ("Pirfenidone 1 mM", s5["pirfenidone_1mM_GSE226249"]["driven_lowered_sig"], ORANGE)]
    for i, (lab, v, col) in enumerate(rows):
        a.barh(i, v * 100, color=col, height=0.55)
        a.text(v * 100 + 1, i, f"{v * 100:.0f}%", va="center", fontsize=6.6)
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([r[0] for r in rows])
    a.set_xlim(0, 60)
    a.set_xlabel("% of 23 TGF-beta-driven bald genes lowered")
    a.set_title("a  Reverses the program at 1 mM")
    vgrid(a)
    items = [("Pirfenidone 0.11 mM", cx.pirf_0p1mM_lfc, cx.pirf_0p1mM_p, GRAY), ("Nintedanib", cx.nint_lfc, cx.nint_p, GRAY),
             ("Pirfenidone 1 mM", cx.pirf_1mM_lfc, cx.pirf_1mM_p, ORANGE)]
    for i, (lab, v, p, col) in enumerate(items):
        b.barh(i, v, color=col, height=0.55)
        b.text(1.03, i, f"{num(v)}  (p {p:.2g})", transform=b.get_yaxis_transform(), va="center", ha="left", fontsize=6.4)
    b.axvline(0, color=MUTED, lw=0.6)
    b.set_yticks(range(len(items)))
    b.set_yticklabels([x[0] for x in items])
    b.set_xlim(-2.5, 1.5)
    b.set_xlabel("CXCL12 in fibroblasts, log2 fold change")
    b.set_title("b  Cuts CXCL12 at 1 mM")
    vgrid(b)
    fig.subplots_adjust(wspace=0.95)
    save(fig, "fig4_pirfenidone_fibroblasts")


def fig5_localization():
    t = pd.read_csv(RESULTS / "08_cluster_scores.csv", index_col=0)
    show = {"TGF-beta-driven bald genes": "TGF-beta-driven bald genes", "Pirfenidone-lowered (1 mM)": "...lowered by pirfenidone",
            "Scarring genes up in bald scalp": "Scarring genes up in bald scalp", "Hair-loss genetics, all strong (>=0.5)": "Strong hair-loss genetics",
            "Bald-down top 200": "Genes lost in bald scalp"}
    imm = t.index[t.compartment == "Immune"]
    keep = [c for c in t.index if t.loc[c, "compartment"] != "Immune"]
    mat = pd.DataFrame({lab: t[f"{k} | score"] for k, lab in show.items()}).T
    qs = pd.DataFrame({lab: t[f"{k} | q"] for k, lab in show.items()}).T
    mat = pd.concat([mat[keep], mat[imm].mean(axis=1).rename("Immune (15, mean)")], axis=1)
    qs = pd.concat([qs[keep], pd.Series(np.nan, index=qs.index, name="Immune (15, mean)")], axis=1)
    groups = [t.loc[c, "compartment"] for c in keep] + ["Immune"]
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", [BLUE, "#F2F2EF", ORANGE])
    fig, ax = plt.subplots(figsize=(DOUBLE, 2.5))
    lim = float(np.nanmax(np.abs(mat.values)))
    im = ax.imshow(mat.values, cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if qs.iloc[i, j] < 0.05 and mat.iloc[i, j] > 0:
                ax.text(j, i, "*", ha="center", va="center", fontsize=8, color=INK)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=6)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=6.8)
    edges = [j for j in range(1, len(groups)) if groups[j] != groups[j - 1]]
    for j in edges:
        ax.axvline(j - 0.5, color="white", lw=2)
    short = {"Dermal fibroblasts": "Dermal\nfibroblasts", "Follicle dermis": "Follicle\ndermis", "Muscle and pericytes": "Muscle,\npericytes",
             "Follicle epithelium": "Follicle\nepithelium", "Melanocytes": "Melano-\ncytes", "Immune": "Imm."}
    for s0, s1 in zip([0] + edges, edges + [len(groups)]):
        ax.text((s0 + s1 - 1) / 2, -0.7, short.get(groups[s0], groups[s0]), ha="center", va="bottom", fontsize=6, color=MUTED)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("mean z across cell types", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Where each gene set is expressed in human scalp (Ober-Reynolds 2023 atlas; * enriched, q < 0.05)", pad=26, fontsize=8)
    save(fig, "fig5_scalp_localization")


def fig6_faers():
    d = pd.read_csv(RESULTS / "06_faers_disproportionality.csv")
    ac = pd.read_csv(RESULTS / "06_faers_pirfenidone_vs_nintedanib.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 2.3), gridspec_kw={"width_ratios": [1.2, 1]})
    g = d[d.group == "growth_control"].sort_values("ror_growth")
    l = d[d.group == "loss_control"].sort_values("ror_alopecia")
    rows = [(f"{r.drug}: excess hair", r.ror_growth, r.ror_growth_lo, r.ror_growth_hi, BLUE) for r in g.itertuples()]
    rows += [(f"{r.drug}: hair loss", r.ror_alopecia, r.ror_alopecia_lo, r.ror_alopecia_hi, GRAY) for r in l.itertuples()]
    for i, (lab, v, lo, hi, col) in enumerate(rows):
        a.hlines(i, lo, hi, color=col, lw=1.2)
        a.scatter(v, i, color=col, s=14, zorder=3)
    a.axvline(1, color=MUTED, lw=0.8, ls="--")
    a.set_xscale("log")
    plain_log(a)
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([r[0] for r in rows], fontsize=6.4)
    a.set_xlabel("reporting odds ratio (log scale)")
    a.set_title("a  The method finds known hair drugs")
    vgrid(a)
    labels = {"all": "All lung-fibrosis patients", "men": "Men", "women": "Women", "age_60_85": "Age 60 to 85"}
    ac = ac.iloc[::-1].reset_index(drop=True)
    for i, r in ac.iterrows():
        b.hlines(i, r.lo, r.hi, color=ORANGE, lw=1.4)
        b.scatter(r.ratio, i, color=ORANGE, s=18, zorder=3)
        b.text(1.03, i, f"{r.ratio:.2f} [{r.lo:.2f}, {r.hi:.2f}]", transform=b.get_yaxis_transform(), va="center", ha="left", fontsize=6.4)
    b.axvline(1, color=MUTED, lw=0.8, ls="--")
    b.set_xscale("log")
    b.set_xlim(0.1, 3)
    plain_log(b)
    b.set_yticks(range(len(ac)))
    b.set_yticklabels([labels[s] for s in ac.stratum])
    b.set_xlabel("hair-loss reporting ratio, pirfenidone vs nintedanib")
    b.set_title("b  Pirfenidone: about 70% fewer hair-loss reports")
    vgrid(b)
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig6_faers")


def fig7_genetics():
    res = pd.read_csv(RESULTS / "10_coloc_mr.csv")
    t = res[res.testable].copy()
    tis = {"fibroblasts (cultured, skin)": "fibroblasts", "skin, not sun-exposed": "skin", "skin, sun-exposed": "sun skin"}
    t["label"] = t.gene + " (" + t.tissue.map(tis) + ")"
    t = t.sort_values(["panel", "PP.H4"], ascending=[False, True]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(SINGLE, 4.3))
    colors = t.panel.map({"hypothesis": ORANGE, "reference": BLUE})
    ax.hlines(range(len(t)), 0, t["PP.H4"], color=LIGHT, lw=1)
    ax.scatter(t["PP.H4"], range(len(t)), c=colors, s=16, zorder=3, edgecolor="white", linewidth=0.8)
    ax.axvline(0.8, color=MUTED, lw=0.8, ls="--")
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels(t.label, fontsize=6)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("P(shared causal variant), coloc PP.H4")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=ORANGE, label="TGF-beta / scarring genes: 0 of 8 above 0.8"),
                       Line2D([], [], marker="o", ls="", color=BLUE, label="strongest baldness genes: 10 of 23")],
              loc="upper center", bbox_to_anchor=(0.4, -0.115), frameon=False, ncol=1)
    ax.set_title("Genetic natural experiment")
    fig.subplots_adjust(left=0.33, right=0.96)
    vgrid(ax)
    save(fig, "fig7_genetics")


def fig8_agent_audit():
    v2 = json.load(open(RESULTS / "09_v2" / "09_summary.json"))["audit"]["accuracy"]
    v3 = json.load(open(RESULTS / "09_summary.json"))["audit"]
    names = [("pipeline", "Full pipeline"), ("meta-llama/llama-3.3-70b-instruct", "Llama 3.3 70B alone"),
             ("deepseek/deepseek-chat-v3.1", "DeepSeek v3.1 alone"), ("openai/gpt-4o-mini (shadow, no vote)", "GPT-4o-mini alone")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 1.9), gridspec_kw={"width_ratios": [1.3, 1]})
    y = np.arange(len(names))
    a.barh(y + 0.18, [v2[k]["population_weighted"] * 100 for k, _ in names], height=0.34, color=LIGHT, label="v2, 60 unseen pairs")
    a.barh(y - 0.18, [v3["accuracy"][k]["population_weighted"] * 100 for k, _ in names], height=0.34, color=BLUE, label="v3, 40 new unseen pairs")
    for i, (k, _) in enumerate(names):
        a.text(v3["accuracy"][k]["population_weighted"] * 100 + 1, i - 0.18, f"{v3['accuracy'][k]['population_weighted'] * 100:.0f}%", va="center", fontsize=6.4)
    a.set_yticks(y)
    a.set_yticklabels([n for _, n in names])
    a.invert_yaxis()
    a.set_xlim(0, 105)
    a.set_xlabel("accuracy vs blind reference labels (population-weighted)")
    a.legend(loc="upper center", bbox_to_anchor=(0.45, -0.45), frameon=False, ncol=2)
    a.set_title("a  Blind audits of the LLM evidence agent")
    vgrid(a)
    p = v3["accuracy"]["pipeline"]
    vals = [p["agree_stratum"] * 100, p["disagree_stratum"] * 100]
    b.barh([0, 1], vals, color=[BLUE, GRAY], height=0.5)
    for i, v in enumerate(vals):
        b.text(v + 1, i, f"{v:.0f}%", va="center", fontsize=6.6)
    b.set_yticks([0, 1])
    b.set_yticklabels(["two reviewers agree", "reviewers disagree"])
    b.invert_yaxis()
    b.set_xlim(0, 105)
    b.set_xlabel("pipeline accuracy, v3")
    b.set_title("b  Agreement is the confidence signal")
    vgrid(b)
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig8_agent_audit")


def fig9_dose():
    s = json.load(open(RESULTS / "14_summary.json"))
    gel = [r for r in s["skin_gel_dermis_average_estimates"] if "plain" not in r["formulation"]]
    items = [("Pill, peak blood level", min(s["oral_cmax_mM"].values()), max(s["oral_cmax_mM"].values()), BLUE),
             ("Gel, dermis average", min(r["dermis_average_mM_low"] for r in gel), max(r["dermis_average_mM_high"] for r in gel), BLUE),
             ("Cells, no effect", s["cell_doses_mM"]["no effect (GSE178518, 20 ug/mL)"], None, GRAY),
             ("Cells, reversal", s["cell_doses_mM"]["reversal (GSE226249)"], None, ORANGE),
             ("Eye drop, site tissue", min(s["topical_eye_site_mM"].values()), max(s["topical_eye_site_mM"].values()), BLUE)]
    fig, ax = plt.subplots(figsize=(SINGLE, 1.95))
    for i, (lab, lo, hi, col) in enumerate(items):
        if hi is None:
            ax.scatter(lo, i, color=col, s=24, zorder=3, marker="D")
        else:
            ax.hlines(i, lo, hi, color=col, lw=3, alpha=0.85)
    for i, (lab, lo, hi, col) in enumerate(items):
        txt = f"{lo:.2g} mM" if hi is None else f"{lo:.2g} to {hi:.2g} mM"
        if (hi or lo) > 5:
            ax.text(lo / 1.35, i, txt, va="center", ha="right", fontsize=6.2, color=INK)
        else:
            ax.text((hi or lo) * 1.35, i, txt, va="center", fontsize=6.2, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(0.01, 100)
    plain_log(ax)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([x[0] for x in items])
    ax.invert_yaxis()
    ax.set_xlabel("pirfenidone, mM (log scale)")
    ax.set_title("The pill works on lungs below the cell dose", x=-0.62, ha="left")
    fig.subplots_adjust(left=0.36, right=0.97, top=0.86, bottom=0.22)
    vgrid(ax)
    save(fig, "fig9_dose")


for f in (fig1_workflow, fig2_calibration, fig3_scarring, fig4_pirfenidone, fig5_localization, fig6_faers, fig7_genetics, fig8_agent_audit, fig9_dose):
    f()
