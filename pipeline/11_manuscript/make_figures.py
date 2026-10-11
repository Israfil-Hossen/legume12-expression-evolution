#!/usr/bin/env python3
"""Main figures for the 12-legume manuscript, drawn to journal (Nature-style) specifications.

Every plotted value is read from a committed result file under results/ or metadata/;
the data behind each panel are written to manuscript_mpe/figures/source_data/.
Widths 89 mm (single column) or 183 mm (double); Arial; 6-7 pt text; 8 pt bold panel letters.
Usage (repository root): python pipeline/11_manuscript/make_figures.py
"""
import gzip, math, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

RES, FIG = "results", "manuscript_mpe/figures"
SRC = os.path.join(FIG, "source_data")
os.makedirs(SRC, exist_ok=True)
MM = 1 / 25.4
W1, W2 = 89, 183
MATERIAL, NEGLIGIBLE = 0.050, 0.025

plt.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.labelpad": 2, "legend.frameon": False, "legend.handlelength": 1.2,
})

SPECIES = ["Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur"]
LATIN = dict(Gmax="Glycine max", Pvul="Phaseolus vulgaris", Vrad="Vigna radiata", Vung="Vigna unguiculata",
             Ccaj="Cajanus cajan", Mtru="Medicago truncatula", Cari="Cicer arietinum", Psat="Pisum sativum",
             Tpra="Trifolium pratense", Ljap="Lotus japonicus", Lang="Lupinus angustifolius", Adur="Arachis duranensis")
SHORT = dict(Gmax="G. max", Pvul="P. vulgaris", Vrad="V. radiata", Vung="V. unguiculata", Ccaj="C. cajan",
             Mtru="M. truncatula", Cari="C. arietinum", Psat="P. sativum", Tpra="T. pratense", Ljap="L. japonicus",
             Lang="L. angustifolius", Adur="A. duranensis")
CLADE = dict(Gmax="Millettioid", Pvul="Millettioid", Vrad="Millettioid", Vung="Millettioid", Ccaj="Millettioid",
             Mtru="IRLC", Cari="IRLC", Psat="IRLC", Tpra="IRLC", Ljap="Robinioid", Lang="Genistoid", Adur="Dalbergioid")
CL_COL = {"Millettioid": "#2F5D8C", "IRLC": "#3A9278", "Robinioid": "#B0569A", "Genistoid": "#D9A021", "Dalbergioid": "#C8553D"}
FOCAL_COL = dict(Gmax="#2F5D8C", Pvul="#C8553D", Mtru="#3A9278")
ACC = "#2F5D8C"          # primary data colour
GREY = "#8C8C8C"
ZONE = dict(tech="#E3F0EC", part="#F2F2F2", evo="#FBE9DD")


def num(v, d=3):
    """Signed number with a typographic minus sign."""
    return f"{v:+.{d}f}".replace("-", chr(0x2212))


def italic(s):
    return "$\\it{" + s.replace(". ", ".\\,").replace(" ", "\\ ") + "}$"


def tag(ax, letter, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="right")


LAYOUT_ISSUES = []


def check_layout(fig, name, pad=0.5):
    """Report text that overlaps other text, leaves the canvas, or sits on plotted data (points or lines)."""
    from matplotlib.text import Text
    from matplotlib.collections import PathCollection
    from matplotlib.lines import Line2D as L2D
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    hidden = set()                                   # tick labels of ticks outside the view limits are never drawn
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not (lo - 1e-9 <= tk.get_loc() <= hi + 1e-9):
                    hidden.update([id(tk.label1), id(tk.label2)])
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip() and t.figure is fig and id(t) not in hidden]
    boxes = []
    for t in texts:
        try:
            b = t.get_window_extent(r)
        except Exception:
            continue
        if b.width > 0 and b.height > 0:
            boxes.append((t, b.expanded(1, 1).padded(-pad)))
    W, H = fig.bbox.width, fig.bbox.height
    for t, b in boxes:
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5:
            LAYOUT_ISSUES.append(f"{name}: text off canvas: '{t.get_text()[:40]}'")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (t1, b1), (t2, b2) = boxes[i], boxes[j]
            if b1.overlaps(b2):
                LAYOUT_ISSUES.append(f"{name}: text overlap: '{t1.get_text()[:30]}' x '{t2.get_text()[:30]}'")
    # text of one panel running into the plotting area of another panel
    for t, b in boxes:
        own = t.axes
        for ax in fig.axes:
            if ax is own or own is None:
                continue
            if b.overlaps(ax.get_window_extent(r).padded(-1)):
                LAYOUT_ISSUES.append(f"{name}: text enters another panel: '{t.get_text()[:40]}'")
    # data under text: scatter offsets and sampled points along plotted lines, inside each axes
    for ax in fig.axes:
        pts = []
        for col in ax.collections:
            if isinstance(col, PathCollection) and len(col.get_offsets()):
                pts.append(col.get_offset_transform().transform(col.get_offsets()))
        for ln in ax.lines:
            if not isinstance(ln, L2D) or not ln.get_visible():
                continue
            xy = ln.get_transform().transform(ln.get_xydata())
            if len(xy) >= 2:
                seg = np.concatenate([np.linspace(xy[k], xy[k + 1], 25) for k in range(len(xy) - 1)])
                pts.append(seg)
        # bars and histogram rectangles (not full-span shading such as decision zones, not schematic boxes)
        from matplotlib.patches import Rectangle as Rect
        axb = ax.get_window_extent(r)
        bars = []
        for pa in ax.patches:
            if isinstance(pa, Rect) and pa.get_visible():
                bb = pa.get_window_extent(r)
                if bb.width < 0.95 * axb.width and bb.height < 0.95 * axb.height and bb.width > 0 and bb.height > 0:
                    bars.append(bb)
        leg = ax.get_legend()
        legbox = leg.get_window_extent(r).padded(-1) if leg is not None and leg.get_visible() else None
        for t, b in boxes:
            if (t.axes is ax or t in ax.texts) and any(b.overlaps(bb) for bb in bars):
                LAYOUT_ISSUES.append(f"{name}: text on a bar: '{t.get_text()[:40]}'")
        if legbox is not None and any(legbox.overlaps(bb) for bb in bars):
            LAYOUT_ISSUES.append(f"{name}: legend covers a bar")
        if not pts:
            continue
        P = np.concatenate(pts)
        if legbox is not None:
            inl = (P[:, 0] > legbox.x0) & (P[:, 0] < legbox.x1) & (P[:, 1] > legbox.y0) & (P[:, 1] < legbox.y1)
            if inl.any():
                LAYOUT_ISSUES.append(f"{name}: legend covers data")
        for t, b in boxes:
            if t.axes is not ax and t not in ax.texts:
                continue
            inside = (P[:, 0] > b.x0) & (P[:, 0] < b.x1) & (P[:, 1] > b.y0) & (P[:, 1] < b.y1)
            if inside.any():
                LAYOUT_ISSUES.append(f"{name}: text on data: '{t.get_text()[:40]}'")


def save(fig, name):
    check_layout(fig, name)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600}), ("svg", {})):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), **kw)
    plt.close(fig)


def decision_zones(ax, horizontal=True):
    """Shade the pre-registered decision zones on the estimate axis."""
    lo, hi = (ax.get_xlim() if horizontal else ax.get_ylim())
    span = ax.axvspan if horizontal else ax.axhspan
    span(lo, NEGLIGIBLE, color=ZONE["tech"], zorder=0, lw=0)
    span(NEGLIGIBLE, MATERIAL, color=ZONE["part"], zorder=0, lw=0)
    span(MATERIAL, hi, color=ZONE["evo"], zorder=0, lw=0)
    (ax.axvline if horizontal else ax.axhline)(0, color="#555555", lw=0.5, zorder=1)
    (ax.set_xlim if horizontal else ax.set_ylim)(lo, hi)


# ------------------------------------------------------------------ inputs
qc = pd.read_csv(f"{RES}/results_step2/qc_samples.tsv", sep="\t")
qc = qc[qc.qc_pass]
q4 = pd.read_csv(f"{RES}/results_step4/qc_samples_step4.tsv", sep="\t")
q4 = q4[q4.qc_pass]
W = pd.read_csv(f"{RES}/results_step2/rho_within_primary.tsv", sep="\t")
B = pd.read_csv(f"{RES}/results_step2/rho_between_primary.tsv", sep="\t")
TT = pd.read_csv("metadata/step2_divergence_times.tsv", sep="\t")
T2 = pd.read_csv(f"{RES}/results_step2/tests.tsv", sep="\t")
T3C = pd.read_csv(f"{RES}/results_step3C/tests.tsv", sep="\t")
P4 = pd.read_csv(f"{RES}/results_step4/pairs_step4.tsv", sep="\t")
T4 = pd.read_csv(f"{RES}/results_step4/tests_step4.tsv", sep="\t")
A3 = pd.read_csv(f"{RES}/results_step3/A_per_species.tsv", sep="\t")
T3 = pd.read_csv(f"{RES}/results_step3/tests_step3AB.tsv", sep="\t")
r3 = open(f"{RES}/results_step3/report_step3AB.txt").read()


# ------------------------------------------------------------------ tree (Newick, no external library)
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = {"name": "", "len": 0.0, "children": []}
        if s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^:,()]*", s[pos:]); n["name"] = m.group(0); pos += len(m.group(0))
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":([0-9.eE+-]+)", s[pos:]); n["len"] = float(m.group(1)); pos += len(m.group(0))
        return n
    return node()


def tips(n):
    return [n["name"]] if not n["children"] else [t for c in n["children"] for t in tips(c)]


tree = parse_newick(open(f"{RES}/server_logs/concat.fa.treefile").read())
# IQ-TREE writes the unrooted tree as a basal trichotomy; rooting on A. thaliana and pruning it
# leaves the legume root at the basal node.
tree["children"] = [c for c in tree["children"] if c["name"] != "Atha"]
assert sorted(tips(tree)) == sorted(SPECIES)


def ladder(n):
    for c in n["children"]:
        ladder(c)
    n["children"].sort(key=lambda c: len(tips(c)))


ladder(tree)


def layout(n, x0, coords, order):
    x = x0 + n["len"]
    if not n["children"]:
        y = len(order); order.append(n["name"]); coords[id(n)] = (x, y); return y
    ys = [layout(c, x, coords, order) for c in n["children"]]
    y = (min(ys) + max(ys)) / 2; coords[id(n)] = (x, y); return y


# ======================================================================= Fig. 1
def fig1():
    fig = plt.figure(figsize=(W2 * MM, 74 * MM))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.30, 0.50, 0.95, 0.62], wspace=0.10,
                          left=0.015, right=0.985, top=0.90, bottom=0.17)
    ax = fig.add_subplot(gs[0]); axn = fig.add_subplot(gs[1], sharey=ax)
    axc = fig.add_subplot(gs[2]); axp = fig.add_subplot(gs[3])
    coords, order = {}, []
    layout(tree, 0.0, coords, order)
    tipx = {}
    SUP = pd.read_csv(f"{RES}/tree_support/support_summary.tsv", sep="\t")
    CF = {frozenset(t.split(",")): (g, s_) for t, g, s_ in zip(SUP.tips, SUP.gCF, SUP.sCF)}
    SUP.to_csv(f"{SRC}/Fig1a_support.csv", index=False)

    segs, cflab = [], []

    def draw(n):
        x, y = coords[id(n)]
        if not n["children"]:
            tipx[n["name"]] = x
            return
        ys = [coords[id(c)][1] for c in n["children"]]
        ax.plot([x, x], [min(ys), max(ys)], color="#333333", lw=0.7, solid_capstyle="butt"); segs.append((x, min(ys), x, max(ys)))
        for c in n["children"]:
            cx, cy = coords[id(c)]
            ax.plot([x, cx], [cy, cy], color="#333333", lw=0.7, solid_capstyle="butt"); segs.append((x, cy, cx, cy))
            if c["children"] and frozenset(tips(c)) in CF:          # gene / site concordance factor of the clade defined by c
                g, s_ = CF[frozenset(tips(c))]
                cflab.append((x, cx, cy, f"{g:.0f}/{s_:.0f}"))
            draw(c)
    draw(tree)
    xmax = max(tipx.values())
    for i, sp in enumerate(order):
        ax.plot([tipx[sp], xmax + 0.004], [i, i], color="#BBBBBB", lw=0.4, ls=(0, (1, 1.5)))   # leader to aligned label
        ax.text(xmax + 0.008, i, italic(LATIN[sp]), va="center", ha="left", fontsize=6.5, color=CL_COL[CLADE[sp]])
    ybar = len(order) + 0.15
    ax.plot([0.0, 0.05], [ybar, ybar], color="#333333", lw=0.7, solid_capstyle="butt")
    ax.text(0.056, ybar, "0.05 substitutions/site", ha="left", va="center", fontsize=6)
    ax.set_xlim(-0.005, xmax + 0.2); ax.set_ylim(len(order) + 0.75, -0.8); ax.axis("off")
    # concordance labels: first candidate position that touches no branch and no other label
    rend = fig.canvas.get_renderer(); inv = ax.transData.inverted(); placed = []
    def clear(bb):
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        x0, x1, y0, y1 = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
        for a0, b0, a1, b1 in segs:
            if b0 == b1 and y0 <= b0 <= y1 and max(min(a0, a1), x0) < min(max(a0, a1), x1):
                return False
            if a0 == a1 and x0 <= a0 <= x1 and max(min(b0, b1), y0) < min(max(b0, b1), y1):
                return False
        return all(x1 < q0 or x0 > q1 or y1 < r0 or y0 > r1 for q0, q1, r0, r1 in placed) and (placed.append((x0, x1, y0, y1)) or True)
    for px, cx, cy, lab in cflab:
        for tx, ty, ha, va in [(cx + 0.003, cy, "left", "center"), (cx - 0.0015, cy - 0.08, "right", "bottom"),
                               (cx + 0.003, cy - 0.36, "left", "center"), (cx + 0.003, cy + 0.36, "left", "center"),
                               (cx - 0.0015, cy + 0.08, "right", "top"), ((px + cx) / 2, cy - 0.08, "center", "bottom")]:
            t = ax.text(tx, ty, lab, ha=ha, va=va, fontsize=4.8, color="#555555")
            if clear(t.get_window_extent(rend).expanded(1.05, 1.15)):
                break
            t.remove()
        else:
            raise RuntimeError(f"no free position for concordance label {lab}")
    handles = [Line2D([], [], marker="s", ls="", color=c, markersize=4, label=k) for k, c in CL_COL.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.01, -0.02), ncol=5, fontsize=6,
              columnspacing=0.7, handletextpad=0.15)
    tag(ax, "a", x=0.02)
    # biological replicates per species: study 1 (Step 2) and study 2 (Step 4)
    n1 = qc.groupby("species").size(); n2 = q4.groupby("species").size()
    rows = []
    for i, sp in enumerate(order):
        a, b = int(n1.get(sp, 0)), int(n2.get(sp, 0))
        axn.barh(i, a, height=0.62, color=CL_COL[CLADE[sp]], lw=0)
        if b:
            axn.barh(i, b, left=a, height=0.62, color="white", edgecolor=CL_COL[CLADE[sp]], lw=0.6, hatch="//////")
        axn.text(a + b + 0.6, i, str(a) if not b else f"{a} + {b}", va="center", fontsize=6)
        rows.append(dict(code=sp, species=LATIN[sp], clade=CLADE[sp], samples_step2=a, samples_step4_study2=b))
    pd.DataFrame(rows).to_csv(f"{SRC}/Fig1a_species.csv", index=False)
    axn.set_xlim(0, 33); axn.set_xticks([0, 10, 20]); axn.tick_params(axis="y", left=False, labelleft=False)
    axn.spines["left"].set_visible(False); axn.set_xlabel("Biological replicates")
    axn.legend(handles=[Rectangle((0, 0), 1, 1, color="#777777", lw=0, label="study 1"),
                        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="#777777", hatch="//////", lw=0.6,
                                  label="study 2 (Step 4)")],
               loc="upper right", bbox_to_anchor=(1.06, 1.0), ncol=1, fontsize=6, handlelength=1.0,
               handletextpad=0.3, borderaxespad=0.1)
    # panel b: the three levels of comparison
    axc.set_xlim(0, 10); axc.set_ylim(0, 10); axc.axis("off"); tag(axc, "b", x=0.03)
    rung = [("Within one study", "same species, same study", "#7A7A7A"),
            ("Between two studies", "same species, two studies", "#D9A021"),
            ("Between two species", "9.7–56 Myr apart", ACC)]
    for k, (t1, t2, col) in enumerate(rung):
        y = 8.7 - k * 2.6
        axc.add_patch(FancyBboxPatch((0.3, y - 0.95), 9.4, 1.9, boxstyle="round,pad=0.02,rounding_size=0.3",
                                     fc=matplotlib.colors.to_rgba(col, 0.10), ec=col, lw=0.7))
        axc.text(0.75, y + 0.35, t1, fontsize=6.8, fontweight="bold", va="center")
        axc.text(0.75, y - 0.4, t2, fontsize=6.2, va="center", color="#444444")
    axc.text(0.3, 0.95, "Statistic at each level: Spearman ρ between\ncopy number and expression divergence,\n2 vs 2 samples",
             fontsize=6.0, va="center", color="#333333")
    # panel c: schematic predictions (levels are arbitrary; only the shape matters)
    xs = [0, 1, 2]
    axp.plot(xs, [0.42, 0.42, 1.0], color="#C8553D", lw=1.1, marker="o", ms=2.8, label="evolution")
    axp.plot(xs, [0.38, 0.38, 0.38], color="#3A9278", lw=1.1, marker="o", ms=2.8, ls=(0, (3, 1.5)), label="non-evolutionary")
    axp.set_xticks(xs, ["within\nstudy", "between\nstudies", "between\nspecies"], fontsize=6)
    axp.set_yticks([]); axp.set_ylim(0, 1.3); axp.set_xlim(-0.3, 2.3)
    axp.set_ylabel("Expected ρ (schematic)", fontsize=6.5)
    axp.legend(loc="upper left", fontsize=6, handlelength=1.6, borderaxespad=0.1)
    tag(axp, "c", x=-0.05)
    save(fig, "Fig1")


# ======================================================================= Step 5 inputs
S5 = f"{RES}/results_step5"
T5 = pd.read_csv(f"{S5}/tests_step5.tsv", sep="\t").set_index("analysis")
P1A = pd.read_csv(f"{S5}/1A_pairs.tsv", sep="\t")
P1B = pd.read_csv(f"{S5}/1B_partial.tsv", sep="\t")
P1C = pd.read_csv(f"{S5}/1C_ambiguity.tsv", sep="\t")
P2A = pd.read_csv(f"{S5}/2A_summary.tsv", sep="\t")
P2B = pd.read_csv(f"{S5}/2B_species_pairs.tsv", sep="\t")
CE = pd.read_csv(f"{S5}/2B_ceilings.tsv", sep="\t")
MRAT = pd.read_csv(f"{S5}/2B_mean_ratios.tsv", sep="\t").set_index("OG")
DIRN = pd.read_csv(f"{S5}/2B_direction.tsv", sep="\t")
D3 = {k: pd.read_csv(f"{S5}/3_dosage_Gmax_{k}.tsv.gz", sep="\t") for k in ("study1", "study2")}
FOCAL = ["Gmax", "Pvul", "Mtru"]
TREE_ORDER = ["Adur", "Lang", "Ljap", "Cari", "Psat", "Mtru", "Tpra", "Ccaj", "Gmax", "Pvul", "Vrad", "Vung"]


def forest(ax, rows, labels, xlim, highlight_first=True):
    ax.set_xlim(*xlim); decision_zones(ax)
    for i, (est, lo, hi) in enumerate(rows):
        y = len(rows) - 1 - i
        ax.plot([lo, hi], [y, y], color="#222222", lw=0.8, solid_capstyle="butt")
        ax.plot(est, y, "o", color=ACC if (i == 0 and highlight_first) else "#222222", ms=4 if i == 0 else 3.2, mec="white", mew=0.4)
    ax.set_yticks(range(len(rows)), labels[::-1]); ax.set_ylim(-0.7, len(rows) - 0.3); ax.tick_params(axis="y", length=0)


def zone_labels(ax, ytop, xs=(-0.0275, 0.0375, 0.065)):
    for x, t in zip(xs, ("technical", "partial", "evolution")):
        ax.text(x, ytop, t, ha="center", va="bottom", fontsize=6, color="#555555")


# ======================================================================= Fig. 2
def fig2():
    """The copy-number gradient is present without evolution and does not grow with divergence time."""
    fig = plt.figure(figsize=(W2 * MM, 66 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.0, 1.05], wspace=0.62, left=0.13, right=0.99, top=0.88, bottom=0.17)
    a, b, c = (fig.add_subplot(gs[i]) for i in range(3))
    w = W.sort_values("rho").reset_index(drop=True); w.to_csv(f"{SRC}/Fig2a_within.csv", index=False)
    for i, r in w.iterrows():
        col = CL_COL[CLADE[r.species]]
        a.plot([0, r.rho], [i, i], color=col, lw=1.0, solid_capstyle="butt"); a.plot(r.rho, i, "o", color=col, ms=3.4, mec="white", mew=0.4)
    a.set_yticks(range(len(w)), [SHORT[s] for s in w.species], fontstyle="italic"); a.tick_params(axis="y", length=0)
    a.set_xlim(0, 0.21); a.set_ylim(-0.7, len(w) - 0.3); a.set_xlabel("ρ within one study\n(replicate halves)")
    a.text(0.205, -0.3, f"{(w.rho > 0).sum()} of {len(w)} > 0", ha="right", fontsize=6.2)
    tag(a, "a", x=-0.5)
    tt = pd.concat([TT[["sp_a", "sp_b", "median_mya"]], TT.rename(columns={"sp_a": "sp_b", "sp_b": "sp_a"})[["sp_a", "sp_b", "median_mya"]]])
    d = B.merge(tt, on=["sp_a", "sp_b"]); d["cajanus"] = (d.sp_a == "Ccaj") | (d.sp_b == "Ccaj")
    d.to_csv(f"{SRC}/Fig2b_rho_between_time.csv", index=False)
    p = d[~d.cajanus]; slope, icpt = np.polyfit(p.median_mya, p.rho, 1)
    b.scatter(p.median_mya, p.rho, s=9, color=ACC, alpha=0.8, lw=0, label=f"primary test ({len(p)} pairs)")
    b.scatter(d[d.cajanus].median_mya, d[d.cajanus].rho, s=9, facecolor="none", edgecolor=GREY, lw=0.6, label=f"{italic('Cajanus')} pairs (D1)")
    xx = np.array([p.median_mya.min(), p.median_mya.max()]); b.plot(xx, icpt + slope * xx, color="#222222", lw=0.9)
    pr = T2[(T2.variant == "primary") & T2.test.str.startswith("primary")].iloc[0]
    b.text(0.03, 0.04, f"change over {xx[1] - xx[0]:.1f} Myr: {num(pr.estimate)}\n95% CI {num(pr.lo)} to {num(pr.hi)}", transform=b.transAxes, fontsize=6.2)
    b.set_xlabel("Divergence time (Myr)"); b.set_ylabel("ρ between species"); b.set_ylim(0.02, 0.15)
    b.legend(loc="upper right", fontsize=6.2, handletextpad=0.1)
    tag(b, "b", x=-0.2)
    lab = {"primary": "Primary", "S1_per_copy": "S1 per-copy", "S2_three_per_side": "S2 3 per side", "S3_no_Pisum": f"S3 no {italic('Pisum')}",
           "S4_thinned_5M": "S4 thinned 5 M", "S5_patristic": "S5 patristic", "C": "C HOGs"}
    f = pd.concat([T2[~T2.test.str.startswith("secondary")][["variant", "estimate", "lo", "hi"]],
                   T3C[T3C.test.str.startswith("primary")][["estimate", "lo", "hi"]].assign(variant="C")]).set_index("variant").loc[list(lab)].reset_index()
    f.to_csv(f"{SRC}/Fig2c_tests.csv", index=False)
    forest(c, list(zip(f.estimate, f.lo, f.hi)), [lab[v] for v in f.variant], (-0.08, 0.1))
    zone_labels(c, len(f) - 0.25, xs=(-0.0275, 0.0375, 0.075))
    c.set_xticks([-0.05, 0, 0.05]); c.set_xlabel("Change in ρ across\nthe divergence range")
    tag(c, "c", x=-0.42)
    save(fig, "Fig2")


# ======================================================================= Fig. 3
def fig3():
    """Study versus species on identical orthogroups: copy-number gradient (a, b) and expression distance (c)."""
    fig = plt.figure(figsize=(W2 * MM, 112 * MM))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.55], height_ratios=[1, 0.8], wspace=0.42, hspace=0.48,
                          left=0.08, right=0.98, top=0.94, bottom=0.1)
    P4.to_csv(f"{SRC}/Fig3a_pairs.csv", index=False); T4.to_csv(f"{SRC}/Fig3b_delta.csv", index=False)
    for k, sp in enumerate(FOCAL):
        ax = fig.add_subplot(gs[0, k]); d = P4[P4.focal == sp]
        for _, r in d.iterrows():
            ax.plot([0, 1], [r.rho_study, r.rho_species], color=FOCAL_COL[sp], alpha=0.45, lw=0.6)
        ax.plot([0, 1], [d.rho_study.mean(), d.rho_species.mean()], color=FOCAL_COL[sp], lw=1.8, marker="o", ms=3.5)
        ax.set_xlim(-0.25, 1.25); ax.set_ylim(0.03, 0.19); ax.set_xticks([0, 1], ["two\nstudies", "two\nspecies"])
        ax.set_title(italic(LATIN[sp]), fontsize=6.6, pad=3)
        if k == 0:
            ax.set_ylabel("Copy-number gradient ρ\n(identical orthogroups)"); tag(ax, "a", x=-0.38)
        else:
            ax.tick_params(axis="y", labelleft=False)
    b = fig.add_subplot(gs[0, 3]); order = FOCAL + ["all"]
    t = T4.set_index("level").loc[order].reset_index(); b.set_xlim(-0.08, 0.08); decision_zones(b)
    for i, r in t.iterrows():
        y = len(t) - 1 - i; col = "#222222" if r.level == "all" else FOCAL_COL[r.level]
        b.plot([r.lo, r.hi], [y, y], color=col, lw=0.9, solid_capstyle="butt")
        b.plot(r.delta, y, "D" if r.level == "all" else "o", color=col, ms=3.6, mec="white", mew=0.4)
    b.set_yticks(range(len(t)), [f"all {len(P4)} pairs" if v == "all" else SHORT[v] for v in order][::-1])
    for lbl, v in zip(b.get_yticklabels(), order[::-1]):
        lbl.set_fontstyle("normal" if v == "all" else "italic")
    b.set_ylim(-0.7, len(t) - 0.3); b.tick_params(axis="y", length=0)
    b.set_xlabel("Δ = ρ$_{species}$ − ρ$_{study}$")
    allr = t[t.level == "all"].iloc[0]
    b.text(0.98, 0.97, f"all pairs\nΔ = {num(allr.delta)}\n({num(allr.lo)} to\n{num(allr.hi)})", transform=b.transAxes, ha="right", va="top", fontsize=6.0)
    tag(b, "b", x=-0.25)
    # c: expression distance (median over orthogroups) within study, between studies, between species (1A);
    # every bar is the median over the 11 orthogroup sets (one per other species) with the range across sets
    c = fig.add_subplot(gs[1, :3]); P1A.to_csv(f"{SRC}/Fig3c_distance.csv", index=False)
    cols = (("d_within", "#BDBDBD", "within one study"), ("d_study", "#D9A021", "between two studies (same species)"),
            ("d_species", ACC, "between two species"))
    for k, sp in enumerate(FOCAL):
        d = P1A[P1A.focal == sp]
        for j_, (col, colour, _) in enumerate(cols):
            x = k * 4 + j_
            c.bar(x, d[col].median(), width=0.8, color=colour, lw=0)
            c.plot([x, x], [d[col].min(), d[col].max()], color="#222222", lw=0.7)
    c.set_xticks([k * 4 + 1 for k in range(3)], [f"{SHORT[s_]}\nF = {P1A[P1A.focal == s_].F.mean():.2f}" for s_ in FOCAL])
    for lbl in c.get_xticklabels():
        lbl.set_fontstyle("italic")
    c.set_ylabel("Expression distance\n(median over orthogroups, log2)"); c.set_ylim(0, 0.85); c.set_xlim(-0.8, 11)
    c.legend(handles=[Rectangle((0, 0), 1, 1, color=colour, lw=0, label=l) for _, colour, l in cols],
             loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=6.2)
    fa = T5.loc["1A_F"]
    c.text(1.02, 0.35, f"study share F,\nmean of 3 species = {fa.estimate:.2f}\n(orthogroup 95% CI\n{fa.lo:.2f} to {fa.hi:.2f})",
           transform=c.transAxes, ha="left", va="top", fontsize=6.2)
    tag(c, "c", x=-0.12, y=1.06)
    save(fig, "Fig3")


# ======================================================================= Fig. 4
def fig4():
    """The copy-number gradient is not explained by expression level or read ambiguity."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2 * MM, 70 * MM), gridspec_kw=dict(width_ratios=[1.5, 1], wspace=0.55))
    fig.subplots_adjust(top=0.85, bottom=0.14, left=0.11, right=0.98)
    w = P1B[P1B.level == "within"].merge(P1C, left_on="label", right_on="species")
    w = w.set_index("label").loc[[s for s in TREE_ORDER if s in set(w.label)]].reset_index()
    w.to_csv(f"{SRC}/Fig4_partial.csv", index=False)
    y = np.arange(len(w))[::-1]
    for yy, r in zip(y, w.itertuples()):
        a.plot([r.rho_raw, r.rho_partial_expr, r.rho_partial_amb_expr], [yy] * 3, color="#CCCCCC", lw=0.8, zorder=1)
    a.scatter(w.rho_raw, y, s=14, color="#BDBDBD", label="raw", zorder=2, lw=0)
    a.scatter(w.rho_partial_expr, y, s=14, color="#D9A021", label="| expression level", zorder=3, lw=0)
    a.scatter(w.rho_partial_amb_expr, y, s=14, color=ACC, label="| expression level + read ambiguity", zorder=4, lw=0)
    a.axvline(0, color="#555555", lw=0.5); a.axvspan(-0.05, NEGLIGIBLE, color=ZONE["tech"], zorder=0, lw=0)
    a.set_yticks(y, [SHORT[s] for s in w.label], fontstyle="italic"); a.tick_params(axis="y", length=0)
    a.set_xlim(-0.05, 0.23); a.set_xlabel("Copy-number gradient within one study (Spearman ρ, partial)")
    a.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=6.2, handletextpad=0.2, columnspacing=0.8)
    tag(a, "a", x=-0.2, y=1.1)
    b.barh(y, w.rho_copies_ambiguity, color="#9DB4CC", height=0.6, lw=0)
    b.set_yticks(y, [SHORT[s] for s in w.label], fontstyle="italic"); b.tick_params(axis="y", length=0)
    b.axvline(0, color="#555555", lw=0.5); b.set_xlabel("ρ (copy number, read ambiguity)")
    tag(b, "b", x=-0.3, y=1.1)
    save(fig, "Fig4")


# ======================================================================= Fig. 5
def fig5():
    """Ancient papilionoid WGD: ohnolog expression asymmetry is conserved, shared by close relatives,
    and erodes with divergence time; within-sample ratios are not protected from laboratory effects."""
    fig = plt.figure(figsize=(W2 * MM, 150 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.08], height_ratios=[1.1, 1], wspace=0.45, hspace=0.5,
                          left=0.07, right=0.93, top=0.95, bottom=0.08)
    a, b, c, d = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    # a: per-species mean ratio of each ohnolog pair
    m = MRAT[TREE_ORDER]; m = m[m.notna().sum(axis=1) >= 10]
    m = m.loc[m.mean(axis=1).sort_values().index]
    m.to_csv(f"{SRC}/Fig5a_mean_ratios.csv")
    cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad("#E6E6E6")
    im = a.imshow(np.clip(m.values, -4, 4), aspect="auto", cmap=cmap, vmin=-4, vmax=4, interpolation="nearest")
    a.set_xticks(range(len(TREE_ORDER)), [SHORT[s] for s in TREE_ORDER], rotation=90, fontstyle="italic", fontsize=5.8)
    a.set_yticks([]); a.set_ylabel(f"{len(m)} ohnolog pairs (ordered by mean)")
    for s in a.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=a, fraction=0.05, pad=0.03, ticks=[-4, -2, 0, 2, 4])
    cb.set_label("log2 (copy A / copy B)", fontsize=6); cb.ax.tick_params(labelsize=5.8, width=0.4, length=2); cb.outline.set_linewidth(0.4)
    tag(a, "a", x=-0.06)
    # b: conservation between every pair of species, ordered by the species tree
    M = pd.DataFrame(np.nan, index=TREE_ORDER, columns=TREE_ORDER)
    for r in P2B.itertuples():
        M.loc[r.sp_a, r.sp_b] = M.loc[r.sp_b, r.sp_a] = r.rho_r
    M.to_csv(f"{SRC}/Fig5b_conservation_matrix.csv")
    cm2 = matplotlib.colormaps["viridis"].copy(); cm2.set_bad("#FFFFFF")
    im2 = b.imshow(M.values, cmap=cm2, vmin=0.3, vmax=0.85, interpolation="nearest")
    b.set_xticks(range(12), [SHORT[s] for s in TREE_ORDER], rotation=90, fontstyle="italic", fontsize=5.8)
    b.set_yticks(range(12), [SHORT[s] for s in TREE_ORDER], fontstyle="italic", fontsize=5.8)
    b.tick_params(length=0)
    for s in b.spines.values():
        s.set_visible(False)
    # clade boundaries in tree order: Adur | Lang | Ljap | IRLC (Cari..Tpra) | Millettioids (Ccaj..Vung)
    for k in (0.5, 1.5, 2.5, 6.5):
        b.axhline(k, color="white", lw=1.0); b.axvline(k, color="white", lw=1.0)
    cb2 = fig.colorbar(im2, ax=b, fraction=0.046, pad=0.03)
    cb2.set_label("Spearman ρ", fontsize=6); cb2.ax.tick_params(labelsize=5.8, width=0.4, length=2)
    cb2.outline.set_linewidth(0.4)
    tag(b, "b", x=-0.25)
    # c: decline with divergence time, against the replicate and laboratory references
    P2B.to_csv(f"{SRC}/Fig5c_conservation_time.csv", index=False); CE.to_csv(f"{SRC}/Fig5c_ceilings.csv", index=False)
    rep = CE[CE.ceiling == "replicate"].rho.median(); labc = CE[CE.ceiling == "laboratory"].rho
    c.axhspan(labc.min(), labc.max(), color="#F3E5C2", lw=0, zorder=0)
    c.axhline(rep, color="#777777", lw=0.6, ls=(0, (3, 2)))
    c.text(9.5, rep - 0.006, "replicate halves (median)", fontsize=5.8, color="#555555", va="top")
    c.text(9.5, labc.max() - 0.006, "two studies, same species", fontsize=5.8, color="#8A6D1F", va="top")
    caj = (P2B.sp_a == "Ccaj") | (P2B.sp_b == "Ccaj"); p = P2B[~caj]
    c.scatter(p.time, p.rho_r, s=10, color=ACC, lw=0, alpha=0.85, label=f"species pairs ({len(p)})")
    c.scatter(P2B[caj].time, P2B[caj].rho_r, s=10, facecolor="none", edgecolor=GREY, lw=0.6, label=f"{italic('Cajanus')} pairs (D1)")
    slope, icpt = np.polyfit(p.time, p.rho_r, 1); xx = np.array([p.time.min(), p.time.max()])
    c.plot(xx, icpt + slope * xx, color="#222222", lw=0.9)
    tr = T5.loc["2B_trend"]
    c.text(0.03, 0.04, f"change over {xx[1] - xx[0]:.1f} Myr: {num(tr.estimate)}\n95% CI {num(tr.lo)} to {num(tr.hi)}", transform=c.transAxes, fontsize=6.2)
    c.set_xlabel("Divergence time (Myr)"); c.set_ylabel("Conservation of ohnolog asymmetry\n(ρ between species)")
    c.set_ylim(0.2, 1.02); c.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=6.2, handletextpad=0.1, columnspacing=0.8)
    tag(c, "c", x=-0.16, y=1.12)
    # d: laboratory excess for ratios versus absolute expression
    P2A.to_csv(f"{SRC}/Fig5d_lab_excess.csv", index=False)
    x = np.arange(len(P2A))
    d.bar(x - 0.19, P2A.X_abs, width=0.36, color="#BDBDBD", lw=0, label="absolute expression")
    d.bar(x + 0.19, P2A.X_ratio, width=0.36, color=ACC, lw=0, label="within-sample ratio")
    d.set_xticks(x, [SHORT[s_] for s_ in P2A.species], fontstyle="italic", fontsize=6.2)
    d.set_ylabel("Study excess\n(between − within studies)"); d.set_ylim(0, 0.8)
    e = T5.loc["2A_E"]
    d.text(0.5, 0.97, f"E = {e.estimate:.2f} (95% CI {e.lo:.2f} to {e.hi:.2f})", transform=d.transAxes, ha="center", va="top", fontsize=6.2)
    d.legend(loc="upper left", bbox_to_anchor=(0, 0.88), fontsize=6.2)
    tag(d, "d", x=-0.2)
    save(fig, "Fig5")


# ======================================================================= Fig. 6
def fig6():
    """Recent Glycine WGD: relative family expression rises less than copy number (partial dosage compensation)."""
    fig, axes = plt.subplots(1, 2, figsize=(W2 * MM, 66 * MM), sharey=True, gridspec_kw=dict(wspace=0.12))
    fig.subplots_adjust(bottom=0.22, top=0.8, left=0.08, right=0.98)
    for ax, (k, lab), letter in zip(axes, (("study1", "study 1"), ("study2", "study 2 (replication)")), "ab"):
        d = D3[k]; d.to_csv(f"{SRC}/Fig6{letter}_dosage_{k}.csv", index=False)
        wgd, cal = d[d.set == "glycine_wgd"].D, d[d.set == "calibration"].D
        bins = np.arange(-3, 4.01, 0.2)
        ax.hist(cal, bins=bins, density=True, color="#BDBDBD", alpha=0.9, lw=0, label=f"single-copy in all 12 ({len(cal)})")
        ax.hist(wgd, bins=bins, density=True, histtype="step", color=ACC, lw=1.0, label=f"Glycine WGD pairs ({len(wgd)})")
        c0 = np.median(cal)
        for xv, col, ls in ((c0, "#555555", "-"), (c0 + 1, "#C8553D", (0, (3, 2))), (np.median(wgd), ACC, "-")):
            ax.axvline(xv, color=col, lw=0.8, ls=ls)
        key = "3_dosage_Gmax" if k == "study1" else "3_dosage_Gmax_study2"
        r = T5.loc[key]
        ax.text(0.98, 0.97, f"{lab}\nrelative dosage {num(r.estimate)}\n95% CI {num(r.lo)} to {num(r.hi)}", transform=ax.transAxes,
                ha="right", va="top", fontsize=6.2)
        ax.set_xlabel("Relative family expression of $\\it{G.\\,max}$ (log2,\nminus mean of single-copy species)")
        ax.set_xlim(-3, 4)
        tag(ax, letter, x=-0.03 if letter == "b" else -0.1, y=1.03)
    axes[0].set_ylabel("Density")
    h, l = axes[0].get_legend_handles_labels()
    h += [Line2D([], [], color="#555555", lw=0.8), Line2D([], [], color=ACC, lw=0.8), Line2D([], [], color="#C8553D", lw=0.8, ls=(0, (3, 2)))]
    l += ["median, single-copy families", "median, Glycine WGD families", "expected with full dosage (+1)"]
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.53, 1.0), ncol=3, fontsize=6, handlelength=1.4, columnspacing=1.2)
    save(fig, "Fig6")


for f in (fig1, fig2, fig3, fig4, fig5, fig6):
    f()
print("figures written to", FIG)


# ======================================================================= Supplementary figures
A3P = pd.read_csv(f"{RES}/results_step3/A_per_species.tsv", sep="\t")
T3AB = pd.read_csv(f"{RES}/results_step3/tests_step3AB.tsv", sep="\t").set_index("analysis")
T1A = pd.read_csv(f"{S5}/1A_time.tsv", sep="\t")
ROB = pd.read_csv(f"{RES}/exploratory_step5/2B_robustness.tsv", sep="\t")
QC1all = pd.read_csv(f"{RES}/results_step2/qc_samples.tsv", sep="\t")
QC2all = pd.read_csv(f"{RES}/results_step4/qc_samples_step4.tsv", sep="\t")


def figS1():
    fig, ax = plt.subplots(figsize=(W1 * MM, 62 * MM)); fig.subplots_adjust(left=0.38, right=0.97, top=0.88, bottom=0.2)
    lab = {"primary": "Primary", "S1_per_copy": "S1 per-copy", "S2_three_per_side": "S2 3 per side", "S3_no_Pisum": f"S3 no {italic('Pisum')}",
           "S4_thinned_5M": "S4 thinned 5 M", "C": "C HOGs"}
    f = pd.concat([T2[T2.test.str.startswith("secondary")][["variant", "estimate", "lo", "hi"]],
                   T3C[T3C.test.str.startswith("secondary")][["estimate", "lo", "hi"]].assign(variant="C")]).set_index("variant").loc[list(lab)].reset_index()
    f.to_csv(f"{SRC}/FigS3_secondary.csv", index=False)
    forest(ax, list(zip(f.estimate, f.lo, f.hi)), [lab[v] for v in f.variant], (-0.09, 0.11))
    zone_labels(ax, len(f) - 0.25, xs=(-0.032, 0.0375, 0.08)); ax.set_xticks([-0.05, 0, 0.05, 0.1])
    ax.set_xlabel("Mean ρ$_{between}$ − mean ρ$_{within}$")
    save(fig, "FigS3")


def figS2():
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2 * MM, 64 * MM), gridspec_kw=dict(width_ratios=[1.4, 1], wspace=0.5))
    fig.subplots_adjust(left=0.12, right=0.98, top=0.78, bottom=0.2)
    d = A3P[np.isfinite(A3P.r_P_ks)].set_index("species").loc[[s for s in TREE_ORDER if s in set(A3P[np.isfinite(A3P.r_P_ks)].species)]].reset_index()
    d.to_csv(f"{SRC}/FigS6a_paralogs.csv", index=False)
    y = np.arange(len(d))[::-1]
    a.axvline(0, color="#555555", lw=0.5); a.axvspan(-NEGLIGIBLE, NEGLIGIBLE, color=ZONE["tech"], lw=0, zorder=0)
    a.scatter(d.r_P_ks, y, s=14, color=ACC, label="A3: Pearson r (noise-corrected divergence, Ks)", zorder=3, lw=0)
    a.scatter(d.rho_N_ks, y, s=14, color="#D9A021", label="A2: Spearman ρ (noise floor N, Ks)", zorder=3, lw=0, marker="s")
    a.set_yticks(y, [SHORT[s] for s in d.species], fontstyle="italic"); a.tick_params(axis="y", length=0)
    a.set_xlabel("Correlation with Ks"); a.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=6, ncol=1)
    tag(a, "a", x=-0.22, y=1.13)
    b1, bme, b2 = T3AB.loc["B1", "estimate"], float(re.search(r"rho\(sigma2_me, +copies\) ([+-][0-9.]+)", r3).group(1)), T3AB.loc["B2"]
    pd.DataFrame(dict(model=["species means only", "with within-species variance"], rho=[b1, bme])).to_csv(f"{SRC}/FigS6b_bm.csv", index=False)
    b.bar([0, 1], [b1, bme], width=0.55, color=["#9DB4CC", ACC], lw=0)
    for x, v in ((0, b1), (1, bme)):
        b.text(x, v + 0.002, num(v), ha="center", va="bottom", fontsize=6.2)
    b.set_xticks([0, 1], ["species means\nonly", "with within-species\nvariance"]); b.set_ylim(0, 0.09)
    b.set_ylabel("ρ (Brownian-motion rate,\nmean log2 copy number)")
    b.text(0.5, 0.083, f"Δρ = {num(b2.estimate)} (95% CI {num(b2.lo)} to {num(b2.hi, 5)})", ha="center", fontsize=6.2)
    tag(b, "b", x=-0.3, y=1.13)
    save(fig, "FigS6")


def figS3():
    fig, ax = plt.subplots(figsize=(W1 * MM, 66 * MM)); fig.subplots_adjust(left=0.2, right=0.97, top=0.93, bottom=0.18)
    T1A.to_csv(f"{SRC}/FigS4_distance_time.csv", index=False)
    ax.scatter(T1A.time, T1A.d_corrected, s=10, color=ACC, lw=0, alpha=0.85)
    sl, ic = np.polyfit(T1A.time, T1A.d_corrected, 1); xx = np.array([T1A.time.min(), T1A.time.max()])
    ax.plot(xx, ic + sl * xx, color="#222222", lw=0.9)
    r = T5.loc["1A_time_secondary"]
    ax.text(0.03, 0.96, f"change over range {num(r.estimate)}\n95% CI {num(r.lo)} to {num(r.hi)}", transform=ax.transAxes, va="top", fontsize=6.2)
    ax.set_xlabel("Divergence time (Myr)"); ax.set_ylabel("Noise-corrected expression distance\n(median, log2)")
    save(fig, "FigS4")


def figS4():
    fig, ax = plt.subplots(figsize=(W1 * MM, 56 * MM)); fig.subplots_adjust(left=0.48, right=0.97, top=0.93, bottom=0.2)
    reg = T5.loc["2B_trend"]
    rows = [(reg.estimate, reg.lo, reg.hi, "Registered estimate")] + \
           [(r.estimate, r.lo, r.hi, lab) for r, lab in zip(ROB.itertuples(), [f"without {italic('G. max')}, {italic('L. angustifolius')}",
                                                                               "≥ 100 shared pairs", "corrected for reliability", "reliability as covariate"])]
    pd.DataFrame(rows, columns=["estimate", "lo", "hi", "check"]).to_csv(f"{SRC}/FigS8_robustness.csv", index=False)
    for i, (e, lo, hi, _) in enumerate(rows):
        y = len(rows) - 1 - i
        if np.isfinite(lo):
            ax.plot([lo, hi], [y, y], color="#222222", lw=0.8)
        ax.plot(e, y, "o", color=ACC if i == 0 else "#222222", ms=3.6, mec="white", mew=0.4)
    ax.axvline(0, color="#555555", lw=0.5); ax.axvline(-MATERIAL, color="#999999", lw=0.5, ls=(0, (3, 2)))
    ax.set_yticks(range(len(rows)), [r[3] for r in rows][::-1]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.35, 0.05); ax.set_xlabel("Change in asymmetry conservation\nover the divergence range")
    save(fig, "FigS8")


def figS5():
    fig, ax = plt.subplots(figsize=(W1 * MM, 70 * MM)); fig.subplots_adjust(left=0.17, right=0.97, top=0.85, bottom=0.17)
    q = pd.concat([QC1all.assign(study="between-species set"), QC2all.assign(study="second studies")])
    q.to_csv(f"{SRC}/FigS1_qc.csv", index=False)
    for st, mk in (("between-species set", "o"), ("second studies", "^")):
        d = q[(q.study == st) & q.qc_pass]
        ax.scatter(d.reads_after / 1e6, d.pct_mapped, s=9, marker=mk, c=[CL_COL[CLADE[s]] for s in d.species], lw=0, alpha=0.85)
        ax.scatter([], [], s=9, marker=mk, color="#777777", label=st)
    x = q[~q.qc_pass]
    ax.scatter(x.reads_after / 1e6, x.pct_mapped, s=22, facecolor="none", edgecolor="#C00000", lw=0.8, label="excluded")
    ax.axvline(5, color="#999999", lw=0.5, ls=(0, (3, 2))); ax.axhline(70, color="#999999", lw=0.5, ls=(0, (3, 2)))
    ax.set_xscale("log"); ax.set_xlabel("Reads after trimming (millions, log scale)"); ax.set_ylabel("salmon mapping rate (%)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=6, handletextpad=0.2, columnspacing=0.8)
    save(fig, "FigS1")


for f in (figS1, figS2, figS3, figS4, figS5):
    f()
print("supplementary figures written")


# ======================================================================= Graphical abstract (13 x 5 cm, Elsevier minimum 1328 x 531 px)
def graphical_abstract():
    fig = plt.figure(figsize=(130 * MM, 50 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.9, 1.0], wspace=0.55, left=0.06, right=0.98, top=0.82, bottom=0.2)
    a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
    for k, (t1, col) in enumerate((("Within one study", "#7A7A7A"), ("Two studies, same species", "#D9A021"), (f"Two species, {TT.median_mya.min():.1f}–{TT.median_mya.max():.0f} Myr", ACC))):
        y = 8.2 - k * 3.1
        a.add_patch(FancyBboxPatch((0.2, y - 1.0), 9.6, 2.0, boxstyle="round,pad=0.02,rounding_size=0.3",
                                   fc=matplotlib.colors.to_rgba(col, 0.12), ec=col, lw=0.7))
        a.text(5, y, t1, ha="center", va="center", fontsize=6.0, fontweight="bold")
    a.set_title("12 legumes, public leaf RNA-seq", fontsize=6.8, pad=4)
    b = fig.add_subplot(gs[1]); t = T4.set_index("level").loc[["Gmax", "Pvul", "Mtru", "all"]].reset_index()
    b.set_xlim(-0.08, 0.06)
    b.axvline(0, color="#555555", lw=0.5)
    for i, r in t.iterrows():
        y = 3 - i; col = "#222222" if r.level == "all" else FOCAL_COL[r.level]
        b.plot([r.lo, r.hi], [y, y], color=col, lw=0.9); b.plot(r.delta, y, "D" if r.level == "all" else "o", color=col, ms=3.4, mec="white", mew=0.4)
    b.set_yticks(range(4), ["all pairs", "M. truncatula", "P. vulgaris", "G. max"])
    for lbl in b.get_yticklabels()[1:]:
        lbl.set_fontstyle("italic")
    b.tick_params(axis="y", length=0, labelsize=6); b.tick_params(axis="x", labelsize=6)
    b.set_xlabel("ρ$_{species}$ − ρ$_{study}$", fontsize=6.5)
    b.set_title("Copy-number gradient:\nlaboratory > species", fontsize=6.8, pad=4)
    c = fig.add_subplot(gs[2]); caj = (P2B.sp_a == "Ccaj") | (P2B.sp_b == "Ccaj"); p = P2B[~caj]
    c.scatter(p.time, p.rho_r, s=6, color=ACC, lw=0, alpha=0.85)
    sl, ic = np.polyfit(p.time, p.rho_r, 1); xx = np.array([p.time.min(), p.time.max()]); c.plot(xx, ic + sl * xx, color="#222222", lw=0.9)
    labc = CE[CE.ceiling == "laboratory"].rho; c.axhspan(labc.min(), labc.max(), color="#F3E5C2", lw=0, zorder=0)
    c.text(p.time.max(), labc.max() + 0.012, "two studies, same species", fontsize=5.5, color="#8A6D1F", ha="right", va="bottom")
    c.set_xlabel("Divergence time (Myr)", fontsize=6.5); c.set_ylabel("Ohnolog asymmetry\nconservation (ρ)", fontsize=6.5)
    c.tick_params(labelsize=6); c.set_ylim(0.25, 0.95)
    c.set_title("Ancient WGD asymmetry\nerodes with time", fontsize=6.8, pad=4)
    check_layout(fig, "GraphicalAbstract")
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {}), ("tif", {"dpi": 600})):
        fig.savefig(os.path.join(FIG, f"GraphicalAbstract.{ext}"), **kw)
    plt.close(fig)


# superseded by the editable PowerPoint version: pipeline/11_manuscript/make_graphical_abstract.js
# graphical_abstract()


# ======================================================================= Fig. S2, S5, S7 (added after the figure audit)
def copy_class(c):
    return np.where(c >= 5, "≥5", c.astype(int).astype(str))


def figS6():
    """Shape of the gradient: median divergence by copy-number class within one study, between two studies
    and between species, for the three species with two studies (display binning of committed per-orthogroup values)."""
    fig, axes = plt.subplots(1, 3, figsize=(W2 * MM, 62 * MM), sharey=True, gridspec_kw=dict(wspace=0.12))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.8, bottom=0.18)
    order = ["1", "2", "3", "4", "≥5"]; rows = []
    for ax, A in zip(axes, FOCAL):
        w = pd.read_csv(f"{RES}/results_step2/div_within_primary_{A}.tsv.gz", sep="\t")
        st_ = pd.read_csv(f"{RES}/results_step4/div_study_{A}.tsv.gz", sep="\t")
        bs = pd.concat([pd.read_csv(f"{RES}/results_step2/div_between_primary_{'_'.join(sorted([A, B], key=SPECIES.index))}.tsv.gz", sep="\t")
                        for B in SPECIES if B != A])
        for lab, df, col, colour in (("within one study", w, "div", "#7A7A7A"), ("two studies, same species", st_, "div_study", "#D9A021"),
                                     ("two species (11 pairs pooled)", bs, "div", ACC)):
            g = df.assign(cls=copy_class(df.copies.values)).groupby("cls")[col].agg(["median", "size"]).reindex(order)
            ax.plot(range(5), g["median"].values, marker="o", ms=3, lw=1.0, color=colour, label=lab)
            for k, cl in enumerate(order):
                rows.append(dict(species=A, comparison=lab, copy_class=cl, median_divergence=g.loc[cl, "median"], n=g.loc[cl, "size"]))
        ax.set_xticks(range(5), order); ax.set_xlabel("Copies in the orthogroup")
        ax.set_title(italic(LATIN[A]), fontsize=6.6, pad=3)
    pd.DataFrame(rows).to_csv(f"{SRC}/FigS2_shape.csv", index=False)
    axes[0].set_ylabel("Median expression divergence (log2)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.53, 1.0), ncol=3, fontsize=6.2, columnspacing=1.2)
    save(fig, "FigS2")


def figS7():
    """EXPLORATORY (E3): within-study gradient against median read depth per species."""
    Z = pd.read_csv(f"{RES}/exploratory_step2/within_rho_vs_depth.tsv", sep="\t")
    Z.to_csv(f"{SRC}/FigS5_depth.csv", index=False)
    fig, ax = plt.subplots(figsize=(W1 * MM, 66 * MM)); fig.subplots_adjust(left=0.18, right=0.97, top=0.95, bottom=0.18)
    for r in Z.itertuples():
        ax.scatter(r.median_reads_M, r.rho, s=14 + 2.5 * r.n, color=CL_COL[CLADE[r.species]], lw=0, alpha=0.9)
        off = {"Vrad": (4, -9), "Ccaj": (-4, -9), "Tpra": (4, 3)}.get(r.species, (3, 2))
        ax.annotate(SHORT[r.species], (r.median_reads_M, r.rho), xytext=off, textcoords="offset points", fontsize=5.6, fontstyle="italic",
                    ha="right" if off[0] < 0 else "left")
    rs = stats_spearman(Z.median_reads_M, Z.rho)
    ax.text(0.97, 0.04, f"Spearman ρ = {num(rs)} (exploratory)", transform=ax.transAxes, ha="right", fontsize=6.2)
    ax.set_xscale("log"); ax.set_xlim(4.5, 140); ax.set_xticks([5, 10, 20, 50, 100], ["5", "10", "20", "50", "100"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Median reads after trimming (millions, log scale)"); ax.set_ylabel("ρ within one study")
    save(fig, "FigS5")


def figS8():
    """EXPLORATORY (pre-declared in the Step 3 specification, no decision rule): median D - N by Ks bin for collinear
    (WGD/segmental) and other paralog pairs, pooled over the species with a half-split."""
    C3 = pd.read_csv(f"{RES}/results_step3/A_exploratory_class_by_ks.tsv", sep="\t"); C3.to_csv(f"{SRC}/FigS7_class_by_ks.csv", index=False)
    bins = list(dict.fromkeys(C3.ks_bin))
    fig, ax = plt.subplots(figsize=(W1 * MM, 62 * MM)); fig.subplots_adjust(left=0.18, right=0.97, top=0.85, bottom=0.2)
    for cls, colour, lab in (("WGD/segmental", "#C8553D", "collinear (WGD/segmental)"), ("other", "#555555", "other (tandem, proximal, dispersed)")):
        d = C3[C3["class"] == cls].set_index("ks_bin").reindex(bins)
        ax.plot(range(len(bins)), d.median_D_minus_N.values, marker="o", ms=3.2, lw=1.0, color=colour, label=lab)
    ax.set_xticks(range(len(bins)), bins, fontsize=6); ax.set_xlabel("Ks bin"); ax.set_ylabel("Median D − N (log2)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=6.2, ncol=1)
    save(fig, "FigS7")


def stats_spearman(x, y):
    from scipy.stats import spearmanr
    return spearmanr(x, y).correlation


for f in (figS6, figS7, figS8):
    f()
print("figures S6-S8 written")


print(f"layout check: {len(LAYOUT_ISSUES)} issue(s)")
for s_ in dict.fromkeys(LAYOUT_ISSUES):
    print("  " + s_)
