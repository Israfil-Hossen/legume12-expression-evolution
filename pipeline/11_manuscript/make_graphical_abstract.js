// Graphical abstract (Elsevier proportions 1328 x 531 px; readable at 13 x 5 cm) as an editable PowerPoint slide.
// Story, left to right: the usual inference -> three reference levels -> verdict -> evolution still visible -> lesson.
// Every number is read from the ledger; the small tree and bars in panel 3 are a labelled schematic.
// Usage (repository root): NODE_PATH=<dir with pptxgenjs, react-icons, react, react-dom, sharp> node pipeline/11_manuscript/make_graphical_abstract.js
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const { FaLeaf, FaFlask, FaRegClock, FaStar, FaDna } = require("react-icons/fa");

const OUT = "manuscript_mpe/figures/GraphicalAbstract.pptx";
const L = {};
for (const line of fs.readFileSync("results/manuscript_numbers/ledger.tsv", "utf8").split(/\r?\n/).slice(1)) {
  const f = line.split("\t");
  if (f[0]) L[f[0]] = f[1];
}
const need = (k) => { if (!(k in L)) throw new Error("ledger key missing: " + k); return L[k]; };
const n = (k) => need(k).replace(/^-/, "−");

const C = { ink: "1F2A26", green: "2F6B4F", greenLt: "E7F0EA", study: "C8553D", studyLt: "FBE9E4", evo: "2E6E8E", evoLt: "E4EEF3",
            grey: "7C8782", line: "C9D2CD", card: "F5F7F6", white: "FFFFFF", gold: "D9A021" };
const FONT = "Arial";
const pres = new pptxgen();
pres.defineLayout({ name: "GA", width: 13.28, height: 5.31 });
pres.layout = "GA";
pres.title = "Graphical abstract";
const s = pres.addSlide();
s.background = { color: C.white };

async function icon(Comp, color) {
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(Comp, { color: "#" + color, size: 256 }));
  return "image/png;base64," + (await sharp(Buffer.from(svg)).png().toBuffer()).toString("base64");
}
const T = (text, o) => s.addText(text, { fontFace: FONT, color: C.ink, margin: 0, isTextBox: true, valign: "top", ...o });
const box = (x, y, w, h, fill, line, r = 0.1) =>
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: line || fill, width: line ? 1 : 0 }, rectRadius: r });
const rect = (x, y, w, h, fill) => s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: fill, width: 0 } });
const hline = (x, y, w, color, width = 1.25, dash) => s.addShape(pres.shapes.LINE, { x, y, w, h: 0, line: { color, width, ...(dash ? { dashType: dash } : {}) } });
const vline = (x, y, h, color, width = 1.25) => s.addShape(pres.shapes.LINE, { x, y, w: 0, h, line: { color, width } });
function heading(x, y, num, text, color) {
  s.addShape(pres.shapes.OVAL, { x, y, w: 0.4, h: 0.4, fill: { color }, line: { color, width: 0 } });
  T(String(num), { x, y, w: 0.4, h: 0.4, fontSize: 15, bold: true, color: C.white, align: "center", valign: "middle" });
  T(text, { x: x + 0.5, y: y - 0.02, w: 3.6, h: 0.44, fontSize: 15, bold: true, color, valign: "middle" });
}
const arrow = (x, y) => s.addShape(pres.shapes.RIGHT_ARROW, { x, y, w: 0.26, h: 0.4, fill: { color: "B5C1BB" }, line: { color: "B5C1BB" } });

(async () => {
  const leafG = await icon(FaLeaf, C.green), leafG2 = await icon(FaLeaf, "7FAF8F"), leafB = await icon(FaLeaf, C.evo);
  const flaskS = await icon(FaFlask, C.study), flaskG = await icon(FaFlask, C.gold), clock = await icon(FaRegClock, C.evo);
  const star = await icon(FaStar, C.gold), dna = await icon(FaDna, C.green);

  // headline question
  T("Do larger gene families diverge more in expression because of evolution?",
    { x: 0.3, y: 0.16, w: 12.68, h: 0.5, fontSize: 19, bold: true, color: C.green, valign: "middle" });

  // ===== 1. the usual inference
  const XA = 0.3, WA = 3.0, Y0 = 0.82, H0 = 3.78;
  box(XA, Y0, WA, H0, C.card);
  heading(XA + 0.15, Y0 + 0.14, 1, "The usual inference", C.green);
  T("copies             difference between species", { x: XA + 0.2, y: 1.5, w: WA - 0.3, h: 0.26, fontSize: 10, color: C.grey, italic: true });
  const rows = [["single-copy gene", 1, 0.55], ["gene family", 4, 1.35]];
  for (let i = 0; i < rows.length; i++) {
    const y = 1.86 + i * 0.78;
    for (let k = 0; k < rows[i][1]; k++) rect(XA + 0.22 + k * 0.17, y + 0.06, 0.12, 0.34, C.green);   // gene copies
    T(rows[i][0], { x: XA + 0.2, y: y + 0.43, w: 1.4, h: 0.24, fontSize: 10.5 });
    hline(XA + 1.2, y + 0.23, rows[i][2], C.evo, 2.25);                                         // divergence bar
    s.addShape(pres.shapes.OVAL, { x: XA + 1.2 + rows[i][2] - 0.05, y: y + 0.18, w: 0.1, h: 0.1, fill: { color: C.evo }, line: { color: C.evo } });
  }
  T("larger families differ more", { x: XA + 0.2, y: 3.4, w: WA - 0.4, h: 0.26, fontSize: 11, color: C.evo, bold: true, align: "center" });
  box(XA + 0.22, 3.78, WA - 0.44, 0.62, C.white, C.line, 0.08);
  T([{ text: "read as relaxed constraint,", options: { breakLine: true } }, { text: "yet each species is another study", options: { bold: true, color: C.study } }],
    { x: XA + 0.26, y: 3.8, w: WA - 0.52, h: 0.58, fontSize: 10, valign: "middle", align: "center" });
  arrow(XA + WA + 0.04, 2.5);

  // ===== 2. three references
  const XB = 3.62, WB = 5.42;
  box(XB, Y0, WB, H0, C.card);
  heading(XB + 0.15, Y0 + 0.14, 2, "Three references, one pipeline", C.study);
  s.addImage({ data: dna, x: XB + WB - 0.42, y: Y0 + 0.17, w: 0.3, h: 0.3 });
  T(`12 legumes · ${need("scale_samples_total")} libraries · ${need("scale_fragments_billion")} billion reads`,
    { x: XB + 0.2, y: 1.48, w: WB - 0.4, h: 0.26, fontSize: 10.5, color: C.grey, italic: true });
  const lev = [
    { leaves: [leafG, leafG], flasks: [flaskS], name: "Replicates, one study", sub: "nothing evolves",
      res: [`gradient in ${need("S2_n_within")}/${need("S2_n_within")} species`, `\u03c1 ${n("S2_rho_within_min")} to ${n("S2_rho_within_max")}`] },
    { leaves: [leafG, leafG], flasks: [flaskS, flaskG], name: "Two studies, one species", sub: "nothing evolves",
      res: [`gradient even larger`, `than between species (${need("S4_n_delta_negative")}/${need("S4_n_pairs")})`] },
    { leaves: [leafG, leafB], flasks: [], name: "Two species", sub: `${need("time_min_Myr_noCcaj")}–${need("time_max_Myr_noCcaj")} Myr apart`,
      res: [`no rise with divergence time`, `(change ${n("S2_primary_prim_est")})`] },
  ];
  for (let i = 0; i < lev.length; i++) {
    const y = 1.86 + i * 0.74, r = lev[i];
    s.addImage({ data: r.leaves[0], x: XB + 0.22, y: y + 0.02, w: 0.3, h: 0.3 });
    s.addImage({ data: r.leaves[1], x: XB + 0.58, y: y + 0.02, w: 0.3, h: 0.3 });
    if (r.flasks.length === 1) s.addImage({ data: r.flasks[0], x: XB + 0.42, y: y + 0.34, w: 0.24, h: 0.24 });
    if (r.flasks.length === 2) { s.addImage({ data: r.flasks[0], x: XB + 0.25, y: y + 0.34, w: 0.24, h: 0.24 }); s.addImage({ data: r.flasks[1], x: XB + 0.61, y: y + 0.34, w: 0.24, h: 0.24 }); }
    if (i === 2) s.addImage({ data: clock, x: XB + 0.43, y: y + 0.34, w: 0.24, h: 0.24 });
    T([{ text: r.name, options: { bold: true, breakLine: true } }, { text: r.sub, options: { color: C.grey } }],
      { x: XB + 1.0, y: y + 0.02, w: 1.95, h: 0.56, fontSize: 10.5, valign: "middle" });
    box(XB + 2.98, y + 0.02, WB - 3.18, 0.56, C.studyLt, null, 0.08);
    T([{ text: r.res[0], options: { bold: true, breakLine: true } }, { text: r.res[1], options: {} }],
      { x: XB + 3.08, y: y + 0.02, w: WB - 3.36, h: 0.56, fontSize: 10.5, valign: "middle" });
  }
  box(XB + 0.22, 4.02, WB - 0.44, 0.48, C.study, null, 0.08);
  T([{ text: "Study effects mimic evolution", options: { bold: true, breakLine: true } },
     { text: `a change of study reproduces ${need("S5_1A_F_species_min")}\u2013${need("S5_1A_F_species_max")} of the species distance`, options: { fontSize: 10 } }],
    { x: XB + 0.3, y: 4.02, w: WB - 0.6, h: 0.48, fontSize: 12.5, color: C.white, valign: "middle", align: "center" });
  arrow(XB + WB + 0.04, 2.5);

  // ===== 3. evolution still visible (schematic tree with ohnolog pairs at the tips)
  const XC = 9.36, WC = 3.62;
  box(XC, Y0, WC, H0, C.card);
  heading(XC + 0.15, Y0 + 0.14, 3, "Evolution still visible", C.evo);
  T(`${need("S5_2_n_ohnolog_OG")} putative ohnolog pairs, shared genome duplication`,
    { x: XC + 0.2, y: 1.42, w: WC - 0.4, h: 0.4, fontSize: 10.5, color: C.grey, italic: true });
  const rx = XC + 0.35, tipx = XC + 1.75, tips = [1.98, 2.36, 2.88, 3.42];
  const n1 = XC + 1.25, n2 = XC + 0.95, n3 = XC + 0.65;                     // internal nodes: (t0,t1) -> n1; +t2 -> n2; +t3 -> n3 = root side
  const ym = (a, b) => (a + b) / 2;
  const y01 = ym(tips[0], tips[1]), y012 = ym(y01, tips[2]), yroot = ym(y012, tips[3]);
  vline(n1, tips[0], tips[1] - tips[0], C.ink); hline(n1, tips[0], tipx - n1, C.ink); hline(n1, tips[1], tipx - n1, C.ink);
  vline(n2, y01, tips[2] - y01, C.ink); hline(n2, y01, n1 - n2, C.ink); hline(n2, tips[2], tipx - n2, C.ink);
  vline(n3, y012, tips[3] - y012, C.ink); hline(n3, y012, n2 - n3, C.ink); hline(n3, tips[3], tipx - n3, C.ink);
  hline(rx, yroot, n3 - rx, C.ink);
  s.addImage({ data: star, x: rx - 0.16, y: yroot - 0.16, w: 0.32, h: 0.32 });
  T("WGD", { x: rx - 0.25, y: yroot + 0.17, w: 0.5, h: 0.2, fontSize: 9, align: "center", color: C.grey });
  // copy A (dark) and copy B (light) at each tip: identical in close relatives, fading with distance
  const bars = [[0.42, 0.14], [0.42, 0.16], [0.36, 0.22], [0.27, 0.27]];
  for (let i = 0; i < 4; i++) {
    const y = tips[i];
    rect(tipx + 0.08, y + 0.09 - bars[i][0] * 0.62, 0.13, bars[i][0] * 0.62, C.evo);
    rect(tipx + 0.24, y + 0.09 - bars[i][1] * 0.62, 0.13, bars[i][1] * 0.62, "9DC3D6");
  }
  T("A  B", { x: tipx + 0.06, y: tips[3] + 0.12, w: 0.36, h: 0.2, fontSize: 9, color: C.grey });
  T([{ text: "close relatives", options: { bold: true, breakLine: true } }, { text: "same copy higher", options: { breakLine: true } }, { text: `\u03c1 = ${n("S5_2B_closest_rho")}`, options: {} }],
    { x: tipx + 0.48, y: tips[0] - 0.16, w: 1.38, h: 0.62, fontSize: 10 });
  T([{ text: "most distant pairs", options: { bold: true, breakLine: true } }, { text: "pattern fades", options: { breakLine: true } }, { text: `median \u03c1 = ${n("S5_2B_farthest_rho_median")}`, options: {} }],
    { x: tipx + 0.48, y: tips[3] - 0.42, w: 1.38, h: 0.62, fontSize: 10 });
  box(XC + 0.22, 4.02, WC - 0.44, 0.48, C.evo, null, 0.08);
  T([{ text: "Asymmetry conserved, erodes with time", options: { bold: true, breakLine: true } },
     { text: `change ${n("S5_2B_trend_est")} over ${need("time_min_Myr_noCcaj")}–${need("time_max_Myr_noCcaj")} Myr`, options: { fontSize: 9.5 } }],
    { x: XC + 0.3, y: 4.02, w: WC - 0.6, h: 0.48, fontSize: 11, color: C.white, valign: "middle", align: "center" });

  // ===== lesson
  T([{ text: "Lesson: ", options: { bold: true, color: C.green } },
     { text: "read cross-species expression differences only against within-study and between-study references", options: { italic: true } }],
    { x: 0.3, y: 4.76, w: 12.68, h: 0.38, fontSize: 13, align: "center", valign: "middle" });

  await pres.writeFile({ fileName: OUT });
  console.log(`wrote ${OUT}`);
})();
