#!/usr/bin/env Rscript
# Main figures for the 12-legume manuscript, redrawn in R (ggplot2 + ggtree + patchwork) to Nature-style
# specifications: 183 mm (double column) width, Arial, 5-7 pt text, 8 pt bold lowercase panel letters,
# 0.25-0.5 pt lines, colour-blind-safe clade palette.
# Every plotted value is read from the same committed result files as make_figures.py.
# Usage (repository root): Rscript pipeline/11_manuscript/make_figures.R [fig1 fig2 ...]

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ape); library(ggtree)
})

RES <- "results"; OUT <- "manuscript_mpe/figures/R_version"
dir.create(file.path(OUT, "source_data"), recursive = TRUE, showWarnings = FALSE)
MM <- 1 / 25.4; W2 <- 183
MATERIAL <- 0.050; NEGLIGIBLE <- 0.025
FONT <- "Arial"

SPECIES <- c("Gmax", "Pvul", "Vrad", "Vung", "Ccaj", "Mtru", "Cari", "Psat", "Tpra", "Ljap", "Lang", "Adur")
LATIN <- c(Gmax = "Glycine max", Pvul = "Phaseolus vulgaris", Vrad = "Vigna radiata", Vung = "Vigna unguiculata",
           Ccaj = "Cajanus cajan", Mtru = "Medicago truncatula", Cari = "Cicer arietinum", Psat = "Pisum sativum",
           Tpra = "Trifolium pratense", Ljap = "Lotus japonicus", Lang = "Lupinus angustifolius", Adur = "Arachis duranensis")
SHORT <- c(Gmax = "G. max", Pvul = "P. vulgaris", Vrad = "V. radiata", Vung = "V. unguiculata", Ccaj = "C. cajan",
           Mtru = "M. truncatula", Cari = "C. arietinum", Psat = "P. sativum", Tpra = "T. pratense", Ljap = "L. japonicus",
           Lang = "L. angustifolius", Adur = "A. duranensis")
CLADE <- c(Gmax = "Millettioid", Pvul = "Millettioid", Vrad = "Millettioid", Vung = "Millettioid", Ccaj = "Millettioid",
           Mtru = "IRLC", Cari = "IRLC", Psat = "IRLC", Tpra = "IRLC", Ljap = "Robinioid", Lang = "Genistoid",
           Adur = "Dalbergioid")
CL_COL <- c(Millettioid = "#2F5D8C", IRLC = "#3A9278", Robinioid = "#B0569A", Genistoid = "#D9A021", Dalbergioid = "#C8553D")
ACC <- "#2F5D8C"; GREY <- "#8C8C8C"; INK <- "#222222"
EVO <- "#C8553D"; NONEVO <- "#3A9278"; STUDY <- "#D9A021"
ZONE <- c(tech = "#E3F0EC", part = "#F2F2F2", evo = "#FBE9DD")

minus <- function(x, d = 3) gsub("-", "−", sprintf(paste0("%+.", d, "f"), x))

theme_nat <- function(base = 6.5) {
  theme_classic(base_size = base, base_family = FONT) %+replace% theme(
    axis.line = element_line(linewidth = 0.25, colour = INK), axis.ticks = element_line(linewidth = 0.25, colour = INK),
    axis.ticks.length = unit(0.8, "mm"), axis.text = element_text(size = 6, colour = INK),
    axis.title = element_text(size = 6.5, colour = INK), legend.text = element_text(size = 6),
    legend.title = element_text(size = 6.5), legend.key.size = unit(2.6, "mm"), legend.background = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "mm"), plot.tag = element_text(size = 8, face = "bold", family = FONT),
    strip.background = element_blank(), strip.text = element_text(size = 6.5, face = "italic"))
}
tagged <- function(p, letter) p + labs(tag = letter) + theme(plot.tag.position = c(0, 1))

save_fig <- function(p, name, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = W2, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".tif")), p, width = W2, height = h, units = "mm", dpi = 600,
         device = "tiff", compression = "lzw", type = "cairo")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = W2, height = h, units = "mm", dpi = 300, type = "cairo")
  message("wrote ", name)
}

# ------------------------------------------------------------------ inputs (same files as make_figures.py)
qc <- fread(file.path(RES, "results_step2/qc_samples.tsv"))[qc_pass == TRUE]
q4 <- fread(file.path(RES, "results_step4/qc_samples_step4.tsv"))[qc_pass == TRUE]
W  <- fread(file.path(RES, "results_step2/rho_within_primary.tsv"))
B  <- fread(file.path(RES, "results_step2/rho_between_primary.tsv"))
TT <- fread("metadata/step2_divergence_times.tsv")
T2 <- fread(file.path(RES, "results_step2/tests.tsv"))
T3C <- fread(file.path(RES, "results_step3C/tests.tsv"))
SUP <- fread(file.path(RES, "tree_support/support_summary.tsv"))

tree <- read.tree(file.path(RES, "server_logs/concat.fa.treefile"))
tree <- drop.tip(tree, "Atha")           # rooted on A. thaliana in the analysis; the outgroup is pruned for display
stopifnot(setequal(tree$tip.label, SPECIES))
tree <- ladderize(tree, right = FALSE)

# =================================================================== Fig. 1
clade_legend <- function() {
  d <- data.frame(cl = factor(names(CL_COL), names(CL_COL)), x = 1:5, y = 1)
  ggplot(d) + geom_point(aes(x, y, colour = cl), shape = 15, size = 2, alpha = 0) +
    scale_colour_manual(values = CL_COL, name = NULL) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2.2))) +
    theme_void(base_family = FONT) +
    theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5), legend.text = element_text(size = 6),
          legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0))
}

fig1 <- function() {
  # ---- a: species tree with concordance factors and per-species data columns
  pt <- ggtree(tree, linewidth = 0.3, colour = INK)
  td <- as.data.table(pt$data)
  tipd <- td[isTip == TRUE, .(code = label, y, x)]
  parts <- prop.part(tree)
  nodes <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)
  cf <- rbindlist(lapply(seq_along(parts), function(k) {
    clade_tips <- sort(attr(parts, "labels")[parts[[k]]])
    keep <- sapply(strsplit(SUP$tips, ","), function(s) setequal(sort(s), clade_tips))
    hit <- SUP[keep]
    if (!nrow(hit)) return(NULL)
    data.table(node = nodes[k], lab = sprintf("%.0f/%.0f", hit$gCF, hit$sCF))
  }))
  cf <- merge(cf, td[, .(node, x, y, branch)], by = "node")
  stopifnot(nrow(cf) == nrow(SUP))
  fwrite(cf[, .(node, lab)], file.path(OUT, "source_data/Fig1a_concordance.csv"))
  xmax <- max(td$x)
  YL <- c(-1.2, 13.2)
  # place each gCF/sCF label at the first candidate position that touches no branch and no other label
  # (box sizes in data units from the panel geometry: ~52 mm for the x range, ~5.6 mm per tip on y)
  xr <- (xmax + 0.26) - (-0.02); ux <- xr / 52; uy <- 1 / 5.6
  cw <- function(s) (nchar(s) - 1) * 0.556 * 4.4 * 0.3528 * ux + 0.278 * 4.4 * 0.3528 * ux
  hh <- 1.7 * uy
  kids <- function(n) td[parent == n & node != n, y]
  segs <- rbind(td[!is.na(parent) & node != parent, .(x0 = td$x[match(parent, td$node)], x1 = x, y0 = y, y1 = y)],
                td[isTip == FALSE, .(x0 = x, x1 = x, y0 = sapply(node, function(n) min(kids(n))),
                                     y1 = sapply(node, function(n) max(kids(n))))])
  placed <- list(); out <- list()
  hit_seg <- function(b) any(segs$x1 >= b[1] & segs$x0 <= b[2] & segs$y1 >= b[3] & segs$y0 <= b[4])
  hit_lab <- function(b) any(vapply(placed, function(q) !(b[2] < q[1] || b[1] > q[2] || b[4] < q[3] || b[3] > q[4]), TRUE))
  for (i in seq_len(nrow(cf))) {
    w <- cw(cf$lab[i]); nx <- cf$x[i]; ny <- cf$y[i]
    px <- td$x[match(td$parent[match(cf$node[i], td$node)], td$node)]
    cands <- list(c(nx - 0.003 - w, nx - 0.003, ny + 0.06, ny + 0.06 + hh), c(nx - 0.003 - w, nx - 0.003, ny - 0.06 - hh, ny - 0.06),
                  c(nx + 0.003, nx + 0.003 + w, ny + 0.12, ny + 0.12 + hh), c(nx + 0.003, nx + 0.003 + w, ny - 0.12 - hh, ny - 0.12),
                  c(px - 0.003 - w, px - 0.003, ny + 0.06, ny + 0.06 + hh), c(px - 0.003 - w, px - 0.003, ny - 0.06 - hh, ny - 0.06))
    ok <- NULL
    for (bx in cands) if (bx[1] > -0.02 && !hit_seg(bx) && !hit_lab(bx)) { ok <- bx; break }
    if (is.null(ok)) stop("no free position for concordance label ", cf$lab[i])
    placed[[length(placed) + 1]] <- ok
    out[[i]] <- data.table(lab = cf$lab[i], xc = mean(ok[1:2]), yc = mean(ok[3:4]))
  }
  cfp <- rbindlist(out)
  pa <- pt +
    geom_text(data = cfp, aes(x = xc, y = yc, label = lab), size = 4.4 / .pt, colour = "#5A5A5A", family = FONT) +
    geom_segment(data = tipd, aes(x = x + 0.002, xend = xmax + 0.006, y = y, yend = y),
                 linewidth = 0.2, colour = "#BDBDBD", linetype = "11") +
    geom_text(data = tipd, aes(x = xmax + 0.009, y = y, label = LATIN[code], colour = CLADE[code]),
              hjust = 0, size = 6.2 / .pt, fontface = "italic", family = FONT) +
    scale_colour_manual(values = CL_COL, guide = "none") +
    annotate("segment", x = 0, xend = 0.05, y = -0.6, yend = -0.6, linewidth = 0.3, colour = INK) +
    annotate("text", x = 0.055, y = -0.6, label = "0.05 substitutions/site", hjust = 0, size = 5.2 / .pt, family = FONT) +
    annotate("text", x = 0.012, y = 12.75, label = "Species tree (gCF/sCF, %)", hjust = 0, size = 6 / .pt,
             fontface = "bold", family = FONT) +
    scale_x_continuous(limits = c(-0.02, xmax + 0.26), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    theme(plot.margin = margin(1, 0, 1, 0, "mm"))

  n1 <- qc[, .N, by = species]; n2 <- q4[, .N, by = species]
  dep <- qc[, .(reads = median(reads_after) / 1e6, mapped = median(pct_mapped)), by = species]
  sp <- merge(tipd[, .(code, y)], merge(merge(n1[, .(code = species, n1 = N)], n2[, .(code = species, n2 = N)], all.x = TRUE),
                                         dep[, .(code = species, reads, mapped)]), by = "code")
  sp[is.na(n2), n2 := 0L]; sp[, clade := CLADE[code]]
  fwrite(sp, file.path(OUT, "source_data/Fig1a_species.csv"))
  col_theme <- theme_nat() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                                   axis.text.y = element_blank(), axis.title.y = element_blank(),
                                   plot.margin = margin(1, 1.2, 1, 1.2, "mm"))
  head_txt <- function(lab, x) annotate("text", x = x, y = 12.75, label = lab, size = 6 / .pt, fontface = "bold",
                                        family = FONT)
  bars <- rbind(sp[, .(code, y, study = "1", n = n1, start = 0, clade)],
                sp[n2 > 0, .(code, y, study = "2", n = n2, start = n1, clade)])
  pb <- ggplot(bars) +
    geom_rect(aes(xmin = start, xmax = start + n, ymin = y - 0.3, ymax = y + 0.3, fill = clade, alpha = study),
              colour = NA) +
    geom_rect(data = bars[study == "2"], aes(xmin = start, xmax = start + n, ymin = y - 0.3, ymax = y + 0.3,
                                              colour = clade), fill = NA, linewidth = 0.25) +
    geom_text(data = sp, aes(x = n1 + n2 + 0.8, y = y, label = ifelse(n2 > 0, paste0(n1, "+", n2), n1)),
              hjust = 0, size = 5.4 / .pt, family = FONT) +
    head_txt("Replicates", 15) +
    scale_fill_manual(values = CL_COL, guide = "none") + scale_colour_manual(values = CL_COL, guide = "none") +
    scale_alpha_manual(values = c("1" = 1, "2" = 0.25), guide = "none") +
    scale_x_continuous(limits = c(0, 37), breaks = c(0, 10, 20), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    labs(x = "n (light: 2nd study)") + col_theme
  pr <- ggplot(sp, aes(y = y)) +
    geom_segment(aes(x = 0, xend = reads, yend = y), linewidth = 0.3, colour = "#CFCFCF") +
    geom_point(aes(x = reads, colour = clade), size = 1.2) +
    head_txt("Depth", 37.5) +
    scale_colour_manual(values = CL_COL, guide = "none") +
    scale_x_continuous(limits = c(0, 75), breaks = c(0, 30, 60), expand = c(0, 0)) +
    scale_y_continuous(limits = YL, expand = c(0, 0)) +
    labs(x = "Reads (M)") + col_theme
  pm <- ggplot(sp, aes(x = 1, y = y, fill = mapped)) +
    geom_tile(width = 0.92, height = 0.74, colour = NA) +
    geom_text(aes(label = sprintf("%.0f", mapped), colour = mapped > 85), size = 5 / .pt, family = FONT) +
    annotate("text", x = 1, y = 12.75, label = "Mapped
(%)", size = 5.6 / .pt, fontface = "bold", family = FONT, lineheight = 0.85, vjust = 0.3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
    scale_fill_gradient(low = "#E8EEF5", high = "#2F5D8C", limits = c(70, 96), guide = "none") +
    scale_y_continuous(limits = YL, expand = c(0, 0)) + scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    labs(x = NULL) + col_theme +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
  a_row <- pa + pb + pr + pm + plot_layout(widths = c(2.3, 0.9, 0.62, 0.36))

  # ---- b: the three levels of comparison (drawn schematic)
  lv <- data.table(cx = c(1.9, 6.5, 11.1),
                   head = c("Within\none study", "Between two studies\nof one species", "Between\ntwo species"),
                   stat = c("rho[within]", "rho[study]", "rho[between]"))
  boxes <- data.table(xmin = c(0.3, 4.5, 6.7, 8.95, 11.25), xmax = c(3.5, 6.3, 8.5, 10.95, 13.25),
                      col = c(GREY, GREY, STUDY, CL_COL[["Millettioid"]], CL_COL[["IRLC"]]),
                      lab = c("one study", "study 1", "study 2", "species A", "species B"))
  boxes[, `:=`(ymin = 2.05, ymax = 4.25)]
  dots <- rbindlist(list(
    data.table(x = seq(0.75, 3.05, length.out = 6), g = rep(c("h1", "h2"), each = 3)),
    data.table(x = c(4.98, 5.4, 5.82), g = "s1"), data.table(x = c(7.18, 7.6, 8.02), g = "s2"),
    data.table(x = c(9.58, 10.0, 10.42), g = "A"), data.table(x = c(11.78, 12.2, 12.62), g = "B")))
  dots[, y := 3.05]
  dcol <- c(h1 = "#4D4D4D", h2 = "#B3B3B3", s1 = "#4D4D4D", s2 = STUDY, A = CL_COL[["Millettioid"]], B = CL_COL[["IRLC"]])
  grid_tab <- data.table(cx = rep(lv$cx, 2), y = rep(c(0.95, 0.3), each = 3), val = c("no", "no", "yes", "no", "yes", "yes"))
  arr <- arrow(length = unit(0.7, "mm"), ends = "both", type = "closed")
  pbx <- ggplot() +
    geom_text(data = lv, aes(x = cx, y = 5.55, label = head), size = 6 / .pt, fontface = "bold", family = FONT,
              lineheight = 0.9) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = I(col)),
              fill = NA, linewidth = 0.35) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = 3.8, label = lab, colour = I(col)), size = 4.5 / .pt,
              family = FONT) +
    geom_point(data = dots, aes(x, y, colour = g), size = 1.7) +
    scale_colour_manual(values = dcol, guide = "none") +
    annotate("segment", x = c(0.6, 1.95), xend = c(1.85, 3.2), y = 2.55, yend = 2.55, linewidth = 0.3, colour = INK) +
    annotate("text", x = 1.9, y = 2.3, label = "replicate halves", size = 4.8 / .pt, family = FONT) +
    annotate("segment", x = 6.3, xend = 6.7, y = 3.05, yend = 3.05, linewidth = 0.3, colour = INK, arrow = arr) +
    annotate("segment", x = 10.95, xend = 11.25, y = 3.05, yend = 3.05, linewidth = 0.3, colour = INK, arrow = arr) +
    annotate("segment", x = c(10.0, 12.2, 10.0), xend = c(10.0, 12.2, 12.2), y = c(4.25, 4.25, 4.6),
             yend = c(4.6, 4.6, 4.6), linewidth = 0.35, colour = INK) +
    annotate("text", x = 11.1, y = 4.78, label = "9.7–56 Myr", size = 4.8 / .pt, family = FONT) +
    geom_text(data = lv, aes(x = cx, y = 1.6, label = stat), parse = TRUE, size = 6.5 / .pt, family = FONT) +
    annotate("text", x = -0.2, y = c(0.95, 0.3), label = c("Evolutionary divergence", "Study difference"),
             hjust = 1, size = 5.4 / .pt, family = FONT, colour = "#444444") +
    geom_label(data = grid_tab, aes(x = cx, y = y, label = val, fill = I(ifelse(val == "yes", "#FBE0D0", "#EDEDED"))),
               size = 5.2 / .pt, family = FONT, linewidth = 0, label.padding = unit(0.5, "mm")) +
    coord_cartesian(xlim = c(-4.6, 13.4), ylim = c(0, 6.2), expand = FALSE, clip = "off") +
    theme_void(base_family = FONT) + theme(plot.margin = margin(1, 1, 0, 1, "mm"))

  # ---- c: what each account predicts (schematic; levels arbitrary)
  LEV <- c("within\nstudy", "between\nstudies", "between\nspecies")
  pred1 <- data.table(level = factor(rep(LEV, 2), levels = LEV), acc = rep(c("evolution", "non-evolutionary"), each = 3),
                      v = c(0.42, 0.42, 1.0, 0.38, 0.38, 0.38))
  pred2 <- data.table(t = rep(c(10, 56), 2), acc = rep(c("evolution", "non-evolutionary"), each = 2), v = c(0.55, 1.05, 0.38, 0.38))
  acc_col <- c(evolution = EVO, `non-evolutionary` = NONEVO); acc_lt <- c(evolution = "solid", `non-evolutionary` = "22")
  pc1 <- ggplot(pred1, aes(level, v, group = acc, colour = acc, linetype = acc)) + geom_line(linewidth = 0.5) +
    geom_point(size = 1.1) + scale_colour_manual(values = acc_col, name = NULL) +
    scale_linetype_manual(values = acc_lt, name = NULL) +
    scale_y_continuous(limits = c(0, 1.2), breaks = NULL) +
    labs(x = NULL, y = "Expected gradient (schematic)") + theme_nat() +
    theme(axis.text.x = element_text(size = 5.4, lineheight = 0.85), legend.position = "inside",
          legend.position.inside = c(0.02, 1.02), legend.justification = c(0, 1), legend.key.width = unit(4, "mm"))
  pc2 <- ggplot(pred2, aes(t, v, colour = acc, linetype = acc)) + geom_line(linewidth = 0.5) +
    scale_colour_manual(values = acc_col, guide = "none") +
    scale_linetype_manual(values = acc_lt, guide = "none") +
    scale_y_continuous(limits = c(0, 1.2), breaks = NULL) + scale_x_continuous(breaks = c(10, 30, 50)) +
    labs(x = "Divergence time (Myr)", y = "Expected ρ between species") + theme_nat()
  pc <- (pc1 + theme(plot.margin = margin(1, 2, 1, 5, "mm")) | pc2)

  p <- (wrap_elements(full = a_row) + labs(tag = "a") |
          ((wrap_elements(full = pbx) + labs(tag = "b")) / (wrap_elements(full = pc) + labs(tag = "c")) +
             plot_layout(heights = c(1, 1.05)))) +
    plot_layout(widths = c(1.2, 1))
  p <- (p / wrap_elements(full = clade_legend())) + plot_layout(heights = c(1, 0.045)) &
    theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(p, "Fig1", 100)
}

# =================================================================== Fig. 2
fig2 <- function() {
  w <- W[order(rho)][, code := species][, lab := SHORT[code]][, clade := CLADE[code]]
  w[, lab := factor(lab, levels = lab)]
  fwrite(w, file.path(OUT, "source_data/Fig2a_within.csv"))
  med_w <- median(w$rho)
  pa <- ggplot(w, aes(y = lab)) +
    geom_vline(xintercept = med_w, linewidth = 0.25, linetype = "22", colour = GREY) +
    geom_segment(aes(x = 0, xend = rho, yend = lab, colour = clade), linewidth = 0.55) +
    geom_point(aes(x = rho, colour = clade), size = 1.6) +
    annotate("text", x = med_w - 0.004, y = 11.85, label = paste0("median ", minus(med_w)), hjust = 1, vjust = 0,
             size = 5.2 / .pt, family = FONT, colour = "#555555") +
    annotate("text", x = 0.215, y = 11.85, label = paste0(sum(w$rho > 0), " of ", nrow(w), " > 0"), hjust = 1,
             vjust = 0, size = 5.4 / .pt, family = FONT) +
    scale_colour_manual(values = CL_COL, guide = "none") +
    scale_x_continuous(limits = c(0, 0.22), breaks = c(0, 0.1, 0.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = expression("Gradient within one study (" * rho[within] * ")"), y = NULL) + theme_nat() +
    theme(axis.text.y = element_text(face = "italic"), axis.ticks.y = element_blank(), plot.margin = margin(5, 2, 2, 2, "mm"))

  tt <- rbind(TT[, .(sp_a, sp_b, median_mya)], TT[, .(sp_a = sp_b, sp_b = sp_a, median_mya)])
  d <- merge(B, tt, by = c("sp_a", "sp_b"))
  d[, cajanus := sp_a == "Ccaj" | sp_b == "Ccaj"]
  d[, kind := fifelse(cajanus, "Cajanus pairs (excluded, D1)",
                      fifelse(CLADE[sp_a] == CLADE[sp_b], "same clade", "different clades"))]
  fwrite(d, file.path(OUT, "source_data/Fig2b_rho_between_time.csv"))
  p <- d[cajanus == FALSE]
  fit <- lm(rho ~ median_mya, data = p)
  pr <- T2[variant == "primary" & startsWith(test, "primary")]
  rng <- range(p$median_mya); span_t <- diff(rng)
  KIND <- c("different clades", "same clade", "Cajanus pairs (excluded, D1)")
  KLAB <- c("different clades", "same clade", expression(italic("Cajanus") ~ "pairs (excluded, D1)"))
  pb <- ggplot(d, aes(median_mya, rho)) +
    geom_hline(yintercept = med_w, linewidth = 0.3, linetype = "22", colour = "#555555") +
    annotate("text", x = 57.5, y = med_w + 0.002, label = "within-study median", hjust = 1,
             vjust = 0, size = 5.2 / .pt, family = FONT, colour = "#555555") +
    annotate("segment", x = rng[1], xend = rng[2], y = predict(fit, data.frame(median_mya = rng[1])),
             yend = predict(fit, data.frame(median_mya = rng[2])), linewidth = 0.5, colour = INK) +
    geom_point(aes(shape = kind, fill = kind, colour = kind), size = 1.35, stroke = 0.35) +
    scale_shape_manual(values = c(21, 24, 21) |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    scale_fill_manual(values = c(ACC, ACC, "white") |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    scale_colour_manual(values = c("white", "white", GREY) |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    annotate("text", x = 9, y = 0.006, hjust = 0, vjust = 0, size = 5.4 / .pt, family = FONT, lineheight = 0.95,
             label = sprintf("change over %.1f Myr: %s\n95%% interval %s to %s", span_t, minus(pr$estimate),
                             minus(pr$lo), minus(pr$hi))) +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05), expand = c(0, 0)) +
    scale_x_continuous(limits = c(8, 58), breaks = c(10, 20, 30, 40, 50), expand = c(0, 0)) +
    labs(x = "Divergence time (Myr)", y = expression("Gradient between species (" * rho[between] * ")")) + theme_nat() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 1.0), legend.justification = c(1, 1),
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 5.6))
  mg <- rbind(w[, .(level = "within", rho)], d[, .(level = "between", rho)])
  mg[, level := factor(level, levels = c("within", "between"))]
  set.seed(1)
  pm <- ggplot(mg, aes(level, rho)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, fill = NA, colour = "#555555") +
    geom_jitter(aes(colour = level), width = 0.13, height = 0, size = 0.7, alpha = 0.9) +
    scale_colour_manual(values = c(within = "#4D4D4D", between = ACC), guide = "none") +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "All values") + theme_nat() +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 5.4, angle = 35, hjust = 1),
          plot.title = element_text(size = 5.8, hjust = 0.5, margin = margin(b = 1)))

  lab <- c(primary = "Primary", S1_per_copy = "S1 per-copy", S2_three_per_side = "S2 three per side",
           S3_no_Pisum = "S3 without Pisum", S4_thinned_5M = "S4 thinned to 5 M", S5_patristic = "S5 patristic",
           C = "C hierarchical OGs")
  labx <- c(primary = "Primary", S1_per_copy = "S1~'per-copy'", S2_three_per_side = "S2~three~per~side",
            S3_no_Pisum = "S3~without~italic('Pisum')", S4_thinned_5M = "S4~thinned~to~5~M", S5_patristic = "S5~patristic",
            C = "C~hierarchical~OGs")
  f <- rbind(T2[!startsWith(test, "secondary"), .(variant, estimate, lo, hi)],
             T3C[startsWith(test, "primary"), .(variant = "C", estimate, lo, hi)])
  f <- f[match(names(lab), variant)][, label := factor(variant, levels = rev(names(lab)))]
  fwrite(f, file.path(OUT, "source_data/Fig2c_tests.csv"))
  zones <- data.table(xmin = c(-0.08, NEGLIGIBLE, MATERIAL), xmax = c(NEGLIGIBLE, MATERIAL, 0.1),
                      z = c("technical", "partial", "evolution"), fill = ZONE[c("tech", "part", "evo")])
  pc <- ggplot(f, aes(y = label)) +
    geom_rect(data = zones, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = I(fill)), inherit.aes = FALSE) +
    geom_text(data = zones, aes(x = (xmin + xmax) / 2, y = 7.85, label = z), size = 5.2 / .pt, family = FONT,
              colour = "#555555", inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#555555") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.45, colour = INK, orientation = "y") +
    geom_point(aes(x = estimate, fill = variant == "primary", size = variant == "primary"), shape = 21,
               colour = "white", stroke = 0.3) +
    scale_fill_manual(values = c(`TRUE` = ACC, `FALSE` = INK), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 1.7), guide = "none") +
    scale_x_continuous(limits = c(-0.08, 0.1), breaks = c(-0.05, 0, 0.05), expand = c(0, 0),
                       labels = function(x) minus(x, 2)) +
    scale_y_discrete(labels = function(v) parse(text = labx[v]), expand = expansion(add = c(0.6, 1.3))) +
    labs(x = "Change in ρ across the divergence range", y = NULL) + theme_nat() +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())

  pbm <- pb + pm + plot_layout(widths = c(1, 0.34))
  p <- (pa + labs(tag = "a") | (wrap_elements(full = pbm) + labs(tag = "b")) | pc + labs(tag = "c")) +
    plot_layout(widths = c(0.8, 1.55, 0.95)) & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(p, "Fig2", 72)
}

# =================================================================== Step 4 and 5 inputs
S5 <- file.path(RES, "results_step5")
P4 <- fread(file.path(RES, "results_step4/pairs_step4.tsv"))
T4 <- fread(file.path(RES, "results_step4/tests_step4.tsv"))
T5 <- fread(file.path(S5, "tests_step5.tsv")); setkey(T5, analysis)
P1A <- fread(file.path(S5, "1A_pairs.tsv")); P1B <- fread(file.path(S5, "1B_partial.tsv"))
P1C <- fread(file.path(S5, "1C_ambiguity.tsv")); P2A <- fread(file.path(S5, "2A_summary.tsv"))
P2B <- fread(file.path(S5, "2B_species_pairs.tsv")); CE <- fread(file.path(S5, "2B_ceilings.tsv"))
MRAT <- fread(file.path(S5, "2B_mean_ratios.tsv"))
D3 <- list(study1 = fread(file.path(S5, "3_dosage_Gmax_study1.tsv.gz")),
           study2 = fread(file.path(S5, "3_dosage_Gmax_study2.tsv.gz")))
FOCAL <- c("Gmax", "Pvul", "Mtru")
FOCAL_COL <- c(Gmax = "#2F5D8C", Pvul = "#C8553D", Mtru = "#3A9278")
TREE_ORDER <- c("Adur", "Lang", "Ljap", "Cari", "Psat", "Mtru", "Tpra", "Ccaj", "Gmax", "Pvul", "Vrad", "Vung")
italic_sp <- function(codes) parse(text = paste0("italic('", SHORT[codes], "')"))
zones_layer <- function(lo, hi, ylab_y = NULL) {
  z <- data.table(xmin = c(lo, NEGLIGIBLE, MATERIAL), xmax = c(NEGLIGIBLE, MATERIAL, hi),
                  fill = ZONE[c("tech", "part", "evo")])
  list(geom_rect(data = z, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = I(fill)), inherit.aes = FALSE),
       geom_vline(xintercept = 0, linewidth = 0.3, colour = "#555555"))
}

# =================================================================== Fig. 3
fig3 <- function() {
  fwrite(P4, file.path(OUT, "source_data/Fig3a_pairs.csv")); fwrite(T4, file.path(OUT, "source_data/Fig3b_delta.csv"))
  sl <- rbind(P4[, .(focal, other, level = "two studies\n(same species)", rho = rho_study)],
              P4[, .(focal, other, level = "two species", rho = rho_species)])
  sl[, level := factor(level, levels = c("two studies\n(same species)", "two species"))]
  sl[, focal := factor(focal, levels = FOCAL)]
  mn <- sl[, .(rho = mean(rho)), by = .(focal, level)]
  pa <- ggplot(sl, aes(level, rho, group = other)) +
    geom_line(aes(colour = focal), linewidth = 0.3, alpha = 0.45) +
    geom_point(aes(colour = focal), size = 0.6, alpha = 0.6) +
    geom_line(data = mn, aes(group = focal, colour = focal), linewidth = 0.9) +
    geom_point(data = mn, aes(group = focal, fill = focal), shape = 21, colour = "white", size = 2.1, stroke = 0.4) +
    facet_wrap(~focal, nrow = 1, labeller = as_labeller(function(v) LATIN[v])) +
    scale_colour_manual(values = FOCAL_COL, guide = "none") + scale_fill_manual(values = FOCAL_COL, guide = "none") +
    scale_y_continuous(limits = c(0.03, 0.19), breaks = seq(0.05, 0.15, 0.05)) +
    scale_x_discrete(expand = expansion(add = 0.3)) +
    labs(x = NULL, y = "Copy-number gradient ρ\n(identical orthogroups)") + theme_nat() +
    theme(axis.text.x = element_text(size = 5.4, lineheight = 0.85), panel.spacing = unit(2.5, "mm"))

  t <- T4[match(c(FOCAL, "all"), level)]
  t[, lab := factor(level, levels = rev(c(FOCAL, "all")))]
  dz <- data.table(xmin = c(-0.08, -NEGLIGIBLE, NEGLIGIBLE, MATERIAL), xmax = c(-NEGLIGIBLE, NEGLIGIBLE, MATERIAL, 0.08),
                   fill = ZONE[c("part", "tech", "part", "evo")])
  pb <- ggplot(t, aes(y = lab)) +
    geom_rect(data = dz, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = I(fill)), inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#555555") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = level), width = 0.25, linewidth = 0.5, orientation = "y") +
    geom_point(aes(x = delta, fill = level, shape = level == "all"), colour = "white", size = 1.7, stroke = 0.3) +
    scale_colour_manual(values = c(FOCAL_COL, all = INK), guide = "none") +
    scale_fill_manual(values = c(FOCAL_COL, all = INK), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21), guide = "none") +
    scale_y_discrete(labels = function(v) parse(text = ifelse(v == "all", "'all 33 pairs'", paste0("italic('", SHORT[v], "')")))) +
    scale_x_continuous(limits = c(-0.08, 0.08), breaks = c(-0.05, 0, 0.05), labels = function(x) minus(x, 2), expand = c(0, 0)) +
    annotate("text", x = c(-0.0525, 0.0125, 0.0375, 0.065), y = 4.75, label = c("partial", "study", "partial", "species"),
             size = 5 / .pt, family = FONT, colour = "#555555") +
    annotate("text", x = 0.077, y = 0.55, hjust = 1, vjust = 0, size = 5.2 / .pt, family = FONT, lineheight = 0.95,
             label = sprintf("all pairs Δ = %s\n(%s to %s)", minus(t[level == "all", delta]), minus(t[level == "all", lo]),
                             minus(t[level == "all", hi]))) +
    labs(x = expression(Delta == rho[species] - rho[study]), y = NULL) + theme_nat() +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
    coord_cartesian(ylim = c(0.5, 5), clip = "off")

  fwrite(P1A, file.path(OUT, "source_data/Fig3c_distance.csv"))
  LEV <- c(d_within = "within one study", d_study = "between two studies", d_species = "between two species")
  dd <- melt(P1A, id.vars = c("focal", "other"), measure.vars = names(LEV), variable.name = "lev", value.name = "dist")
  ds <- dd[, .(med = median(dist), lo = min(dist), hi = max(dist)), by = .(focal, lev)]
  ds[, lev := factor(LEV[as.character(lev)], levels = LEV)]
  dd[, lev := factor(LEV[as.character(lev)], levels = LEV)]
  Fm <- P1A[, .(F = mean(F)), by = focal]
  ds[, focal := factor(focal, levels = FOCAL)]; dd[, focal := factor(focal, levels = FOCAL)]
  Fm[, focal := factor(focal, levels = FOCAL)]
  LCOL <- c("within one study" = "#9E9E9E", "between two studies" = STUDY, "between two species" = ACC)
  set.seed(4)
  pc <- ggplot(ds, aes(lev, med)) +
    geom_jitter(data = dd, aes(lev, dist, colour = lev), width = 0.12, height = 0, size = 0.55, alpha = 0.45) +
    geom_line(aes(group = focal), linewidth = 0.35, colour = "#777777") +
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = lev), width = 0, linewidth = 0.6) +
    geom_point(aes(fill = lev), shape = 21, colour = "white", size = 2.3, stroke = 0.35) +
    geom_text(data = Fm, aes(x = 2, y = 0.84, label = sprintf("study share F = %.2f", F)), inherit.aes = FALSE,
              size = 5.4 / .pt, family = FONT) +
    facet_wrap(~focal, nrow = 1, labeller = as_labeller(function(v) LATIN[v])) +
    scale_colour_manual(values = LCOL, guide = "none") + scale_fill_manual(values = LCOL, name = NULL) +
    scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
    scale_x_discrete(labels = c("within\nstudy", "between\nstudies", "between\nspecies")) +
    labs(x = NULL, y = "Expression distance\n(median over orthogroups, log2)") + theme_nat() +
    theme(axis.text.x = element_text(size = 5.4, lineheight = 0.85), panel.spacing = unit(4, "mm"),
          legend.position = "none")
  fa <- T5["1A_F"]
  pc_note <- ggplot() + theme_void(base_family = FONT) +
    annotate("text", x = 0, y = 0.62, hjust = 0, vjust = 1, size = 5.6 / .pt, family = FONT, lineheight = 1,
             label = sprintf("Study share F: the fraction of the\nbetween-species excess over replicate\nnoise that a change of study reproduces.\n\nMean of the three species: %.2f\n(orthogroup 95%% interval %.2f to %.2f)",
                             fa$estimate, fa$lo, fa$hi)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.7), clip = "off")
  top <- (pa + labs(tag = "a") | pb + labs(tag = "b")) + plot_layout(widths = c(2.6, 1))
  bot <- (pc + labs(tag = "c") | pc_note) + plot_layout(widths = c(2.6, 1))
  p <- (top / bot) + plot_layout(heights = c(1, 1)) & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(p, "Fig3", 118)
}

# =================================================================== Fig. 4
fig4 <- function() {
  w <- merge(P1B[level == "within"], P1C, by.x = "label", by.y = "species")
  w <- w[match(intersect(TREE_ORDER, label), label)]
  fwrite(w, file.path(OUT, "source_data/Fig4_partial.csv"))
  w[, sp := factor(label, levels = rev(label))]
  long <- rbind(w[, .(sp, label, kind = "raw", v = rho_raw)], w[, .(sp, label, kind = "| expression level", v = rho_partial_expr)],
                w[, .(sp, label, kind = "| expression level + read ambiguity", v = rho_partial_amb_expr)])
  KC <- c(raw = "#BDBDBD", "| expression level" = STUDY, "| expression level + read ambiguity" = ACC)
  long[, kind := factor(kind, levels = names(KC))]
  rng <- long[, .(lo = min(v), hi = max(v)), by = sp]
  b1 <- T5["1B_partial_expr"]; c1 <- T5["1C_partial_amb_expr"]
  pa <- ggplot(long, aes(y = sp)) +
    annotate("rect", xmin = -NEGLIGIBLE, xmax = NEGLIGIBLE, ymin = -Inf, ymax = Inf, fill = ZONE[["tech"]]) +
    annotate("text", x = 0, y = 11.75, label = "negligible", hjust = 0.5, size = 5 / .pt, family = FONT,
             colour = "#555555") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#555555") +
    geom_segment(data = rng, aes(x = lo, xend = hi, yend = sp), linewidth = 0.5, colour = "#D5D5D5") +
    geom_point(aes(x = v, fill = kind), shape = 21, colour = "white", size = 2.0, stroke = 0.3) +
    scale_fill_manual(values = KC, name = NULL) +
    scale_y_discrete(labels = function(v) italic_sp(v)) +
    scale_x_continuous(limits = c(-0.05, 0.31), breaks = c(0, 0.1, 0.2, 0.3), expand = c(0, 0)) +
    annotate("text", x = 0.307, y = 0.75, hjust = 1, vjust = 0, size = 5.2 / .pt, family = FONT, lineheight = 0.95,
             label = sprintf("median partial ρ:\n| expression %s (%s to %s)\n| + ambiguity %s (%s to %s)",
                             minus(b1$estimate), minus(b1$lo), minus(b1$hi), minus(c1$estimate), minus(c1$lo), minus(c1$hi))) +
    labs(x = "Copy-number gradient within one study (Spearman ρ, partial)", y = NULL) + theme_nat() +
    theme(axis.ticks.y = element_blank(), legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, -3, 0)) +
    coord_cartesian(ylim = c(0.5, 12), clip = "off")
  pb <- ggplot(w, aes(y = sp)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#555555") +
    geom_segment(aes(x = 0, xend = rho_copies_ambiguity, yend = sp, colour = CLADE[label]), linewidth = 0.55) +
    geom_point(aes(x = rho_copies_ambiguity, colour = CLADE[label]), size = 1.6) +
    scale_colour_manual(values = CL_COL, guide = "none") +
    scale_y_discrete(labels = function(v) italic_sp(v)) +
    scale_x_continuous(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3), expand = c(0, 0)) +
    labs(x = "ρ (copy number, read ambiguity)", y = NULL) + theme_nat() + theme(axis.ticks.y = element_blank())
  p <- (pa + labs(tag = "a") | pb + labs(tag = "b")) + plot_layout(widths = c(1.6, 1)) &
    theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(p, "Fig4", 74)
}

# =================================================================== Fig. 5
fig5 <- function() {
  m <- as.matrix(MRAT[, ..TREE_ORDER]); rownames(m) <- MRAT$OG
  m <- m[rowSums(!is.na(m)) >= 10, ]
  m <- m[order(rowMeans(m, na.rm = TRUE)), ]
  fwrite(data.table(OG = rownames(m), m), file.path(OUT, "source_data/Fig5a_mean_ratios.csv"))
  hm <- data.table(row = rep(seq_len(nrow(m)), ncol(m)), sp = rep(colnames(m), each = nrow(m)), v = as.vector(m))
  hm[, sp := factor(sp, levels = TREE_ORDER)]
  hm[, vc := pmax(pmin(v, 4), -4)]
  pa <- ggplot(hm, aes(sp, row, fill = vc)) + geom_raster() +
    scale_fill_gradientn(colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"), limits = c(-4, 4),
                         na.value = "#E3E3E3", name = "log2 (copy A /\ncopy B)", breaks = c(-4, -2, 0, 2, 4),
                         labels = function(x) gsub("-", "−", x)) +
    scale_x_discrete(labels = function(v) italic_sp(v), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = sprintf("%d putative ohnolog pairs (ordered by mean ratio)", nrow(m))) + theme_nat() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.6, colour = CL_COL[CLADE[TREE_ORDER]]),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          legend.key.height = unit(5, "mm"), legend.key.width = unit(1.8, "mm"), legend.title = element_text(size = 5.6))

  M <- rbind(P2B[, .(a = sp_a, b = sp_b, rho = rho_r)], P2B[, .(a = sp_b, b = sp_a, rho = rho_r)])
  M[, ia := match(a, TREE_ORDER)][, ib := match(b, TREE_ORDER)]
  M <- M[ia > ib]
  fwrite(M, file.path(OUT, "source_data/Fig5b_conservation_matrix.csv"))
  M[, `:=`(a = factor(a, levels = TREE_ORDER), b = factor(b, levels = TREE_ORDER))]
  sep <- c(1.5, 2.5, 3.5, 7.5)
  pb <- ggplot(M, aes(b, a, fill = rho)) + geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", rho), colour = rho > 0.58), size = 4.6 / .pt, family = FONT) +
    scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = "white"), guide = "none") +
    scale_fill_viridis_c(limits = c(0.3, 0.85), name = "Conservation\n(Spearman ρ)", breaks = c(0.3, 0.5, 0.7)) +
    scale_x_discrete(labels = function(v) italic_sp(v), limits = TREE_ORDER[-12]) +
    scale_y_discrete(labels = function(v) italic_sp(v), limits = rev(TREE_ORDER[-1])) +
    labs(x = NULL, y = NULL) + theme_nat() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.6, colour = CL_COL[CLADE[TREE_ORDER[-12]]]),
          axis.text.y = element_text(size = 5.6, colour = rev(CL_COL[CLADE[TREE_ORDER[-1]]])),
          axis.line = element_blank(), axis.ticks = element_blank(), legend.key.height = unit(5, "mm"),
          legend.key.width = unit(1.8, "mm"), legend.title = element_text(size = 5.6),
          legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1)) +
    coord_fixed()

  fwrite(P2B, file.path(OUT, "source_data/Fig5c_conservation_time.csv")); fwrite(CE, file.path(OUT, "source_data/Fig5c_ceilings.csv"))
  rep_med <- median(CE[ceiling == "replicate", rho]); labr <- range(CE[ceiling == "laboratory", rho])
  P2B[, caj := sp_a == "Ccaj" | sp_b == "Ccaj"]
  P2B[, kind := fifelse(caj, "Cajanus pairs (excluded, D1)", fifelse(CLADE[sp_a] == CLADE[sp_b], "same clade", "different clades"))]
  pp <- P2B[caj == FALSE]; fit <- lm(rho_r ~ time, pp); rng <- range(pp$time)
  tr <- T5["2B_trend"]
  KIND <- c("different clades", "same clade", "Cajanus pairs (excluded, D1)")
  KLAB <- c("different clades", "same clade", expression(italic("Cajanus") ~ "pairs (excluded, D1)"))
  NL <- fread(file.path(RES, "exploratory_review/R7_random_pair_null_lines.tsv"))
  NS <- fread(file.path(RES, "exploratory_review/R7_random_pair_null_summary.tsv"))[statistic == "change_over_range"]
  fwrite(NL, file.path(OUT, "source_data/Fig5c_null_lines.csv"))
  pc <- ggplot(P2B, aes(time, rho_r)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = labr[1], ymax = labr[2], fill = "#F6E7C1") +
    geom_ribbon(data = NL, aes(x = time, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "#BDBDBD", alpha = 0.45) +
    geom_line(data = NL, aes(x = time, y = median), inherit.aes = FALSE, colour = "#6E6E6E", linewidth = 0.45,
              linetype = "42") +
    annotate("text", x = 40, y = NL[which.min(abs(time - 40)), hi] + 0.008, label = "matched non-paralogous pairs",
             hjust = 0.5, vjust = 0, size = 5 / .pt, family = FONT, colour = "#555555") +
    geom_hline(yintercept = rep_med, linewidth = 0.3, linetype = "22", colour = "#555555") +
    annotate("text", x = 57.5, y = rep_med - 0.012, label = "replicate halves of one study (median)", hjust = 1,
             vjust = 1, size = 5 / .pt, family = FONT, colour = "#555555") +
    annotate("text", x = 57.5, y = labr[2] - 0.012, label = "two studies of one species (range)", hjust = 1, vjust = 1,
             size = 5 / .pt, family = FONT, colour = "#8A6D1F") +
    annotate("segment", x = rng[1], xend = rng[2], y = predict(fit, data.frame(time = rng[1])),
             yend = predict(fit, data.frame(time = rng[2])), linewidth = 0.5, colour = INK) +
    geom_point(aes(shape = kind, fill = kind, colour = kind), size = 1.35, stroke = 0.35) +
    scale_shape_manual(values = c(21, 24, 21) |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    scale_fill_manual(values = c(ACC, ACC, "white") |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    scale_colour_manual(values = c("white", "white", GREY) |> setNames(KIND), breaks = KIND, labels = KLAB, name = NULL) +
    annotate("text", x = 9, y = 0.215, hjust = 0, vjust = 0, size = 5.4 / .pt, family = FONT, lineheight = 0.95,
             label = sprintf("ohnologs, change over %.1f Myr: %s (%s to %s)\nmatched null: %s (%s to %s)", diff(rng),
                             minus(tr$estimate), minus(tr$lo), minus(tr$hi), minus(NS$null_median), minus(NS$null_lo),
                             minus(NS$null_hi))) +
    scale_y_continuous(limits = c(0.2, 1.0), breaks = seq(0.2, 1, 0.2), expand = c(0, 0)) +
    scale_x_continuous(limits = c(8, 58), breaks = c(10, 20, 30, 40, 50), expand = c(0, 0)) +
    labs(x = "Divergence time (Myr)", y = "Conservation of ohnolog asymmetry\n(ρ between species)") + theme_nat() +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, -3, 0),
          legend.text = element_text(size = 5.6))

  fwrite(P2A, file.path(OUT, "source_data/Fig5d_study_excess.csv"))
  ex <- rbind(P2A[, .(species, kind = "absolute expression", v = X_abs)], P2A[, .(species, kind = "within-sample ratio", v = X_ratio)])
  ex[, species := factor(species, levels = FOCAL)]
  ex[, kind := factor(kind, levels = c("absolute expression", "within-sample ratio"))]
  e <- T5["2A_E"]
  pd <- ggplot(ex, aes(species, v, fill = kind)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    geom_text(aes(label = sprintf("%.2f", v)), position = position_dodge(width = 0.72), vjust = -0.4, size = 5 / .pt,
              family = FONT) +
    annotate("text", x = 2, y = 0.76, label = sprintf("E = %.2f (95%% interval %.2f to %.2f)", e$estimate, e$lo, e$hi),
             size = 5.6 / .pt, family = FONT) +
    scale_fill_manual(values = c("absolute expression" = "#BDBDBD", "within-sample ratio" = ACC), name = NULL) +
    scale_x_discrete(labels = function(v) italic_sp(v)) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
    labs(x = NULL, y = "Study excess\n(between − within studies)") + theme_nat() +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, -3, 0))

  p <- ((pa + labs(tag = "a") | pb + labs(tag = "b")) + plot_layout(widths = c(1, 1.15))) /
    ((pc + labs(tag = "c") | pd + labs(tag = "d")) + plot_layout(widths = c(1.35, 1))) +
    plot_layout(heights = c(1.15, 1)) & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  save_fig(p, "Fig5", 160)
}

# =================================================================== Fig. 6
fig6 <- function() {
  mk <- function(k, lab, key) {
    d <- D3[[k]]; fwrite(d, file.path(OUT, sprintf("source_data/Fig6%s_dosage_%s.csv", if (k == "study1") "a" else "b", k)))
    d[, grp := fifelse(set == "glycine_wgd", sprintf("Glycine WGD families (%d)", sum(set == "glycine_wgd")),
                       sprintf("single-copy in all 12 (%d)", sum(set == "calibration")))]
    c0 <- median(d[set == "calibration", D]); g0 <- median(d[set == "glycine_wgd", D])
    r <- T5[key]
    vl <- data.table(x = c(c0, c0 + 1, g0), what = c("median, single-copy families", "full dosage (single-copy median + 1)",
                                                      "median, duplicated families"))
    ggplot(d, aes(D)) +
      geom_histogram(data = d[set == "calibration"], aes(y = after_stat(density)), binwidth = 0.2, boundary = 0,
                     fill = "#C9C9C9", colour = NA) +
      geom_histogram(data = d[set == "glycine_wgd"], aes(y = after_stat(density)), binwidth = 0.2, boundary = 0,
                     fill = ACC, alpha = 0.28, colour = NA) +
      geom_step(data = d[set == "glycine_wgd"], aes(y = after_stat(density)), stat = "bin", binwidth = 0.2, boundary = 0,
                colour = ACC, linewidth = 0.35, direction = "mid") +
      geom_vline(data = vl, aes(xintercept = x, linetype = what, colour = what), linewidth = 0.45) +
      scale_linetype_manual(values = c("median, single-copy families" = "solid", "full dosage (single-copy median + 1)" = "22",
                                       "median, duplicated families" = "solid"), name = NULL) +
      scale_colour_manual(values = c("median, single-copy families" = "#555555", "full dosage (single-copy median + 1)" = EVO,
                                     "median, duplicated families" = ACC), name = NULL) +
      labs(title = sprintf("%s: relative dosage %s", lab, minus(r$estimate)),
           subtitle = sprintf("95%% interval %s to %s", minus(r$lo), minus(r$hi))) +
      scale_x_continuous(limits = c(-3, 4), breaks = -3:4, labels = function(x) gsub("-", "−", x), expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = expression("Relative family expression of " * italic("G. max") * " (log2)"),
           y = "Density") + theme_nat() +
      theme(plot.title = element_text(size = 6, hjust = 0, face = "bold", margin = margin(b = 0.5)),
            plot.subtitle = element_text(size = 5.4, hjust = 0, margin = margin(b = 1.5)))
  }
  pa <- mk("study1", "Study 1", "3_dosage_Gmax"); pb <- mk("study2", "Study 2", "3_dosage_Gmax_study2")
  # c: relative dosage by quartile of reference expression (exploratory), against full dosage under log2(x + 1)
  q <- fread(file.path(RES, "exploratory_review/R1_dosage_by_expression.tsv"))
  fwrite(q, file.path(OUT, "source_data/Fig6c_dosage_by_quartile.csv"))
  full <- unique(q[, .(bin, full_dosage_expected)])
  dodge <- position_dodge(width = 0.45)
  pc <- ggplot(q, aes(bin, rel_dosage, colour = study)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#777777") +
    geom_line(data = full, aes(bin, full_dosage_expected), inherit.aes = FALSE, colour = EVO, linewidth = 0.45,
              linetype = "22") +
    geom_point(data = full, aes(bin, full_dosage_expected), inherit.aes = FALSE, colour = EVO, shape = 4, size = 1.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45, position = dodge) +
    geom_point(aes(fill = study), shape = 21, colour = "white", size = 2.0, stroke = 0.3, position = dodge) +
    annotate("text", x = 0.75, y = 1.12, hjust = 0, size = 5 / .pt, family = FONT, colour = EVO,
             label = "full dosage (×, dashed)") +
    annotate("text", x = 4.35, y = 0.05, hjust = 1, vjust = 0, size = 5 / .pt, family = FONT, colour = "#555555",
             label = "complete compensation") +
    scale_colour_manual(values = c("study 1" = ACC, "study 2" = STUDY), guide = "none") +
    scale_fill_manual(values = c("study 1" = ACC, "study 2" = STUDY), name = NULL) +
    scale_x_continuous(breaks = 1:4, labels = c("Q1\n<3.1", "Q2\n3.1–9.8", "Q3\n9.8–22.6", "Q4\n>22.6"),
                       limits = c(0.6, 4.4)) +
    scale_y_continuous(limits = c(-0.05, 1.2), breaks = seq(0, 1, 0.25)) +
    labs(x = "Reference expression quartile (scaled TPM)", y = "Relative dosage (log2)") + theme_nat() +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.16), legend.justification = c(1, 0),
          axis.text.x = element_text(size = 5.2, lineheight = 0.85))
  leg <- ggplot() + theme_void()
  p <- (pa + labs(tag = "a") + theme(legend.position = "none") |
          pb + labs(tag = "b") + theme(axis.title.y = element_blank()) |
          pc + labs(tag = "c")) + plot_layout(widths = c(1, 1, 0.82))
  p <- p & theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
  p[[2]] <- p[[2]] + theme(legend.position = "bottom", legend.key.width = unit(4, "mm"), legend.text = element_text(size = 5.4),
                           legend.margin = margin(-2, 0, 0, 0)) + guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
  save_fig(p, "Fig6", 80)

}

args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else paste0("fig", 1:6)
for (f in todo) get(f)()
