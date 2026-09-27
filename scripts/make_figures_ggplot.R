#!/usr/bin/env Rscript
# Manuscript Figures 2-5 in ggplot2 (Figure 1, the study flowchart, stays in matplotlib:
# scripts/make_study_flowchart.py). Colours are ggplot2's own palettes: the default hue scale for the
# five sufficiency regimes (shared by Figures 3 and 4b so one key serves both) and for the two
# signal-to-noise classes of Figure 5, and the built-in viridis scale for the three modalities of
# Figure 4a, so that modality colours never coincide with regime colours.
#
# Inputs (all committed): fixtures/fig3_hist.json, fixtures/fig4_kde.json, fixtures/unified_spectrum.json,
# fixtures/tahoe_direct_curves_full.json, fixtures/emeraldbay_falsification_full.json,
# fixtures/orion_{HCT116,HEK293T}_heldout.json, fixtures/trade_{jurkat,hepg2}_falsification.json.
# Figure 2 is a synthetic schematic (fixed seed). Numbers printed on the figures come from the fixtures.
#
#   Rscript scripts/make_figures_ggplot.R            # writes figures/fig{2,3,4,5}_*.{pdf,png}
# Needs R >= 4.1 with ggplot2, jsonlite, dplyr, scales, patchwork, ggridges.

suppressPackageStartupMessages({
  library(ggplot2); library(jsonlite); library(dplyr); library(scales); library(patchwork); library(ggridges)
})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), ".."))
FX <- file.path(root, "fixtures"); FIG <- file.path(root, "figures")
dir.create(FIG, showWarnings = FALSE)

BASE <- 8.5   # base font size (pt) for print at 6.3 in width
theme_paper <- function(base = BASE) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
          axis.line.x = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks.x = element_line(linewidth = 0.3, colour = "grey30"),
          axis.title = element_text(size = base),
          axis.text = element_text(size = base - 1, colour = "grey20"),
          plot.title = element_text(size = base + 0.5, face = "plain", hjust = 0, margin = margin(b = 4)),
          plot.title.position = "plot",
          legend.title = element_text(size = base - 0.5),
          legend.text = element_text(size = base - 1),
          legend.key.size = unit(9, "pt"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(4, 6, 4, 4))
}
save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  cat("  wrote figures/", name, ".{pdf,png}\n", sep = "")
}

# regime factor: one order everywhere, so the default hue palette gives the same colour in every figure
REGIMES <- c("over-sampled", "treated-depth-limited", "ghost", "control-pool-limited", "not detectable")
REGIME_KEY <- c("over-sampled" = "over-sampled\nn* within the acquired depth",
                "treated-depth-limited" = "treated-depth-limited\nn* above the acquired depth",
                "ghost" = "ghost\nn* > 50,000 treated cells",
                "control-pool-limited" = "control-pool-limited\nn* = ∞, control pool too small",
                "not detectable" = "not detectable\neffect below the sampling floor")
SCREENS <- c("Tahoe-100M", "EmeraldBay", "Orion HCT116", "Orion HEK293T", "TRADE Jurkat", "TRADE HepG2")

# ------------------------------------------------------------------------------------------------
# Figure 2: geometry of the angular error (schematic)
# ------------------------------------------------------------------------------------------------
fig2 <- function() {
  set.seed(7)
  v <- c(3, 0); NOISE <- 2.6
  panel <- function(n, letter, decompose) {
    sd <- NOISE / sqrt(n)
    est <- data.frame(x = v[1] + rnorm(60, 0, sd), y = v[2] + rnorm(60, 0, sd))
    circ <- data.frame(t = seq(0, 2 * pi, length.out = 200)) |>
      transmute(x = v[1] + sd * cos(t), y = v[2] + sd * sin(t))
    p <- ggplot() +
      geom_polygon(data = circ, aes(x, y), fill = "grey85", alpha = 0.5, colour = "grey60", linewidth = 0.3, linetype = "22") +
      geom_point(data = est, aes(x, y), colour = "grey55", size = 0.7, alpha = 0.7) +
      annotate("segment", x = 0, y = 0, xend = v[1], yend = v[2], linewidth = 1.1,
               arrow = arrow(length = unit(5, "pt"), type = "closed"), colour = "black") +
      annotate("point", x = 0, y = 0, size = 1.2) +
      annotate("text", x = 1.15, y = -0.25, label = "italic(v)*'  true effect'", parse = TRUE, size = 2.8, hjust = 0.5, vjust = 1) +
      coord_equal(xlim = c(-0.3, 6.2), ylim = c(-1.75, 1.6), expand = FALSE, clip = "off") +
      labs(title = sprintf("(%s)  %s cells per arm,  n = %d", letter, if (n == 12) "few" else "many", n)) +
      theme_void(base_size = BASE) +
      theme(plot.title = element_text(size = BASE + 0.5, hjust = 0, margin = margin(b = 2)),
            plot.title.position = "plot", legend.position = "none",
            plot.background = element_rect(fill = "white", colour = NA))
    if (decompose) {
      vh <- c(3.5, 0.7); u <- v / sqrt(sum(v^2)); e <- vh - v
      along <- sum(e * u) * u; foot <- v + along
      comp <- data.frame(x = c(v[1], foot[1]), y = c(v[2], foot[2]), xend = c(foot[1], vh[1]), yend = c(foot[2], vh[2]),
                         part = factor(c("along-signal: changes length only", "across-signal: rotates the direction"),
                                       levels = c("along-signal: changes length only", "across-signal: rotates the direction")))
      th <- atan2(vh[2], vh[1]); arc <- data.frame(t = seq(0, th, length.out = 40)) |> transmute(x = 1.2 * cos(t), y = 1.2 * sin(t))
      s <- 0.13; ra <- data.frame(x = c(foot[1] - s, foot[1] - s, foot[1]), y = c(foot[2], foot[2] + s, foot[2] + s))
      p <- p +
        annotate("segment", x = 0, y = 0, xend = vh[1], yend = vh[2], linewidth = 0.6, colour = "grey45",
                 arrow = arrow(length = unit(4, "pt"), type = "closed")) +
        geom_segment(data = comp, aes(x, y, xend = xend, yend = yend, colour = part), linewidth = 1.4, lineend = "round") +
        geom_path(data = ra, aes(x, y), colour = "grey40", linewidth = 0.3) +
        geom_path(data = arc, aes(x, y), colour = "grey40", linewidth = 0.4) +
        annotate("point", x = vh[1], y = vh[2], size = 1.4) +
        annotate("text", x = vh[1] - 0.12, y = vh[2] + 0.12, label = "hat(italic(v))*'  estimate'", parse = TRUE, size = 2.8, hjust = 1, vjust = 0) +
        annotate("text", x = 1.45, y = 0.2, label = "theta", parse = TRUE, size = 3.6) +
        annotate("text", x = v[1] - 0.95, y = -sd - 0.12, label = "'RMS scatter' %prop% 1/sqrt(italic(n))", parse = TRUE, size = 2.6, colour = "grey35", vjust = 1) +
        annotate("segment", x = foot[1] + 0.06, y = foot[2] + along[2] * 0 + e[2] / 2, xend = 4.2, yend = 0.85, linewidth = 0.3, colour = "grey50") +
        annotate("text", x = 4.25, y = 0.85, label = "across-signal\nrotates the direction", size = 2.6, hjust = 0, vjust = 0.5, colour = "grey20") +
        annotate("segment", x = v[1] + along[1] / 2, y = -0.04, xend = 5.05, yend = -0.7, linewidth = 0.3, colour = "grey50") +
        annotate("text", x = 5.05, y = -0.75, label = "along-signal\nchanges length only", size = 2.6, hjust = 0.5, vjust = 1, colour = "grey20")
    } else {
      p <- p +
        annotate("text", x = v[1], y = sd + 0.3, label = "'RMS scatter' %prop% 1/sqrt(italic(n))", parse = TRUE, size = 2.6, colour = "grey35", vjust = 0) +
        annotate("text", x = v[1] + 0.1, y = -0.55, label = "hat(italic(v))*' locks onto '*italic(v)*', '*theta%->%0", parse = TRUE, size = 2.8, vjust = 1)
    }
    p
  }
  p <- panel(12, "a", TRUE) + panel(200, "b", FALSE) + plot_layout(ncol = 2)
  save_fig(p, "fig2_geometry", 7.4, 2.7)
}

# ------------------------------------------------------------------------------------------------
# Figure 3: Tahoe-100M sufficiency spectrum
# ------------------------------------------------------------------------------------------------
fig3 <- function() {
  H <- fromJSON(file.path(FX, "fig3_hist.json"))
  pct <- H$pct
  strip <- data.frame(regime = factor(REGIMES, levels = REGIMES),
                      pct = c(pct$OVER, pct$UNDER, pct$GHOST, pct$`POOL-LIMITED`, pct$`NOT-DETECTABLE`)) |>
    mutate(xmax = cumsum(pct), xmin = xmax - pct, xmid = (xmin + xmax) / 2,
           lab = sprintf("%.1f%%", pct), inside = pct >= 5)
  pa <- ggplot(strip) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = regime), colour = "white", linewidth = 0.6) +
    geom_text(data = filter(strip, inside), aes(xmid, 0.5, label = lab), size = 2.8, fontface = "bold") +
    geom_text(data = filter(strip, !inside), aes(xmid, 1.25, label = lab), size = 2.6) +
    geom_segment(data = filter(strip, !inside), aes(x = xmid, xend = xmid, y = 1.02, yend = 1.15), linewidth = 0.3, colour = "grey50") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.55), expand = c(0, 0)) +
    scale_fill_hue(name = NULL, labels = REGIME_KEY, drop = FALSE) +
    labs(title = sprintf("(a)  All %s conditions by regime", format(H$n_all, big.mark = ","))) +
    theme_void(base_size = BASE) +
    theme(plot.title = element_text(size = BASE + 0.5, hjust = 0, margin = margin(b = 3)), plot.title.position = "plot",
          legend.position = "right", legend.text = element_text(size = BASE - 1.5, lineheight = 0.9),
          legend.key.size = unit(9, "pt"), legend.key.spacing.y = unit(5, "pt"),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 6, 2, 4))

  edges <- H$bin_edges; nb <- length(edges) - 1
  hist <- bind_rows(lapply(c("OVER", "UNDER", "GHOST"), function(k)
    data.frame(xmin = edges[-length(edges)], xmax = edges[-1], count = H$counts[[k]],
               regime = c(OVER = "over-sampled", UNDER = "treated-depth-limited", GHOST = "ghost")[[k]]))) |>
    mutate(regime = factor(regime, levels = REGIMES)) |>
    filter(xmin < 2e6) |>
    group_by(xmin) |> arrange(regime, .by_group = TRUE) |>
    mutate(ymax = cumsum(count), ymin = ymax - count) |> ungroup()
  refs <- data.frame(x = c(H$N0, H$median_n_star, H$ghost), lt = c("22", "solid", "11"),
                     lab = c(sprintf("'median depth '*italic(N)[0]*' = %s'", format(H$N0, big.mark = ",")),
                             sprintf("median ~ italic(n)*'*' == '%s'", format(H$median_n_star, big.mark = ",")),
                             sprintf("ghost ~ italic(n)*'*' >= '%s'", format(H$ghost, big.mark = ","))),
                     hj = c(1, 0, 0), yf = c(1.04, 1.14, 1.04))
  ymax <- max(tapply(hist$count, hist$xmin, sum)) * 1.02
  pb <- ggplot(hist) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = regime), colour = NA) +
    geom_vline(data = refs, aes(xintercept = x, linetype = lt), linewidth = 0.5, colour = "grey15", show.legend = FALSE) +
    geom_text(data = refs, aes(x = x, y = ymax * yf, label = lab, hjust = hj), parse = TRUE, size = 2.6) +
    scale_linetype_identity() +
    scale_x_log10(limits = c(20, 2e6), breaks = 10^(2:6), labels = trans_format("log10", math_format(10^.x)), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.24)), labels = comma) +
    scale_fill_hue(drop = FALSE, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(title = sprintf("(b)  Finite quotas of the %.1f%% detectable conditions", H$pct_detectable),
         x = "required treated cells per condition, two-arm quota n*", y = "number of conditions") +
    theme_paper() + theme(plot.margin = margin(4, 6, 4, 4))
  p <- (pa / pb) + plot_layout(heights = c(1, 5.2), guides = "collect") &
    theme(legend.position = "right", legend.justification = c(0, 1), legend.box.margin = margin(t = 22))
  save_fig(p, "fig3_tahoe_spectrum", 7.0, 4.4)
}

# ------------------------------------------------------------------------------------------------
# Figure 4: cross-modality summary
# ------------------------------------------------------------------------------------------------
fig4 <- function() {
  K <- fromJSON(file.path(FX, "fig4_kde.json"), simplifyVector = FALSE)
  S <- fromJSON(file.path(FX, "unified_spectrum.json"), simplifyVector = FALSE)$headline
  S <- setNames(S, sapply(S, `[[`, "name"))
  grid <- unlist(K$grid_log10m)
  MODLAB <- c(chemical = "chemical, small-molecule (Vevo Mosaic)", `genome-wide` = "genome-wide CRISPRi (X-Atlas/Orion)",
              `essential-gene` = "essential-gene CRISPRi (TRADE)")
  ridge <- bind_rows(lapply(K$screens, function(sc) {
    d <- unlist(sc$density); data.frame(screen = sc$name, modality = MODLAB[[sc$modality]], x = 10^grid, h = d / max(d) * 0.85)
  })) |> mutate(screen = factor(screen, levels = rev(SCREENS)), modality = factor(modality, levels = MODLAB))
  meds <- data.frame(screen = factor(sapply(K$screens, `[[`, "name"), levels = rev(SCREENS)),
                     med = sapply(K$screens, `[[`, "median_m"))
  pa <- ggplot(ridge, aes(x = x, y = screen)) +
    geom_ridgeline(aes(height = h, fill = modality, colour = modality), alpha = 0.35, linewidth = 0.5, scale = 1) +
    geom_segment(data = meds, aes(x = med, xend = med, y = as.numeric(screen), yend = as.numeric(screen) + 0.85),
                 inherit.aes = FALSE, linewidth = 0.7, colour = "black") +
    geom_text(data = meds, aes(x = med, y = as.numeric(screen) + 0.85, label = sprintf("%.2f", med)),
              inherit.aes = FALSE, hjust = -0.15, vjust = 1, size = 2.5) +
    geom_hline(yintercept = c(2.5, 4.5) - 0 + 0.0, colour = "grey85", linewidth = 0.3) +
    scale_x_log10(limits = c(0.03, 30), breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10"), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.15, 1.05))) +
    scale_fill_viridis_d(name = NULL, end = 0.85) + scale_colour_viridis_d(name = NULL, end = 0.85) +
    labs(title = "(a)  Effect magnitude per screen", x = "bias-corrected effect magnitude m (log scale)", y = NULL) +
    theme_paper() +
    theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey88"),
          legend.position = "bottom", legend.direction = "vertical", legend.justification = "left",
          legend.key.spacing.y = unit(1, "pt"), legend.margin = margin(t = 0))

  keys <- c("pct_over", "pct_under", "pct_ghost", "pct_pool_limited", "not_detectable_pct")
  bars <- bind_rows(lapply(SCREENS, function(nm) data.frame(screen = nm, regime = REGIMES, pct = sapply(keys, function(k) S[[nm]][[k]])))) |>
    mutate(screen = factor(screen, levels = SCREENS), regime = factor(regime, levels = REGIMES)) |>
    group_by(screen) |> arrange(desc(regime), .by_group = TRUE) |>
    mutate(top = cumsum(pct), mid = top - pct / 2) |> ungroup()
  overs <- filter(bars, regime == "over-sampled")
  xlabs <- c("Tahoe-\n100M", "Emerald-\nBay", "Orion\nHCT116", "Orion\nHEK293T", "TRADE\nJurkat", "TRADE\nHepG2")
  pb <- ggplot(bars, aes(x = screen)) +
    geom_col(aes(y = pct, fill = regime), width = 0.62, colour = NA) +
    geom_text(data = filter(bars, pct >= 7), aes(y = mid, label = ifelse(regime == "over-sampled", sprintf("%.1f%%", pct), sprintf("%.0f%%", pct)),
                                                colour = regime == "not detectable"), size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    geom_text(data = overs, aes(y = 102, label = sprintf("%.1f%%", pct)), size = 2.6, fontface = "bold", vjust = 0) +
    annotate("text", x = 0.55, y = 104, label = "over-\nsampled", size = 2.2, hjust = 1, vjust = 0, lineheight = 0.85, colour = "grey30") +
    annotate("segment", x = 0.7, xend = 2.3, y = -13, yend = -13, colour = "grey55", linewidth = 0.4) +
    annotate("segment", x = 2.7, xend = 6.3, y = -13, yend = -13, colour = "grey55", linewidth = 0.4) +
    annotate("text", x = 1.5, y = -15, label = "chemical", size = 2.5, vjust = 1, colour = "grey40") +
    annotate("text", x = 4.5, y = -15, label = "genetic (CRISPRi)", size = 2.5, vjust = 1, colour = "grey40") +
    scale_x_discrete(labels = xlabs, expand = expansion(add = c(1.0, 0.5))) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    scale_fill_hue(name = NULL, labels = REGIME_KEY, drop = FALSE, guide = guide_legend(reverse = TRUE)) +
    coord_cartesian(ylim = c(0, 112), clip = "off") +
    labs(title = "(b)  Sufficiency regime per screen", x = NULL, y = "% of conditions") +
    theme_paper() +
    theme(axis.text.x = element_text(size = BASE - 1.5, lineheight = 0.85), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "right", legend.text = element_text(size = BASE - 1.5, lineheight = 0.9),
          legend.key.spacing.y = unit(6, "pt"), plot.margin = margin(4, 6, 18, 4))
  p <- pa + pb + plot_layout(widths = c(1, 1.15))
  save_fig(p, "fig4_crossmodality_summary", 8.0, 3.9)
}

# ------------------------------------------------------------------------------------------------
# Figure 5: held-out downsample-and-measure test, one panel per screen
# ------------------------------------------------------------------------------------------------
fig5 <- function() {
  ld <- function(f) { p <- file.path(FX, f); if (file.exists(p)) fromJSON(p) else NULL }
  rows_of <- function(kind) {
    r <- switch(kind,
      tahoe = ld("tahoe_direct_curves_full.json")$per_group,
      emb = ld("emeraldbay_falsification_full.json")$per_group_slope,
      HCT116 = ld("orion_HCT116_heldout.json")$per_group,
      HEK293T = ld("orion_HEK293T_heldout.json")$per_group,
      jurkat = ld("trade_jurkat_falsification.json")$rows,
      hepg2 = ld("trade_hepg2_falsification.json")$rows)
    if (is.null(r)) return(NULL)
    as.data.frame(r) |> filter(N >= 400) |> select(N, rho2, ratio)
  }
  kinds <- c(tahoe = "Tahoe-100M (chemical)", emb = "EmeraldBay (chemical)", HCT116 = "Orion HCT116 (genetic)",
             HEK293T = "Orion HEK293T (genetic)", jurkat = "TRADE Jurkat (genetic)", hepg2 = "TRADE HepG2 (genetic)")
  dat <- bind_rows(lapply(names(kinds), function(k) { d <- rows_of(k); if (is.null(d)) NULL else mutate(d, kind = k) }))
  dat <- dat |> mutate(x = 1 / pmax(rho2, 1e-3), regime = factor(ifelse(rho2 >= 3, "first-order regime (ρ² ≥ 3)", "below it (ρ² < 3)"),
                                                            levels = c("first-order regime (ρ² ≥ 3)", "below it (ρ² < 3)")))
  summ <- dat |> group_by(kind) |> summarise(
    n = n(), nin = sum(rho2 >= 3), medin = if (any(rho2 >= 3)) median(ratio[rho2 >= 3]) else NA_real_,
    maxrho = max(rho2), sp = suppressWarnings(cor(x, 1 - ratio, method = "spearman")), .groups = "drop") |>
    mutate(label = sprintf("%s, %s groups\n%s, Spearman %.2f", kinds[kind], format(n, big.mark = ",", trim = TRUE),
                           ifelse(is.na(medin), sprintf("all ρ² < 3 (max %.1f)", maxrho),
                                  sprintf("ρ² ≥ 3: median ratio %.2f (n = %s)", medin, format(nin, big.mark = ",", trim = TRUE))), sp))
  dat <- dat |> left_join(select(summ, kind, label), by = "kind") |>
    mutate(label = factor(label, levels = summ$label[match(names(kinds), summ$kind)]))
  p <- ggplot(dat, aes(x, ratio, colour = regime)) +
    geom_hline(yintercept = 1, linetype = "22", linewidth = 0.4, colour = "grey20") +
    geom_vline(xintercept = 1 / 3, linetype = "11", linewidth = 0.4, colour = "grey55") +
    geom_point(size = 0.55, alpha = 0.45, stroke = 0) +
    facet_wrap(~ label, ncol = 3) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_colour_hue(name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    coord_cartesian(ylim = c(-0.05, 1.5)) +
    labs(x = "1 / ρ²  (signal-to-noise decreases to the right)", y = "realized / predicted slope") +
    theme_paper() +
    theme(strip.text = element_text(size = BASE - 2, hjust = 0, lineheight = 1.05, margin = margin(b = 3)),
          panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
          legend.position = "bottom", panel.spacing = unit(10, "pt"))
  save_fig(p, "fig5_heldout_per_dataset", 8.6, 5.4)
}

cat("Generating ggplot2 figures ->", FIG, "\n")
fig2(); fig3(); fig4(); fig5()
cat("done.\n")
