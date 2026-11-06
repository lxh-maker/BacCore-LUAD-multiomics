suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

repo_dir <- if (dir.exists("source_data")) "." else ".."
src_dir <- file.path(repo_dir, "source_data")
fig_dir <- file.path(repo_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(blue = "#0072B2", sky = "#56B4E9", green = "#009E73",
         orange = "#E69F00", vermillion = "#D55E00", purple = "#CC79A7",
         yellow = "#F0E442", black = "#000000", grey = "#999999")

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(text = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
          plot.subtitle = element_text(size = base_size - 0.5, colour = "#333333"),
          legend.title = element_blank(), legend.key = element_blank(),
          plot.margin = margin(6, 7, 6, 7))
}

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
}

fmt_p <- function(p) ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p))

# Figure 1: study design and interpretive scope
blocks <- data.table(
  xmin = c(0.25, 2.15, 4.05, 5.95), xmax = c(1.85, 3.75, 5.65, 7.55),
  ymin = 3.2, ymax = 4.7,
  label = c("TCGA-LUAD\nBulk RNA-seq\nDiscovery",
            "GSE131907\nSingle-cell RNA-seq\nLocalization",
            "CPTAC-LUAD\nProteomics\nValidation",
            "Four cohorts\nOverall survival\nMeta-analysis"),
  fill = c(pal["blue"], pal["green"], pal["orange"], pal["purple"])
)
genes <- data.table(x = seq(0.8, 7.0, length.out = 6), y = 8.1,
                    gene = c("BAX", "BAK1", "VDAC1", "VDAC2", "PHB", "PHB2"))
p1 <- ggplot() +
  annotate("rect", xmin = 1.45, xmax = 6.35, ymin = 6.0, ymax = 7.15,
           fill = "#EAF2F8", colour = pal["blue"], linewidth = 0.7) +
  annotate("text", x = 3.9, y = 6.75, label = "Prespecified equal-weight BacCore score",
           fontface = "bold", size = 4.0, family = "sans") +
  annotate("text", x = 3.9, y = 6.32,
           label = "mean of gene-wise standardized abundance values", size = 3.1, family = "sans") +
  geom_tile(data = genes, aes(x, y), width = 0.92, height = 0.78,
            fill = "white", colour = pal["blue"], linewidth = 0.6) +
  geom_text(data = genes, aes(x, y, label = gene), fontface = "bold", size = 3.4) +
  geom_segment(data = genes, aes(x = x, y = 7.68, xend = 3.9, yend = 7.17),
               arrow = arrow(length = grid::unit(0.10, "inches")), linewidth = 0.35,
               colour = "#555555") +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = "white", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = blocks, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
            colour = "white", fontface = "bold", lineheight = 1.02, size = 2.9) +
  annotate("segment", x = 3.9, xend = 3.9, y = 5.95, yend = 4.83,
           arrow = arrow(length = grid::unit(0.10, "inches")), colour = "#555555") +
  annotate("text", x = 3.9, y = 2.3,
           label = "Integrated interpretation: a mitochondrial-stress-associated\ntranscriptional/protein state with modest prognostic information",
           fontface = "bold", size = 3.3, lineheight = 1.1) +
  annotate("text", x = 3.9, y = 0.49,
           label = "Interpretive boundary: machinery abundance does not prove complex assembly,\nmtDNA efflux, or Baciptosis execution",
           colour = pal["vermillion"], fontface = "bold", size = 2.8) +
  coord_cartesian(xlim = c(0, 7.8), ylim = c(0, 8.75), clip = "off") +
  theme_void(base_family = "sans") + theme(plot.margin = margin(8, 10, 8, 10))
save_figure(p1, "Figure_1_study_design", 7.1, 7.4)

# Figure 2: differential expression and WGCNA
de <- fread(file.path(src_dir, "Figure_2A_BacCore_DE.tsv"))
de[, gene_symbol := factor(gene_symbol, levels = gene_symbol[order(log2FoldChange_shrunk)])]
p2a <- ggplot(de, aes(gene_symbol, log2FoldChange_shrunk)) +
  geom_hline(yintercept = 1, linetype = 2, colour = pal["vermillion"], linewidth = 0.45) +
  geom_segment(aes(xend = gene_symbol, y = 0, yend = log2FoldChange_shrunk),
               colour = "#B8C9D6", linewidth = 1.0) +
  geom_errorbar(aes(ymin = log2FoldChange_shrunk - 1.96 * lfcSE_shrunk,
                    ymax = log2FoldChange_shrunk + 1.96 * lfcSE_shrunk),
                width = 0.16, colour = "#333333", linewidth = 0.45) +
  geom_point(size = 3.0, colour = pal["blue"]) + coord_flip() +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(title = "BacCore gene expression", x = NULL, y = "Shrunken log2 fold change") + theme_pub()

mt <- fread(file.path(src_dir, "Figure_2B_module_trait.tsv"))
mt[, module_label := factor(as.character(module_label), levels = c("turquoise", "yellow", "blue"))]
mt[, trait_label := factor(as.character(trait_label),
                           levels = c("BacCore", "Pathological stage", "Death status", "Tumor purity"))]
p2b <- ggplot(mt, aes(trait_label, module_label, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.7) + geom_text(aes(label = lab), size = 3.0) +
  scale_fill_gradient2(low = pal["blue"], mid = "white", high = pal["vermillion"],
                       midpoint = 0, limits = c(-0.5, 0.5), name = "r") +
  labs(title = "Module-trait correlations", x = NULL, y = "Module") + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

over <- fread(file.path(src_dir, "Figure_2C_hub_DEG_overlap.tsv"))
hc <- over[, .N, by = .(module_color, main_regulation)]
hc[, module_color := factor(module_color, levels = c("blue", "turquoise", "yellow"))]
p2c <- ggplot(hc, aes(module_color, N, fill = main_regulation)) +
  geom_col(width = 0.68, colour = "white") +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), colour = "white",
            fontface = "bold", size = 3.3) +
  scale_fill_manual(values = c(Down = unname(pal["blue"]), Up = unname(pal["vermillion"]))) +
  labs(title = "Hub-DEG overlap", x = "Module", y = "Hub-DEG count") + theme_pub() +
  theme(legend.position = "top")
fig2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1.08, 0.82)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 12)))
save_figure(fig2, "Figure_2_expression_WGCNA", 7.2, 7.4)

# Figure 3: pathway context and apoptosome boundary
pa <- fread(file.path(src_dir, "Figure_3A_TCGA_pathways.tsv"))
pa[, label := factor(label, levels = label[order(adjusted_beta)])]
p3a <- ggplot(pa, aes(adjusted_beta, label, colour = class)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = adjusted_beta, yend = label), linewidth = 1.0, alpha = 0.75) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c("Mitochondrial and DNA stress" = unname(pal["vermillion"]),
                                 "Other prespecified programs" = unname(pal["blue"]))) +
  labs(title = "Adjusted pathway associations", x = "Adjusted standardized beta", y = NULL) +
  theme_pub() + theme(legend.position = "bottom")

tcga_sc <- fread(file.path(src_dir, "Figure_3B_TCGA_independent_mtDNA.tsv"))
p3b <- ggplot(tcga_sc, aes(BacCore_z, mtDNA_independent_z)) +
  geom_point(size = 1.0, alpha = 0.38, colour = pal["blue"]) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, colour = pal["vermillion"], fill = "#F7C9B6") +
  labs(title = "Independent mtDNA stress", x = "BacCore score (z)",
       y = "Independent mtDNA-stress score (z)") + theme_pub()

it <- fread(file.path(src_dir, "Figure_3C_apoptosome_interactions.tsv"))
it[, gate_label := factor(as.character(gate_label),
                          levels = c("APAF1", "APAF1 + CASP9", "APAF1 + CYCS + CASP9"))]
p3c <- ggplot(it, aes(interaction_beta, gate_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#666666") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.55) +
  geom_point(size = 2.8, colour = pal["purple"]) +
  annotate("text", x = mean(range(c(it$lo, it$hi))), y = 3.35,
           label = "No interaction passed FDR correction", hjust = 0.5, size = 2.7, fontface = "bold") +
  coord_cartesian(xlim = c(min(it$lo) - 0.01, max(it$hi) + 0.02), clip = "off") +
  labs(title = "Apoptosome interactions", x = "BacCore x apoptosome interaction beta", y = NULL) +
  theme_pub()
fig3 <- p3a / (p3b | p3c) + plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 12)))
save_figure(fig3, "Figure_3_pathways_mechanistic_boundary", 7.2, 7.4)

# Figure 4: sample-level single-cell localization
pb <- fread(file.path(src_dir, "Figure_4_single_cell_pseudobulk.tsv"))
state_levels <- c("Normal epithelial", "Primary tumor epithelial", "Metastatic malignant")
pb[, epithelial_state := factor(epithelial_state, levels = state_levels)]
state_cols <- c("Normal epithelial" = unname(pal["sky"]),
                "Primary tumor epithelial" = unname(pal["orange"]),
                "Metastatic malignant" = unname(pal["vermillion"]))
p4a <- ggplot(pb, aes(epithelial_state, BACIPTOSIS_CORE, fill = epithelial_state)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.78, linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.8, shape = 21, colour = "white") +
  scale_fill_manual(values = state_cols) +
  labs(title = "BacCore across epithelial states", x = NULL, y = "BacCore score") + theme_pub() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1))
cor_dt <- fread(file.path(src_dir, "Figure_4B_single_cell_correlations.tsv"))
cor_dt[, pathway := factor(pathway, levels = pathway[order(rho)])]
p4b <- ggplot(cor_dt, aes(rho, pathway)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = rho, yend = pathway), linewidth = 1.0, colour = "#9EC7E2") +
  geom_point(size = 2.8, colour = pal["blue"]) +
  geom_text(aes(label = paste0("P=", fmt_p(p))), hjust = -0.08, size = 2.7) +
  coord_cartesian(xlim = c(0, 0.72), clip = "off") +
  labs(title = "Pathway correlations", x = "Spearman rho", y = NULL) + theme_pub()
p4c <- ggplot(pb, aes(epithelial_state, mtDNA_independent, fill = epithelial_state)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.78, linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.8, shape = 21, colour = "white") +
  scale_fill_manual(values = state_cols) +
  labs(title = "mtDNA stress across states", x = NULL, y = "Independent mtDNA-stress score") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1))
fig4 <- p4a / (p4b | p4c) + plot_layout(heights = c(1.0, 0.95)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 12)))
save_figure(fig4, "Figure_4_single_cell_localization", 7.2, 7.4)

# Figure 5: CPTAC protein validation
pair_score <- fread(file.path(src_dir, "Figure_5A_CPTAC_paired_scores.tsv"))
pair_score[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
p5a <- ggplot(pair_score, aes(Tissue, BacCore_protein, group = pair_id)) +
  geom_line(alpha = 0.16, linewidth = 0.35, colour = "#777777") +
  geom_point(aes(colour = Tissue), size = 1.35, alpha = 0.75) +
  geom_boxplot(aes(group = Tissue), width = 0.42, outlier.shape = NA, fill = NA, linewidth = 0.65) +
  scale_colour_manual(values = c(Normal = unname(pal["sky"]), Tumor = unname(pal["vermillion"]))) +
  labs(title = "Paired protein validation", x = NULL, y = "Six-protein BacCore score") +
  theme_pub() + theme(legend.position = "none")
rp <- fread(file.path(src_dir, "Figure_5B_CPTAC_RNA_protein_scores.tsv"))
p5b <- ggplot(rp, aes(RNA_BacCore, Protein_BacCore)) +
  geom_point(size = 1.35, alpha = 0.55, colour = pal["blue"]) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.75, colour = pal["vermillion"], fill = "#F7C9B6") +
  labs(title = "RNA-protein concordance", x = "RNA BacCore score", y = "Protein BacCore score") + theme_pub()
pc <- fread(file.path(src_dir, "Figure_5C_CPTAC_pathways.tsv"))
pc[, label := factor(as.character(label),
                     levels = c("G2M checkpoint", "E2F targets", "cGAS-STING type I IFN", "Independent mtDNA stress"))]
p5c <- ggplot(pc, aes(rho, label)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = rho, yend = label), linewidth = 1.0,
               colour = ifelse(pc$FDR < 0.05, pal["vermillion"], "#B8B8B8")) +
  geom_point(aes(colour = FDR < 0.05), size = 2.8) +
  scale_colour_manual(values = c(`TRUE` = unname(pal["vermillion"]), `FALSE` = "#777777"), guide = "none") +
  labs(title = "Protein pathway associations", x = "Spearman rho", y = NULL) + theme_pub()
fig5 <- p5a / (p5b | p5c) + plot_layout(heights = c(1.0, 0.95)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 12)))
save_figure(fig5, "Figure_5_CPTAC_validation", 7.2, 7.4)

# Figure 6: survival meta-analysis
sv <- fread(file.path(src_dir, "Figure_6_survival_cohorts.tsv"))
meta <- fread(file.path(src_dir, "Figure_6_meta_analysis.tsv"))
sv[, cohort := factor(cohort, levels = rev(c("TCGA-LUAD", "GSE68465", "GSE31210", "GSE72094")))]
unadj <- sv[, .(cohort, HR, lo = CI_low, hi = CI_high, p)]
unadj <- rbind(unadj,
               data.table(cohort = factor("Random effects", levels = c(levels(sv$cohort), "Random effects")),
                          HR = meta$pooled_HR, lo = meta$CI_low, hi = meta$CI_high, p = meta$p), fill = TRUE)
unadj[, cohort_label := factor(as.character(cohort),
                               levels = c("Random effects", "GSE72094", "GSE31210", "GSE68465", "TCGA-LUAD"))]
p6a <- ggplot(unadj, aes(HR, cohort_label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#666666") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14, linewidth = 0.55) +
  geom_point(aes(shape = cohort_label == "Random effects"), size = 3.0, colour = pal["blue"], fill = pal["blue"]) +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
  scale_x_log10(limits = c(0.72, 2.05), breaks = c(0.75, 1, 1.25, 1.5, 2)) +
  labs(title = "Univariable association", x = "Hazard ratio per SD increase in BacCore", y = NULL) + theme_pub()
adj <- sv[, .(cohort, HR = adjusted_HR, lo = adjusted_CI_low, hi = adjusted_CI_high, p = adjusted_p)]
p6b <- ggplot(adj, aes(HR, cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#666666") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14, linewidth = 0.55) +
  geom_point(size = 3.0, colour = pal["vermillion"]) +
  scale_x_log10(limits = c(0.72, 2.05), breaks = c(0.75, 1, 1.25, 1.5, 2)) +
  labs(title = "Adjusted association", x = "Adjusted hazard ratio per SD", y = NULL) + theme_pub()
p6c <- ggplot(sv, aes(C_index, cohort)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#777777") +
  geom_segment(aes(x = 0.5, xend = C_index, yend = cohort), colour = "#9EC7E2", linewidth = 1.0) +
  geom_point(size = 3.0, colour = pal["blue"]) +
  geom_text(aes(label = sprintf("%.3f", C_index)), hjust = -0.25, size = 2.8) +
  coord_cartesian(xlim = c(0.49, 0.66), clip = "off") +
  labs(title = "C index", x = "Harrell C index", y = NULL) + theme_pub()
fig6 <- (p6a | p6b) / p6c + plot_layout(heights = c(1.08, 0.82)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 12)))
save_figure(fig6, "Figure_6_survival_meta_analysis", 7.2, 6.6)

message("Created 12 files (PNG and PDF) in: ", fig_dir)

