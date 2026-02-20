#!/usr/bin/env Rscript
# =============================================================================
# generate_purine_pyrimidine_panels.R
#
# Comprehensive GSEA visualization panel:
#   - Enrichment curves (focused + KEGG canonical)
#   - Dot plots, Ridge plots, Upset plots
#   - emap plots (pathway network)
#   - cnet plots (gene-concept network)
#   - Tree plots (hierarchical clustering of pathways)
#   - Combined purine & pyrimidine panels
#
# Databases: KEGG, GO-BP, Reactome + Custom focused gene sets
# Comparison: NAD-treated vs Anti-CD3/28 control T cells
# =============================================================================

user_lib <- path.expand("~/R/library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# Install any missing packages
bioc_pkgs <- c("clusterProfiler", "enrichplot", "fgsea", "msigdbr",
               "org.Hs.eg.db", "ReactomePA", "DOSE")
for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, lib = user_lib, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(patchwork)
  library(ggplot2)
  library(openxlsx)
})

# Resolve namespace conflicts
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# ── Paths ─────────────────────────────────────────────────────────────────── #
BASE <- "./"
DATA <- file.path(BASE, "data")
RES  <- file.path(BASE, "results")
FIGS <- file.path(RES,  "figures")
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# ── Load DE results ───────────────────────────────────────────────────────── #
message("Loading DE results...")
de_raw <- read.delim(file.path(DATA, "DE_results.txt"),
                     stringsAsFactors = FALSE, check.names = FALSE)

de <- de_raw %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue), Symbol != "") %>%
  mutate(
    pvalue_floor = pmax(pvalue, 1e-300),
    rank_score   = sign(log2FoldChange) * (-log10(pvalue_floor)) +
                   runif(n(), -1e-6, 1e-6),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up (NAD)",
      padj < 0.05 & log2FoldChange < -1 ~ "Down (NAD)",
      TRUE                               ~ "NS"
    )
  )

# Symbol-based ranked list (for fgsea / custom GSEA)
ranked_sym <- de %>%
  arrange(desc(rank_score)) %>%
  { setNames(.$rank_score, .$Symbol) }
ranked_sym <- ranked_sym[!duplicated(names(ranked_sym))]

# Entrez-based ranked list (for clusterProfiler KEGG / GO / Reactome)
ranked_entrez <- de %>%
  filter(!is.na(GeneID), GeneID != "") %>%
  arrange(desc(rank_score)) %>%
  { setNames(.$rank_score, as.character(.$GeneID)) }
ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]

message(sprintf("Ranked: %d symbols | %d Entrez IDs",
                length(ranked_sym), length(ranked_entrez)))

# ── Focused gene sets ─────────────────────────────────────────────────────── #
purine_denovo <- unique(c(
  "PPAT", "GART", "PFAS", "PAICS", "ADSS", "ADSL", "ATIC",
  "IMPDH1", "IMPDH2", "GMPS", "ADSS1", "ADSS2", "PRPS1", "PRPS2"
))
purine_salvage <- unique(c(
  "HPRT1", "APRT", "ADA", "ADK", "PNP", "DGUOK",
  "RRM1", "RRM2", "RRM2B",
  "NT5E", "NT5C1A", "NT5C2", "NT5C3A",
  "ENTPD1", "ENTPD2", "ENTPD3",
  "AK1", "AK2", "AK3", "AK4", "GUK1"
))
pyrim_denovo <- unique(c(
  "CAD", "DHODH", "UMPS", "CTPS1", "CTPS2"
))
pyrim_salvage <- unique(c(
  "TK1", "TK2", "TYMS", "TYMP",
  "UCK1", "UCK2", "CMPK1", "CMPK2",
  "DCTD", "UPP1", "UPP2", "DPYS", "DPYD",
  "NME1", "NME2", "UNG", "SMUG1", "MBD4"
))
glut_genes <- unique(c(
  "GLS", "GLS2", "GLUD1", "GLUD2", "GLUL",
  "GOT1", "GOT2", "GPT", "GPT2", "PSAT1", "PSPH",
  "SLC1A1","SLC1A2","SLC1A3","SLC1A4","SLC1A5",
  "SLC38A1","SLC38A2","SLC38A5",
  "IDH1","IDH2","IDH3A","IDH3B","IDH3G",
  "OGDH","OGDHL","SUCLA2","SUCLG1","SUCLG2",
  "ASNS","ASS1","ASL","PYCR1","PYCR2","PYCR3","ALDH18A1","CPS1","PPAT"
))
nad_genes <- unique(c(
  "NAMPT","NMNAT1","NMNAT2","NMNAT3","NADK","NADK2","NAPRT","NADSYN1",
  "IDO1","IDO2","TDO2","KYNU","HAAO","ACMSD","QPRT",
  "CD38","BST1","PARP1","PARP2","PARP3",
  "SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7",
  "NNMT","NUDT12","NUDT5"
))

focused_sets <- list(
  Purine_De_Novo       = purine_denovo,
  Purine_Salvage       = purine_salvage,
  Pyrimidine_De_Novo   = pyrim_denovo,
  Pyrimidine_Salvage   = pyrim_salvage,
  Glutamate_Glutamine  = glut_genes,
  NAD_Metabolism       = nad_genes
)

# Build TERM2GENE for clusterProfiler::GSEA()
term2gene_focused <- purrr::imap_dfr(focused_sets,
  ~ data.frame(term = .y, gene = .x))

# ── Helper: safe ggsave ───────────────────────────────────────────────────── #
safe_save <- function(p, file, width = 12, height = 9, ...) {
  tryCatch(
    ggsave(file, p, width = width, height = height, ...),
    error = function(e) message("  [WARN] Could not save ", basename(file), ": ", e$message)
  )
  message("  Saved: ", basename(file))
}

# ── Helper: enrichment curve with stats overlay ───────────────────────────── #
make_curve <- function(set_name, sets, ranked, stats_df = NULL) {
  if (!set_name %in% names(sets)) return(NULL)
  genes <- sets[[set_name]]
  n_matched <- sum(genes %in% names(ranked))
  sub_str <- sprintf("Genes in set: %d | Matched: %d", length(genes), n_matched)
  if (!is.null(stats_df)) {
    row <- stats_df %>% filter(pathway == set_name)
    if (nrow(row) == 1) {
      stars <- case_when(row$padj < 0.001 ~ "***", row$padj < 0.01 ~ "**",
                         row$padj < 0.05  ~ "*",   TRUE ~ "ns")
      sub_str <- paste0(sub_str,
        sprintf("  |  NES = %.2f  |  padj = %.3g  %s", row$NES, row$padj, stars))
    }
  }
  plotEnrichment(genes, ranked) +
    labs(title    = str_replace_all(set_name, "_", " "),
         subtitle = sub_str,
         x = "Gene rank  (left = activated in NAD)",
         y = "Enrichment score") +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}

# ── Helper: save all enrichplot visualizations for a result object ─────────── #
save_enrichplot_suite <- function(res_obj, prefix, label,
                                  showCategory = 20, font_size = 9) {
  if (is.null(res_obj)) { message("  NULL object, skipping: ", label); return(invisible()) }

  result_df <- tryCatch(as.data.frame(res_obj), error = function(e) NULL)
  if (is.null(result_df) || nrow(result_df) == 0) {
    message("  No results for: ", label); return(invisible())
  }

  n_sig <- sum(result_df$p.adjust < 0.05, na.rm = TRUE)
  message(sprintf("  %s: %d significant pathways (p.adjust<0.05)", label, n_sig))

  # 1. Dot plot
  tryCatch({
    p <- dotplot(res_obj, showCategory = showCategory) +
      ggtitle(paste(label, "– Dot plot")) +
      theme_bw(base_size = font_size)
    safe_save(p, file.path(FIGS, paste0(prefix, "_dotplot.pdf")), 11, 9)
  }, error = function(e) message("  [WARN] dotplot failed: ", e$message))

  # 2. Ridge plot (requires log2 fold change column — only for GSEA objects)
  tryCatch({
    p <- ridgeplot(res_obj, showCategory = min(showCategory, 20)) +
      ggtitle(paste(label, "– Ridge plot")) +
      theme_bw(base_size = font_size)
    safe_save(p, file.path(FIGS, paste0(prefix, "_ridgeplot.pdf")), 12, 10)
  }, error = function(e) message("  [WARN] ridgeplot failed: ", e$message))

  # 3. Upset plot
  tryCatch({
    pdf(file.path(FIGS, paste0(prefix, "_upsetplot.pdf")), width = 12, height = 7)
    print(upsetplot(res_obj, n = min(showCategory, 10)))
    dev.off()
    message("  Saved: ", paste0(prefix, "_upsetplot.pdf"))
  }, error = function(e) {
    dev.off()
    message("  [WARN] upsetplot failed: ", e$message)
  })

  # 4. emap plot (needs pairwise similarity)
  tryCatch({
    res_sim <- pairwise_termsim(res_obj)
    p <- emapplot(res_sim, showCategory = min(showCategory, 30)) +
      ggtitle(paste(label, "– Pathway network (emap)"))
    safe_save(p, file.path(FIGS, paste0(prefix, "_emapplot.pdf")), 12, 10)
  }, error = function(e) message("  [WARN] emapplot failed: ", e$message))

  # 5. cnet plot (gene-concept network)
  tryCatch({
    # Build fold-change vector with gene symbols for cnetplot
    fc_sym <- setNames(de$log2FoldChange, de$Symbol)
    p <- cnetplot(res_obj,
                  showCategory = min(showCategory, 15),
                  foldChange   = fc_sym,
                  circular     = FALSE,
                  colorEdge    = TRUE) +
      ggtitle(paste(label, "– Gene-concept network (cnet)")) +
      scale_colour_gradient2(low = "#377EB8", mid = "grey90", high = "#E41A1C",
                             name = "log2FC")
    safe_save(p, file.path(FIGS, paste0(prefix, "_cnetplot.pdf")), 14, 12)
  }, error = function(e) message("  [WARN] cnetplot failed: ", e$message))

  # 6. Tree plot
  tryCatch({
    res_sim <- pairwise_termsim(res_obj)
    p <- treeplot(res_sim, showCategory = min(showCategory, 30),
                  fontsize = font_size - 1) +
      ggtitle(paste(label, "– Pathway hierarchy (tree)"))
    safe_save(p, file.path(FIGS, paste0(prefix, "_treeplot.pdf")), 14, 10)
  }, error = function(e) message("  [WARN] treeplot failed: ", e$message))

  invisible()
}

# ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 1: Focused gene-set GSEA (custom sets via clusterProfiler::GSEA)
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n── Section 1: Focused GSEA (custom gene sets) ──")

gsea_focused_cp <- GSEA(
  geneList     = ranked_sym,
  TERM2GENE    = term2gene_focused,
  minGSSize    = 5,
  maxGSSize    = 600,
  pvalueCutoff = 1.0,       # keep all; we'll filter visually
  eps          = 0,
  nPermSimple  = 100000,
  verbose      = FALSE
)

message("\nFocused GSEA results:")
print(as.data.frame(gsea_focused_cp) %>%
      select(ID, setSize, NES, pvalue, p.adjust) %>%
      arrange(p.adjust))

# fgsea version for plotEnrichment (needs list + stats)
gsea_focused_fgsea <- fgsea(
  pathways    = focused_sets,
  stats       = ranked_sym,
  minSize     = 5, maxSize = 600,
  eps = 0, nPermSimple = 100000
) %>% as_tibble()

save_enrichplot_suite(gsea_focused_cp, "23_focused", "Focused Pathways GSEA",
                      showCategory = 6)

# ── Enrichment curves for all focused sets ────────────────────────────────── #
message("  Generating individual enrichment curves (focused sets)...")
curve_plots <- lapply(names(focused_sets), function(nm) {
  make_curve(nm, focused_sets, ranked_sym, gsea_focused_fgsea)
})
names(curve_plots) <- names(focused_sets)

# Combined panel: all 6 focused sets (2 rows × 3 cols)
panel_all6 <- wrap_plots(curve_plots, nrow = 2, ncol = 3) +
  plot_annotation(
    title    = "Focused Pathway GSEA – Enrichment Curves",
    subtitle = "NAD-treated vs Anti-CD3/28 control T cells",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
safe_save(panel_all6, file.path(FIGS, "24_GSEA_curves_all_focused.pdf"), 18, 11)

# ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 2: KEGG GSEA (clusterProfiler gseKEGG)
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n── Section 2: KEGG GSEA (gseKEGG) ──")

gsea_kegg_cp <- gseKEGG(
  geneList     = ranked_entrez,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  nPermSimple  = 10000,
  verbose      = FALSE
)

message("KEGG GSEA: ", nrow(gsea_kegg_cp@result), " results")
save_enrichplot_suite(gsea_kegg_cp, "25_KEGG", "KEGG GSEA")

# ── KEGG curves for purine & pyrimidine canonical pathways ─────────────────── #
message("  Generating KEGG purine/pyrimidine enrichment curves...")
kegg_sets_msig <- msigdbr(species = "Homo sapiens", collection = "C2",
                           subcollection = "CP:KEGG_LEGACY") %>%
  split(.$gs_name) %>% lapply(`[[`, "gene_symbol")

kegg_pur_fgsea <- fgsea(
  pathways    = kegg_sets_msig[c("KEGG_PURINE_METABOLISM", "KEGG_PYRIMIDINE_METABOLISM")],
  stats       = ranked_sym,
  minSize = 10, maxSize = 500, eps = 0, nPermSimple = 50000
) %>% as_tibble()

p_kegg_pur <- make_curve("KEGG_PURINE_METABOLISM",    kegg_sets_msig, ranked_sym, kegg_pur_fgsea)
p_kegg_pyr <- make_curve("KEGG_PYRIMIDINE_METABOLISM", kegg_sets_msig, ranked_sym, kegg_pur_fgsea)

# ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 3: GO Biological Process GSEA
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n── Section 3: GO BP GSEA (gseGO) ──")

gsea_go_cp <- gseGO(
  geneList     = ranked_entrez,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  nPermSimple  = 10000,
  verbose      = FALSE
)
gsea_go_cp <- setReadable(gsea_go_cp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

message("GO-BP GSEA: ", nrow(gsea_go_cp@result), " results")
save_enrichplot_suite(gsea_go_cp, "26_GO_BP", "GO BP GSEA", showCategory = 25)

# ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 4: Reactome GSEA (ReactomePA)
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n── Section 4: Reactome GSEA (gsePathway) ──")

gsea_react_cp <- gsePathway(
  geneList     = ranked_entrez,
  organism     = "human",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  verbose      = FALSE
)
gsea_react_cp <- setReadable(gsea_react_cp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

message("Reactome GSEA: ", nrow(gsea_react_cp@result), " results")
save_enrichplot_suite(gsea_react_cp, "27_Reactome", "Reactome GSEA", showCategory = 25)

# ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 5: Purine & Pyrimidine combined panels
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n── Section 5: Combined Purine & Pyrimidine panels ──")

# Enrichment curves (focused + KEGG canonical)
p_pur_dn  <- curve_plots[["Purine_De_Novo"]]
p_pur_sal <- curve_plots[["Purine_Salvage"]]
p_pyr_dn  <- curve_plots[["Pyrimidine_De_Novo"]]
p_pyr_sal <- curve_plots[["Pyrimidine_Salvage"]]

# Panel A: Purine (3-panel with KEGG canonical)
if (!is.null(p_kegg_pur)) {
  panel_purine <- (p_pur_dn | p_pur_sal) / p_kegg_pur +
    plot_annotation(
      title    = "Purine Metabolism – GSEA Enrichment",
      subtitle = "NAD-treated vs Anti-CD3/28 control T cells  |  Rank metric: sign(LFC) × −log10(p)",
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, colour = "grey30"))
    )
  safe_save(panel_purine, file.path(FIGS, "19_GSEA_panel_Purine_enrichment.pdf"), 14, 10)
}

# Panel B: Pyrimidine (3-panel with KEGG canonical)
if (!is.null(p_kegg_pyr)) {
  panel_pyrimidine <- (p_pyr_dn | p_pyr_sal) / p_kegg_pyr +
    plot_annotation(
      title    = "Pyrimidine Metabolism – GSEA Enrichment",
      subtitle = "NAD-treated vs Anti-CD3/28 control T cells  |  Rank metric: sign(LFC) × −log10(p)",
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, colour = "grey30"))
    )
  safe_save(panel_pyrimidine, file.path(FIGS, "20_GSEA_panel_Pyrimidine_enrichment.pdf"), 14, 10)
}

# Panel C: 2×2 all four focused sets
panel_2x2 <- (p_pur_dn | p_pur_sal) / (p_pyr_dn | p_pyr_sal) +
  plot_annotation(
    title    = "Purine & Pyrimidine Metabolism – GSEA Enrichment",
    subtitle = "NAD-treated vs Anti-CD3/28 control T cells",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, colour = "grey30"))
  )
safe_save(panel_2x2, file.path(FIGS, "21_GSEA_panel_Purine_Pyrimidine_combined.pdf"), 14, 10)

# ── Lollipop NES summary (focused + KEGG canonical) ──────────────────────── #
message("  Building NES lollipop summary...")

summary_df <- bind_rows(
  gsea_focused_fgsea %>%
    select(pathway, NES, pval, padj) %>%
    mutate(source = "Custom gene set"),
  kegg_pur_fgsea %>%
    select(pathway, NES, pval, padj) %>%
    mutate(source = "KEGG canonical")
) %>%
  mutate(
    label     = str_replace_all(pathway, "_", " ") %>% str_wrap(30),
    direction = ifelse(NES > 0, "Activated (NAD)", "Suppressed (NAD)"),
    sig_tag   = case_when(
      padj < 0.001 ~ "***", padj < 0.01 ~ "**",
      padj < 0.05  ~ "*",   TRUE         ~ "ns"
    ),
    group = case_when(
      str_detect(pathway, "Purine|PURINE")         ~ "Purine",
      str_detect(pathway, "Pyrimidine|PYRIMIDINE") ~ "Pyrimidine",
      str_detect(pathway, "Glutamate|GLUTAMATE")   ~ "Glutamate/Glutamine",
      str_detect(pathway, "NAD")                   ~ "NAD+ Metabolism",
      TRUE                                          ~ "Other"
    )
  ) %>%
  arrange(group, NES)

p_lollipop <- ggplot(summary_df,
    aes(NES, reorder(label, NES),
        colour = direction, size = -log10(pval), shape = source)) +
  geom_segment(aes(x = 0, xend = NES, yend = label),
               linewidth = 0.8, colour = "grey70") +
  geom_point() +
  geom_text(aes(label = sig_tag),
            hjust = ifelse(summary_df$NES > 0, -0.5, 1.5),
            size = 4.5, colour = "black", fontface = "bold") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(
    values = c("Activated (NAD)" = "#E41A1C", "Suppressed (NAD)" = "#377EB8"),
    name = NULL
  ) +
  scale_size_continuous(range = c(4, 12), name = "-log10(p)") +
  scale_shape_manual(values = c("Custom gene set" = 16, "KEGG canonical" = 17),
                     name = "Gene set source") +
  labs(
    title    = "Purine & Pyrimidine Metabolism: GSEA Summary",
    subtitle = "NAD-treated vs Control T cells  |  * p<0.05  ** p<0.01  *** p<0.001",
    x = "Normalized Enrichment Score (NES)  [+ = Activated in NAD]",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey92"),
    plot.title       = element_text(face = "bold", size = 14)
  )

safe_save(p_lollipop, file.path(FIGS, "22_GSEA_lollipop_Purine_Pyrimidine.pdf"), 11, 9)

# ═══════════════════════════════════════════════════════════════════════════ #
#  Done
# ═══════════════════════════════════════════════════════════════════════════ #
message("\n==============================================")
message("  All figures saved to: ", FIGS)
message("  Section 1 – Focused (custom) gene sets:")
message("    23_focused_dotplot.pdf")
message("    23_focused_ridgeplot.pdf")
message("    23_focused_emapplot.pdf")
message("    23_focused_cnetplot.pdf")
message("    23_focused_treeplot.pdf")
message("    24_GSEA_curves_all_focused.pdf")
message("  Section 2 – KEGG GSEA:")
message("    25_KEGG_dotplot / ridgeplot / upsetplot / emapplot / cnetplot / treeplot")
message("  Section 3 – GO BP GSEA:")
message("    26_GO_BP_dotplot / ridgeplot / upsetplot / emapplot / cnetplot / treeplot")
message("  Section 4 – Reactome GSEA:")
message("    27_Reactome_dotplot / ridgeplot / upsetplot / emapplot / cnetplot / treeplot")
message("  Section 5 – Combined Purine & Pyrimidine panels:")
message("    19_GSEA_panel_Purine_enrichment.pdf")
message("    20_GSEA_panel_Pyrimidine_enrichment.pdf")
message("    21_GSEA_panel_Purine_Pyrimidine_combined.pdf")
message("    22_GSEA_lollipop_Purine_Pyrimidine.pdf")
message("==============================================\n")
