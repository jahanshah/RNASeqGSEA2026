#!/usr/bin/env Rscript
# =============================================================================
# RNASeqGSEA2026 – Full GSEA & Pathway Analysis
# Comparison : Anti-CD3/28 (control, stimulated T cells)
#              vs Anti-CD3/28 + NAD (treatment)
# Rank metric: sign(log2FC) × -log10(p-value)
# Focus      : Glutamate/Glutamine metabolism,
#              Purine  (de novo + salvage),
#              Pyrimidine (de novo + salvage),
#              NAD+ metabolism
#
# NOTE: DARs (ATAC-seq) – add a DAR results file to data/ and see
#       Section 17 for the integration stub.
# =============================================================================

# ── 0. PERSONAL LIBRARY (avoids system-level permission errors) ───────────── #
user_lib <- path.expand("~/R/library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))   # personal lib takes priority

# ── 1. PACKAGES ──────────────────────────────────────────────────────────── #
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = user_lib)

bioc_pkgs <- c("clusterProfiler", "enrichplot", "fgsea", "msigdbr",
               "org.Hs.eg.db", "ReactomePA", "DOSE", "pathview")
cran_pkgs <- c("tidyverse", "ggrepel", "pheatmap", "RColorBrewer",
               "viridis", "patchwork", "openxlsx")

for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, lib = user_lib, ask = FALSE, update = FALSE)
for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = user_lib)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(pathview)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(openxlsx)
  library(ggrepel)
})

# ── 2. PATHS ─────────────────────────────────────────────────────────────── #
BASE <- "/Users/jahanshah/Documents/Consultant-NGS/Ahmed/Projects2026/Project-RNAseq/RNASeqGSEA2026"
DATA <- file.path(BASE, "data")
RES  <- file.path(BASE, "results")
FIGS <- file.path(RES,  "figures")
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# Sample metadata
CONTROL_COLS <- c("KO-13_S31_R1", "KO-14_S32_R2", "KO-15_S33_R1")
NAD_COLS     <- c("KO-16_S34_R1", "KO-17_S35_R1", "KO-18_S36_R1")
COUNT_COLS   <- c(CONTROL_COLS, NAD_COLS)

# ── 3. LOAD & PROCESS DE RESULTS ─────────────────────────────────────────── #
message("Loading DE results...")
de_raw <- read.delim(file.path(DATA, "DE_results.txt"), stringsAsFactors = FALSE)

de <- de_raw %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue), Symbol != "") %>%
  mutate(
    pvalue_floor = pmax(pvalue, 1e-300),           # floor for -log10
    rank_score   = sign(log2FoldChange) * (-log10(pvalue_floor)),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up (NAD)",
      padj < 0.05 & log2FoldChange < -1 ~ "Down (NAD)",
      TRUE                               ~ "NS"
    )
  )

message(sprintf(
  "Total genes: %d | DEGs (padj<0.05, |LFC|>1): Up=%d, Down=%d",
  nrow(de), sum(de$sig == "Up (NAD)"), sum(de$sig == "Down (NAD)")
))

# ── 4. VOLCANO PLOT ───────────────────────────────────────────────────────── #
message("Plotting volcano...")
top_labels <- de %>%
  filter(sig != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 35)

p_volcano <- ggplot(de, aes(log2FoldChange, -log10(pvalue_floor), colour = sig)) +
  geom_point(data = filter(de, sig == "NS"), alpha = 0.3, size = 0.6) +
  geom_point(data = filter(de, sig != "NS"), alpha = 0.8, size = 1.0) +
  geom_text_repel(
    data = top_labels, aes(label = Symbol),
    size = 2.5, max.overlaps = 25, segment.size = 0.2, segment.alpha = 0.5
  ) +
  scale_colour_manual(
    values = c("Up (NAD)" = "#E41A1C", "Down (NAD)" = "#377EB8", NS = "grey75")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  labs(
    title    = "Volcano: Anti-CD3/28 vs Anti-CD3/28 + NAD",
    subtitle = sprintf(
      "Up in NAD: %d | Down in NAD: %d  (padj < 0.05, |LFC| > 1)",
      sum(de$sig == "Up (NAD)"), sum(de$sig == "Down (NAD)")
    ),
    x      = "log\u2082 Fold Change  (NAD / Control)",
    y      = "-log\u2081\u2080(p-value)",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(FIGS, "01_volcano_plot.pdf"), p_volcano, width = 8, height = 7)
message("  Saved: 01_volcano_plot.pdf")

# ── 5. GENE RANKING VECTORS ───────────────────────────────────────────────── #
# Symbol-based (for MSigDB / fgsea)
ranked_sym <- de %>%
  arrange(desc(rank_score)) %>%
  { setNames(.$rank_score, .$Symbol) }
ranked_sym <- ranked_sym[!duplicated(names(ranked_sym))]

# Entrez-based (for clusterProfiler GO/KEGG / ReactomePA)
ranked_entrez <- de %>%
  filter(!is.na(GeneID), GeneID != "") %>%
  arrange(desc(rank_score)) %>%
  { setNames(.$rank_score, as.character(.$GeneID)) }
ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]

# DEG subsets for ORA
up_entrez   <- de %>% filter(sig == "Up (NAD)")   %>% pull(GeneID) %>% as.character()
down_entrez <- de %>% filter(sig == "Down (NAD)") %>% pull(GeneID) %>% as.character()
universe    <- de %>% filter(!is.na(GeneID))      %>% pull(GeneID) %>% as.character()

message(sprintf("Ranked list: %d symbols | %d Entrez IDs", length(ranked_sym), length(ranked_entrez)))

# ── 6. LOAD MSigDB GENE SETS ─────────────────────────────────────────────── #
message("Loading MSigDB gene sets...")

to_list <- function(df) split(df$gene_symbol, df$gs_name)

h_sets     <- msigdbr(species = "Homo sapiens", category = "H") %>% to_list()
kegg_sets  <- msigdbr(species = "Homo sapiens", category = "C2",
                      subcategory = "CP:KEGG_LEGACY") %>% to_list()
react_sets <- msigdbr(species = "Homo sapiens", category = "C2",
                      subcategory = "CP:REACTOME") %>% to_list()

# ── 7. FULL PRERANKED GSEA (fgsea) ───────────────────────────────────────── #
run_fgsea <- function(pathways, ranked, label, nPerm = 10000) {
  message(sprintf("  fgsea: %s ...", label))
  fgsea(pathways    = pathways,
        stats       = ranked,
        minSize     = 10,
        maxSize     = 500,
        nPermSimple = nPerm) %>%
    as_tibble() %>%
    mutate(
      direction = ifelse(NES > 0, "Activated (NAD)", "Suppressed (NAD)"),
      database  = label
    ) %>%
    arrange(padj, pval)
}

message("Running GSEA (Hallmark, KEGG, Reactome)...")
gsea_h     <- run_fgsea(h_sets,     ranked_sym, "Hallmark")
gsea_kegg  <- run_fgsea(kegg_sets,  ranked_sym, "KEGG")
gsea_react <- run_fgsea(react_sets, ranked_sym, "Reactome")

message("  fgsea: GO Biological Process (clusterProfiler)...")
gsea_go <- gseGO(
  geneList     = ranked_entrez,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  nPermSimple  = 10000,
  verbose      = FALSE
)

message("  GSEA: Reactome (ReactomePA)...")
gsea_reactome_pa <- gsePathway(
  geneList     = ranked_entrez,
  organism     = "human",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# ── 8. ORA ───────────────────────────────────────────────────────────────── #
message("Running ORA (KEGG + GO) for Up and Down DEGs...")

ora_up_kegg <- enrichKEGG(gene = up_entrez,   universe = universe,
                           organism = "hsa",   pvalueCutoff = 0.05)
ora_dn_kegg <- enrichKEGG(gene = down_entrez, universe = universe,
                           organism = "hsa",   pvalueCutoff = 0.05)
ora_up_react <- enrichPathway(gene = up_entrez,   organism = "human",
                               universe = universe, pvalueCutoff = 0.05,
                               readable = TRUE)
ora_dn_react <- enrichPathway(gene = down_entrez, organism = "human",
                               universe = universe, pvalueCutoff = 0.05,
                               readable = TRUE)
ora_up_go <- enrichGO(gene = up_entrez,   OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID", ont = "BP",
                       universe = universe,  pvalueCutoff = 0.05,
                       readable = TRUE)
ora_dn_go <- enrichGO(gene = down_entrez, OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID", ont = "BP",
                       universe = universe,  pvalueCutoff = 0.05,
                       readable = TRUE)

# ── 9. SAVE ALL RESULT TABLES ────────────────────────────────────────────── #
message("Saving result tables to Excel...")

tidy_fgsea <- function(df)
  df %>% mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>% arrange(padj)

wb <- createWorkbook()
add_sheet <- function(wb, sheet, data) {
  addWorksheet(wb, sheet); writeData(wb, sheet, as.data.frame(data))
}

add_sheet(wb, "GSEA_Hallmark",      tidy_fgsea(gsea_h))
add_sheet(wb, "GSEA_KEGG",          tidy_fgsea(gsea_kegg))
add_sheet(wb, "GSEA_Reactome_fgsea",tidy_fgsea(gsea_react))
add_sheet(wb, "GSEA_GO_BP",         gsea_go@result)
add_sheet(wb, "GSEA_Reactome_PA",   gsea_reactome_pa@result)
add_sheet(wb, "ORA_Up_KEGG",        ora_up_kegg@result)
add_sheet(wb, "ORA_Down_KEGG",      ora_dn_kegg@result)
add_sheet(wb, "ORA_Up_Reactome",    ora_up_react@result)
add_sheet(wb, "ORA_Down_Reactome",  ora_dn_react@result)
add_sheet(wb, "ORA_Up_GO_BP",       ora_up_go@result)
add_sheet(wb, "ORA_Down_GO_BP",     ora_dn_go@result)

saveWorkbook(wb, file.path(RES, "01_GSEA_ORA_all_results.xlsx"), overwrite = TRUE)
message("  Saved: 01_GSEA_ORA_all_results.xlsx")

# ── 10. GLOBAL GSEA DOT PLOTS ────────────────────────────────────────────── #
message("Plotting global GSEA summaries...")

plot_nes_dot <- function(df, title, n_top = 20, strip_prefix = NULL) {
  top <- df %>%
    filter(padj < 0.05) %>%
    group_by(direction) %>%
    slice_min(padj, n = n_top %/% 2, with_ties = FALSE) %>%
    ungroup()
  if (!is.null(strip_prefix))
    top <- mutate(top, pathway = str_remove(pathway, strip_prefix))
  top <- mutate(top, pathway = str_wrap(pathway, 45))
  if (nrow(top) == 0) { message("  No sig pathways: ", title); return(NULL) }
  ggplot(top, aes(NES, reorder(pathway, NES), colour = padj, size = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_viridis_c(option = "plasma", direction = -1, name = "padj") +
    scale_size_continuous(range = c(2, 8), name = "Set size") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
}

p_h    <- plot_nes_dot(gsea_h,     "Hallmark GSEA",  strip_prefix = "HALLMARK_")
p_kegg <- plot_nes_dot(gsea_kegg,  "KEGG GSEA",      strip_prefix = "KEGG_")
p_react<- plot_nes_dot(gsea_react, "Reactome GSEA",  strip_prefix = "REACTOME_")

if (!is.null(p_h))     ggsave(file.path(FIGS, "02_GSEA_Hallmark.pdf"),  p_h,     width = 10, height = 8)
if (!is.null(p_kegg))  ggsave(file.path(FIGS, "03_GSEA_KEGG.pdf"),      p_kegg,  width = 11, height = 9)
if (!is.null(p_react)) ggsave(file.path(FIGS, "04_GSEA_Reactome.pdf"),  p_react, width = 12, height = 10)
message("  Saved: 02-04 GSEA dot plots")

# GO dotplot (clusterProfiler native)
if (nrow(gsea_go@result) > 0) {
  p_go <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GO Biological Process GSEA") +
    theme_bw(base_size = 10)
  ggsave(file.path(FIGS, "05_GSEA_GO_BP.pdf"), p_go, width = 14, height = 10)
  message("  Saved: 05_GSEA_GO_BP.pdf")
}

# ORA bar charts
plot_ora_bar <- function(obj, title) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0) { message("  No ORA hits: ", title); return(NULL) }
  df %>%
    slice_min(p.adjust, n = 20, with_ties = FALSE) %>%
    mutate(Description = str_wrap(Description, 45)) %>%
    ggplot(aes(reorder(Description, -log10(p.adjust)),
               -log10(p.adjust), fill = -log10(p.adjust))) +
    geom_col() +
    scale_fill_viridis_c(option = "magma") +
    coord_flip() +
    labs(title = title, x = NULL, y = "-log\u2081\u2080(adj. p-value)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
}

for (lst in list(
  list(ora_up_kegg,  "ORA KEGG: Up-regulated in NAD",   "06_ORA_KEGG_up.pdf"),
  list(ora_dn_kegg,  "ORA KEGG: Down-regulated in NAD",  "07_ORA_KEGG_down.pdf"),
  list(ora_up_react, "ORA Reactome: Up-regulated in NAD","08_ORA_Reactome_up.pdf"),
  list(ora_dn_react, "ORA Reactome: Down in NAD",        "09_ORA_Reactome_down.pdf"),
  list(ora_up_go,    "ORA GO-BP: Up-regulated in NAD",   "10_ORA_GO_up.pdf"),
  list(ora_dn_go,    "ORA GO-BP: Down-regulated in NAD", "11_ORA_GO_down.pdf")
)) {
  p <- plot_ora_bar(lst[[1]], lst[[2]])
  if (!is.null(p)) ggsave(file.path(FIGS, lst[[3]]), p, width = 10, height = 7)
}
message("  Saved: 06-11 ORA plots")

# ── 11. FOCUSED GENE SETS ────────────────────────────────────────────────── #
# ---------------------------------------------------------------------------
# GLUTAMATE / GLUTAMINE METABOLISM
# (Glutamine is the primary nitrogen donor for de novo purine & pyrimidine
#  synthesis; NAD+ directly influences glutamine-dependent PPAT and CAD)
# ---------------------------------------------------------------------------
glut_genes <- c(
  # Glutaminases (Gln → Glu)
  "GLS", "GLS2",
  # Glutamate dehydrogenases
  "GLUD1", "GLUD2",
  # Glutamine synthetase
  "GLUL",
  # Aminotransferases
  "GOT1", "GOT2", "GPT", "GPT2", "PSAT1", "PSPH",
  # SLC transporters (glutamate / glutamine uptake)
  "SLC1A1","SLC1A2","SLC1A3","SLC1A4","SLC1A5","SLC1A6","SLC1A7",
  "SLC38A1","SLC38A2","SLC38A5",
  # TCA anaplerosis / alpha-ketoglutarate
  "IDH1","IDH2","IDH3A","IDH3B","IDH3G",
  "OGDH","OGDHL","SUCLA2","SUCLG1","SUCLG2",
  # Asparagine / aspartate interconversion
  "ASNS", "ASS1", "ASL",
  # Proline synthesis (from glutamate)
  "PYCR1","PYCR2","PYCR3","ALDH18A1",
  # CPS1 / CPS2 (urea cycle & pyrimidine de novo, Gln-dependent)
  "CPS1",
  # Nitrogen metabolism
  "PPAT"  # also in purine de novo – Gln is substrate
)

# ---------------------------------------------------------------------------
# PURINE DE NOVO SYNTHESIS  (10-step IMP synthesis from PRPP + Gln)
# ---------------------------------------------------------------------------
purine_denovo <- c(
  "PPAT",          # Step  1: PRPP + Gln → PRA        (rate-limiting)
  "GART",          # Steps 2/3/5: trifunctional enzyme
  "PFAS",          # Step  4: FGAM synthase (Gln-dependent)
  "PAICS",         # Steps 6/7: bifunctional
  "ADSS","ADSL",   # Step  8: adenylosuccinate synthase/lyase
  "ATIC",          # Steps 9/10: AICAR transformylase / IMP cyclohydrolase
  "IMPDH1","IMPDH2",# IMP → XMP
  "GMPS",          # XMP → GMP (Gln-dependent)
  "ADSS1","ADSS2", # IMP → AMP branch
  "PFAS",
  "PRPS1","PRPS2"  # PRPP synthetase (makes substrate for step 1)
)

# ---------------------------------------------------------------------------
# PURINE SALVAGE PATHWAY
# ---------------------------------------------------------------------------
purine_salvage <- c(
  "HPRT1",                    # Hypoxanthine/guanine + PRPP → IMP/GMP
  "APRT",                     # Adenine + PRPP → AMP
  "ADA",                      # Adenosine deaminase
  "ADK",                      # Adenosine kinase
  "PNP",                      # Purine nucleoside phosphorylase
  "DGUOK",                    # Deoxyguanosine kinase
  "RRM1","RRM2","RRM2B",      # Ribonucleotide reductase
  "NT5E",                     # CD73: AMP → adenosine (extracellular)
  "NT5C1A","NT5C2","NT5C3A",  # 5'-nucleotidases
  "ENTPD1","ENTPD2","ENTPD3", # CD39 family NTPDases
  "AK1","AK2","AK3","AK4",    # Adenylate kinases
  "GUK1"                      # Guanylate kinase
)

# ---------------------------------------------------------------------------
# PYRIMIDINE DE NOVO SYNTHESIS
# ---------------------------------------------------------------------------
pyrim_denovo <- c(
  "CAD",    # Trifunctional: CPS2 (Gln-dep.) + ATCase + DHOase (steps 1-3)
  "DHODH",  # Dihydroorotate dehydrogenase (step 4, mitochondrial)
  "UMPS",   # Bifunctional OPRT + ODCase (steps 5-6) → UMP
  "CTPS1","CTPS2"  # UTP → CTP synthase (step 7, Gln-dependent)
)

# ---------------------------------------------------------------------------
# PYRIMIDINE SALVAGE PATHWAY
# ---------------------------------------------------------------------------
pyrim_salvage <- c(
  "TK1","TK2",           # Thymidine kinase (cytosolic / mitochondrial)
  "TYMS",                # Thymidylate synthase: dUMP → dTMP
  "TYMP",                # Thymidine phosphorylase
  "UCK1","UCK2",         # Uridine-cytidine kinase
  "CMPK1","CMPK2",       # UMP-CMP kinase
  "DCTD",                # dCMP deaminase
  "UPP1","UPP2",         # Uridine phosphorylase
  "DPYS","DPYD",         # Dihydropyrimidinase / dihydropyrimidine dehydrogenase
  "NME1","NME2",         # Nucleoside diphosphate kinase
  "UNG","SMUG1","MBD4"   # Uracil-DNA glycosylase (uracil removal from DNA)
)

# ---------------------------------------------------------------------------
# NAD+ METABOLISM  (relevant because treatment IS NAD+)
# ---------------------------------------------------------------------------
nad_genes <- c(
  # Biosynthesis (salvage from nicotinamide — main route in T cells)
  "NAMPT",                         # Nicotinamide → NMN (rate-limiting)
  "NMNAT1","NMNAT2","NMNAT3",      # NMN → NAD+
  "NADK","NADK2",                  # NAD+ → NADP+
  # Biosynthesis from nicotinic acid (Preiss-Handler)
  "NAPRT","NADSYN1",
  # De novo from tryptophan (kynurenine pathway)
  "IDO1","IDO2","TDO2","KYNU","HAAO","ACMSD","QPRT",
  # NAD+ consumers / ADP-ribosylation
  "CD38","BST1",                   # NAD+ hydrolases / cyclases
  "PARP1","PARP2","PARP3",         # Poly-ADP-ribose polymerases
  "SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7",  # Sirtuins
  "NNMT",                          # Nicotinamide N-methyltransferase
  "NUDT12","NUDT5"                 # NAD+ hydrolases
)

focused_sets <- list(
  Glutamate_Glutamine   = unique(glut_genes),
  Purine_De_Novo        = unique(purine_denovo),
  Purine_Salvage        = unique(purine_salvage),
  Pyrimidine_De_Novo    = unique(pyrim_denovo),
  Pyrimidine_Salvage    = unique(pyrim_salvage),
  NAD_Metabolism        = unique(nad_genes)
)

# ── 12. FOCUSED GSEA ─────────────────────────────────────────────────────── #
message("Running focused pathway GSEA...")
gsea_focused <- fgsea(
  pathways    = focused_sets,
  stats       = ranked_sym,
  minSize     = 5,
  maxSize     = 600,
  nPermSimple = 100000
) %>%
  as_tibble() %>%
  mutate(
    direction       = ifelse(NES > 0, "Activated (NAD)", "Suppressed (NAD)"),
    leadingEdge_str = sapply(leadingEdge, paste, collapse = ";")
  ) %>%
  arrange(padj)

message("\nFocused pathway GSEA results:")
print(gsea_focused %>% select(pathway, size, NES, pval, padj, direction))

# Merge DE stats for each focused set gene
focused_gene_table <- purrr::map2_dfr(focused_sets, names(focused_sets),
  function(genes, nm) {
    de %>%
      filter(Symbol %in% genes) %>%
      mutate(pathway = nm) %>%
      select(pathway, Symbol, GeneID, log2FoldChange, pvalue, padj, sig,
             all_of(COUNT_COLS))
  }
)

# Save focused results
wb2 <- createWorkbook()
addWorksheet(wb2, "GSEA_Focused_Summary")
writeData(wb2, "GSEA_Focused_Summary",
          gsea_focused %>% select(-leadingEdge))

for (nm in names(focused_sets)) {
  dat <- focused_gene_table %>% filter(pathway == nm)
  if (nrow(dat) > 0) {
    addWorksheet(wb2, nm)
    writeData(wb2, nm, dat)
  }
}
saveWorkbook(wb2, file.path(RES, "02_Focused_pathway_results.xlsx"), overwrite = TRUE)
message("  Saved: 02_Focused_pathway_results.xlsx")

# ── 13. GSEA ENRICHMENT CURVES ───────────────────────────────────────────── #
message("Plotting GSEA enrichment curves...")

plot_enrichment_curve <- function(set_name, sets, ranked, filename) {
  if (!set_name %in% names(sets)) {
    message("  Not found: ", set_name); return(invisible(NULL))
  }
  p <- plotEnrichment(sets[[set_name]], ranked) +
    labs(
      title    = str_replace_all(set_name, "_", " "),
      subtitle = sprintf("Genes in set: %d  |  Genes with stats: %d",
                         length(sets[[set_name]]),
                         sum(sets[[set_name]] %in% names(ranked))),
      x = "Gene rank (high = activated in NAD)",
      y = "Enrichment score"
    ) +
    theme_bw(base_size = 12)
  ggsave(filename, p, width = 8, height = 5)
  message("  Saved: ", basename(filename))
}

# Focused sets
for (nm in names(focused_sets)) {
  plot_enrichment_curve(
    nm, focused_sets, ranked_sym,
    file.path(FIGS, paste0("12_GSEA_curve_", nm, ".pdf"))
  )
}

# KEGG canonical pathways (from MSigDB)
kegg_focused <- c(
  "KEGG_PURINE_METABOLISM",
  "KEGG_PYRIMIDINE_METABOLISM",
  "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM",
  "KEGG_NICOTINATE_AND_NICOTINAMIDE_METABOLISM",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_TCA_CYCLE"
)
for (kp in kegg_focused) {
  plot_enrichment_curve(
    kp, kegg_sets, ranked_sym,
    file.path(FIGS, paste0("13_GSEA_", kp, ".pdf"))
  )
}

# ── 14. HEATMAPS FOR FOCUSED PATHWAYS ────────────────────────────────────── #
message("Plotting heatmaps...")

meta_annot <- data.frame(
  row.names = COUNT_COLS,
  Treatment = c(rep("Control", 3), rep("NAD", 3))
)
annot_col_colors <- list(
  Treatment = c(Control = "#377EB8", NAD = "#E41A1C")
)

plot_focused_heatmap <- function(gene_list, de_df, count_cols, meta,
                                  annot_colors, title, filename) {
  sub <- de_df %>%
    filter(Symbol %in% gene_list, !is.na(Symbol),
           rowSums(across(all_of(count_cols))) > 0) %>%
    distinct(Symbol, .keep_all = TRUE)

  if (nrow(sub) < 3) {
    message("  Too few genes for heatmap: ", title); return(invisible(NULL))
  }

  mat <- as.matrix(sub[, count_cols])
  rownames(mat) <- sub$Symbol

  # Z-score per row, clamp to ±2.5
  mat_z <- t(scale(t(mat)))
  mat_z[is.nan(mat_z)] <- 0
  mat_z <- pmax(pmin(mat_z, 2.5), -2.5)

  # Row annotation: direction
  row_annot <- data.frame(
    row.names = sub$Symbol,
    Direction = factor(sub$sig, levels = c("Up (NAD)", "Down (NAD)", "NS"))
  )
  annot_row_colors <- list(
    Direction = c("Up (NAD)" = "#E41A1C", "Down (NAD)" = "#377EB8", NS = "grey85")
  )

  pdf(filename, width = 9, height = max(6, nrow(mat_z) * 0.28 + 3))
  pheatmap(
    mat_z,
    annotation_col    = meta,
    annotation_row    = row_annot,
    annotation_colors = c(annot_colors, annot_row_colors),
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 9,
    border_color      = NA,
    main              = title
  )
  dev.off()
  message("  Saved: ", basename(filename))
}

plot_focused_heatmap(glut_genes, de, COUNT_COLS, meta_annot, annot_col_colors,
  "Glutamate / Glutamine Metabolism",
  file.path(FIGS, "14_heatmap_glutamate_glutamine.pdf"))

plot_focused_heatmap(c(purine_denovo, purine_salvage), de, COUNT_COLS, meta_annot, annot_col_colors,
  "Purine Metabolism: De Novo + Salvage",
  file.path(FIGS, "15_heatmap_purine.pdf"))

plot_focused_heatmap(c(pyrim_denovo, pyrim_salvage), de, COUNT_COLS, meta_annot, annot_col_colors,
  "Pyrimidine Metabolism: De Novo + Salvage",
  file.path(FIGS, "16_heatmap_pyrimidine.pdf"))

plot_focused_heatmap(nad_genes, de, COUNT_COLS, meta_annot, annot_col_colors,
  "NAD+ Metabolism",
  file.path(FIGS, "17_heatmap_NAD.pdf"))

# ── 15. FOCUSED LOLLIPOP SUMMARY ─────────────────────────────────────────── #
message("Plotting focused pathway lollipop summary...")

p_lollipop <- gsea_focused %>%
  mutate(
    label   = str_replace_all(pathway, "_", " "),
    sig_tag = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                        padj < 0.05  ~ "*",   TRUE         ~ "")
  ) %>%
  ggplot(aes(NES, reorder(label, NES),
             colour = direction, size = -log10(pval))) +
  geom_segment(aes(x = 0, xend = NES, yend = label),
               linewidth = 0.9, colour = "grey70") +
  geom_point() +
  geom_text(aes(label = sig_tag),
            hjust = ifelse(gsea_focused$NES > 0, -0.4, 1.4),
            size = 5, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(
    values = c("Activated (NAD)" = "#E41A1C", "Suppressed (NAD)" = "#377EB8")
  ) +
  scale_size_continuous(range = c(4, 11), name = "-log\u2081\u2080(p)") +
  labs(
    title    = "Focused Pathway GSEA: NAD-treated vs Control T cells",
    subtitle = "Significance: * p<0.05  ** p<0.01  *** p<0.001",
    x = "Normalized Enrichment Score (NES)  \u2192  Activated in NAD",
    y = NULL, colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 11))

ggsave(file.path(FIGS, "18_focused_lollipop_GSEA.pdf"),
       p_lollipop, width = 10, height = 6)
message("  Saved: 18_focused_lollipop_GSEA.pdf")

# ── 16. PATHVIEW (KEGG PATHWAY MAPS WITH LFC OVERLAY) ───────────────────── #
message("Generating Pathview KEGG maps...")

lfc_entrez <- setNames(de$log2FoldChange, as.character(de$GeneID))
lfc_entrez <- lfc_entrez[!is.na(names(lfc_entrez)) & names(lfc_entrez) != ""]

kegg_map_ids <- c(
  "hsa00230",  # Purine metabolism
  "hsa00240",  # Pyrimidine metabolism
  "hsa00250",  # Alanine, aspartate & glutamate metabolism
  "hsa00760",  # Nicotinate and nicotinamide metabolism (NAD+)
  "hsa00900",  # Terpenoid backbone biosynthesis (for context)
  "hsa04064",  # NF-kB signaling (immune context)
  "hsa04660"   # T cell receptor signaling pathway
)

old_wd <- getwd()
setwd(FIGS)

for (kid in kegg_map_ids) {
  tryCatch({
    pathview(
      gene.data   = lfc_entrez,
      pathway.id  = kid,
      species     = "hsa",
      gene.idtype = "entrez",
      limit       = list(gene = 2, cpd = 1),
      bins        = list(gene = 10, cpd = 10),
      low         = list(gene = "blue"),
      high        = list(gene = "red"),
      out.suffix  = "NAD_vs_Control",
      kegg.native = TRUE
    )
    message("  Pathview saved: ", kid)
  }, error = function(e)
    message("  Pathview error for ", kid, ": ", e$message))
}

setwd(old_wd)

# ── 17. DAR INTEGRATION STUB ─────────────────────────────────────────────── #
# When ATAC-seq DAR results are available, place the file in data/ and
# uncomment + adapt the block below:
#
# dar <- read.delim(file.path(DATA, "DAR_results.txt"))
#   # Expected columns: peak_id, nearest_gene, log2FC_ATAC, padj_ATAC
#
# # Overlap DEGs with DARs at nearest gene level
# deg_dar <- de %>%
#   inner_join(dar %>% select(nearest_gene, log2FC_ATAC, padj_ATAC),
#              by = c("Symbol" = "nearest_gene")) %>%
#   mutate(
#     concordant = sign(log2FoldChange) == sign(log2FC_ATAC),
#     category   = case_when(
#       sig != "NS" & padj_ATAC < 0.05 & concordant  ~ "DEG+DAR concordant",
#       sig != "NS" & padj_ATAC < 0.05 & !concordant ~ "DEG+DAR discordant",
#       sig != "NS"                                   ~ "DEG only",
#       padj_ATAC < 0.05                              ~ "DAR only",
#       TRUE                                          ~ "NS"
#     )
#   )
#
# # Focused pathway overlap
# focused_dar <- deg_dar %>%
#   filter(Symbol %in% unlist(focused_sets)) %>%
#   mutate(pathway = purrr::map_chr(Symbol, function(g) {
#     hit <- purrr::keep(focused_sets, ~ g %in% .)
#     if (length(hit) > 0) names(hit)[1] else NA_character_
#   }))

# ── 18. SESSION INFO ─────────────────────────────────────────────────────── #
message("\n========================================")
message("  Analysis complete!")
message("  Results : ", RES)
message("  Figures : ", FIGS)
message("========================================\n")

sink(file.path(RES, "session_info.txt"))
print(sessionInfo())
sink()
