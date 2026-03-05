################################################################################
##  09_functional_enrichment.R                                                 ##
##                                                                             ##
##  Input:  DEGs_SingleCell.RData  (da sezione 07)                            ##
##  Output: outdir/functional/ con Excel, CSV e PNG                           ##
##                                                                             ##
##  Analisi:                                                                   ##
##    - GSEA (fgsea): GO:BP, Hallmark, Reactome — per ogni CT × contrasto    ##
##    - ORA  (clusterProfiler): GO:BP, KEGG                                   ##
##    - Heatmap NES cross-CT per contrasto                                    ##
##    - Heatmap cross-contrasto (CTX_90 vs CTX_24 vs THAL_90 vs THAL_24)    ##
##    - Dotplot e enrichment plot per pathway di interesse LSD                ##
################################################################################

source("code/00_config.R")
outdir<-"results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "DEGs_SingleCell.RData"))

# ==============================================================================
# SETUP OUTPUT
# ==============================================================================

func_dir <- file.path(outdir, "functional")
dir.create(func_dir, recursive = TRUE, showWarnings = FALSE)

savefunc <- function(p, filename, width = 3500, height = 2800, res = 300) {
  png(file.path(func_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# ==============================================================================
# GENE SETS
# ==============================================================================

message("Caricamento gene sets MSigDB (Mus musculus)...")

gmt_gobp <- msigdbr(species = "Mus musculus", category = "C5",
                    subcategory = "GO:BP")
gs_gobp  <- split(gmt_gobp$gene_symbol, gmt_gobp$gs_name)

gmt_hall <- msigdbr(species = "Mus musculus", category = "H")
gs_hall  <- split(gmt_hall$gene_symbol, gmt_hall$gs_name)

gmt_reac <- msigdbr(species = "Mus musculus", category = "C2",
                    subcategory = "CP:REACTOME")
gs_reac  <- split(gmt_reac$gene_symbol, gmt_reac$gs_name)

message(sprintf("  GO:BP: %d | Hallmark: %d | Reactome: %d",
                length(gs_gobp), length(gs_hall), length(gs_reac)))

# Conversione symbol → Entrez (per KEGG e ORA)
symbol_to_entrez <- function(genes) {
  ids <- suppressMessages(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(genes),
                          column = "ENTREZID", keytype = "SYMBOL",
                          multiVals = "first")
  )
  ids[!is.na(ids)]
}

# ==============================================================================
# INPUT DEG
# ==============================================================================

deg_inputs <- list(
  CTX_90min  = DEGs_SC_CTX_90,
  CTX_24h    = DEGs_SC_CTX_24
)

# ==============================================================================
# 09.1 — GSEA (fgsea) — per ogni CT × contrasto × database
# ==============================================================================

message("=== 09.1 GSEA ===")

all_gsea_gobp <- list()
all_gsea_hall <- list()
all_gsea_reac <- list()

for (contrast_nm in names(deg_inputs)) {
  degs_list <- deg_inputs[[contrast_nm]]
  dir.create(file.path(func_dir, contrast_nm), showWarnings = FALSE)

  for (ct in names(degs_list)) {
    df  <- degs_list[[ct]]
    lbl <- paste(contrast_nm, ct, sep = "_")

    df_valid <- df[!is.na(df$log2FoldChange), ]
    if (nrow(df_valid) < 15) next

    # Ranked list: LFC shrinkato
    ranks <- setNames(df_valid$log2FoldChange, df_valid$gene)
    ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
    if (length(ranks) < 10) next

    run_gsea_db <- function(gene_sets, db_name) {
      res <- tryCatch(
        fgsea(pathways    = gene_sets,
              stats       = ranks,
              minSize     = 5,
              maxSize     = 500,
              eps         = 0,
              nPermSimple = 10000),
        error = function(e) {
          message(sprintf("  GSEA errore %s/%s: %s", lbl, db_name, e$message))
          NULL
        }
      )
      if (is.null(res) || nrow(res) == 0) return(NULL)

      # Soglia adattiva
      n_sig <- sum(!is.na(res$padj) & res$padj < 0.05)
      p_use <- 0.05
      if (n_sig < 3) {
        for (pb in c(0.10, 0.20)) {
          if (sum(!is.na(res$pval) & res$pval < pb) >= 3) {
            p_use <- pb; break
          }
        }
      }
      res$p_used  <- p_use
      res$is_sig  <- !is.na(res$pval) & res$pval < p_use
      res$label   <- lbl
      res
    }

    all_gsea_gobp[[lbl]] <- run_gsea_db(gs_gobp, "GO:BP")
    all_gsea_hall[[lbl]] <- run_gsea_db(gs_hall, "Hallmark")
    all_gsea_reac[[lbl]] <- run_gsea_db(gs_reac, "Reactome")

    n_gobp <- if (!is.null(all_gsea_gobp[[lbl]]))
      sum(all_gsea_gobp[[lbl]]$is_sig) else 0
    n_hall <- if (!is.null(all_gsea_hall[[lbl]]))
      sum(all_gsea_hall[[lbl]]$is_sig) else 0
    message(sprintf("  %s | %-20s GO:BP sig=%d | Hall sig=%d",
                    contrast_nm, ct, n_gobp, n_hall))
  }
}

# ==============================================================================
# 09.2 — ORA (clusterProfiler) — GO:BP + KEGG
# ==============================================================================

message("=== 09.2 ORA ===")

all_ora <- list()

for (contrast_nm in names(deg_inputs)) {
  degs_list <- deg_inputs[[contrast_nm]]

  for (ct in names(degs_list)) {
    df  <- degs_list[[ct]]
    lbl <- paste(contrast_nm, ct, sep = "_")

    df_valid <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    sig_genes <- df_valid$gene[df_valid$padj < 0.05]
    universe  <- df_valid$gene

    if (length(sig_genes) < 5) next

    entrez_sig  <- symbol_to_entrez(sig_genes)
    entrez_univ <- symbol_to_entrez(universe)
    if (length(entrez_sig) < 5) next

    ora_gobp <- tryCatch(
      enrichGO(gene          = entrez_sig,
               universe      = entrez_univ,
               OrgDb         = org.Mm.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.1,
               qvalueCutoff  = 0.2,
               readable      = TRUE),
      error = function(e) NULL
    )

    ora_kegg <- tryCatch(
      enrichKEGG(gene          = entrez_sig,
                 universe      = entrez_univ,
                 organism      = "mmu",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.1),
      error = function(e) NULL
    )

    all_ora[[lbl]] <- list(gobp = ora_gobp, kegg = ora_kegg, label = lbl)

    n_bp   <- if (!is.null(ora_gobp)) nrow(ora_gobp) else 0
    n_kegg <- if (!is.null(ora_kegg)) nrow(ora_kegg) else 0
    message(sprintf("  %s | %-20s GO:BP=%d KEGG=%d",
                    contrast_nm, ct, n_bp, n_kegg))
  }
}

# Salva oggetti raw
save(all_gsea_gobp, all_gsea_hall, all_gsea_reac, all_ora,
     file = file.path(func_dir, "functional_results_raw.RData"))
message("Salvato: functional_results_raw.RData")

# ==============================================================================
# 09.3 — EXPORT EXCEL GSEA
# ==============================================================================

message("=== 09.3 Export Excel ===")

flatten_fgsea <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  df <- as.data.frame(res_df)
  if ("leadingEdge" %in% colnames(df))
    df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = ";")
  df
}

for (contrast_nm in names(deg_inputs)) {
  wb <- createWorkbook()

  for (ct in names(deg_inputs[[contrast_nm]])) {
    lbl      <- paste(contrast_nm, ct, sep = "_")
    sheet_nm <- substr(ct, 1, 31)

    rows <- do.call(rbind, list(
      if (!is.null(all_gsea_gobp[[lbl]]))
        cbind(db = "GO:BP",    flatten_fgsea(all_gsea_gobp[[lbl]])),
      if (!is.null(all_gsea_hall[[lbl]]))
        cbind(db = "Hallmark", flatten_fgsea(all_gsea_hall[[lbl]])),
      if (!is.null(all_gsea_reac[[lbl]]))
        cbind(db = "Reactome", flatten_fgsea(all_gsea_reac[[lbl]]))
    ))
    if (is.null(rows) || nrow(rows) == 0) next

    rows <- rows[order(rows$pval), ]
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, rows)

    # Evidenzia righe significative (padj < 0.05)
    sig_rows <- which(!is.na(rows$padj) & rows$padj < 0.05) + 1
    if (length(sig_rows) > 0)
      addStyle(wb, sheet_nm,
               style      = createStyle(fgFill = "#FFF2CC"),
               rows       = sig_rows,
               cols       = seq_len(ncol(rows)),
               gridExpand = TRUE)
  }

  saveWorkbook(wb,
               file.path(func_dir, paste0("GSEA_", contrast_nm, ".xlsx")),
               overwrite = TRUE)
}
message("Excel GSEA salvati.")

# ==============================================================================
# 09.4 — HEATMAP NES cross-CT (Hallmark + GO:BP)
# ==============================================================================

message("=== 09.4 Heatmap NES per contrasto ===")

plot_nes_heatmap <- function(gsea_list, contrast_nm, db_label,
                             p_thr = 0.20, top_n = 30) {

  ct_names <- names(deg_inputs[[contrast_nm]])
  res_list <- lapply(ct_names, function(ct) {
    lbl <- paste(contrast_nm, ct, sep = "_")
    res <- gsea_list[[lbl]]
    if (is.null(res) || nrow(res) == 0) return(NULL)
    df <- as.data.frame(res)
    df$cell_type <- ct
    df
  })
  df_all <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(df_all) || nrow(df_all) == 0) return(invisible(NULL))

  # Pathway significativi in almeno un CT
  sig_paths <- df_all[!is.na(df_all$pval) & df_all$pval < p_thr, "pathway"]
  if (length(sig_paths) == 0) return(invisible(NULL))

  # Top per |NES| massimo
  path_rank <- df_all[df_all$pathway %in% sig_paths, ] %>%
    group_by(pathway) %>%
    summarise(max_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_nes)) %>%
    slice_head(n = top_n)
  sig_paths <- path_rank$pathway

  # Matrice NES (pathway × CT)
  nes_mat <- df_all[df_all$pathway %in% sig_paths, ] %>%
    select(pathway, cell_type, NES) %>%
    pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  # Matrice p-value per asterischi
  pval_mat <- df_all[df_all$pathway %in% sig_paths, ] %>%
    select(pathway, cell_type, pval) %>%
    pivot_wider(names_from = cell_type, values_from = pval, values_fill = 1) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  pval_mat <- pval_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]

  sig_lab <- matrix(nrow=nrow(pval_mat), ncol=ncol(pval_mat),
    ifelse(pval_mat < 0.01,  "**",
    ifelse(pval_mat < 0.05,  "*",
    ifelse(pval_mat < 0.20,  "·", ""))),
    dimnames = dimnames(nes_mat)
  )

  rownames(nes_mat) <- str_wrap(
    gsub("HALLMARK_|GOBP_|REACTOME_", "", rownames(nes_mat)), 45
  )
  rownames(sig_lab) <- rownames(nes_mat)

  lim     <- max(abs(nes_mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))

  ht <- Heatmap(
    nes_mat,
    name              = "NES",
    col               = col_fun,
    cell_fun          = function(j, i, x, y, w, h, fill)
      grid.text(sig_lab[i, j], x, y, gp = gpar(fontsize = 9)),
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 8),
    column_names_rot  = 45,
    column_title      = paste("GSEA NES —", db_label, "|", contrast_nm),
    heatmap_legend_param = list(title = "NES")
  )

  fn <- paste0("GSEA_NES_heatmap_", db_label, "_", contrast_nm, ".png")
  png(file.path(func_dir, fn),
      width  = max(3000, ncol(nes_mat) * 250 + 1500),
      height = max(3000, nrow(nes_mat) * 60  + 800),
      res    = 300)
  draw(ht)
  dev.off()
  message(sprintf("  Salvato: %s", fn))
}

for (contrast_nm in names(deg_inputs)) {
  plot_nes_heatmap(all_gsea_hall, contrast_nm, "Hallmark", top_n = 30)
  plot_nes_heatmap(all_gsea_gobp, contrast_nm, "GOBP",     top_n = 40)
  plot_nes_heatmap(all_gsea_reac, contrast_nm, "Reactome", top_n = 30)
}

# ==============================================================================
# 09.5 — HEATMAP CROSS-CONTRASTO
# ==============================================================================

message("=== 09.5 Heatmap cross-contrasto ===")

build_cross_contrast_heatmap <- function(gsea_list, db_label,
                                          p_thr = 0.20, top_n = 40) {

  df_all <- do.call(rbind, lapply(names(deg_inputs), function(cnt) {
    do.call(rbind, lapply(names(deg_inputs[[cnt]]), function(ct) {
      lbl <- paste(cnt, ct, sep = "_")
      res <- gsea_list[[lbl]]
      if (is.null(res) || nrow(res) == 0) return(NULL)
      df <- as.data.frame(res)
      df$contrast  <- cnt
      df$cell_type <- ct
      df$group     <- paste(cnt, ct, sep = "|")
      df
    }))
  }))
  if (is.null(df_all) || nrow(df_all) == 0) return(invisible(NULL))

  top_paths <- df_all[!is.na(df_all$pval) & df_all$pval < p_thr, ] %>%
    group_by(pathway) %>%
    summarise(score = max(abs(NES), na.rm=TRUE) * n(), .groups="drop") %>%
    arrange(desc(score)) %>%
    slice_head(n = top_n) %>%
    pull(pathway)
  if (length(top_paths) == 0) return(invisible(NULL))

  nes_mat <- df_all[df_all$pathway %in% top_paths, ] %>%
    select(pathway, group, NES) %>%
    pivot_wider(names_from = group, values_from = NES, values_fill = 0) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  rownames(nes_mat) <- str_wrap(
    gsub("HALLMARK_|GOBP_|REACTOME_", "", rownames(nes_mat)), 40
  )

  lim     <- max(abs(nes_mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))

  # Annotazione colonne per contrasto
  col_groups <- data.frame(
    Contrast = sub("\\|.*$", "", colnames(nes_mat)),
    row.names = colnames(nes_mat)
  )
  pal_contrasts <- setNames(
    brewer.pal(max(3, length(unique(col_groups$Contrast))), "Set2"),
    unique(col_groups$Contrast)
  )
  col_ha <- HeatmapAnnotation(
    df  = col_groups,
    col = list(Contrast = pal_contrasts[unique(col_groups$Contrast)])
  )

  ht <- Heatmap(
    nes_mat, name = "NES", col = col_fun,
    top_annotation   = col_ha,
    cluster_rows     = TRUE,
    cluster_columns  = TRUE,
    row_names_gp     = gpar(fontsize = 6.5),
    column_names_gp  = gpar(fontsize = 7),
    column_names_rot = 45,
    column_title     = paste("Cross-contrasto —", db_label)
  )

  fn <- paste0("GSEA_crosscontrast_", db_label, ".png")
  png(file.path(func_dir, fn),
      width  = max(4500, ncol(nes_mat) * 160 + 1500),
      height = max(3500, nrow(nes_mat) * 55  + 800),
      res    = 300)
  draw(ht)
  dev.off()
  message(sprintf("  Salvato: %s", fn))
}

build_cross_contrast_heatmap(all_gsea_hall, "Hallmark")
build_cross_contrast_heatmap(all_gsea_gobp, "GOBP", top_n = 50)

# ==============================================================================
# 09.6 — ENRICHMENT PLOT per pathway di interesse biologico
# ==============================================================================

message("=== 09.6 Enrichment plot pathway LSD ===")

pathways_focus <- c(
  "GOBP_IMMEDIATE_EARLY_GENE_ACTIVATION",
  "GOBP_SEROTONIN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_SYNAPTIC_PLASTICITY",
  "GOBP_LONG_TERM_POTENTIATION",
  "GOBP_REGULATION_OF_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_THALAMUS_DEVELOPMENT",
  "GOBP_NEURON_PROJECTION_DEVELOPMENT"
)

gs_all <- c(gs_gobp, gs_hall, gs_reac)

for (contrast_nm in names(deg_inputs)) {
  for (ct in names(deg_inputs[[contrast_nm]])) {
    df_ct  <- deg_inputs[[contrast_nm]][[ct]]
    df_v   <- df_ct[!is.na(df_ct$log2FoldChange), ]
    if (nrow(df_v) < 10) next

    ranks <- setNames(df_v$log2FoldChange, df_v$gene)
    ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)

    for (pw in pathways_focus) {
      if (!pw %in% names(gs_all)) next
      if (length(intersect(gs_all[[pw]], names(ranks))) < 5) next

      tryCatch({
        p_enr <- plotEnrichment(gs_all[[pw]], stats = ranks) +
          labs(title    = paste(gsub("GOBP_|HALLMARK_|REACTOME_","",pw),
                                "\n", contrast_nm, "|", ct),
               x = "Rank", y = "ES") +
          theme_bw(base_size = 9)
        fn <- file.path(func_dir, contrast_nm,
                        paste0("enrichplot_", ct, "_",
                               gsub("[^A-Za-z0-9]","_",pw), ".png"))
        png(fn, width = 3000, height = 2000, res = 300)
        print(p_enr)
        dev.off()
      }, error = function(e) NULL)
    }
  }
}

message(sprintf("=== 09 completato — output in: %s ===", func_dir))


# ==============================================================================
# 09b — Frequenza pathway GO:BP significativi per timepoint
# ==============================================================================
# Input:  all_gsea_gobp  (lista fgsea results, nomi = "CTX_90min_CT" ecc.)
# Output: 09b_gobp_freq_90min.png
#         09b_gobp_freq_24h.png
#
# Per ogni timepoint:
#   - conta in quanti CT un pathway e' significativo (pval < P_THR)
#   - prende i top TOP_N per frequenza (in caso di parita', ordina per NES medio)
#   - rappresenta come barplot orizzontale:
#       barra = n CT in cui e' significativo
#       colore = NES medio (blu negativo, rosso positivo)
#       testo  = n CT su totale testato
# ==============================================================================

P_THR  <- 0.05   # soglia p-value per "significativo"
TOP_N  <- 30     # pathway da mostrare per timepoint

# ------------------------------------------------------------------------------
# 1. Appiattisci all_gsea_gobp in un unico data.frame con colonne
#    contrast, cell_type, pathway, NES, pval, padj
# ------------------------------------------------------------------------------

gsea_flat <- do.call(rbind, lapply(names(all_gsea_gobp), function(lbl) {
  res <- all_gsea_gobp[[lbl]]
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # Estrai contrasto e CT dal nome (formato: "CTX_90min_ExL23_1" ecc.)
  # Il contrasto e' sempre il primo token con "_90min" o "_24h"
  contrast <- regmatches(lbl, regexpr("^[^_]+_(90min|24h)", lbl))
  ct       <- sub(paste0("^", contrast, "_"), "", lbl)
  
  df           <- as.data.frame(res)
  df$contrast  <- contrast
  df$cell_type <- ct
  df$timepoint <- ifelse(grepl("90min", contrast), "90min", "24h")
  df
}))

if (is.null(gsea_flat) || nrow(gsea_flat) == 0)
  stop("all_gsea_gobp e' vuoto o mal formato")

# Pulisci nome pathway per la visualizzazione
gsea_flat$pathway_label <- gsub("^GOBP_", "", gsea_flat$pathway)
gsea_flat$pathway_label <- gsub("_", " ", gsea_flat$pathway_label)
gsea_flat$pathway_label <- tolower(gsea_flat$pathway_label)
gsea_flat$pathway_label <- paste0(
  toupper(substr(gsea_flat$pathway_label, 1, 1)),
  substr(gsea_flat$pathway_label, 2, nchar(gsea_flat$pathway_label))
)

# ------------------------------------------------------------------------------
# 2. Funzione plot per singolo timepoint
# ------------------------------------------------------------------------------

plot_gobp_freq <- function(tp, filename) {
  
  df_tp <- gsea_flat[gsea_flat$timepoint == tp, ]
  if (nrow(df_tp) == 0) {
    message(sprintf("  Nessun dato per timepoint %s", tp))
    return(invisible(NULL))
  }
  
  # N CT testati in questo timepoint (denominatore per la label)
  n_ct_total <- length(unique(df_tp$cell_type))
  
  # Per ogni pathway: n CT significativi + NES medio tra i CT sig
  df_sig <- df_tp[!is.na(df_tp$pval) & df_tp$pval < P_THR, ]
  
  if (nrow(df_sig) == 0) {
    message(sprintf("  Nessun pathway significativo a p<%g per %s", P_THR, tp))
    return(invisible(NULL))
  }
  
  freq <- do.call(rbind, lapply(
    split(df_sig, df_sig$pathway),
    function(d) {
      data.frame(
        pathway       = d$pathway[1],
        pathway_label = d$pathway_label[1],
        n_sig         = nrow(d),
        nes_mean      = mean(d$NES, na.rm = TRUE),
        nes_sd        = sd(d$NES,   na.rm = TRUE),
        ct_list       = paste(sort(d$cell_type), collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  # Ordina: prima per n_sig decrescente, poi per |nes_mean| per i pareggi
  freq <- freq[order(-freq$n_sig, -abs(freq$nes_mean)), ]
  freq <- freq[seq_len(min(TOP_N, nrow(freq))), ]
  
  # Ordine asse y: dal meno frequente (in alto) al piu' frequente (in basso)
  # cosi' le barre piu' lunghe sono in fondo (convenzione lollipop/barplot)
  freq$pathway_label <- factor(freq$pathway_label,
                               levels = rev(freq$pathway_label))
  
  # Scala colore NES: simmetrica attorno a 0
  nes_lim <- max(abs(freq$nes_mean), na.rm = TRUE)
  nes_lim <- ceiling(nes_lim * 10) / 10   # arrotonda al decimo sopra
  
  p <- ggplot(freq,
              aes(x    = n_sig,
                  y    = pathway_label,
                  fill = nes_mean)) +
    
    # Barra
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    
    # Linea verticale: soglia "presente in tutti i CT"
    geom_vline(xintercept = n_ct_total,
               linetype  = "dotted",
               color     = "grey40",
               linewidth = 0.5) +
    
    # Etichetta n_sig / n_ct_total a destra della barra
    geom_text(aes(label = sprintf("%d/%d", n_sig, n_ct_total)),
              hjust  = -0.15,
              size   = 2.6,
              color  = "grey25") +
    
    scale_fill_gradient2(
      low      = "#2471A3",
      mid      = "grey92",
      high     = "#CB4335",
      midpoint = 0,
      limits   = c(-nes_lim, nes_lim),
      name     = "NES medio"
    ) +
    
    scale_x_continuous(
      breaks = seq(0, n_ct_total, by = max(1, n_ct_total %/% 5)),
      limits = c(0, n_ct_total * 1.22),   # spazio per le etichette
      expand = c(0, 0)
    ) +
    
    theme_bw(base_size = 10) +
    theme(
      axis.text.y      = element_text(size = 8, lineheight = 0.9),
      axis.text.x      = element_text(size = 8),
      axis.title.y     = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position  = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.title     = element_text(size = 8),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    ) +
    labs(
      title    = sprintf("GO:BP — pathway piu' frequentemente significativi | %s", tp),
      subtitle = sprintf(
        "p < %.2f (fgsea) | top %d pathway per frequenza tra %d CT | colore = NES medio tra CT sig.",
        P_THR, TOP_N, n_ct_total),
      x = sprintf("N cell type con pathway significativo (su %d)", n_ct_total)
    )
  
  savepng(p, filename,
          width  = 4200,
          height = max(2200, nrow(freq) * 75 + 700))
  message(sprintf("  Salvato: %s (%d pathway, %d CT)", filename, nrow(freq), n_ct_total))
  invisible(p)
}

# ------------------------------------------------------------------------------
# 3. Esegui per i due timepoint
# ------------------------------------------------------------------------------

message("=== 09b: frequenza GO:BP per timepoint ===")

plot_gobp_freq("90min", "09b_gobp_freq_90min.png")
plot_gobp_freq("24h",   "09b_gobp_freq_24h.png")

message("=== 09b completato ===")
