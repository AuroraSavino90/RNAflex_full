################################################################################
##  09_functional_enrichment_THAL.R                                            ##
##                                                                             ##
##  Input:  THAL/DEGs_SingleCell_THAL.RData  (da 07_deseq2_THAL)             ##
##  Output: THAL/functional/ con Excel, CSV e PNG                             ##
##                                                                             ##
##  Analisi:                                                                   ##
##    - GSEA (fgsea): GO:BP, Hallmark, Reactome — per ogni CT × contrasto    ##
##    - ORA  (clusterProfiler): GO:BP, KEGG                                   ##
##    - Heatmap NES cross-CT per contrasto                                    ##
##    - Heatmap cross-contrasto (THAL_90 vs THAL_24)                         ##
##    - Dotplot e enrichment plot per pathway di interesse LSD / talamo       ##
##    - Barplot frequenza GO:BP per timepoint (09b)                           ##
##                                                                             ##
##  Pathway di interesse biologico aggiuntivi (talamo):                       ##
##    GOBP_THALAMUS_DEVELOPMENT                                               ##
##    GOBP_THALAMOCORTICAL_PROJECTION_NEURON_DIFFERENTIATION                 ##
##    GOBP_REGULATION_OF_THALAMIC_GENE_EXPRESSION                            ##
##    GOBP_SENSORY_THALAMUS_DEVELOPMENT                                       ##
##    GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY                              ##
##    GOBP_GABA_RECEPTOR_SIGNALING (TRN)                                      ##
################################################################################

source("code/00_config.R")

outdir   <- "results/20260224"
thal_dir <- file.path(outdir, "THAL")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(thal_dir, "DEGs_SingleCell_THAL.RData"))
# Oggetti: DEGs_SC_THAL_90, DEGs_SC_THAL_24

# ==============================================================================
# SETUP OUTPUT
# ==============================================================================

func_dir <- file.path(thal_dir, "functional")
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

gmt_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
gs_gobp  <- split(gmt_gobp$gene_symbol, gmt_gobp$gs_name)

gmt_hall <- msigdbr(species = "Mus musculus", category = "H")
gs_hall  <- split(gmt_hall$gene_symbol, gmt_hall$gs_name)

gmt_reac <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
gs_reac  <- split(gmt_reac$gene_symbol, gmt_reac$gs_name)

message(sprintf("  GO:BP: %d | Hallmark: %d | Reactome: %d",
                length(gs_gobp), length(gs_hall), length(gs_reac)))

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
  THAL_90min = DEGs_SC_THAL_90,
  THAL_24h   = DEGs_SC_THAL_24
)

# ==============================================================================
# 09.1 — GSEA (fgsea) — per ogni CT × contrasto × database
# ==============================================================================

message("=== 09.1 GSEA THAL ===")

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

      n_sig <- sum(!is.na(res$padj) & res$padj < 0.05)
      p_use <- 0.05
      if (n_sig < 3) {
        for (pb in c(0.10, 0.20)) {
          if (sum(!is.na(res$pval) & res$pval < pb) >= 3) {
            p_use <- pb; break
          }
        }
      }
      res$p_used <- p_use
      res$is_sig <- !is.na(res$pval) & res$pval < p_use
      res$label  <- lbl
      res
    }

    all_gsea_gobp[[lbl]] <- run_gsea_db(gs_gobp, "GO:BP")
    all_gsea_hall[[lbl]] <- run_gsea_db(gs_hall, "Hallmark")
    all_gsea_reac[[lbl]] <- run_gsea_db(gs_reac, "Reactome")

    n_gobp <- if (!is.null(all_gsea_gobp[[lbl]])) sum(all_gsea_gobp[[lbl]]$is_sig) else 0
    n_hall <- if (!is.null(all_gsea_hall[[lbl]])) sum(all_gsea_hall[[lbl]]$is_sig) else 0
    message(sprintf("  %s | %-20s GO:BP sig=%d | Hall sig=%d",
                    contrast_nm, ct, n_gobp, n_hall))
  }
}

# ==============================================================================
# 09.2 — ORA (clusterProfiler) — GO:BP + KEGG
# ==============================================================================

message("=== 09.2 ORA THAL ===")

all_ora <- list()

for (contrast_nm in names(deg_inputs)) {
  degs_list <- deg_inputs[[contrast_nm]]

  for (ct in names(degs_list)) {
    df  <- degs_list[[ct]]
    lbl <- paste(contrast_nm, ct, sep = "_")

    df_valid  <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
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
    message(sprintf("  %s | %-20s GO:BP=%d KEGG=%d", contrast_nm, ct, n_bp, n_kegg))
  }
}

save(all_gsea_gobp, all_gsea_hall, all_gsea_reac, all_ora,
     file = file.path(func_dir, "functional_results_raw_THAL.RData"))
message("Salvato: functional_results_raw_THAL.RData")

# ==============================================================================
# 09.3 — EXPORT EXCEL GSEA
# ==============================================================================

message("=== 09.3 Export Excel THAL ===")

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
      if (!is.null(all_gsea_gobp[[lbl]])) cbind(db="GO:BP",    flatten_fgsea(all_gsea_gobp[[lbl]])),
      if (!is.null(all_gsea_hall[[lbl]])) cbind(db="Hallmark", flatten_fgsea(all_gsea_hall[[lbl]])),
      if (!is.null(all_gsea_reac[[lbl]])) cbind(db="Reactome", flatten_fgsea(all_gsea_reac[[lbl]]))
    ))
    if (is.null(rows) || nrow(rows) == 0) next

    rows <- rows[order(rows$pval), ]
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, rows)

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

message("=== 09.4 Heatmap NES per contrasto THAL ===")

plot_nes_heatmap <- function(gsea_list, contrast_nm, db_label, p_thr = 0.20, top_n = 30) {

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

  sig_paths <- df_all[!is.na(df_all$pval) & df_all$pval < p_thr, "pathway"]
  if (length(sig_paths) == 0) return(invisible(NULL))

  path_rank <- df_all[df_all$pathway %in% sig_paths, ] %>%
    group_by(pathway) %>%
    summarise(max_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_nes)) %>%
    slice_head(n = top_n)
  sig_paths <- path_rank$pathway

  nes_mat <- df_all[df_all$pathway %in% sig_paths, ] %>%
    select(pathway, cell_type, NES) %>%
    pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  pval_mat <- df_all[df_all$pathway %in% sig_paths, ] %>%
    select(pathway, cell_type, pval) %>%
    pivot_wider(names_from = cell_type, values_from = pval, values_fill = 1) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  pval_mat <- pval_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]

  sig_lab <- matrix(
    nrow = nrow(pval_mat), ncol = ncol(pval_mat),
    ifelse(pval_mat < 0.01, "**",
    ifelse(pval_mat < 0.05, "*",
    ifelse(pval_mat < 0.20, "\u00b7", ""))),
    dimnames = dimnames(nes_mat)
  )

  rownames(nes_mat) <- str_wrap(
    gsub("HALLMARK_|GOBP_|REACTOME_", "", rownames(nes_mat)), 45
  )
  rownames(sig_lab) <- rownames(nes_mat)

  lim     <- max(abs(nes_mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))

  ht <- Heatmap(
    nes_mat, name = "NES", col = col_fun,
    cell_fun = function(j, i, x, y, w, h, fill)
      grid.text(sig_lab[i, j], x, y, gp = gpar(fontsize = 9)),
    cluster_rows     = TRUE, cluster_columns = TRUE,
    row_names_gp     = gpar(fontsize = 7),
    column_names_gp  = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title     = paste("GSEA NES \u2014", db_label, "|", contrast_nm),
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
# 09.5 — HEATMAP CROSS-CONTRASTO (THAL_90 vs THAL_24)
# ==============================================================================

message("=== 09.5 Heatmap cross-contrasto THAL ===")

build_cross_contrast_heatmap <- function(gsea_list, db_label, p_thr = 0.20, top_n = 40) {

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
    top_annotation  = col_ha,
    cluster_rows    = TRUE, cluster_columns = TRUE,
    row_names_gp    = gpar(fontsize = 6.5),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    column_title    = paste("Cross-contrasto THAL \u2014", db_label)
  )

  fn <- paste0("GSEA_crosscontrast_THAL_", db_label, ".png")
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
# 09.6 — ENRICHMENT PLOT per pathway di interesse biologico (talamo)
# ==============================================================================

message("=== 09.6 Enrichment plot pathway LSD — THAL ===")

# Pathway di interesse biologico specifici per il talamo:
#   - IEG / plasticità sinaptica (come nel CTX)
#   - Sviluppo talamica / proiezioni talamocorticali
#   - Segnalazione glutamatergica (TC relay) e GABAergica (TRN)
#   - LSD: serotonina, mTOR, infiammazione
pathways_focus_thal <- c(
  # IEG / attivazione precoce
  "GOBP_IMMEDIATE_EARLY_GENE_ACTIVATION",
  "GOBP_REGULATION_OF_GENE_EXPRESSION",
  # Plasticità sinaptica e proiezioni
  "GOBP_SYNAPTIC_PLASTICITY",
  "GOBP_LONG_TERM_POTENTIATION",
  "GOBP_THALAMOCORTICAL_AXON_GUIDANCE",
  "GOBP_THALAMUS_DEVELOPMENT",
  # Segnalazione dei recettori talamici
  "GOBP_REGULATION_OF_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_GABA_RECEPTOR_SIGNALING",
  "GOBP_IONOTROPIC_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY",
  # Serotonina (target LSD)
  "GOBP_SEROTONIN_RECEPTOR_SIGNALING_PATHWAY",
  # Pathways LSD generali
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  # Sviluppo neuronale
  "GOBP_NEURON_PROJECTION_DEVELOPMENT",
  "GOBP_AXON_GUIDANCE"
)

gs_all <- c(gs_gobp, gs_hall, gs_reac)

for (contrast_nm in names(deg_inputs)) {
  for (ct in names(deg_inputs[[contrast_nm]])) {
    df_ct  <- deg_inputs[[contrast_nm]][[ct]]
    df_v   <- df_ct[!is.na(df_ct$log2FoldChange), ]
    if (nrow(df_v) < 10) next

    ranks <- setNames(df_v$log2FoldChange, df_v$gene)
    ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)

    for (pw in pathways_focus_thal) {
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

# ==============================================================================
# 09b — Frequenza pathway GO:BP significativi per timepoint
# ==============================================================================

message("=== 09b: frequenza GO:BP per timepoint THAL ===")

P_THR <- 0.05
TOP_N <- 30

gsea_flat <- do.call(rbind, lapply(names(all_gsea_gobp), function(lbl) {
  res <- all_gsea_gobp[[lbl]]
  if (is.null(res) || nrow(res) == 0) return(NULL)

  contrast <- regmatches(lbl, regexpr("^[^_]+_(90min|24h)", lbl))
  ct       <- sub(paste0("^", contrast, "_"), "", lbl)

  df           <- as.data.frame(res)
  df$contrast  <- contrast
  df$cell_type <- ct
  df$timepoint <- ifelse(grepl("90min", contrast), "90min", "24h")
  df
}))

if (is.null(gsea_flat) || nrow(gsea_flat) == 0) {
  message("  all_gsea_gobp vuoto — salto 09b")
} else {

  gsea_flat$pathway_label <- gsub("^GOBP_", "",   gsea_flat$pathway)
  gsea_flat$pathway_label <- gsub("_",      " ",   gsea_flat$pathway_label)
  gsea_flat$pathway_label <- tolower(gsea_flat$pathway_label)
  gsea_flat$pathway_label <- paste0(
    toupper(substr(gsea_flat$pathway_label, 1, 1)),
    substr(gsea_flat$pathway_label, 2, nchar(gsea_flat$pathway_label))
  )

  plot_gobp_freq <- function(tp, filename) {

    df_tp <- gsea_flat[gsea_flat$timepoint == tp, ]
    if (nrow(df_tp) == 0) {
      message(sprintf("  Nessun dato per timepoint %s", tp))
      return(invisible(NULL))
    }

    n_ct_total <- length(unique(df_tp$cell_type))
    df_sig     <- df_tp[!is.na(df_tp$pval) & df_tp$pval < P_THR, ]

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

    freq <- freq[order(-freq$n_sig, -abs(freq$nes_mean)), ]
    freq <- freq[seq_len(min(TOP_N, nrow(freq))), ]

    freq$pathway_label <- factor(freq$pathway_label,
                                 levels = rev(freq$pathway_label))
    nes_lim <- ceiling(max(abs(freq$nes_mean), na.rm = TRUE) * 10) / 10

    p <- ggplot(freq, aes(x = n_sig, y = pathway_label, fill = nes_mean)) +
      geom_col(width = 0.7, color = "white", linewidth = 0.2) +
      geom_vline(xintercept = n_ct_total,
                 linetype = "dotted", color = "grey40", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%d/%d", n_sig, n_ct_total)),
                hjust = -0.15, size = 2.6, color = "grey25") +
      scale_fill_gradient2(
        low = "#2471A3", mid = "grey92", high = "#CB4335",
        midpoint = 0, limits = c(-nes_lim, nes_lim), name = "NES medio"
      ) +
      scale_x_continuous(
        breaks = seq(0, n_ct_total, by = max(1, n_ct_total %/% 5)),
        limits = c(0, n_ct_total * 1.22),
        expand = c(0, 0)
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.y        = element_text(size = 8, lineheight = 0.9),
        axis.text.x        = element_text(size = 8),
        axis.title.y       = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right",
        legend.key.height  = unit(1.2, "cm"),
        legend.title       = element_text(size = 8),
        plot.title         = element_text(size = 11, face = "bold"),
        plot.subtitle      = element_text(size = 8, color = "grey40")
      ) +
      labs(
        title    = sprintf("GO:BP \u2014 pathway pi\u00f9 frequenti | THAL | %s", tp),
        subtitle = sprintf(
          "p < %.2f (fgsea) | top %d pathway per frequenza tra %d CT | colore = NES medio",
          P_THR, TOP_N, n_ct_total),
        x = sprintf("N cell type con pathway significativo (su %d)", n_ct_total)
      )

    savepng(p, file.path(thal_dir, filename),
            width  = 4200,
            height = max(2200, nrow(freq) * 75 + 700))
    message(sprintf("  Salvato: %s (%d pathway, %d CT)", filename, nrow(freq), n_ct_total))
    invisible(p)
  }

  plot_gobp_freq("90min", "09b_gobp_freq_THAL_90min.png")
  plot_gobp_freq("24h",   "09b_gobp_freq_THAL_24h.png")
}

message(sprintf("=== 09 THAL completato — output in: %s ===", func_dir))
message("=== 09b THAL completato ===")



# ==============================================================================
# GENE SETS CUSTOM — Sistemi neurotrasmettitoriali e funzioni talamiche
# ==============================================================================
# Fonti: Allen Brain Atlas, IUPHAR-DB, letteratura serotonina/LSD
# Geni validati per Mus musculus

custom_gene_sets <- list(
  
  # --- SEROTONINA (target primario LSD) ---
  # Recettori 5-HT espressi nel talamo (relay e TRN)
  SEROTONIN_RECEPTORS = c(
    "Htr1a","Htr1b","Htr1d","Htr1f",          # 5-HT1: Gi, inibitori
    "Htr2a","Htr2b","Htr2c",                    # 5-HT2: Gq — TARGET PRIMARIO LSD
    "Htr3a","Htr3b",                            # 5-HT3: ionotropici
    "Htr4","Htr5a","Htr6","Htr7"               # altri sottotipi
  ),
  # Macchina serotoninergica completa (sintesi, trasporto, degradazione)
  SEROTONIN_SYSTEM = c(
    "Tph1","Tph2",                             # triptofano idrossilasi (sintesi)
    "Ddc",                                     # DOPA decarbossilasi
    "Slc6a4",                                  # SERT (trasportatore reuptake)
    "Slc18a2",                                 # VMAT2 (storage vescicolare)
    "Maoa","Maob",                             # MAO (degradazione)
    "Htr2a","Htr2c","Htr1a","Htr1b","Htr7"   # recettori chiave
  ),
  
  # --- DOPAMINA ---
  DOPAMINE_SYSTEM = c(
    "Th","Ddc",                                # sintesi (TH è rate-limiting)
    "Slc6a3",                                  # DAT (reuptake)
    "Slc18a2",                                 # VMAT2
    "Comt","Maoa","Maob","Drd",               # degradazione
    "Drd1","Drd2","Drd3","Drd4","Drd5"       # recettori D1-D5
  ),
  DOPAMINE_RECEPTORS_D1_CLASS = c("Drd1","Drd5"),     # Gs, attivatori AC
  DOPAMINE_RECEPTORS_D2_CLASS = c("Drd2","Drd3","Drd4"), # Gi, inibitori
  
  # --- GLUTAMMATO (relay talamocorticale) ---
  # Il talamo usa principalmente glutammato come neurotrasmettitore principale
  GLUTAMATE_RECEPTORS_IONOTROPIC = c(
    "Gria1","Gria2","Gria3","Gria4",          # AMPA (GluA1-4)
    "Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b", # NMDA
    "Grik1","Grik2","Grik3","Grik4","Grik5"  # Kainato
  ),
  GLUTAMATE_RECEPTORS_METABOTROPIC = c(
    "Grm1","Grm5",                            # Gruppo I: Gq, postsinaptici
    "Grm2","Grm3",                            # Gruppo II: Gi, presinaptici
    "Grm4","Grm6","Grm7","Grm8"              # Gruppo III: Gi
  ),
  GLUTAMATE_SYNTHESIS_TRANSPORT = c(
    "Gls","Gls2",                             # glutaminasi
    "Glul",                                   # glutammina sintetasi (glia)
    "Slc17a5","Slc17a6","Slc17a7","Slc17a8", # VGLUT1-4 (caricamento vescicolare)
    "Slc1a1","Slc1a2","Slc1a3","Slc1a6","Slc1a7" # EAAT1-5 (reuptake)
  ),
  
  # --- GABA (TRN — Thalamic Reticular Nucleus) ---
  # Il TRN è puramente GABAergico e fondamentale per il gating talamico
  GABA_SYSTEM = c(
    "Gad1","Gad2",                            # GAD65/GAD67 — sintesi GABA
    "Slc32a1",                                # VGAT (caricamento vescicolare)
    "Slc6a1",                                 # GAT-1 (reuptake)
    "Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6", # GABA-A alpha
    "Gabrb1","Gabrb2","Gabrb3",              # GABA-A beta
    "Gabrg1","Gabrg2","Gabrg3",              # GABA-A gamma
    "Gabrd","Gabre","Gabrp","Gabrq",         # GABA-A delta/epsilon/pi/theta
    "Gabbr1","Gabbr2"                        # GABA-B (metabotropico)
  ),
  GABA_RECEPTORS_A = c(
    "Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6",
    "Gabrb1","Gabrb2","Gabrb3",
    "Gabrg1","Gabrg2","Gabrg3","Gabrd"
  ),
  
  # --- ACETILCOLINA (modulazione talamica — input da nucleo basale) ---
  ACETYLCHOLINE_SYSTEM = c(
    "Chat",                                   # ChAT — sintesi ACh
    "Slc18a3",                                # VAChT — caricamento vescicolare
    "Ache","Bche",                            # acetilcolinesterasi (degradazione)
    "Slc5a7",                                 # CHT1 (reuptake colina)
    "Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7", # nAChR alpha
    "Chrnb2","Chrnb4",                        # nAChR beta
    "Chrm1","Chrm2","Chrm3","Chrm4","Chrm5" # mAChR M1-M5
  ),
  MUSCARINIC_RECEPTORS = c("Chrm1","Chrm2","Chrm3","Chrm4","Chrm5"),
  NICOTINIC_RECEPTORS  = c(
    "Chrna2","Chrna4","Chrna7","Chrnb2","Chrnb4" # principali nel SNC
  ),
  
  # --- NORADRENALINA (locus coeruleus → talamo) ---
  NORADRENALINE_SYSTEM = c(
    "Th","Ddc","Dbh",                         # sintesi (TH→DOPA→NE)
    "Slc6a2",                                 # NET (reuptake)
    "Slc18a2",                                # VMAT2
    "Maoa","Maob","Comt",                     # degradazione
    "Adra1a","Adra1b","Adra1d",              # alpha-1 (Gq)
    "Adra2a","Adra2b","Adra2c",              # alpha-2 (Gi, presinaptici)
    "Adrb1","Adrb2","Adrb3"                  # beta (Gs)
  ),
  
  # --- PEPTIDI NEUROMODULATORI ---
  NEUROPEPTIDE_SIGNALING_THAL = c(
    # Orexina/ipocretina (arousal, gating talamico)
    "Hcrt","Hcrtr1","Hcrtr2",
    # Neuropeptide Y
    "Npy","Npy1r","Npy2r","Npy5r",
    # Sostanza P / Tachichinine
    "Tac1","Tacr1","Tacr2","Tacr3",
    # Somatostatina (interneuroni inibitori)
    "Sst","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5",
    # CRF/CRH (stress)
    "Crh","Crhr1","Crhr2",
    # VIP (ritmo circadiano, modulazione SCN→talamo)
    "Vip","Vipr1","Vipr2"
  ),
  
  # --- IEG / PLASTICITA' SINAPTICA (risposta acuta LSD 90min) ---
  IMMEDIATE_EARLY_GENES = c(
    "Fos","Fosb","Fosl1","Fosl2",
    "Jun","Junb","Jund",
    "Egr1","Egr2","Egr3","Egr4",
    "Arc","Homer1","Nptx2",
    "Npas4","Nr4a1","Nr4a2","Nr4a3",
    "Bdnf","Vgf","Pcsk1"
  ),
  
  # --- mTOR / TRASDUZIONE DEL SEGNALE (LSD attiva mTOR via 5-HT2A) ---
  MTOR_SIGNALING_CORE = c(
    "Mtor","Rptor","Mlst8",                   # mTORC1 core
    "Rictor","Mapkap1",                       # mTORC2 core
    "Rps6kb1","Rps6kb2","Eif4ebp1",          # effettori mTORC1
    "Akt1","Akt2","Akt3",                     # AKT upstream
    "Pik3ca","Pik3cb","Pik3cd","Pik3r1",     # PI3K
    "Tsc1","Tsc2","Rheb",                     # regolatori
    "Pten","Pdk1"
  ),
  
  # --- PLASTICITA' SINAPTICA STRUTTURALE ---
  SYNAPTIC_REMODELING = c(
    "Bdnf","Ntrk2",                           # BDNF-TrkB
    "Arc","Homer1","Homer2","Homer3",
    "Shank1","Shank2","Shank3",
    "Dlg1","Dlg2","Dlg3","Dlg4",              # PSD-95 family
    "Syngap1","Kalrn",
    "Limk1","Limk2","Cofilin1",               # rimodellamento actina
    "Cdc42","Rac1","RhoA"
  ),
  
  # --- GATING TALAMOCORTICALE (funzione specifica del relay talamica) ---
  # Canali e correnti che controllano le modalità firing (tonic vs burst)
  THALAMIC_GATING_ION_CHANNELS = c(
    "Cacna1g","Cacna1h","Cacna1i",            # Cav3 (T-type Ca2+) — burst firing
    "Cacna1a","Cacna1b","Cacna1c","Cacna1d",  # HVA Ca2+
    "Hcn1","Hcn2","Hcn3","Hcn4",            # HCN (Ih corrente — pacemaker)
    "Kcna1","Kcna2","Kcna4",                  # Kv1 (AHP)
    "Kcnq2","Kcnq3","Kcnq5",                 # Kv7 (M-current)
    "Scn1a","Scn2a","Scn3a","Scn8a",         # Nav — spike initiation
    "Kcnip1","Kcnip2","Kcnip4"               # Kv4 modulatori (A-current)
  ),
  
  # --- NEUROINFIAMMAZIONE / MICROGLIA (risposta a LSD) ---
  NEUROINFLAMMATION = c(
    "Tnf","Il1b","Il6","Il18",
    "Nlrp3","Casp1",
    "Cx3cr1","P2ry12","Tmem119",             # markers microglia omeostasi
    "Cd68","Iba1",                           # attivazione microglia (Aif1)
    "Tgfb1","Tgfb2",
    "Nfkb1","Rela","Ikbkb",
    "Ccl2","Ccl5","Cxcl10"
  ),
  
  # --- CICLO CIRCADIANO (rilevante per effetti temporali LSD 90min vs 24h) ---
  CIRCADIAN_CLOCK_CORE = c(
    "Clock","Bmal1",                          # Arntl
    "Per1","Per2","Per3",
    "Cry1","Cry2",
    "Nr1d1","Nr1d2",                          # REV-ERB alpha/beta
    "Rora","Rorb","Rorc",
    "Ciart","Nfil3","Dbp","Tef","Hlf"
  )
)

message(sprintf("  Custom gene sets caricati: %d categorie", length(custom_gene_sets)))


# ==============================================================================
# 09.1b — GSEA su gene sets custom (neurotrasmettitori)
# ==============================================================================

message("=== 09.1b GSEA custom neurotransmitter sets ===")

all_gsea_custom <- list()

for (contrast_nm in names(deg_inputs)) {
  degs_list <- deg_inputs[[contrast_nm]]
  
  for (ct in names(degs_list)) {
    df    <- degs_list[[ct]]
    lbl   <- paste(contrast_nm, ct, sep = "_")
    df_v  <- df[!is.na(df$log2FoldChange), ]
    if (nrow(df_v) < 15) next
    
    ranks <- setNames(df_v$log2FoldChange, df_v$gene)
    ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
    
    res <- tryCatch(
      fgsea(pathways    = custom_gene_sets,
            stats       = ranks,
            minSize     = 4,   # set piccoli: abbassa soglia
            maxSize     = 200,
            eps         = 0,
            nPermSimple = 10000),
      error = function(e) { message(sprintf("  Custom GSEA errore %s: %s", lbl, e$message)); NULL }
    )
    if (is.null(res) || nrow(res) == 0) next
    res$label      <- lbl
    res$contrast   <- contrast_nm
    res$cell_type  <- ct
    all_gsea_custom[[lbl]] <- res
    
    n_sig <- sum(!is.na(res$padj) & res$padj < 0.05)
    message(sprintf("  %s | %-20s custom sig=%d", contrast_nm, ct, n_sig))
  }
}


# ==============================================================================
# 09.7 — Heatmap NES per sistema neurotrasmettitoriale × CT × contrasto
# ==============================================================================

message("=== 09.7 Heatmap custom neurotransmitter systems ===")

if (length(all_gsea_custom) > 0) {
  
  df_custom_flat <- do.call(rbind, lapply(all_gsea_custom, as.data.frame))
  
  for (contrast_nm in names(deg_inputs)) {
    
    df_cnt <- df_custom_flat[df_custom_flat$contrast == contrast_nm, ]
    if (nrow(df_cnt) == 0) next
    
    # Matrice NES: righe = gene set, colonne = cell type
    nes_mat <- df_cnt %>%
      select(pathway, cell_type, NES) %>%
      pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
      column_to_rownames("pathway") %>%
      as.matrix()
    
    pval_mat <- df_cnt %>%
      select(pathway, cell_type, pval) %>%
      pivot_wider(names_from = cell_type, values_from = pval, values_fill = 1) %>%
      column_to_rownames("pathway") %>%
      as.matrix()
    pval_mat <- pval_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]
    
    sig_lab <- matrix(
      ifelse(pval_mat < 0.01, "**",
             ifelse(pval_mat < 0.05, "*",
                    ifelse(pval_mat < 0.10, "·", ""))),
      nrow = nrow(pval_mat), ncol = ncol(pval_mat),
      dimnames = dimnames(nes_mat)
    )
    
    # Raggruppa le righe per sistema (per row_split in ComplexHeatmap)
    system_group <- case_when(
      grepl("SEROTONIN",    rownames(nes_mat)) ~ "Serotonin",
      grepl("DOPAMINE",     rownames(nes_mat)) ~ "Dopamine",
      grepl("GLUTAMATE",    rownames(nes_mat)) ~ "Glutamate",
      grepl("GABA",         rownames(nes_mat)) ~ "GABA",
      grepl("ACETYLCHOLINE|MUSCARINIC|NICOTINIC", rownames(nes_mat)) ~ "ACh",
      grepl("NORADRENALINE",rownames(nes_mat)) ~ "Noradrenaline",
      grepl("NEUROPEPTIDE", rownames(nes_mat)) ~ "Neuropeptides",
      grepl("IEG|IMMEDIATE",rownames(nes_mat)) ~ "IEG",
      grepl("MTOR",         rownames(nes_mat)) ~ "mTOR",
      grepl("SYNAPTIC",     rownames(nes_mat)) ~ "Synaptic plasticity",
      grepl("THALAMIC_GATING",rownames(nes_mat))~ "Thalamic gating",
      grepl("NEUROINFLAM",  rownames(nes_mat)) ~ "Neuroinflammation",
      grepl("CIRCADIAN",    rownames(nes_mat)) ~ "Circadian",
      TRUE ~ "Other"
    )
    
    lim     <- max(abs(nes_mat), na.rm = TRUE)
    col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))
    
    ht <- Heatmap(
      nes_mat, name = "NES", col = col_fun,
      row_split       = system_group,
      row_gap         = unit(2, "mm"),
      row_title_gp    = gpar(fontsize = 8, fontface = "bold"),
      cell_fun = function(j, i, x, y, w, h, fill)
        grid.text(sig_lab[i, j], x, y, gp = gpar(fontsize = 8)),
      cluster_rows     = TRUE, cluster_columns = TRUE,
      cluster_row_slices = FALSE,  # mantieni l'ordine dei sistemi
      row_names_gp     = gpar(fontsize = 7.5),
      column_names_gp  = gpar(fontsize = 8),
      column_names_rot = 45,
      column_title     = paste("Neurotransmitter systems — THAL |", contrast_nm),
      heatmap_legend_param = list(title = "NES")
    )
    
    fn <- paste0("custom_NT_heatmap_", contrast_nm, ".png")
    png(file.path(func_dir, fn),
        width  = max(3500, ncol(nes_mat) * 280 + 1800),
        height = max(3000, nrow(nes_mat) * 65 + 1000),
        res    = 300)
    draw(ht)
    dev.off()
    message(sprintf("  Salvato: %s", fn))
  }
}

# ==============================================================================
# 09.8 — UP/DOWN per gene set custom — Bubble plot migliorato
# ==============================================================================
# - Scala fissa con 0 visibile in tutti i pannelli (facet per contrasto)
# - Etichette solo sui punti con n_up + n_down >= LABEL_THR
# - Diagonale di riferimento (net score = 0)
# - Colore = net score normalizzato [-1, +1]
# - Size = pct overlap (n_up+n_down) / gene_set_size

message("=== 09.8 UP/DOWN bubble plot custom sets ===")

LABEL_THR <- 3   # etichetta solo se almeno 3 geni DEG nel gene set

updown_results <- list()

for (contrast_nm in names(deg_inputs)) {
  for (ct in names(deg_inputs[[contrast_nm]])) {
    df  <- deg_inputs[[contrast_nm]][[ct]]
    lbl <- paste(contrast_nm, ct, sep = "_")
    
    df_sig <- df[!is.na(df$padj) & df$padj < 0.05 & !is.na(df$log2FoldChange), ]
    if (nrow(df_sig) < 5) next
    
    genes_up   <- df_sig$gene[df_sig$log2FoldChange > 0]
    genes_down <- df_sig$gene[df_sig$log2FoldChange < 0]
    
    rows <- lapply(names(custom_gene_sets), function(gs_nm) {
      gs   <- custom_gene_sets[[gs_nm]]
      n_up <- sum(genes_up   %in% gs)
      n_dn <- sum(genes_down %in% gs)
      n_tot <- length(gs)
      data.frame(
        gene_set    = gs_nm,
        contrast    = contrast_nm,
        cell_type   = ct,
        n_up        = n_up,
        n_down      = n_dn,
        n_total_gs  = n_tot,
        pct_overlap = (n_up + n_dn) / max(n_tot, 1),
        net_score   = ifelse((n_up + n_dn) > 0,
                             (n_up - n_dn) / (n_up + n_dn),
                             NA_real_),
        stringsAsFactors = FALSE
      )
    })
    updown_results[[lbl]] <- do.call(rbind, rows)
  }
}

ud_flat <- do.call(rbind, updown_results)

# Rimuovi righe senza alcun hit in nessun CT/contrasto
gs_any_hit <- ud_flat %>%
  group_by(gene_set) %>%
  summarise(total_hits = sum(n_up + n_down, na.rm = TRUE), .groups = "drop") %>%
  filter(total_hits > 0) %>%
  pull(gene_set)
ud_plot <- ud_flat[ud_flat$gene_set %in% gs_any_hit, ]

# Sistema neurotrasmettitoriale per colore shape/shape alternativo
ud_plot$system <- dplyr::case_when(
  grepl("SEROTONIN",       ud_plot$gene_set) ~ "Serotonin",
  grepl("DOPAMINE",        ud_plot$gene_set) ~ "Dopamine",
  grepl("GLUTAMATE",       ud_plot$gene_set) ~ "Glutamate",
  grepl("GABA",            ud_plot$gene_set) ~ "GABA",
  grepl("ACETYLCHOLINE|MUSCARINIC|NICOTINIC", ud_plot$gene_set) ~ "ACh",
  grepl("NORADRENALINE",   ud_plot$gene_set) ~ "Noradrenaline",
  grepl("NEUROPEPTIDE",    ud_plot$gene_set) ~ "Neuropeptides",
  grepl("IEG|IMMEDIATE",   ud_plot$gene_set) ~ "IEG",
  grepl("MTOR",            ud_plot$gene_set) ~ "mTOR",
  grepl("SYNAPTIC",        ud_plot$gene_set) ~ "Synaptic plasticity",
  grepl("THALAMIC_GATING", ud_plot$gene_set) ~ "Thalamic gating",
  grepl("NEUROINFLAM",     ud_plot$gene_set) ~ "Neuroinflammation",
  grepl("CIRCADIAN",       ud_plot$gene_set) ~ "Circadian",
  TRUE ~ "Other"
)

# Etichetta pulita (rimuove prefisso, sostituisce _ con spazio)
ud_plot$label_clean <- gsub("_", " ", tolower(ud_plot$gene_set))

# Asse fisso: calcola il max globale su tutti i contrasti
axis_max <- max(c(ud_plot$n_up, ud_plot$n_down), na.rm = TRUE)
axis_lim <- ceiling(axis_max * 1.08)  # margine per etichette

# Un plot per CT (facet per contrasto), così ogni pannello ha stessa scala
# Un unico plot: facet_grid CT (righe) × contrasto (colonne)
df_plot_all <- ud_plot[(ud_plot$n_up + ud_plot$n_down) > 0, ]

if (nrow(df_plot_all) == 0) {
  message("  Nessun hit — salto 09.8")
} else {
  
  df_label_all <- df_plot_all[(df_plot_all$n_up + df_plot_all$n_down) >= LABEL_THR, ]
  
  p <- ggplot(df_plot_all, aes(x = n_up, y = n_down)) +
    
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    
    geom_point(aes(size = pct_overlap, fill = net_score),
               shape = 21, color = "white", stroke = 0.3, alpha = 0.88) +
    
    ggrepel::geom_text_repel(
      data              = df_label_all,
      aes(label         = label_clean),
      size              = 2.2,
      color             = "grey20",
      segment.color     = "grey65",
      segment.size      = 0.3,
      box.padding       = 0.35,
      point.padding     = 0.2,
      max.overlaps      = 15,
      min.segment.length = 0.2
    ) +
    
    scale_fill_gradient2(
      low      = "#2471A3", mid = "grey92", high = "#CB4335",
      midpoint = 0, limits = c(-1, 1), na.value = "grey75",
      name     = "Net score\n(UP\u2212DOWN)/total",
      breaks   = c(-1, -0.5, 0, 0.5, 1),
      labels   = c("\u22121", "\u22120.5", "0", "+0.5", "+1")
    ) +
    
    scale_size_continuous(
      range  = c(1.5, 7),
      name   = "% gene set\ncon DEGs",
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    scale_x_continuous(
      limits = c(-0.5, axis_lim),
      breaks = pretty(c(0, axis_lim), n = 4),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(-0.5, axis_lim),
      breaks = pretty(c(0, axis_lim), n = 4),
      expand = c(0, 0)
    ) +
    
    # Righe = CT, Colonne = contrasto — scala fissa uguale in tutti i pannelli
    facet_grid(cell_type ~ contrast) +
    
    theme_bw(base_size = 8) +
    theme(
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text.x       = element_text(size = 8, face = "bold"),
      strip.text.y       = element_text(size = 7, face = "bold", angle = 0),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
      legend.position    = "right",
      legend.key.height  = unit(1.0, "cm"),
      legend.title       = element_text(size = 7.5),
      legend.text        = element_text(size = 7),
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 7.5, color = "grey40"),
      axis.title         = element_text(size = 7.5),
      axis.text          = element_text(size = 6.5)
    ) +
    labs(
      title    = "UP vs DOWN DEGs in neurotransmitter gene sets | THAL",
      subtitle = sprintf(
        "Etichette se \u2265%d geni DEG (padj<0.05) nel gene set | diagonale = net score 0 | scala fissa",
        LABEL_THR),
      x = "N geni UP-regolati",
      y = "N geni DOWN-regolati"
    )
  
  n_ct  <- length(unique(df_plot_all$cell_type))
  n_cnt <- length(unique(df_plot_all$contrast))
  
  savefunc(p, "updown_bubble_allCT.png",
           width  = max(2000, n_cnt * 50 + 100),
           height = max(2500, n_ct  * 300 + 800))
  message("  Salvato: updown_bubble_allCT.png")
}

#TODO 
# WikiPathways (ottimo per pathway di segnalazione neurotrasmettitoriale)
gmt_wiki <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS")
gs_wiki  <- split(gmt_wiki$gene_symbol, gmt_wiki$gs_name)
# Pathway rilevanti: WP_SEROTONIN_RECEPTOR_4_6_7_AND_NR3C_SIGNALING,
# WP_DOPAMINE_METABOLISM, WP_GABA_RECEPTOR_SIGNALING, ecc.

# Cell type signatures (utile per verificare identità cellulare)
gmt_cell <- msigdbr(species = "Mus musculus", category = "C8")
gs_cell  <- split(gmt_cell$gene_symbol, gmt_cell$gs_name)

# SynGO (sinapsi) — non è in MSigDB, scaricabile da syngoportal.org
# Ottimale per analisi di compartimenti sinaptici pre/postsinaptici