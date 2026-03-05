################################################################################
##  11_htr2a.R                                                                 ##
##                                                                             ##
##  Input:  CTX_clean.RData, THAL_clean.RData  (da sezione 05)               ##
##  Output: CTX_htr2a.RData, THAL_htr2a.RData  (con score in @meta.data)     ##
##                                                                             ##
##  Analisi:                                                                   ##
##    1. FeaturePlot Htr2a + IEG per area × trattamento                      ##
##    2. AddModuleScore: IEG score + Gq signaling score                       ##
##    3. VlnPlot score per CT × trattamento                                   ##
##    4. Scatter Htr2a medio vs IEG score (per CT × trattamento)             ##
##    5. Confronto CTX vs THAL: IEG score boxplot                             ##
##    6. Heatmap firma completa Htr2a                                         ##
################################################################################

source("00_config.R")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_clean.RData"))   # CTX_clean
load(file.path(outdir, "THAL_clean.RData"))  # THAL_clean

# ==============================================================================
# 11.1 — FEATUREPLOT HTR2A E IEG (pre-score, ispezione visiva)
# ==============================================================================

message("=== 11.1 FeaturePlot Htr2a e IEG ===")

genes_show <- unique(c("Htr2a","Fos","Arc","Egr1","Nr4a1","Fosb","Junb"))

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  # Tutti i geni split per trattamento (Ctrl | LSD)
  genes_present <- filter_genes(genes_show, obj)
  if (length(genes_present) == 0) next

  p_feat <- FeaturePlot(
    obj,
    features  = genes_present,
    split.by  = "treatment",
    reduction = "umap",
    ncol      = 2,
    pt.size   = 0.2
  ) & scale_color_viridis_c(option = "inferno") &
    theme(legend.position = "right",
          plot.title = element_text(size = 8))

  savepng(p_feat,
          paste0("11_FeaturePlot_Htr2a_IEG_", label, ".png"),
          width  = 5500,
          height = max(2000, ceiling(length(genes_present) / 2) * 2200))

  # Htr2a da solo, per campione
  if ("Htr2a" %in% rownames(obj)) {
    p_htr2a_samp <- FeaturePlot(
      obj, features = "Htr2a", split.by = "orig.ident",
      reduction = "umap", pt.size = 0.2
    ) & scale_color_viridis_c(option = "plasma") &
      theme(plot.title = element_text(size = 7))
    savepng(p_htr2a_samp,
            paste0("11_Htr2a_by_sample_", label, ".png"),
            width = 8000, height = 2500)
  }
}

# ==============================================================================
# 11.2 — MODULE SCORE: IEG + Gq signaling
# ==============================================================================

message("=== 11.2 Module score IEG + Gq ===")

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  ieg_genes <- filter_genes(
    c("Fos","Arc","Egr1","Nr4a1","Npas4","Fosb","Junb","Dusp1","Dusp5"),
    obj
  )
  gq_genes  <- filter_genes(
    c("Gnaq","Plcb1","Prkca","Camk2a","Camk4"),
    obj
  )

  if (length(ieg_genes) < 3) {
    message(sprintf("  %s: troppo pochi geni IEG (%d) — skip score",
                    label, length(ieg_genes)))
    next
  }

  # AddModuleScore: aggiunge colonne Score1, Score2 a @meta.data
  # Poi le rinominiamo per chiarezza
  obj <- AddModuleScore(
    obj,
    features = list(ieg_genes, gq_genes),
    name     = "ModScore_"
  )
  obj$IEG_score <- obj$ModScore_1
  obj$Gq_score  <- obj$ModScore_2
  obj$ModScore_1 <- NULL
  obj$ModScore_2 <- NULL

  message(sprintf("  %s: IEG score calcolato su %d geni | Gq su %d geni",
                  label, length(ieg_genes), length(gq_genes)))

  # Riassegna l'oggetto aggiornato
  assign(obj_nm, obj)
}

# ==============================================================================
# 11.3 — FEATUREPLOT SCORE SUL UMAP
# ==============================================================================

message("=== 11.3 FeaturePlot module score ===")

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  score_cols <- intersect(c("IEG_score","Gq_score"), colnames(obj@meta.data))
  if (length(score_cols) == 0) next

  plots <- lapply(score_cols, function(sc) {
    FeaturePlot(obj, features = sc, reduction = "umap", pt.size = 0.3) +
      scale_color_gradient2(low = "#2471A3", mid = "white", high = "#E74C3C",
                            midpoint = 0, name = sc) +
      ggtitle(paste(label, "—", sc))
  })

  savepng(
    wrap_plots(plots, ncol = 2) +
      plot_annotation(title = paste("Module score —", label)),
    paste0("11_ModuleScore_UMAP_", label, ".png"),
    width = 5000, height = 2500
  )

  # FeaturePlot split per trattamento
  for (sc in score_cols) {
    p_split <- FeaturePlot(
      obj, features = sc, split.by = "treatment",
      reduction = "umap", pt.size = 0.3
    ) + scale_color_gradient2(low = "#2471A3", mid = "white",
                               high = "#E74C3C", midpoint = 0)
    savepng(p_split,
            paste0("11_ModScore_", sc, "_split_", label, ".png"),
            width = 5000, height = 2500)
  }
}

# ==============================================================================
# 11.4 — VIOLINPLOT SCORE PER CT × TRATTAMENTO
# ==============================================================================

message("=== 11.4 VlnPlot score per CT ===")

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  score_cols <- intersect(c("IEG_score","Gq_score"), colnames(obj@meta.data))
  if (length(score_cols) == 0) next

  n_ct <- length(unique(obj$cell_type))

  p_vln <- VlnPlot(
    obj,
    features  = score_cols,
    group.by  = "cell_type",
    split.by  = "treatment",
    pt.size   = 0,
    ncol      = length(score_cols)
  ) & scale_fill_manual(values = pal_trt) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "top")

  savepng(p_vln,
          paste0("11_Vln_score_", label, ".png"),
          width  = max(3000, n_ct * 200 + 600) * length(score_cols),
          height = 2500)

  # VlnPlot anche per Htr2a
  if ("Htr2a" %in% rownames(obj)) {
    p_htr2a_vln <- VlnPlot(
      obj, features = "Htr2a",
      group.by = "cell_type", split.by = "treatment",
      pt.size  = 0, cols = pal_trt
    ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste("Htr2a —", label))
    savepng(p_htr2a_vln,
            paste0("11_Vln_Htr2a_", label, ".png"),
            width = max(3000, n_ct * 200 + 600), height = 2500)
  }
}

# ==============================================================================
# 11.5 — SCATTER: Htr2a medio vs IEG score per CT × trattamento
# ==============================================================================

message("=== 11.5 Scatter Htr2a vs IEG score ===")

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  if (!"IEG_score" %in% colnames(obj@meta.data)) next
  if (!"Htr2a"    %in% rownames(obj))            next

  # Media per CT × trattamento
  scatter_df <- obj@meta.data %>%
    group_by(cell_type, treatment) %>%
    summarise(
      IEG_mean = mean(IEG_score, na.rm = TRUE),
      .groups  = "drop"
    )

  # Espressione media Htr2a per CT × trattamento
  # In Seurat v5 usa FetchData che gestisce i layer correttamente
  htr2a_df <- FetchData(obj, vars = c("Htr2a","cell_type","treatment")) %>%
    group_by(cell_type, treatment) %>%
    summarise(Htr2a_mean = mean(Htr2a, na.rm = TRUE), .groups = "drop")

  scatter_df <- left_join(scatter_df, htr2a_df, by = c("cell_type","treatment"))

  p_sc <- ggplot(scatter_df,
                 aes(x = Htr2a_mean, y = IEG_mean,
                     color = treatment, label = cell_type)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = pal_trt) +
    theme_bw() +
    labs(title    = paste("Htr2a vs IEG score —", label),
         subtitle = "Un punto per CT × trattamento",
         x = "Htr2a (media log-norm)",
         y = "IEG module score")
  savepng(p_sc, paste0("11_Scatter_Htr2a_IEG_", label, ".png"),
          width = 3000, height = 2500)
}

# ==============================================================================
# 11.6 — CONFRONTO CTX vs THAL: IEG score boxplot
# ==============================================================================

message("=== 11.6 CTX vs THAL — IEG score ===")

has_ieg_ctx  <- "IEG_score" %in% colnames(CTX_clean@meta.data)
has_ieg_thal <- "IEG_score" %in% colnames(THAL_clean@meta.data)

if (has_ieg_ctx && has_ieg_thal) {

  df_cmp <- bind_rows(
    CTX_clean@meta.data  %>%
      select(cell_type, treatment, time, IEG_score) %>%
      mutate(area = "Cortex"),
    THAL_clean@meta.data %>%
      select(cell_type, treatment, time, IEG_score) %>%
      mutate(area = "Thalamus")
  )

  p_cmp <- ggplot(df_cmp,
                  aes(x = cell_type, y = IEG_score, fill = treatment)) +
    geom_boxplot(outlier.size = 0.2, position = position_dodge(0.8),
                 linewidth = 0.4) +
    facet_grid(area ~ time, scales = "free_x", space = "free") +
    scale_fill_manual(values = pal_trt) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text   = element_text(size = 8)) +
    labs(title    = "IEG score — CTX vs THAL | LSD vs Ctrl",
         subtitle = "Facet: area × timepoint",
         x = "", y = "IEG module score")

  n_ct_max <- max(length(unique(CTX_clean$cell_type)),
                  length(unique(THAL_clean$cell_type)))
  savepng(p_cmp, "11_IEG_score_CTX_vs_THAL.png",
          width  = max(5000, n_ct_max * 250 + 1000),
          height = 4500)
}

# ==============================================================================
# 11.7 — HEATMAP FIRMA COMPLETA Htr2a
# ==============================================================================

message("=== 11.7 Heatmap firma Htr2a ===")

for (obj_nm in c("CTX_clean","THAL_clean")) {
  obj   <- get(obj_nm)
  label <- sub("_clean","", obj_nm)

  sig_genes <- filter_genes(markers_htr2a_signature, obj)
  if (length(sig_genes) < 5) next

  # Espressione media per CT × trattamento
  # In Seurat v5 DefaultAssay deve essere RNA o SCT, non matrice sparse raw
  DefaultAssay(obj) <- "SCT"   # usa normalizzazione SCT per la visualizzazione

  avg <- AverageExpression(
    obj,
    features = sig_genes,
    group.by = c("cell_type","treatment"),
    assays   = "SCT",
    layer    = "data"
  )$SCT

  # Pulisci nomi colonne (formato "CT_treatment")
  colnames(avg) <- gsub("_Ctrl$","_C", gsub("_LSD$","_L", colnames(avg)))

  mat_z <- t(scale(t(as.matrix(avg))))
  mat_z[is.nan(mat_z)] <- 0

  lim     <- min(3, max(abs(mat_z), na.rm = TRUE))
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))

  # Annotazione colonne: Ctrl vs LSD
  col_ann_df <- data.frame(
    Treatment = ifelse(grepl("_C$", colnames(mat_z)), "Ctrl", "LSD"),
    row.names = colnames(mat_z)
  )
  col_ha <- HeatmapAnnotation(
    df  = col_ann_df,
    col = list(Treatment = pal_trt),
    show_legend = TRUE
  )

  ht <- Heatmap(
    mat_z,
    name              = "z-score",
    col               = col_fun,
    top_annotation    = col_ha,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    row_names_gp      = gpar(fontsize = 8),
    column_names_gp   = gpar(fontsize = 7),
    column_names_rot  = 60,
    column_title      = paste("Firma Htr2a (z-score SCT) —", label),
    heatmap_legend_param = list(title = "z-score")
  )

  png(file.path(outdir, paste0("11_Heatmap_Htr2a_signature_", label, ".png")),
      width  = max(3500, ncol(mat_z) * 90 + 1000),
      height = max(2500, nrow(mat_z) * 60 + 600),
      res    = 300)
  draw(ht)
  dev.off()
  message(sprintf("  Salvato: 11_Heatmap_Htr2a_signature_%s.png", label))

  DefaultAssay(obj) <- "RNA"  # ripristina default
  assign(obj_nm, obj)
}

# ==============================================================================
# 11.8 — SALVATAGGIO
# ==============================================================================

message("=== 11.8 Salvataggio ===")

CTX_clean  <- CTX_clean
THAL_clean <- THAL_clean

save(CTX_clean,  file = file.path(outdir, "CTX_htr2a.RData"))
save(THAL_clean, file = file.path(outdir, "THAL_htr2a.RData"))

message(sprintf("Salvati: CTX_htr2a.RData, THAL_htr2a.RData"))
message("=== 11 completato ===")
