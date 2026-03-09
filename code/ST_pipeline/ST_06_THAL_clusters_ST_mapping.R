################################################################################
##  ST_06_THAL_clusters_ST_mapping.R
##
##  Mappa i cluster talamici dal subclustering (10b) sui dati spatial Visium.
##
##  Usa DUE livelli di annotazione come riferimento scRNA:
##    1. cell_type  (da 05_cell_type_annotation_THAL — annotazione biologica)
##    2. cluster_full (da 10b — subclustering su tutti i geni)
##
##  Per ciascun livello:
##    A. Label transfer scRNA → ST (score per cluster per spot)
##    B. Plot spaziali: score per cluster, facettato per campione
##    C. Plot UMAP ST colorato per area talamiche predette
##    D. Heatmap score medio per cluster × area anatomica ST
##    E. Dotplot: score medio TH spot vs non-TH spot per cluster
##    F. Test Wilcoxon: arricchimento nei TH spot vs resto
##    G. Confronto cell_type vs cluster_full
##
##  Input:
##    ST_03_experiment_merged_ABA.RData   (experiment.merged, df con area.spot)
##    ST_03_df_toplot_ABA.RData
##    THAL/THAL_celltype.RData            (THAL con cell_type)
##    THAL_neurons.RData                  (neu con cluster_full, label_marker)
##
##  Output:
##    ST_06_predictions_celltype.RData
##    ST_06_predictions_cluster_full.RData
##    ST_06_df_toplot.RData
##    ST_06_spatial_celltype/     — plot spaziali per cell_type
##    ST_06_spatial_cluster_full/ — plot spaziali per cluster_full
##    ST_06_heatmap_celltype.png
##    ST_06_heatmap_cluster_full.png
##    ST_06_wilcox_results.csv
##    ST_06_comparison_ct_vs_full.png
################################################################################

source("code/ST_pipeline/ST_00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
})

# ==============================================================================
# PARAMETRI
# ==============================================================================

N_PCS         <- 30     # dims per TransferData
SCORE_THRESH  <- 0.05   # soglia minima score per considerare uno spot "positivo"
TH_AREA_LABEL <- "Thalamus"   # label di area.spot che corrisponde al talamo

# ==============================================================================
# CARICAMENTO
# ==============================================================================

message("=== Caricamento dati ===")

load(file.path(outdir, "ST_03_experiment_merged_ABA.RData"))   # experiment.merged
load(file.path(outdir, "ST_03_df_toplot_ABA.RData"))           # df

load("D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/Psychedelics/4-RNAflex/2_Full_experiment/4_Data_Analysis/RNAflex_full/results/20260224/THAL/THAL_celltype.RData")         # THAL

# Verifica campi necessari  
stopifnot("cell_type"    %in% colnames(THAL@meta.data))

message(sprintf("  THAL: %d cellule | %d cell_type",
                ncol(THAL), length(unique(THAL$cell_type))))
message(sprintf("  ST:   %d spot",    ncol(experiment.merged)))

# Directory output
dir_ct   <- file.path(outdir, "ST_06_spatial_celltype")
dir_full <- file.path(outdir, "ST_06_spatial_cluster_full")
dir.create(dir_ct,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_full, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HELPER
# ==============================================================================

# Assicura SCT su reference e query
prep_sct <- function(obj, assay_name = "SCT") {
  if (!assay_name %in% names(obj@assays)) {
    message("  SCTransform su ", deparse(substitute(obj)), "...")
    obj <- SCTransform(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
  }
  DefaultAssay(obj) <- assay_name
  obj
}

# Label transfer generico: reference → query ST
do_transfer <- function(ref, query, ref_labels, n_pcs = N_PCS,
                        label_name = "transfer") {
  ref   <- prep_sct(ref)
  query <- prep_sct(query)

  tf <- intersect(VariableFeatures(ref), rownames(query))
  message(sprintf("  [%s] Features comuni: %d", label_name, length(tf)))
  if (length(tf) < 100)
    stop(sprintf("Troppo poche features comuni per %s (%d)", label_name, length(tf)))

  anchors <- FindTransferAnchors(
    reference            = ref,
    query                = query,
    normalization.method = "SCT",
    features             = tf,
    verbose              = FALSE
  )
  message(sprintf("  [%s] Anchors trovati: %d", label_name, nrow(anchors@anchors)))

  preds <- TransferData(
    anchorset        = anchors,
    refdata          = ref_labels,
    prediction.assay = TRUE,
    weight.reduction = query[["pca"]],
    dims             = 1:n_pcs,
    verbose          = FALSE
  )
  preds
}

# Estrai score matrix da prediction assay e aggiungi a df
add_scores_to_df <- function(df, preds, prefix) {
  score_mat           <- t(preds@data[rownames(preds@data) != "max", , drop = FALSE])
  colnames(score_mat) <- paste0(prefix, colnames(score_mat))
  shared              <- intersect(rownames(df), rownames(score_mat))
  df[shared, colnames(score_mat)] <- score_mat[shared, ]
  df
}

# Plot spaziale singolo cluster
plot_spatial_score <- function(df, score_col, title, subtitle = "",
                               color_high = "#E63946") {
  scores <- df[[score_col]]
  scores[is.na(scores)] <- 0
  ggplot(df, aes(x = x, y = -y, colour = !!sym(score_col))) +
    geom_point(size = 0.6, alpha = 0.85) +
    facet_wrap(~ area, ncol = 2) +
    scale_colour_gradientn(
      colors = c("grey93", scales::muted(color_high, l = 70), color_high),
      name   = "Score", limits = c(0, NA), na.value = "grey93"
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey92", color = NA),
          strip.text       = element_text(size = 7))
}

# ==============================================================================
# PARTE A — LABEL TRANSFER: cell_type
# ==============================================================================

message("=== PARTE A: Label transfer cell_type ===")

preds_ct_path <- file.path(outdir, "ST_06_predictions_celltype.RData")

if (file.exists(preds_ct_path)) {
  message("  Carico predictions_ct salvate...")
  load(preds_ct_path)
} else {
  labels_ct <- as.factor(THAL$cell_type)
  names(labels_ct) <- colnames(THAL)

  predictions_ct <- do_transfer(
    ref        = THAL,
    query      = experiment.merged,
    ref_labels = labels_ct,
    label_name = "cell_type"
  )
  save(predictions_ct, file = preds_ct_path)
}

experiment.merged[["predictions_ct"]] <- predictions_ct
df <- add_scores_to_df(df, predictions_ct, prefix = "ct_")

ct_labels <- setdiff(rownames(predictions_ct@data), "max")
message(sprintf("  Cell type trasferiti: %d", length(ct_labels)))

# ==============================================================================
# PARTE B — LABEL TRANSFER: cluster_full
# ==============================================================================

message("=== PARTE B: Label transfer cluster_full ===")

preds_full_path <- file.path(outdir, "ST_06_predictions_cluster_full.RData")

if (file.exists(preds_full_path)) {
  message("  Carico predictions_full salvate...")
  load(preds_full_path)
} else {
  labels_full <- as.factor(neu$cluster_full)
  names(labels_full) <- colnames(neu)

  predictions_full <- do_transfer(
    ref        = neu,
    query      = experiment.merged,
    ref_labels = labels_full,
    label_name = "cluster_full"
  )
  save(predictions_full, file = preds_full_path)
}

experiment.merged[["predictions_full"]] <- predictions_full
df <- add_scores_to_df(df, predictions_full, prefix = "full_")

full_labels <- setdiff(rownames(predictions_full@data), "max")
message(sprintf("  Cluster_full trasferiti: %d", length(full_labels)))

# ==============================================================================
# PARTE C — PLOT SPAZIALI: cell_type
# ==============================================================================

message("=== PARTE C: Plot spaziali cell_type ===")

for (ct in ct_labels) {
  col <- paste0("ct_", ct)
  if (!col %in% colnames(df)) next

  p <- plot_spatial_score(
    df         = df,
    score_col  = col,
    title      = sprintf("cell_type: %s", ct),
    subtitle   = "Score label transfer — tutti i campioni",
    color_high = "#E63946"
  )
  ggsave(
    file.path(dir_ct, sprintf("ct_%s.png", gsub("[/ ]", "_", ct))),
    p, width = 14, height = 16, dpi = 200
  )
}
message(sprintf("  Salvati %d plot in ST_06_spatial_celltype/", length(ct_labels)))

# ==============================================================================
# PARTE D — PLOT SPAZIALI: cluster_full
# ==============================================================================

message("=== PARTE D: Plot spaziali cluster_full ===")

# Aggiungi label_marker al cluster_full per titoli più informativi
cl_to_label <- setNames(
  neu$label_marker[match(full_labels, neu$cluster_full)],
  full_labels
)

for (cl in full_labels) {
  col <- paste0("full_", cl)
  if (!col %in% colnames(df)) next

  lbl <- cl_to_label[cl]
  lbl <- if (!is.na(lbl)) lbl else cl

  p <- plot_spatial_score(
    df         = df,
    score_col  = col,
    title      = sprintf("cluster_full %s — %s", cl, lbl),
    subtitle   = "Score label transfer — tutti i campioni",
    color_high = "#457B9D"
  )
  ggsave(
    file.path(dir_full, sprintf("full_%s.png", cl)),
    p, width = 14, height = 16, dpi = 200
  )
}
message(sprintf("  Salvati %d plot in ST_06_spatial_cluster_full/", length(full_labels)))

# ==============================================================================
# PARTE E — HEATMAP: score medio per cluster × area anatomica ST
# ==============================================================================

message("=== PARTE E: Heatmap cluster × area ===")

make_area_heatmap <- function(df, score_prefix, labels, title, filepath,
                              anno_df = NULL) {
  score_cols <- paste0(score_prefix, labels)
  score_cols <- score_cols[score_cols %in% colnames(df)]

  mat <- df %>%
    filter(!is.na(area.spot)) %>%
    select(area.spot, all_of(score_cols)) %>%
    group_by(area.spot) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    tibble::column_to_rownames("area.spot") %>%
    as.matrix()

  colnames(mat) <- sub(score_prefix, "", colnames(mat))

  # scala per colonna (Z-score) per rendere comparabili i pattern
  mat_z <- scale(mat)
  mat_z[is.nan(mat_z)] <- 0

  h <- max(2000, nrow(mat) * 140 + 600)
  w <- max(3000, ncol(mat) * 120 + 800)
  png(filepath, width = w, height = h, res = 300)
  pheatmap(
    mat_z,
    annotation_col  = anno_df,
    color           = colorRampPalette(c("#1D3557","grey97","#E63946"))(80),
    cluster_rows    = TRUE, cluster_cols = TRUE,
    fontsize        = 8, fontsize_row = 8, fontsize_col = 7,
    angle_col       = "45", border_color = "white",
    main            = title
  )
  dev.off()
  message("  Salvato: ", basename(filepath))

  # Salva anche la versione con score raw (non z-scored) per interpretazione
  filepath_raw <- sub("\\.png$", "_raw.png", filepath)
  png(filepath_raw, width = w, height = h, res = 300)
  pheatmap(
    mat,
    annotation_col  = anno_df,
    color           = colorRampPalette(c("white","#F4A261","#E63946"))(80),
    cluster_rows    = TRUE, cluster_cols = TRUE,
    fontsize        = 8, fontsize_row = 8, fontsize_col = 7,
    angle_col       = "45", border_color = "white",
    main            = paste(title, "(score raw)")
  )
  dev.off()
  message("  Salvato: ", basename(filepath_raw))

  invisible(mat)
}

# Annotazione colonne per cell_type: nessuna (la label è già il nome)
mat_ct <- make_area_heatmap(
  df           = df,
  score_prefix = "ct_",
  labels       = ct_labels,
  title        = "Score cell_type talamico × area anatomica ST",
  filepath     = file.path(outdir, "ST_06_heatmap_celltype.png")
)

# Annotazione colonne per cluster_full: aggiungi label_marker
anno_full <- data.frame(
  label_marker = as.character(cl_to_label[full_labels]),
  row.names    = full_labels
)
anno_full$label_marker[is.na(anno_full$label_marker)] <- "unknown"

mat_full <- make_area_heatmap(
  df           = df,
  score_prefix = "full_",
  labels       = full_labels,
  title        = "Score cluster_full × area anatomica ST",
  filepath     = file.path(outdir, "ST_06_heatmap_cluster_full.png"),
  anno_df      = anno_full
)

# ==============================================================================
# PARTE F — DOTPLOT: score TH spot vs non-TH spot
# ==============================================================================

message("=== PARTE F: Dotplot TH vs non-TH spot ===")

make_th_dotplot <- function(df, score_prefix, labels, title, filepath) {
  df_th    <- df[!is.na(df$area.spot) & df$area.spot == TH_AREA_LABEL, ]
  df_nonth <- df[!is.na(df$area.spot) & df$area.spot != TH_AREA_LABEL, ]

  res <- lapply(labels, function(lb) {
    col <- paste0(score_prefix, lb)
    if (!col %in% colnames(df)) return(NULL)
    s_th    <- df_th[[col]];    s_th[is.na(s_th)]       <- 0
    s_nonth <- df_nonth[[col]]; s_nonth[is.na(s_nonth)] <- 0
    data.frame(
      label        = lb,
      mean_TH      = mean(s_th),
      mean_nonTH   = mean(s_nonth),
      ratio        = mean(s_th) / (mean(s_nonth) + 1e-6),
      stringsAsFactors = FALSE
    )
  })
  res_df <- do.call(rbind, Filter(Negate(is.null), res))

  p <- ggplot(res_df, aes(x = mean_nonTH, y = mean_TH,
                           label = label, colour = log2(ratio + 0.01))) +
    geom_point(size = 3, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_text_repel(size = 2.5, max.overlaps = 30) +
    scale_colour_gradientn(
      colors = c("#1D3557","grey80","#E63946"),
      name   = "log2(ratio\nTH/nonTH)"
    ) +
    labs(title = title,
         x     = "Score medio — spot non-TH",
         y     = "Score medio — spot TH") +
    theme_classic(base_size = 11)

  ggsave(filepath, p, width = 12, height = 10, dpi = 250)
  message("  Salvato: ", basename(filepath))
  invisible(res_df)
}

ratio_ct <- make_th_dotplot(
  df           = df,
  score_prefix = "ct_",
  labels       = ct_labels,
  title        = "Score TH vs non-TH per cell_type talamico",
  filepath     = file.path(outdir, "ST_06_dotplot_TH_ct.png")
)

ratio_full <- make_th_dotplot(
  df           = df,
  score_prefix = "full_",
  labels       = full_labels,
  title        = "Score TH vs non-TH per cluster_full",
  filepath     = file.path(outdir, "ST_06_dotplot_TH_full.png")
)

# ==============================================================================
# PARTE G — TEST WILCOXON: TH spot vs non-TH spot per cluster
# ==============================================================================

message("=== PARTE G: Test Wilcoxon TH vs non-TH ===")

wilcox_test_clusters <- function(df, score_prefix, labels, label_map = NULL) {
  df_th    <- df[!is.na(df$area.spot) & df$area.spot == TH_AREA_LABEL, ]
  df_nonth <- df[!is.na(df$area.spot) & df$area.spot != TH_AREA_LABEL, ]

  res <- lapply(labels, function(lb) {
    col <- paste0(score_prefix, lb)
    if (!col %in% colnames(df)) return(NULL)
    s_th    <- df_th[[col]];    s_th[is.na(s_th)]       <- 0
    s_nonth <- df_nonth[[col]]; s_nonth[is.na(s_nonth)] <- 0
    if (length(s_th) < 3 || length(s_nonth) < 3) return(NULL)

    wt   <- wilcox.test(s_th, s_nonth, alternative = "greater")
    r_rb <- 1 - (2 * wt$statistic) / (length(s_th) * length(s_nonth))

    data.frame(
      label          = lb,
      annotation     = if (!is.null(label_map)) label_map[lb] else lb,
      mean_TH        = round(mean(s_th),    5),
      mean_nonTH     = round(mean(s_nonth), 5),
      ratio          = round(mean(s_th) / (mean(s_nonth) + 1e-6), 3),
      wilcox_p       = wt$p.value,
      rank_biserial  = round(as.numeric(r_rb), 4),
      n_spots_TH     = length(s_th),
      n_spots_nonTH  = length(s_nonth),
      stringsAsFactors = FALSE
    )
  })
  res_df <- do.call(rbind, Filter(Negate(is.null), res))
  res_df$wilcox_padj <- p.adjust(res_df$wilcox_p, method = "fdr")
  res_df[order(res_df$wilcox_padj), ]
}

wlx_ct <- wilcox_test_clusters(
  df           = df,
  score_prefix = "ct_",
  labels       = ct_labels
)

wlx_full <- wilcox_test_clusters(
  df           = df,
  score_prefix = "full_",
  labels       = full_labels,
  label_map    = cl_to_label
)

# Combina e salva
wlx_ct$type   <- "cell_type"
wlx_full$type <- "cluster_full"
wlx_all       <- bind_rows(wlx_ct, wlx_full)
write.csv(wlx_all,
          file.path(outdir, "ST_06_wilcox_results.csv"),
          row.names = FALSE)

message("  Wilcoxon cell_type — significativi (FDR<5%):")
print(wlx_ct[wlx_ct$wilcox_padj < 0.05, c("label","mean_TH","mean_nonTH",
                                             "rank_biserial","wilcox_padj")])
message("  Wilcoxon cluster_full — significativi (FDR<5%):")
print(wlx_full[wlx_full$wilcox_padj < 0.05, c("label","annotation",
                                                "rank_biserial","wilcox_padj")])

# Volcano plot
make_volcano <- function(res_df, title, filepath) {
  res_df$sig <- res_df$wilcox_padj < 0.05
  p <- ggplot(res_df,
              aes(x     = rank_biserial,
                  y     = -log10(wilcox_padj + 1e-10),
                  colour = sig,
                  label  = paste0(label, "\n", annotation))) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(data = filter(res_df, sig),
                    size = 2.5, max.overlaps = 30) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", colour = "grey50") +
    scale_colour_manual(
      values = c("TRUE" = "#E63946", "FALSE" = "grey60"),
      name   = "FDR < 5%"
    ) +
    labs(title = title,
         x     = "Effect size (rank-biserial r)",
         y     = "-log10(p_adj)") +
    theme_classic(base_size = 11)
  ggsave(filepath, p, width = 10, height = 8, dpi = 250)
  message("  Salvato: ", basename(filepath))
}

make_volcano(wlx_ct,
             "Arricchimento TH spot — cell_type",
             file.path(outdir, "ST_06_volcano_ct.png"))

make_volcano(wlx_full,
             "Arricchimento TH spot — cluster_full",
             file.path(outdir, "ST_06_volcano_full.png"))

# ==============================================================================
# PARTE H — CONFRONTO cell_type vs cluster_full per area TH
# ==============================================================================

message("=== PARTE H: Confronto cell_type vs cluster_full ===")

# Per ogni spot TH: argmax cell_type e argmax cluster_full
df_th <- df[!is.na(df$area.spot) & df$area.spot == TH_AREA_LABEL, ]

cols_ct_ok   <- paste0("ct_",   ct_labels)[paste0("ct_",   ct_labels)   %in% colnames(df_th)]
cols_full_ok <- paste0("full_", full_labels)[paste0("full_", full_labels) %in% colnames(df_th)]

df_th$pred_ct   <- ct_labels[apply(df_th[, cols_ct_ok,   drop = FALSE], 1, which.max)]
df_th$pred_full <- full_labels[apply(df_th[, cols_full_ok, drop = FALSE], 1, which.max)]
df_th$pred_full_label <- cl_to_label[df_th$pred_full]

# Heatmap di contingenza (numero spot)
ct_full_tab <- table(ct = df_th$pred_ct, full = df_th$pred_full_label)
ct_full_tab <- ct_full_tab[rowSums(ct_full_tab) > 0, colSums(ct_full_tab) > 0]

png(file.path(outdir, "ST_06_comparison_ct_vs_full.png"),
    width  = max(3000, ncol(ct_full_tab) * 150 + 800),
    height = max(2000, nrow(ct_full_tab) * 140 + 600),
    res    = 300)
pheatmap(
  log1p(ct_full_tab),
  color        = colorRampPalette(c("white","#457B9D","#1D3557"))(80),
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize     = 8, fontsize_row = 8, fontsize_col = 7,
  angle_col    = "45", border_color = "white",
  main         = "Spot TH: cell_type vs cluster_full (log1p N spot)"
)
dev.off()
message("  Salvato: ST_06_comparison_ct_vs_full.png")

# UMAP ST colorato per label predetta
if ("umap" %in% names(experiment.merged@reductions)) {
  pred_ct_disc   <- ct_labels[apply(
    t(predictions_ct@data[predictions_ct@data %>% rownames() %>%
                            setdiff("max"), , drop = FALSE]),
    1, which.max)]
  pred_full_disc <- full_labels[apply(
    t(predictions_full@data[predictions_full@data %>% rownames() %>%
                              setdiff("max"), , drop = FALSE]),
    1, which.max)]

  experiment.merged$pred_celltype    <- pred_ct_disc
  experiment.merged$pred_cluster_full <- pred_full_disc
  experiment.merged$pred_full_label  <-
    cl_to_label[experiment.merged$pred_cluster_full]

  p_umap_ct <- DimPlot(experiment.merged, group.by = "pred_celltype",
                       label = TRUE, label.size = 2.5, repel = TRUE,
                       pt.size = 0.3) +
    ggtitle("UMAP ST — cell_type talamico predetto") + NoLegend()

  p_umap_full <- DimPlot(experiment.merged, group.by = "pred_full_label",
                         label = TRUE, label.size = 2.2, repel = TRUE,
                         pt.size = 0.3) +
    ggtitle("UMAP ST — cluster_full (ThalamoSeq label) predetto") + NoLegend()

  p_umap_area <- DimPlot(experiment.merged, group.by = "area.spot",
                         label = TRUE, label.size = 2.5, repel = TRUE,
                         pt.size = 0.3) +
    ggtitle("UMAP ST — area anatomica") + NoLegend()

  ggsave(file.path(outdir, "ST_06_UMAP_ST_predictions.png"),
         p_umap_area | p_umap_ct | p_umap_full,
         width = 18, height = 6, dpi = 250)
  message("  Salvato: ST_06_UMAP_ST_predictions.png")
}

# ==============================================================================
# PARTE I — PLOT SPAZIALE RIASSUNTIVO: label predetta discretizzata per TH spot
# ==============================================================================

message("=== PARTE I: Plot spaziale riassuntivo ===")

# Aggiungi pred discreta a df completo
spots_common <- intersect(rownames(df), names(pred_ct_disc))
df$pred_celltype    <- NA_character_
df$pred_full_label  <- NA_character_
df[spots_common, "pred_celltype"]   <- pred_ct_disc[spots_common]
df[spots_common, "pred_full_label"] <- cl_to_label[pred_full_disc[spots_common]]

# Solo spot TH: colora per pred, resto grigio
df$pred_ct_TH <- ifelse(
  !is.na(df$area.spot) & df$area.spot == TH_AREA_LABEL,
  df$pred_celltype, "non-TH"
)
df$pred_full_TH <- ifelse(
  !is.na(df$area.spot) & df$area.spot == TH_AREA_LABEL,
  df$pred_full_label, "non-TH"
)

# Palette cell_type
ct_unique   <- setdiff(unique(df$pred_ct_TH), "non-TH")
pal_ct_pred <- setNames(
  colorRampPalette(c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
                     "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85",
                     "#FFDC91","#6DB6FF"))(length(ct_unique)),
  ct_unique
)
pal_ct_pred["non-TH"] <- "grey90"

p_ct_space <- ggplot(df, aes(x = x, y = -y, colour = pred_ct_TH)) +
  geom_point(data = filter(df, pred_ct_TH == "non-TH"),
             size = 0.4, alpha = 0.4) +
  geom_point(data = filter(df, pred_ct_TH != "non-TH"),
             size = 0.7, alpha = 0.9) +
  facet_wrap(~ area, ncol = 2) +
  scale_colour_manual(values = pal_ct_pred, name = "cell_type") +
  labs(title    = "Cell type talamico predetto — spot TH",
       subtitle = "Spot non-TH in grigio") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(size = 7),
        legend.text      = element_text(size = 7))

ggsave(file.path(outdir, "ST_06_spatial_ct_discrete.png"),
       p_ct_space, width = 14, height = 16, dpi = 250)
message("  Salvato: ST_06_spatial_ct_discrete.png")

# Palette cluster_full / label
full_unique   <- setdiff(unique(df$pred_full_TH), c("non-TH", NA))
pal_full_pred <- setNames(
  colorRampPalette(c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
                     "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85",
                     "#FFDC91","#6DB6FF","#FF6DB6","#490092"))(length(full_unique)),
  full_unique
)
pal_full_pred["non-TH"] <- "grey90"

p_full_space <- ggplot(df, aes(x = x, y = -y, colour = pred_full_TH)) +
  geom_point(data = filter(df, pred_full_TH == "non-TH"),
             size = 0.4, alpha = 0.4) +
  geom_point(data = filter(df, pred_full_TH != "non-TH"),
             size = 0.7, alpha = 0.9) +
  facet_wrap(~ area, ncol = 2) +
  scale_colour_manual(values = pal_full_pred, name = "cluster\n(ThalamoSeq)") +
  labs(title    = "Cluster_full (label ThalamoSeq) — spot TH",
       subtitle = "Spot non-TH in grigio") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(size = 7),
        legend.text      = element_text(size = 7))

ggsave(file.path(outdir, "ST_06_spatial_full_discrete.png"),
       p_full_space, width = 14, height = 16, dpi = 250)
message("  Salvato: ST_06_spatial_full_discrete.png")

# ==============================================================================
# SALVATAGGIO
# ==============================================================================

message("=== Salvataggio ===")
save(df, file = file.path(outdir, "ST_06_df_toplot.RData"))
message("  Salvato: ST_06_df_toplot.RData")
message("=== ST_06 completato ===")
