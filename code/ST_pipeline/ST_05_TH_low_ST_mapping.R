################################################################################
##  ST_05_TH_low_ST_mapping.R
##
##  Input:  ST_03_experiment_merged_ABA.RData
##          ST_03_df_toplot_ABA.RData
##          THAL_annotated_MMC.RData  (da pipeline scRNA MapMyCells)
##            contiene: THAL, cluster_TH_stats, low_aliases,
##                      cluster_roi_norm, clust_meta
##  Output: ST_05_predictions_TH_low.RData
##          ST_05_df_toplot_THlow.RData
##          ST_05_TH_low_ST_test.csv
##          ST_05_spatial_TH_low/    — plot spaziali per cluster TH_low
##          ST_05_heatmap_TH_low_ST_areas.png
##          ST_05_dotplot_TH_low_test.png
##
##  Flusso:
##    1. Transfer THAL TH_low -> ST (score per cluster TH_low per spot)
##    2. Visualizzazione spaziale per ogni cluster TH_low
##    3. Test Wilcoxon: area ST dominante vs ROI "scorretta" da MapMyCells
##    4. Heatmap e dotplot riassuntivi
################################################################################

source("code/ST_pipeline/ST_00_config.R")

load(file.path(outdir, "ST_03_experiment_merged_ABA.RData"))
load(file.path(outdir, "ST_03_df_toplot_ABA.RData"))

# Carica output pipeline scRNA MapMyCells
MMC_PATH <- file.path("results/20260224", "THAL_annotated_MMC.RData")
load(MMC_PATH)
# Oggetti attesi: THAL, cluster_TH_stats, low_aliases, clust_meta

# ==============================================================================
# 1. Label transfer THAL TH_low -> ST
# ==============================================================================

message("=== ST_05A. Label transfer THAL TH_low -> ST ===")

THAL_low <- THAL[, THAL$TH_group == "TH_low"]
cell_labels_low <- as.factor(THAL_low$mmc_cluster_alias)
names(cell_labels_low) <- colnames(THAL_low)

cat(sprintf("Cellule TH_low: %d | Cluster distinti: %d\n",
            ncol(THAL_low), nlevels(cell_labels_low)))

if (!"SCT" %in% names(THAL_low@assays)) {
  THAL_low <- SCTransform(THAL_low, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
}

if (!"SCT" %in% names(experiment.merged@assays)) {
  experiment.merged <- SCTransform(experiment.merged, assay = "Spatial",
                                   verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
}

DefaultAssay(experiment.merged) <- "SCT"
DefaultAssay(THAL_low)        <- "SCT"

transfer_features <- intersect(
  VariableFeatures(THAL_low),
  rownames(experiment.merged)
)
cat(sprintf("Variable features comuni: %d\n", length(transfer_features)))

anchors_THlow <- FindTransferAnchors(
  reference            = THAL_low,
  query                = experiment.merged,
  normalization.method = "SCT",
  features             = transfer_features
)

predictions_TH_low <- TransferData(
  anchorset        = anchors_THlow,
  refdata          = cell_labels_low,
  prediction.assay = TRUE,
  weight.reduction = experiment.merged[["pca"]],
  dims             = 1:N_PCS
)

experiment.merged[["predictions_TH_low"]] <- predictions_TH_low
save(predictions_TH_low,
     file = file.path(outdir, "ST_05_predictions_TH_low.RData"))

# Aggiungi score a df
score_mat <- t(predictions_TH_low@data[
  rownames(predictions_TH_low@data) != "max", , drop = FALSE
])
colnames(score_mat) <- paste0("TH_low_cl_", colnames(score_mat))
shared <- intersect(rownames(df), rownames(score_mat))
df[shared, colnames(score_mat)] <- score_mat[shared, ]

message("=== ST_05A completato ===")


# ==============================================================================
# 2. Visualizzazione spaziale per cluster TH_low
# ==============================================================================

message("=== ST_05B. Plot spaziali per ROI (aggregato) ===")

dir.create(file.path(outdir, "ST_05_spatial_TH_low"), showWarnings = FALSE)
# Se structure_df è in memoria (dalla pipeline 04_mapmycells):
# Scarica ontologia CCFv3 e crea roi_name_map
library(httr)

ontology_resp <- GET("http://api.brain-map.org/api/v2/structure_graph_download/1.json")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

parse_structure_tree <- function(node) {
  row <- data.frame(
    id      = node$id      %||% NA_integer_,
    acronym = node$acronym %||% NA_character_,
    name    = if (!is.null(node$`safe-name`) && length(node$`safe-name`) > 0)
      node$`safe-name`
    else if (!is.null(node$name) && length(node$name) > 0)
      node$name
    else NA_character_,
    stringsAsFactors = FALSE
  )
  children <- node$children
  if (length(children) > 0) {
    child_rows <- bind_rows(lapply(children, parse_structure_tree))
    row <- bind_rows(row, child_rows)
  }
  row
}

ont_data     <- content(ontology_resp, as = "parsed")
structure_df <- parse_structure_tree(ont_data$msg[[1]])
roi_name_map <- setNames(structure_df$name, structure_df$acronym)

cat("Strutture caricate:", nrow(structure_df), "\n")

cluster_top_roi_map <- cluster_TH_stats %>%
  filter(cluster_alias %in% low_aliases) %>%
  select(cluster_alias, top_roi, pct_TH_norm)

# Raggruppa cluster per top_roi
roi_groups <- split(
  as.character(cluster_top_roi_map$cluster_alias),
  cluster_top_roi_map$top_roi
)

cat("ROI con cluster TH_low assegnati:\n")
print(sapply(roi_groups, length))

# Per ogni ROI: score medio dei cluster assegnati
for (roi in names(roi_groups)) {
  aliases_in_roi <- roi_groups[[roi]]
  cols_in_roi    <- paste0("TH_low_cl_", aliases_in_roi)
  cols_in_roi    <- cols_in_roi[cols_in_roi %in% colnames(df)]
  
  if (length(cols_in_roi) == 0) next
  
  # Score aggregato: media dei cluster assegnati a questa ROI
  if (length(cols_in_roi) == 1) {
    df$score_roi_tmp <- df[[cols_in_roi]]
  } else {
    df$score_roi_tmp <- rowMeans(df[, cols_in_roi], na.rm = TRUE)
  }
  
  n_clusters <- length(cols_in_roi)
  roi_name   <- roi_name_map[roi] %||% roi  # usa nomi estesi CCFv3 se disponibile
  
  p <- ggplot(df, aes(x = x, y = -y, colour = score_roi_tmp)) +
    geom_point(size = 0.7, alpha = 0.8) +
    facet_wrap(~ area, ncol = 2) +
    scale_colour_gradientn(
      colors = c("grey93", "#F4A261", "#E63946", "#6B0000"),
      name   = "Score\nmedio", limits = c(0, NA)
    ) +
    labs(
      title    = sprintf("TH_low → %s", roi_name),
      subtitle = sprintf("ROI: %s | %d cluster aggregati", roi, n_clusters)
    ) +
    theme_classic(base_size = 9)
  
  savepng(
    p,
    file.path(outdir, "ST_05_spatial_TH_low",
              sprintf("roi_%s.png", gsub(" |/|\\(|\\)", "_", roi))),
    width = 4500, height = 5500
  )
  
  message(sprintf("  Plottato: %s (%d cluster)", roi, n_clusters))
}

df$score_roi_tmp <- NULL  # pulizia colonna temporanea
message("=== ST_05B completato ===")


# ==============================================================================
# ST_05B2. Plot spaziale cellule TH_high (score aggregato cluster TH_high)
# ==============================================================================

message("=== ST_05B2. Plot spaziale TH_high ===")

dir.create(file.path(outdir, "ST_05_spatial_TH_high"), showWarnings = FALSE)

# Cluster TH_high: cluster_alias presenti in THAL con TH_group == "TH_high"
high_aliases <- unique(THAL$mmc_cluster_alias[THAL$TH_group == "TH_high"])
high_aliases <- high_aliases[!is.na(high_aliases)]

cat(sprintf("Cluster TH_high: %d\n", length(high_aliases)))

# Mappa top_roi per i cluster TH_high
cluster_top_roi_high <- cluster_TH_stats %>%
  filter(cluster_alias %in% high_aliases) %>%
  select(cluster_alias, top_roi, pct_TH_norm)

# Per i cluster TH_high il transfer è già stato fatto su tutto THAL —
# le colonne TH_low_cl_* esistono solo per i TH_low.
# Usa predictions_THAL se disponibile, altrimenti fai transfer su TH_high.

if (!"predictions_TH_high" %in% ls()) {
  
  message("  Eseguo label transfer THAL_high -> ST...")
  
  THAL_high <- THAL[, THAL$TH_group == "TH_high"]
  cell_labels_high <- as.factor(THAL_high$mmc_cluster_alias)
  names(cell_labels_high) <- colnames(THAL_high)
  
  if (!"SCT" %in% names(THAL_high@assays)) {
    THAL_high <- SCTransform(THAL_high, verbose = FALSE) %>%
      RunPCA(verbose = FALSE)
  }
  
  DefaultAssay(THAL_high)        <- "SCT"
  DefaultAssay(experiment.merged) <- "SCT"
  
  transfer_features_high <- intersect(
    VariableFeatures(THAL_high),
    rownames(experiment.merged)
  )
  cat(sprintf("  Variable features comuni TH_high: %d\n",
              length(transfer_features_high)))
  
  anchors_THhigh <- FindTransferAnchors(
    reference            = THAL_high,
    query                = experiment.merged,
    normalization.method = "SCT",
    features             = transfer_features_high
  )
  
  predictions_TH_high <- TransferData(
    anchorset        = anchors_THhigh,
    refdata          = cell_labels_high,
    prediction.assay = TRUE,
    weight.reduction = experiment.merged[["pca"]],
    dims             = 1:N_PCS
  )
  
  save(predictions_TH_high,
       file = file.path(outdir, "ST_05_predictions_TH_high.RData"))
  
  # Aggiungi score a df
  score_mat_high <- t(predictions_TH_high@data[
    rownames(predictions_TH_high@data) != "max", , drop = FALSE
  ])
  colnames(score_mat_high) <- paste0("TH_high_cl_", colnames(score_mat_high))
  shared_high <- intersect(rownames(df), rownames(score_mat_high))
  df[shared_high, colnames(score_mat_high)] <- score_mat_high[shared_high, ]
}

# Raggruppa cluster TH_high per top_roi e plotta score medio
roi_groups_high <- split(
  as.character(cluster_top_roi_high$cluster_alias),
  cluster_top_roi_high$top_roi
)

cat("ROI con cluster TH_high assegnati:\n")
print(sapply(roi_groups_high, length))

for (roi in names(roi_groups_high)) {
  aliases_in_roi <- roi_groups_high[[roi]]
  cols_in_roi    <- paste0("TH_high_cl_", aliases_in_roi)
  cols_in_roi    <- cols_in_roi[cols_in_roi %in% colnames(df)]
  
  if (length(cols_in_roi) == 0) next
  
  df$score_roi_tmp <- if (length(cols_in_roi) == 1) {
    df[[cols_in_roi]]
  } else {
    rowMeans(df[, cols_in_roi], na.rm = TRUE)
  }
  
  roi_name  <- roi_name_map[roi] %||% roi
  n_clusters <- length(cols_in_roi)
  
  p <- ggplot(df, aes(x = x, y = -y, colour = score_roi_tmp)) +
    geom_point(size = 0.7, alpha = 0.8) +
    facet_wrap(~ area, ncol = 2) +
    scale_colour_gradientn(
      colors = c("grey93", "#A8DADC", "#457B9D", "#1D3557"),
      name   = "Score\nmedio", limits = c(0, NA)
    ) +
    labs(
      title    = sprintf("TH_high → %s", roi_name),
      subtitle = sprintf("ROI: %s | %d cluster aggregati", roi, n_clusters)
    ) +
    theme_classic(base_size = 9)
  
  savepng(
    p,
    file.path(outdir, "ST_05_spatial_TH_high",
              sprintf("roi_%s.png", gsub(" |/|\\(|\\)", "_", roi))),
    width = 4500, height = 5500
  )
  
  message(sprintf("  Plottato TH_high: %s (%d cluster)", roi, n_clusters))
}

df$score_roi_tmp <- NULL
message("=== ST_05B2 completato ===")


# Score aggregato unico TH_high — media di tutti i cluster TH_high
cols_TH_high <- grep("^TH_high_cl_", colnames(df), value = TRUE)

df$score_TH_high <- if (length(cols_TH_high) == 1) {
  df[[cols_TH_high]]
} else {
  rowMeans(df[, cols_TH_high], na.rm = TRUE)
}

p_high <- ggplot(df, aes(x = x, y = -y, colour = score_TH_high)) +
  geom_point(size = 0.7, alpha = 0.8) +
  facet_wrap(~ area, ncol = 2) +
  scale_colour_gradientn(
    colors = c("grey93", "#A8DADC", "#457B9D", "#1D3557"),
    name   = "Score\nmedio", limits = c(0, NA)
  ) +
  labs(
    title    = "TH_high — score aggregato",
    subtitle = sprintf("%d cluster TH_high aggregati", length(cols_TH_high))
  ) +
  theme_classic(base_size = 9)

savepng(p_high,
        file.path(outdir, "ST_05_spatial_TH_high", "TH_high_aggregated.png"),
        width = 4500, height = 5500)

df$score_roi_tmp <- NULL
message("=== ST_05B2 completato ===")

# ==============================================================================
# 3. Test Wilcoxon: area ST dominante vs ROI attesa da MapMyCells
# ==============================================================================

message("=== ST_05C. Test area ST vs ROI MapMyCells ===")

# <<< Aggiorna questo mapping dopo aver visto unique(df$area.spot) >>>
roi_to_areaspot <- c(
  "TH"               = "Thalamus",
  "HIP"              = "Hippocampal formation",
  "RHP"              = "Retrohippocampal region",
  "STRd"             = "Striatum",
  "STRv"             = "Striatum",
  "HY"               = "Hypothalamus",
  "MB"               = "Midbrain",
  "MY"               = "Medulla",
  "PAL"              = "Pallidum",
  "LSX"              = "Lateral septal complex",
  "sAMY"             = "Striatum",
  "OLF"              = "Olfactory areas",
  "CTXsp"            = "Cortical subplate",
  "Isocortex"        = "Isocortex",
  "P"                = "Pons",
  "CB"               = "Cerebellum",
  "AUD-TEa-PERI-ECT" = "Isocortex",
  "MO-FRP"           = "Isocortex",
  "ACA"              = "Isocortex",
  "RSP"              = "Isocortex",
  "AI"               = "Isocortex",
  "VIS"              = "Isocortex",
  "VIS-PTLp"         = "Isocortex",
  "SSp"              = "Isocortex",
  "PL-ILA-ORB"       = "Isocortex"
)

cat("\nLabel area.spot in df:\n")
print(sort(table(df$area.spot), decreasing = TRUE))

test_cluster_vs_ST <- function(cl, df, cluster_top_roi_map, roi_to_areaspot) {
  col_name <- paste0("TH_low_cl_", cl)
  if (!col_name %in% colnames(df)) return(NULL)

  scores     <- df[[col_name]]; scores[is.na(scores)] <- 0
  top_roi_cl <- cluster_top_roi_map$top_roi[
    cluster_top_roi_map$cluster_alias == as.integer(cl)]
  if (!length(top_roi_cl) || is.na(top_roi_cl)) return(NULL)

  expected_area <- roi_to_areaspot[top_roi_cl]
  if (is.na(expected_area)) {
    message(sprintf("  WARN cluster %s: no mapping per ROI '%s'", cl, top_roi_cl))
    expected_area <- top_roi_cl
  }

  group_A <- scores[!is.na(df$area.spot) & df$area.spot == expected_area]
  group_B <- scores[!is.na(df$area.spot) & df$area.spot != expected_area]
  if (length(group_A) < 3 || length(group_B) < 3) return(NULL)

  wt   <- wilcox.test(group_A, group_B, alternative = "greater")
  r_rb <- 1 - (2 * wt$statistic) / (length(group_A) * length(group_B))

  area_means <- df %>%
    filter(!is.na(area.spot)) %>%
    group_by(area.spot) %>%
    summarise(mean_score = mean(.data[[col_name]], na.rm = TRUE),
              .groups    = "drop") %>%
    arrange(desc(mean_score))

  data.frame(
    cluster_alias       = cl,
    top_roi_mmc         = top_roi_cl,
    expected_area_spot  = expected_area,
    n_spots_expected    = length(group_A),
    n_spots_other       = length(group_B),
    mean_score_expected = round(mean(group_A), 4),
    mean_score_other    = round(mean(group_B), 4),
    wilcox_p            = wt$p.value,
    rank_biserial_r     = round(as.numeric(r_rb), 3),
    top_area_ST         = area_means$area.spot[1],
    top_area_ST_score   = round(area_means$mean_score[1], 4),
    match_roi_ST        = (area_means$area.spot[1] == expected_area),
    stringsAsFactors    = FALSE
  )
}

results_list <- lapply(as.character(low_aliases), function(cl) {
  tryCatch(
    test_cluster_vs_ST(cl, df, cluster_top_roi_map, roi_to_areaspot),
    error = function(e) { message("  ERROR cl ", cl, ": ", e$message); NULL }
  )
})

results_ST <- bind_rows(results_list) %>%
  mutate(wilcox_p_adj = p.adjust(wilcox_p, method = "fdr")) %>%
  arrange(wilcox_p_adj)

# Aggiungi nome annotazione cluster
if ("cluster_annotation_term_name" %in% colnames(clust_meta)) {
  results_ST <- results_ST %>%
    left_join(
      clust_meta %>%
        mutate(cluster_alias = as.character(cluster_alias)) %>%
        select(cluster_alias, cluster_annotation_term_name),
      by = "cluster_alias"
    )
}

cat("\n=== Risultati test ===\n")
print(as.data.frame(results_ST), n = 40)

n_match     <- sum(results_ST$match_roi_ST,  na.rm = TRUE)
n_sig       <- sum(results_ST$wilcox_p_adj < 0.05, na.rm = TRUE)
n_sig_match <- sum(results_ST$match_roi_ST & results_ST$wilcox_p_adj < 0.05,
                   na.rm = TRUE)
cat(sprintf("\nCluster testati:             %d\n",   nrow(results_ST)))
cat(sprintf("Match ROI/area_ST:           %d (%.0f%%)\n",
            n_match, 100 * n_match / nrow(results_ST)))
cat(sprintf("Significativi (FDR<5%%):     %d\n",    n_sig))
cat(sprintf("Significativi + match:       %d\n",    n_sig_match))

write.csv(results_ST,
          file.path(outdir, "ST_05_TH_low_ST_test.csv"),
          row.names = FALSE)


# ==============================================================================
# 4. Plot riassuntivi
# ==============================================================================

message("=== ST_05D. Plot riassuntivi ===")

# -- Heatmap score medio per cluster TH_low x area ST -------------------------
score_cols <- paste0("TH_low_cl_", as.character(low_aliases))
score_cols <- score_cols[score_cols %in% colnames(df)]

area_mean_mat <- df %>%
  filter(!is.na(area.spot)) %>%
  select(area.spot, all_of(score_cols)) %>%
  group_by(area.spot) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  tibble::column_to_rownames("area.spot") %>%
  as.matrix()

colnames(area_mean_mat) <- sub("TH_low_cl_", "", colnames(area_mean_mat))

# Annota colonne con top_roi
col_anno <- data.frame(
  top_ROI_MMC = cluster_top_roi_map$top_roi[
    match(colnames(area_mean_mat),
          as.character(cluster_top_roi_map$cluster_alias))
  ],
  row.names = colnames(area_mean_mat)
)

png(file.path(outdir, "ST_05_heatmap_TH_low_ST_areas.png"),
    width  = max(3000, ncol(area_mean_mat) * 120 + 1500),
    height = max(2000, nrow(area_mean_mat) * 160 + 800),
    res    = 300)
pheatmap(
  area_mean_mat,
  annotation_col = col_anno,
  color          = colorRampPalette(c("white", "#F4A261", "#E63946"))(50),
  fontsize       = 8,
  fontsize_col   = 6,
  main           = "Score medio transfer per area ST x cluster TH_low"
)
dev.off()

# -- Dotplot effect size vs significatività ------------------------------------
p_dot <- ggplot(results_ST,
                aes(x     = rank_biserial_r,
                    y     = -log10(wilcox_p_adj + 1e-10),
                    colour = match_roi_ST,
                    label  = paste0(cluster_alias, "\n", top_roi_mmc))) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey50") +
  scale_colour_manual(
    values = c("TRUE" = "#2D6A4F", "FALSE" = "#E63946"),
    name   = "Area ST\ncorrisponde\na ROI MMC"
  ) +
  labs(
    title    = "Cluster TH_low: arricchimento nell'area ST attesa",
    subtitle = "Linea tratteggiata = FDR 5%",
    x        = "Effect size (rank-biserial r)",
    y        = "-log10(p_adj)"
  ) +
  theme_classic(base_size = 12)
savepng(p_dot, "ST_05_dotplot_TH_low_test.png", width = 4000, height = 3000)

# -- Barplot match per ROI -----------------------------------------------------
roi_match_df <- results_ST %>%
  group_by(top_roi_mmc) %>%
  summarise(
    n_total = n(),
    n_match = sum(match_roi_ST, na.rm = TRUE),
    n_sig   = sum(wilcox_p_adj < 0.05, na.rm = TRUE),
    n_sig_match = sum(match_roi_ST & wilcox_p_adj < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pct_match = 100 * n_match / n_total) %>%
  arrange(desc(pct_match))

p_bar_roi <- ggplot(roi_match_df,
                    aes(x = reorder(top_roi_mmc, pct_match),
                        y = pct_match, fill = n_total)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%d/%d", n_match, n_total)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_viridis_c(name = "N cluster") +
  labs(title = "% cluster TH_low con match ROI/area_ST per ROI",
       x = "ROI MapMyCells", y = "% match") +
  theme_classic(base_size = 11) +
  expand_limits(y = 110)
savepng(p_bar_roi, "ST_05_barplot_match_per_ROI.png",
        width = 3500, height = 2500)


# ==============================================================================
# 5. Salvataggio
# ==============================================================================

save(df, file = file.path(outdir, "ST_05_df_toplot_THlow.RData"))

message("=== ST_05 completato ===")
message("  Risultati test: ", file.path(outdir, "ST_05_TH_low_ST_test.csv"))



# Evidenzia cellule TH_low mappate su Medulla nella UMAP
THAL$highlight_MY <- case_when(
  THAL$TH_group != "TH_low"                        ~ "TH_high / no_mapping",
  THAL$cluster_top_roi == "MY"                      ~ "TH_low → Medulla",
  THAL$TH_group == "TH_low"                         ~ "TH_low → altra ROI"
)

p_MY <- DimPlot(
  THAL,
  group.by = "highlight_MY",
  cols     = c(
    "TH_low → Medulla"     = "#E63946",
    "TH_low → altra ROI"   = "#F4A261",
    "TH_high / no_mapping" = "grey85"
  ),
  pt.size = 0.3,
  order   = c("TH_low → Medulla", "TH_low → altra ROI")
) +
  ggtitle("Cellule TH_low mappate su Medulla (MY)") +
  theme_classic(base_size = 12)

savepng(p_MY, "04_UMAP_TH_low_Medulla.png", width = 3500, height = 3000)

cat("Cellule TH_low → Medulla:", sum(THAL$highlight_MY == "TH_low → Medulla", na.rm = TRUE), "\n")
cat("% su totale TH_low:", 
    round(100 * mean(THAL$highlight_MY[THAL$TH_group == "TH_low"] == "TH_low → Medulla", na.rm = TRUE), 1), "%\n")
