################################################################################
##  ST_03_label_transfer_ABA.R
##
##  Input:  ST_02_experiment_merged_clustered.RData
##          ST_02_df_toplot.RData
##          data/Reference/expr_raw_counts_table.tsv
##          data/Reference/meta_table.tsv
##  Output: ST_03_experiment_merged_ABA.RData
##          ST_03_df_toplot_ABA.RData
##
##  Flusso:
##    1. Carica e prepara reference ABA (spot ISH atlas)
##    2. Label transfer ABA -> experiment.merged (area anatomica per spot)
##    3. Assegna area.spot (label discreto) per ogni spot
##    4. Plot spaziali per area anatomica
##    5. Salvataggio
################################################################################

source("code/ST_pipeline/ST_00_config.R")

load(file.path(outdir, "ST_02_experiment_merged_clustered.RData"))
load(file.path(outdir, "ST_02_df_toplot.RData"))

# ==============================================================================
# 1. Prepara reference ABA
# ==============================================================================

message("=== ST_03A. Caricamento reference ABA ===")

reference_data <- read.csv(ABA_REF_COUNTS, sep = "\t", row.names = 1)
reference_meta <- read.csv(ABA_REF_META,   sep = "\t", row.names = 1)
reference_data <- t(reference_data)

# Allinea meta a data
reference_meta <- reference_meta[colnames(reference_data), ]

# Filtra sezioni vicine all'esperimento
reference_data <- reference_data[,
  reference_meta$section_index %in% ABA_SECTIONS
]
reference_meta <- reference_meta[
  reference_meta$section_index %in% ABA_SECTIONS, ,
  drop = FALSE
]

cat(sprintf("Reference ABA: %d spot, %d geni (%d sezioni)\n",
            ncol(reference_data), nrow(reference_data),
            length(unique(reference_meta$section_index))))

# Rimuovi geni contaminanti
reference_data <- reference_data[,
  !colnames(reference_data) %in% CONTAMINATED_GENES
]

# Crea oggetto Seurat reference
areas_reference <- CreateSeuratObject(reference_data,
                                      project = "ABA_Reference",
                                      assay   = "RNA")

# QC e filtraggio reference
selected_c_ref <- colnames(areas_reference)[
  Matrix::colSums(GetAssayData(areas_reference, assay = "RNA", layer = "counts")) >= MIN_UMI_REF
]
selected_f_ref <- rownames(areas_reference)[
  Matrix::rowSums(GetAssayData(areas_reference, assay = "RNA", layer = "counts")) >= MIN_EXPR_REF
]

# Rimuovi cluster ABA con < 3 spot (TransferData richiede minimo per classe)
clusters_toremove <- names(which(
  table(reference_meta[
    intersect(colnames(areas_reference), selected_c_ref),
    "ABA_parent"
  ]) < 3
))
spots_toremove <- colnames(areas_reference)[
  reference_meta[colnames(areas_reference), "ABA_parent"] %in% clusters_toremove
]
selected_c_ref <- setdiff(selected_c_ref, spots_toremove)

areas_reference <- subset(areas_reference,
                           features = selected_f_ref,
                           cells    = selected_c_ref)
areas_reference$cluster_name <- reference_meta[
  colnames(areas_reference), "ABA_parent"
]

cat(sprintf("Reference filtrato: %d spot, %d geni, %d aree ABA\n",
            ncol(areas_reference), nrow(areas_reference),
            nlevels(factor(areas_reference$cluster_name))))

# SCTransform + PCA reference
areas_reference <- SCTransform(areas_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

p_umap_ref <- DimPlot(areas_reference, group.by = "cluster_name",
                      label = TRUE, label.size = 2.5, pt.size = 0.5) +
  ggtitle("Reference ABA — aree anatomiche") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")
savepng(p_umap_ref, "ST_03_ABA_reference_UMAP.png", width = 3500, height = 3000)


# ==============================================================================
# 2. Label transfer ABA -> ST
# ==============================================================================

message("=== ST_03B. Label transfer ABA -> ST ===")

# SCTransform ST se non già fatto (normalizzazione per transfer)
if (!"SCT" %in% names(experiment.merged@assays)) {
  experiment.merged <- SCTransform(experiment.merged, assay = "Spatial",
                                   verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
}

anchors_ABA <- FindTransferAnchors(
  reference            = areas_reference,
  query                = experiment.merged,
  normalization.method = "SCT"
)

predictions_ABA <- TransferData(
  anchorset        = anchors_ABA,
  refdata          = areas_reference$cluster_name,
  prediction.assay = TRUE,
  weight.reduction = experiment.merged[["pca"]],
  dims             = 1:N_PCS
)

experiment.merged[["cluster_name"]] <- predictions_ABA
DefaultAssay(experiment.merged)     <- "cluster_name"

cat(sprintf("\nAree ABA nel prediction assay: %d\n",
            nrow(predictions_ABA) - 1))  # -1 per riga "max"


# ==============================================================================
# 3. Assegna area.spot discreto per spot
# ==============================================================================

message("=== ST_03C. Assegna area.spot discreto ===")

# Prendi l'area con score massimo; se max <= soglia -> "Ambiguous"
area_spot_idx <- apply(
  predictions_ABA@data[-nrow(predictions_ABA@data), , drop = FALSE],
  2, which.max
)

# Spot con confidence bassa -> Ambiguous
ambiguous_mask <- predictions_ABA@data["max", ] <= LABEL_TRANSFER_MIN_SCORE
area_spot_idx[ambiguous_mask] <- nrow(predictions_ABA@data)

experiment.merged$area.spot <- rownames(predictions_ABA@data)[area_spot_idx]
experiment.merged$area.spot[experiment.merged$area.spot == "max"] <- "Ambiguous"

cat("\nDistribuzione area.spot:\n")
print(sort(table(experiment.merged$area.spot), decreasing = TRUE))

# Proporzione ambiguous per campione
ambig_by_sample <- experiment.merged@meta.data %>%
  group_by(capture.area) %>%
  summarise(
    pct_ambiguous = round(100 * mean(area.spot == "Ambiguous"), 1),
    .groups = "drop"
  )
cat("\nAmbiguous per campione:\n")
print(ambig_by_sample)


# ==============================================================================
# 4. Aggiorna df e plot spaziali
# ==============================================================================

message("=== ST_03D. Plot spaziali area.spot ===")

df$area.spot <- experiment.merged$area.spot[match(rownames(df),
                                                   Cells(experiment.merged))]

# Aggiungi score continui per area al df
score_mat_ABA <- t(predictions_ABA@data[
  rownames(predictions_ABA@data) != "max", , drop = FALSE
])
shared_spots <- intersect(rownames(df), rownames(score_mat_ABA))
for (col in colnames(score_mat_ABA)) {
  safe_col <- paste0("ABA_", gsub(" |/|-", "_", col))
  df[shared_spots, safe_col] <- score_mat_ABA[shared_spots, col]
}

# Plot area.spot discreta nello spazio
area_levels <- sort(unique(na.omit(df$area.spot)))
area_levels <- c(area_levels[area_levels != "Ambiguous"], "Ambiguous")
n_areas     <- length(area_levels)

pal_area <- setNames(
  c(scales::hue_pal()(n_areas - 1), "grey70"),
  area_levels
)

p_area_space <- ggplot(df, aes(x = x, y = -y,
                               colour = area.spot)) +
  geom_point(size = 0.6, alpha = 0.8) +
  facet_wrap(~ area, ncol = 2) +
  scale_colour_manual(values = pal_area, na.value = "grey90") +
  labs(title = "Aree anatomiche ABA (label transfer)",
       colour = "Area") +
  theme_classic(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
savepng(p_area_space, "ST_03_area_spot_space.png",
        width = 2500, height = 2500)

# Plot score continuo per ogni area ABA (solo 90min per brevità)
aba_cols <- grep("^ABA_", colnames(df), value = TRUE)
df_90min <- subset(df, area %in% c("LSD90min_1", "CTRL90min_1",
                                    "LSD90min_2", "CTRL90min_2"))

dir.create(file.path(outdir, "ST_03_ABA_scores"), showWarnings = FALSE)
for (acol in aba_cols) {
  p_sc <- ggplot(df_90min, aes(x = x, y = -y,
                               colour = .data[[acol]])) +
    geom_point(size = 0.6) +
    facet_wrap(~ area, ncol = 2) +
    scale_colour_viridis_c(option = "D", na.value = "grey90") +
    labs(title = sub("^ABA_", "", acol), colour = "Score") +
    theme_classic(base_size = 8)
  savepng(p_sc,
          file.path(outdir, "ST_03_ABA_scores",
                    paste0(acol, "_90min.png")),
          width = 2300, height = 2000)
}

# Barplot composizione anatomica per campione
comp_df <- df %>%
  filter(!is.na(area.spot)) %>%
  count(area, area.spot) %>%
  group_by(area) %>%
  mutate(pct = 100 * n / sum(n))

p_comp <- ggplot(comp_df, aes(x = area, y = pct, fill = area.spot)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_area) +
  labs(title = "Composizione anatomica per campione",
       x = "", y = "% spot", fill = "Area") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
savepng(p_comp, "ST_03_area_composition_barplot.png",
        width = 3500, height = 2500)


# ==============================================================================
# 5. Salvataggio
# ==============================================================================

save(experiment.merged, areas_reference, predictions_ABA,
     file = file.path(outdir, "ST_03_experiment_merged_ABA.RData"))

save(df, file = file.path(outdir, "ST_03_df_toplot_ABA.RData"))

message("=== ST_03 completato ===")
message("  area.spot assegnato a ",
        sum(!is.na(experiment.merged$area.spot)), " spot")
