################################################################################
##  ST_07_THAL_atlas_mapping.R
##
##  Proietta i cluster THAL sul reference atlas Visium ABA usando label transfer
##  (identico a ST_06 ma con areas_reference come query invece dei dati propri).
##
##  Flusso:
##    A. Ricostruisce areas_reference da tutte le sezioni disponibili
##       (senza filtro ABA_SECTIONS — vogliamo tutte le sezioni)
##    B. SCTransform THAL se non già presente
##    C. FindTransferAnchors + TransferData: THAL → areas_reference
##    D. Plot score continui per cell_type sulle sezioni del reference
##    E. Plot discreto (argmax per spot) su tutte le sezioni
##    F. Approccio B: purity talamica (ABA_parent degli spot con score alto)
##
##  Input:
##    - THAL_celltype.RData
##    - data/Reference/expr_raw_counts_table.tsv
##    - data/Reference/meta_table.tsv
##
##  Output:
##    - ST_07_atlas_mapping/ST_07_<ct>_atlas.png   (score continuo per cell_type)
##    - ST_07_atlas_discrete.png                    (argmax per spot)
##    - ST_07_thalamic_purity.csv / .png            (purity talamica)
##    - ST_07_synthesis.csv
################################################################################

source("code/ST_pipeline/ST_00_config.R")

library(dplyr)
library(ggplot2)

# --------------------------------------------------------------------------
# PERCORSI
# --------------------------------------------------------------------------
load("D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/Psychedelics/4-RNAflex/2_Full_experiment/4_Data_Analysis/RNAflex_full/results/20260224/THAL/THAL_celltype.RData")

atlas_dir <- file.path(outdir, "ST_07_atlas_mapping")
dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)

CACHE_FILE <- file.path(outdir, "ST_07_transfer_cache.RData")

ct_colors <- c(
  "TC_FO"  = "#1F77B4", "TC_FO2" = "#FF7F0E", "TC_FO3" = "#2CA02C",
  "TC_MTX" = "#D62728", "TRN"    = "#9467BD", "TRN2"   = "#8C564B",
  "TRN3"   = "#E377C2", "NA"     = "#7F7F7F", "NA2"    = "#BCBD22",
  "NA3"    = "#17BECF", "Astro"  = "#AEC7E8", "Astro2" = "#FFBB78",
  "Astro3" = "#98DF8A", "ChP"    = "#FF9896", "Oligo"  = "#C5B0D5",
  "OPC"    = "#C49C94", "Micro"  = "#F7B6D2", "Endo"   = "#DBD8A0"
)

# ==============================================================================
# A. Ricostruisce areas_reference (TUTTE le sezioni, non solo ABA_SECTIONS)
# ==============================================================================

message("=== A. Costruzione areas_reference (tutte le sezioni) ===")

reference_data <- read.csv(ABA_REF_COUNTS, sep = "\t", row.names = 1)
reference_meta <- read.csv(ABA_REF_META,   sep = "\t", row.names = 1)
reference_data <- t(reference_data)

# Allinea meta a data
reference_meta <- reference_meta[colnames(reference_data), ]

# NON filtriamo per ABA_SECTIONS: vogliamo tutte le sezioni
message(sprintf("  Sezioni disponibili: %s",
                paste(sort(unique(reference_meta$section_index)), collapse = ", ")))

# Rimuovi geni contaminanti
reference_data <- reference_data[
  !rownames(reference_data) %in% CONTAMINATED_GENES, ]

# Crea oggetto Seurat
areas_reference <- CreateSeuratObject(reference_data,
                                      project = "ABA_Reference",
                                      assay   = "RNA")

# QC
selected_c_ref <- colnames(areas_reference)[
  Matrix::colSums(GetAssayData(areas_reference, assay = "RNA",
                               layer = "counts")) >= MIN_UMI_REF
]
selected_f_ref <- rownames(areas_reference)[
  Matrix::rowSums(GetAssayData(areas_reference, assay = "RNA",
                               layer = "counts")) >= MIN_EXPR_REF
]

# Rimuovi aree con < 3 spot
clusters_toremove <- names(which(
  table(reference_meta[
    intersect(colnames(areas_reference), selected_c_ref), "ABA_parent"
  ]) < 3
))
spots_toremove <- colnames(areas_reference)[
  reference_meta[colnames(areas_reference), "ABA_parent"] %in% clusters_toremove
]
selected_c_ref <- setdiff(selected_c_ref, spots_toremove)

areas_reference <- subset(areas_reference,
                           features = selected_f_ref,
                           cells    = selected_c_ref)

# Aggiungi metadati utili per i plot
areas_reference$ABA_parent    <- reference_meta[colnames(areas_reference), "ABA_parent"]
areas_reference$ABA_acronym   <- reference_meta[colnames(areas_reference), "ABA_acronym"]
areas_reference$section_index <- reference_meta[colnames(areas_reference), "section_index"]
areas_reference$stereo_ML     <- reference_meta[colnames(areas_reference), "stereo_ML"]
areas_reference$stereo_DV     <- reference_meta[colnames(areas_reference), "stereo_DV"]

message(sprintf("  areas_reference: %d spot, %d geni, %d sezioni, %d aree ABA",
                ncol(areas_reference), nrow(areas_reference),
                n_distinct(areas_reference$section_index),
                n_distinct(areas_reference$ABA_parent)))

# SCTransform + PCA per il transfer
options(future.globals.maxSize = 4 * 1024^3)
areas_reference <- SCTransform(areas_reference, ncells = 3000,
                                verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# ==============================================================================
# B. Prepara THAL come reference per il transfer
# ==============================================================================

message("=== B. Preparazione THAL come reference ===")

DefaultAssay(THAL) <- "SCT"
Idents(THAL)       <- "cell_type"

# SCTransform se non presente
if (!"SCT" %in% names(THAL@assays)) {
  message("  SCTransform THAL...")
  options(future.globals.maxSize = 4 * 1024^3)
  THAL <- SCTransform(THAL, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
}

message(sprintf("  THAL: %d cellule, %d cell_type: %s",
                ncol(THAL),
                n_distinct(THAL$cell_type),
                paste(sort(unique(THAL$cell_type)), collapse = ", ")))

# ==============================================================================
# C. Label transfer THAL → areas_reference
# ==============================================================================

if (file.exists(CACHE_FILE)) {
  message("=== C. Carico transfer da cache ===")
  load(CACHE_FILE)   # → predictions_thal
} else {
  message("=== C. FindTransferAnchors + TransferData ===")

  # Feature set condiviso
  transfer_features <- intersect(
    VariableFeatures(THAL),
    rownames(areas_reference[["SCT"]])
  )
  message(sprintf("  Features per transfer: %d", length(transfer_features)))

  DefaultAssay(THAL)            <- "SCT"
  DefaultAssay(areas_reference) <- "SCT"

  anchors_thal <- FindTransferAnchors(
    reference            = THAL,
    query                = areas_reference,
    normalization.method = "SCT",
    features             = transfer_features,
    dims                 = 1:N_PCS,
    verbose              = TRUE
  )

  predictions_thal <- TransferData(
    anchorset        = anchors_thal,
    refdata          = THAL$cell_type,
    prediction.assay = TRUE,
    weight.reduction = areas_reference[["pca"]],
    dims             = 1:N_PCS
  )

  save(predictions_thal, file = CACHE_FILE)
  message("  Transfer salvato in cache: ", CACHE_FILE)
}

# Aggiungi scores a areas_reference
areas_reference[["cell_type_pred"]] <- predictions_thal
DefaultAssay(areas_reference)       <- "cell_type_pred"

# Estrai score in data.frame per i plot
score_mat <- t(as.matrix(
  predictions_thal@data[rownames(predictions_thal@data) != "max", , drop = FALSE]
))

df_atlas <- areas_reference@meta.data %>%
  select(section_index, stereo_ML, stereo_DV, ABA_parent, ABA_acronym) %>%
  cbind(score_mat) %>%
  mutate(
    pred_max   = predictions_thal@data["max", rownames(.)],
    pred_label = rownames(predictions_thal@data)[
      apply(predictions_thal@data[rownames(predictions_thal@data) != "max", ],
            2, which.max)
    ]
  )

# Argmax label: NA se score < soglia
df_atlas$pred_label[df_atlas$pred_max < LABEL_TRANSFER_MIN_SCORE] <- NA

message(sprintf("  Spot con predizione >= %.2f: %d / %d",
                LABEL_TRANSFER_MIN_SCORE,
                sum(!is.na(df_atlas$pred_label)),
                nrow(df_atlas)))

# ==============================================================================
# D. Plot score continui per cell_type su tutte le sezioni
# ==============================================================================

message("=== D. Plot score continui per cell_type ===")

ct_types <- setdiff(rownames(predictions_thal@data), "max")
all_sections <- sort(unique(df_atlas$section_index))
n_rows <- ceiling(length(all_sections) / 4)

colnames(df_atlas) <- gsub("-", "_", colnames(df_atlas))
ct_types           <- gsub("-", "_", ct_types)
df_atlas$pred_label <- gsub("-", "_", df_atlas$pred_label)

# Verifica
colnames(df_atlas)
summary(df_atlas[["TC_FO2"]])

for (ct in ct_types) {

  if (!ct %in% colnames(df_atlas)) next

  col <- ct_colors[ct] %||% "#E63946"

  p <- ggplot(df_atlas,
              aes(x = stereo_ML, y = -stereo_DV,
                  color = .data[[ct]])) +
    geom_point(size = 0.5, alpha = 0.8, shape = 16) +
    scale_color_gradientn(
      colours = c("grey92", "grey70", col),
      values  = c(0, 0.1, 1),
      limits  = c(0, 1),
      name    = "Score",
      oob     = scales::squish
    ) +
    scale_y_reverse() +
    facet_wrap(~ section_index, ncol = 4) +
    coord_equal() +
    labs(
      title    = sprintf("cell_type: %s", ct),
      subtitle = "Score label transfer THAL → reference Visium ABA (tutte le sezioni)"
    ) +
    theme_bw(base_size = 7) +
    theme(
      strip.text      = element_text(size = 6, face = "bold"),
      legend.position = "right",
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", size = 10, color = col),
      plot.subtitle   = element_text(size = 7, color = "grey40")
    )

  fname <- sprintf("ST_07_%s_atlas.png", gsub("[^A-Za-z0-9]", "_", ct))
  savepng(p, file.path(atlas_dir, fname),
          width  = 3000,
          height = 300 * n_rows + 500)
  message(sprintf("  Salvato: %s", fname))
}

# ==============================================================================
# E. Plot discreto argmax su tutte le sezioni
# ==============================================================================

message("=== E. Plot discreto argmax ===")

df_discrete <- df_atlas %>%
  filter(!is.na(pred_label))

df_bg <- df_atlas %>%
  filter(is.na(pred_label))

p_disc <- ggplot() +
  geom_point(data = df_bg,
             aes(x = stereo_ML, y = stereo_DV),
             color = "grey92", size = 0.3, shape = 16) +
  geom_point(data = df_discrete,
             aes(x = stereo_ML, y = stereo_DV, color = pred_label),
             size = 0.7, alpha = 0.8, shape = 16) +
  scale_color_manual(values = ct_colors, na.value = "grey80",
                     name = "cell_type") +
  scale_y_reverse() +
  facet_wrap(~ section_index, ncol = 4) +
  coord_equal() +
  labs(
    title    = "Predizione cell_type THAL — reference Visium ABA (argmax)",
    subtitle = sprintf("Score >= %.2f | Grigio = sotto soglia",
                       LABEL_TRANSFER_MIN_SCORE)
  ) +
  theme_bw(base_size = 7) +
  theme(
    strip.text      = element_text(size = 6, face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 7),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1),
                               nrow = 3))

savepng(p_disc, "ST_07_atlas_discrete.png",
        width  = 5500,
        height = 1200 * n_rows + 700)
message("Salvato: ST_07_atlas_discrete.png")

# ==============================================================================
# F. Purity talamica: per ogni cell_type, % spot con ABA_parent talamico
# ==============================================================================

message("=== F. Purity talamica ===")

# Strutture talamiche da ABA_parent del reference
# (usiamo direttamente ABA_parent, più affidabile di CCFv3 per questo reference)
thal_aba_parents <- unique(
  areas_reference$ABA_parent[grepl(
    "thalamus|Thalamus|TH|thal",
    areas_reference$ABA_parent,
    ignore.case = TRUE
  )]
)
message(sprintf("  Aree ABA talamiche: %s",
                paste(sort(thal_aba_parents), collapse = ", ")))

# Per ogni spot: score massimo e cell_type predetto
# Per ogni cell_type: quali spot hanno score > soglia, e qual è la loro ABA_parent
SCORE_THRESH <- 0.3   # soglia più bassa per purity (non solo argmax)

ct_purity <- lapply(ct_types, function(ct) {
  if (!ct %in% colnames(df_atlas)) return(NULL)

  spots_ct <- df_atlas %>%
    filter(.data[[ct]] >= SCORE_THRESH) %>%
    summarise(
      n_spots    = n(),
      n_thalamic = sum(ABA_parent %in% thal_aba_parents),
      purity_th  = n_thalamic / max(n_spots, 1),
      top_rois   = paste(
        names(sort(table(ABA_parent), decreasing = TRUE)[1:3]),
        collapse = " | "
      )
    ) %>%
    mutate(cell_type = ct)

  spots_ct
}) %>%
  bind_rows() %>%
  arrange(purity_th)

message("Purity talamica per cell_type:")
print(ct_purity %>% mutate(across(starts_with("purity"), ~ round(.x * 100, 1))))

write.csv(ct_purity, file.path(outdir, "ST_07_thalamic_purity.csv"),
          row.names = FALSE)

# Plot purity
p_purity <- ggplot(ct_purity,
                   aes(x = reorder(cell_type, purity_th),
                       y = purity_th, fill = cell_type)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_text(aes(label = sprintf("n=%d\n%s", n_spots, top_rois)),
            hjust = -0.05, size = 2.3, lineheight = 0.85) +
  scale_fill_manual(values = ct_colors, na.value = "grey70") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.4),
                     expand = c(0, 0)) +
  coord_flip() +
  labs(
    title    = "Purity talamica (reference Visium ABA) per cell_type THAL",
    subtitle = sprintf("Spot con score >= %.2f | Linea rossa = 50%%\nTop 3 ROI anatomiche per cell_type",
                       SCORE_THRESH),
    x = NULL, y = "% spot in ROI talamica"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40"))

savepng(p_purity, "ST_07_thalamic_purity.png", width = 5000, height = 3500)
message("Salvato: ST_07_thalamic_purity.png")

# ==============================================================================
# SINTESI
# ==============================================================================

synthesis <- ct_purity %>%
  mutate(
    flag = case_when(
      purity_th >= 0.8 ~ "OK",
      purity_th >= 0.5 ~ "verificare",
      TRUE             ~ "DUBBIO"
    ),
    purity_th = round(purity_th * 100, 1)
  ) %>%
  arrange(purity_th)

message("\n=== SINTESI PURITY ===")
print(synthesis)

write.csv(synthesis, file.path(outdir, "ST_07_synthesis.csv"),
          row.names = FALSE)

message("\n=== ST_07 completato ===")
message("Output in: ", outdir)
message("  ST_07_atlas_mapping/ST_07_<ct>_atlas.png")
message("  ST_07_atlas_discrete.png")
message("  ST_07_thalamic_purity.csv / .png")
message("  ST_07_synthesis.csv")
