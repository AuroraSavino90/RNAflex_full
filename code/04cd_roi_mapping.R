################################################################################
##  04cd_roi_mapping.R                                                         ##
##                                                                             ##
##  Input:  CTX_celltype.RData     (da sezione 03, con cell_type assegnato)   ##
##          CTX_annotated.RData    (da 04ab, per estrarre mmc_cluster)        ##
##          cell_metadata.csv      (ABCA S3, WMB-10X)                         ##
##                                                                             ##
##  Output: CTX_celltype con colonne aggiuntive:                              ##
##            mmc_cluster_alias, mmc_roi, mmc_roi_purity,                     ##
##            mmc_macro_area, mmc_cluster_entropy                              ##
##          04cd_roi_assignments.csv                                           ##
##          Plot validation checks (04D_check*.png / pdf)                     ##
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)

# ==============================================================================
# CARICAMENTO
# ==============================================================================

# CTX_celltype: oggetto con cell_type assegnato (da sezione 03)
load(file.path(outdir, "CTX_celltype.RData"))    # -> CTX

# Trasferisci mmc_cluster da CTX_annotated a CTX_celltype
# (le stesse cellule, ma CTX_annotated ha i risultati MapMyCells)
load(file.path(outdir, "CTX_annotated.RData"))   # -> sovrascrive CTX
CTX_annotated <- CTX
load(file.path(outdir, "CTX_celltype.RData"))    # -> ripristina CTX con cell_type

shared <- intersect(colnames(CTX), colnames(CTX_annotated))
cat("Cellule condivise tra CTX_celltype e CTX_annotated:", length(shared), "\n")

CTX$mmc_cluster    <- NA_character_
CTX$mmc_subclass   <- NA_character_
CTX$mmc_class      <- NA_character_
idx <- match(shared, colnames(CTX))
CTX$mmc_cluster[ idx] <- CTX_annotated$mmc_cluster[ match(shared, colnames(CTX_annotated))]
CTX$mmc_subclass[idx] <- CTX_annotated$mmc_subclass[match(shared, colnames(CTX_annotated))]
CTX$mmc_class[   idx] <- CTX_annotated$mmc_class[   match(shared, colnames(CTX_annotated))]

rm(CTX_annotated)
cat("mmc_cluster trasferito. NA:", sum(is.na(CTX$mmc_cluster)), "\n")


# ==============================================================================
# DEFINIZIONI COMUNI
# ==============================================================================

macro_areas <- list(
  Visual        = c("VIS", "VIS-PTLp"),
  Somatosensory = c("SSp", "SS-GU-VISC"),
  Motor         = c("MOp", "MO-FRP"),
  Auditory      = c("AUD-TEa-PERI-ECT"),
  Prefrontal    = c("PL-ILA-ORB", "ACA"),
  Retrosplenial = c("RSP"),
  Temporal      = c("AUD-TEa-PERI-ECT"),  # stessa ROI; distinto da mmc_subclass
  Insular       = c("AI")
)

area_colors <- c(
  Visual        = "#E63946",
  Somatosensory = "#457B9D",
  Motor         = "#2D6A4F",
  Auditory      = "#F4A261",
  Prefrontal    = "#9B5DE5",
  Retrosplenial = "#00B4D8",
  Temporal      = "#F77F00",
  Insular       = "#588157",
  Unassigned    = "#CCCCCC"
)

area_signatures <- list(
  Visual        = c("Nr4a2", "Pcp4", "Nefm", "Calb1", "Crhr2", "Batf3", "Tbr1", "Ccbe1"),
  Somatosensory = c("Etv1", "Lgr5", "Ldb2", "Gpr26", "Krt80", "Tafa1", "Kcnk9", "Nxph1"),
  Motor         = c("Rbp4", "Crym", "Tshz2", "Hpse", "Bcl6", "Slco2a1", "Syt6", "Colgalt2"),
  Auditory      = c("Fam19a1", "Reln", "Asap2", "Gabrg1", "Kcnmb2", "Sdk2", "Cntn5"),
  Prefrontal    = c("Foxp2", "Zbtb20", "Bcl11b", "Cdh12", "Ptpro", "Nxph3", "Sema3e", "Fezf2"),
  Retrosplenial = c("Nrp1", "Syt6", "Ctgf", "Col6a1", "Adcyap1", "Nts", "Fn1"),
  Temporal      = c("Cdh13", "Cbln2", "Lhx2", "Sema5a", "Nrg1", "Kcnh5", "Prss12"),
  Insular       = c("Htr2a", "Tac1", "Npy", "Vip", "Nos1", "Chrm2", "Pcp4l1")
)

neuronal_types <- c("Ex_L23_1", "Ex_L23_2", "Ex_L5PT_1", "Ex_L5PT_2",
                    "Ex_L6CT_1", "Ex_L6CT_2", "Ex_L6CT_3",
                    "Ex1", "Ex2", "Ex3",
                    "Inh_PV", "Inh_SST", "Inh_VIP")

CELL_META_CSV <- "data/cell_metadata.csv"


# ==============================================================================
# PARTE C — ASSEGNAZIONE ROI DA MAPMYCELLS CLUSTER
# ==============================================================================

message("=== 04C. Assegnazione ROI da cluster MapMyCells ===")

# ── 1. Carica cell_metadata e filtra isocortex ────────────────────────────────
meta_full <- read.csv(CELL_META_CSV, row.names = 1)
meta_iso  <- meta_full[meta_full$feature_matrix_label %in%
                         c("WMB-10Xv3-Isocortex-1", "WMB-10Xv3-Isocortex-2"), ]
cat("Cellule isocortex in cell_metadata:", nrow(meta_iso), "\n")

# ── 2. Costruisci mapping cluster_alias -> ROI maggioritaria ──────────────────
cluster_roi_map <- meta_iso %>%
  filter(!is.na(region_of_interest_acronym), !is.na(cluster_alias)) %>%
  group_by(cluster_alias, region_of_interest_acronym) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster_alias) %>%
  mutate(total        = sum(n),
         pct          = round(100 * n / total, 1),
         roi_majority = region_of_interest_acronym[which.max(n)]) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(cluster_alias, roi_majority, pct_majority = pct, n_cells_in_roi = n) %>%
  ungroup()

cat("Cluster totali con ROI assegnata:", nrow(cluster_roi_map), "\n")
cat("Purezza mapping (% cellule nella ROI maggioritaria):\n")
print(summary(cluster_roi_map$pct_majority))
cat("Cluster con purezza >80%:", sum(cluster_roi_map$pct_majority > 80), "\n")
cat("Cluster con purezza >50%:", sum(cluster_roi_map$pct_majority > 50), "\n")

# ── 3. Estrai cluster_alias da mmc_cluster ────────────────────────────────────
# Verifica se cluster_alias è direttamente disponibile nel CSV MapMyCells
MAPPING_CSV <- file.path(outdir, "CTX_from_mapmycells/CTX_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1772014289213.csv")
mapping_raw <- read.csv(MAPPING_CSV, comment.char = "#")

if ("cluster_alias" %in% colnames(mapping_raw)) {
  # CASO A: cluster_alias diretto nel CSV
  cat("Usando colonna cluster_alias diretta dal CSV MapMyCells\n")
  alias_vec <- setNames(as.character(mapping_raw$cluster_alias), mapping_raw$cell_id)
  CTX$mmc_cluster_alias <- alias_vec[colnames(CTX)]
} else {
  # CASO B: estrai numero iniziale dal nome cluster (es. "0042 L2/3 IT CTX..." -> "42")
  cat("Estraendo cluster_alias dal numero iniziale del cluster_name\n")
  CTX$mmc_cluster_alias <- as.character(
    as.integer(sub("^(\\d+).*", "\\1", CTX$mmc_cluster))
  )
}

cat("Esempi mmc_cluster_alias:\n")
print(head(unique(CTX$mmc_cluster_alias), 10))

# ── 4. Join: cluster_alias -> ROI ─────────────────────────────────────────────
CTX$mmc_roi         <- cluster_roi_map$roi_majority[
  match(CTX$mmc_cluster_alias, cluster_roi_map$cluster_alias)]
CTX$mmc_roi_purity  <- cluster_roi_map$pct_majority[
  match(CTX$mmc_cluster_alias, cluster_roi_map$cluster_alias)]

n_mapped   <- sum(!is.na(CTX$mmc_roi))
n_unmapped <- sum( is.na(CTX$mmc_roi))
cat(sprintf("\nCellule con ROI assegnata: %d (%.1f%%)\n",
            n_mapped, 100 * n_mapped / ncol(CTX)))
cat(sprintf("Cellule senza ROI:         %d (%.1f%%)\n",
            n_unmapped, 100 * n_unmapped / ncol(CTX)))
cat("\nDistribuzione ROI fine:\n")
print(table(CTX$mmc_roi, useNA = "always"))

# ── 5. Aggrega ROI -> macro-area ──────────────────────────────────────────────
roi_to_macro <- stack(macro_areas)
colnames(roi_to_macro) <- c("roi", "macro_area")
roi_to_macro$roi       <- as.character(roi_to_macro$roi)
roi_to_macro$macro_area <- as.character(roi_to_macro$macro_area)
# Rimuovi duplicati (AUD-TEa-PERI-ECT compare sia in Auditory che Temporal)
# -> primo match = Auditory; Temporal va distinto via mmc_subclass
roi_to_macro <- roi_to_macro[!duplicated(roi_to_macro$roi), ]

CTX$mmc_macro_area <- roi_to_macro$macro_area[match(CTX$mmc_roi, roi_to_macro$roi)]
CTX$mmc_macro_area <- as.character(CTX$mmc_macro_area)
CTX$mmc_macro_area[is.na(CTX$mmc_macro_area)] <- "Unassigned"

cat("\nDistribuzione macro-aree:\n")
print(table(CTX$mmc_macro_area))

# ── 6. Plot diagnostici parte C ───────────────────────────────────────────────
p_roi <- DimPlot(CTX, group.by = "mmc_macro_area",
                 cols = area_colors, label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Area Corticale (MapMyCells cluster -> ROI)") +
  theme_classic(base_size = 12)
savepng(p_roi, "04C_UMAP_mmc_macro_area.png", width = 4000, height = 3500)

roi_levels <- sort(unique(na.omit(CTX$mmc_roi)))
pal_roi    <- setNames(scales::hue_pal()(length(roi_levels)), roi_levels)
p_roi_fine <- DimPlot(CTX, group.by = "mmc_roi", cols = pal_roi,
                      label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 2.5) +
  ggtitle("ROI ABCA fine (MapMyCells cluster)") +
  theme_classic(base_size = 12)
savepng(p_roi_fine, "04C_UMAP_mmc_roi_fine.png", width = 5000, height = 3500)

if ("cell_type" %in% colnames(CTX@meta.data)) {
  comp_df <- CTX@meta.data %>%
    filter(!is.na(mmc_roi)) %>%
    count(cell_type, mmc_roi) %>%
    group_by(cell_type) %>%
    mutate(pct = n / sum(n) * 100)
  p_bar <- ggplot(comp_df, aes(x = cell_type, y = pct, fill = mmc_roi)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Composizione ROI per cell type", x = "", y = "% cellule", fill = "ROI") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  savepng(p_bar, "04C_barplot_roi_by_celltype.png", width = 4000, height = 2500)
}

# ── 7. Salvataggio intermedio -------------------------------------------------
save(CTX, cluster_roi_map, meta_iso,
     file = file.path(outdir, "CTX_celltype_roi.RData"))

write.csv(
  CTX@meta.data %>%
    select(any_of(c("cell_type", "mmc_class", "mmc_subclass", "mmc_cluster",
                    "mmc_cluster_alias", "mmc_roi", "mmc_roi_purity", "mmc_macro_area"))),
  file.path(outdir, "04C_roi_assignments.csv"),
  quote = FALSE
)

message(sprintf("=== 04C completato: ROI assegnata a %d cellule (%.1f%%) ===",
                n_mapped, 100 * n_mapped / ncol(CTX)))


# ==============================================================================
# PARTE D — VALIDATION CHECKS
# ==============================================================================

message("=== 04D. Validation checks ROI MapMyCells ===")

# Normalizza se necessario
  CTX <- NormalizeData(CTX, assay = "RNA", verbose = FALSE)
  CTX <- JoinLayers(CTX, assay = "RNA")


# ── CHECK 1: marker areali -> espressione coerente con mmc_macro_area ─────────
cat("=== CHECK 1: espressione marker per mmc_macro_area ===\n")

marker_genes <- c("Rbp4", "Crym", "Nr4a2", "Etv1", "Foxp2", "Nrp1")
marker_genes <- marker_genes[marker_genes %in% rownames(CTX)]

expr_df <- as.data.frame(
  t(as.matrix(GetAssayData(CTX, assay = "RNA", layer = "data")[marker_genes, ]))
)
expr_df$mmc_macro_area <- CTX$mmc_macro_area
expr_df$cell_type      <- CTX$cell_type

expr_neurons_df <- expr_df[
  expr_df$cell_type %in% neuronal_types &
    !is.na(expr_df$mmc_macro_area) &
    expr_df$mmc_macro_area != "Unassigned", ]

marker_expected <- c(Rbp4 = "Motor", Crym = "Motor", Nr4a2 = "Visual",
                     Etv1 = "Somatosensory", Foxp2 = "Prefrontal", Nrp1 = "Retrosplenial")

plots_check1 <- lapply(marker_genes, function(g) {
  ggplot(expr_neurons_df,
         aes_string(x = "mmc_macro_area", y = g, fill = "mmc_macro_area")) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = area_colors, guide = "none") +
    labs(title = paste0(g, "  [atteso: ", marker_expected[g], "]"),
         x = "", y = "log-norm") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
})
savepng(wrap_plots(plots_check1, ncol = 3),
        "04D_check1_markers_by_roi.png", width = 5000, height = 3500)

cat("Media Rbp4 per macro-area (solo neuroni):\n")
print(round(tapply(expr_neurons_df$Rbp4, expr_neurons_df$mmc_macro_area, mean), 3))
cat("Media Nr4a2 per macro-area (solo neuroni):\n")
print(round(tapply(expr_neurons_df$Nr4a2, expr_neurons_df$mmc_macro_area, mean), 3))

if ("mmc_roi_purity" %in% colnames(CTX@meta.data)) {
  motor_cells  <- expr_df$mmc_macro_area == "Motor"  & !is.na(expr_df$mmc_macro_area)
  visual_cells <- expr_df$mmc_macro_area == "Visual" & !is.na(expr_df$mmc_macro_area)
  cat("\nCor(Rbp4,  mmc_roi_purity | Motor cells): ",
      round(cor(expr_df$Rbp4[motor_cells],
                CTX$mmc_roi_purity[motor_cells],  use = "complete"), 3), "\n")
  cat("Cor(Nr4a2, mmc_roi_purity | Visual cells):",
      round(cor(expr_df$Nr4a2[visual_cells],
                CTX$mmc_roi_purity[visual_cells], use = "complete"), 3), "\n")
}

# ── CHECK 2: non-neuronali -> Unassigned o distribuzione piatta ───────────────
cat("\n=== CHECK 2: cellule non-neuronali -> Unassigned ===\n")

CTX$is_nonneuronal <- !(CTX$cell_type %in% neuronal_types)

nn_summary <- CTX@meta.data %>%
  filter(is_nonneuronal) %>%
  summarise(n_tot          = n(),
            n_unassigned   = sum(mmc_macro_area == "Unassigned" | is.na(mmc_macro_area)),
            pct_unassigned = round(100 * n_unassigned / n_tot, 1))
cat("Non-neuronali totali:    ", nn_summary$n_tot, "\n")
cat("Non-neuronali Unassigned:", nn_summary$n_unassigned,
    sprintf("(%.1f%%)\n", nn_summary$pct_unassigned))
# Atteso: pct_unassigned > 50%

p_nn <- ggplot(
  CTX@meta.data %>% filter(is_nonneuronal, !is.na(mmc_macro_area)),
  aes(x = cell_type, fill = mmc_macro_area)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = area_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title    = "Distribuzione ROI per cellule non-neuronali",
       subtitle = "Atteso: distribuzione piatta o dominanza Unassigned",
       x = "", y = "Proporzione", fill = "ROI") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
savepng(p_nn, "04D_check2_nonneuronal_roi.png", width = 3500, height = 2500)

# ── CHECK 3: purezza mapping cluster -> ROI ───────────────────────────────────
cat("\n=== CHECK 3: purezza mapping cluster -> ROI ===\n")

if ("mmc_roi_purity" %in% colnames(CTX@meta.data)) {
  cat("Distribuzione purezza per cellula:\n")
  print(summary(CTX$mmc_roi_purity))
  cat("% cellule con purezza >80%:",
      round(100 * mean(CTX$mmc_roi_purity > 80, na.rm = TRUE), 1), "%\n")
  cat("% cellule con purezza >50%:",
      round(100 * mean(CTX$mmc_roi_purity > 50, na.rm = TRUE), 1), "%\n")
  
  p_purity <- ggplot(
    CTX@meta.data %>%
      filter(!is.na(mmc_roi_purity), !is.na(mmc_macro_area), mmc_macro_area != "Unassigned"),
    aes(x = mmc_macro_area, y = mmc_roi_purity, fill = mmc_macro_area)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "orange") +
    scale_fill_manual(values = area_colors, guide = "none") +
    labs(title    = "Purezza mapping cluster -> ROI per macro-area",
         subtitle = "Rosso = 80% | Arancio = 50%",
         x = "", y = "% cellule nella ROI maggioritaria") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  savepng(p_purity, "04D_check3_roi_purity.png", width = 3500, height = 2500)
  
  low_purity <- cluster_roi_map %>% filter(pct_majority < 50) %>% arrange(pct_majority)
  cat("\nCluster con purezza < 50% (ROI inaffidabile):\n")
  print(low_purity)
}

# ── CHECK 4: entropia Shannon cluster -> ROI ──────────────────────────────────
cat("\n=== CHECK 4: entropia Shannon cluster -> ROI ===\n")

cluster_roi_entropy <- meta_iso %>%
  filter(!is.na(region_of_interest_acronym), !is.na(cluster_alias)) %>%
  group_by(cluster_alias, region_of_interest_acronym) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster_alias) %>%
  mutate(p = n / sum(n)) %>%
  summarise(entropy = -sum(p * log(p)), n_roi = n_distinct(region_of_interest_acronym),
            .groups = "drop")

CTX$mmc_cluster_entropy <- cluster_roi_entropy$entropy[
  match(CTX$mmc_cluster_alias, cluster_roi_entropy$cluster_alias)]

cat("Entropia media per cellula:", round(mean(CTX$mmc_cluster_entropy, na.rm = TRUE), 3), "\n")
cat("% cellule in cluster con H < 0.5:",
    round(100 * mean(CTX$mmc_cluster_entropy < 0.5, na.rm = TRUE), 1), "%\n")
# Atteso: la maggior parte delle cellule in cluster con entropia bassa

p_entropy <- ggplot(
  CTX@meta.data %>% filter(!is.na(mmc_cluster_entropy), mmc_macro_area != "Unassigned"),
  aes(x = mmc_cluster_entropy, fill = mmc_macro_area)) +
  geom_histogram(bins = 60, alpha = 0.8, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = area_colors) +
  labs(title    = "Entropia Shannon cluster -> ROI per cellula",
       subtitle = paste0("Media = ",
                         round(mean(CTX$mmc_cluster_entropy, na.rm = TRUE), 2),
                         " | Soglia accettabile = 0.5"),
       x = "Entropia di Shannon", y = "N cellule", fill = "Macro-area") +
  theme_classic(base_size = 12)
savepng(p_entropy, "04D_check4_cluster_roi_entropy.png", width = 4000, height = 2500)

# ── CHECK 5: confusion matrix cell_type vs mmc_macro_area ────────────────────
cat("\n=== CHECK 5: confusion matrix cell_type vs mmc_macro_area ===\n")

if ("cell_type" %in% colnames(CTX@meta.data)) {
  conf_tab  <- table(CellType = CTX$cell_type, ROI = CTX$mmc_macro_area)
  conf_norm <- conf_tab / rowSums(conf_tab)
  conf_mat  <- as.matrix(conf_norm)
  area_order <- c(names(macro_areas), "Unassigned")
  conf_mat   <- conf_mat[, area_order[area_order %in% colnames(conf_mat)], drop = FALSE]
  
  pdf(file.path(outdir, "04D_check5_confusion_matrix.pdf"), width = 14, height = 10)
  draw(Heatmap(
    conf_mat,
    name           = "Proporzione",
    top_annotation = HeatmapAnnotation(
      Area = colnames(conf_mat),
      col  = list(Area = area_colors[colnames(conf_mat)])),
    col               = colorRamp2(c(0, 0.5, 1), c("white", "#457B9D", "#1A3A5C")),
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 9),
    column_names_gp   = gpar(fontsize = 10, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- conf_mat[i, j]
      if (val > 0.15)
        grid.text(sprintf("%.2f", val), x, y,
                  gp = gpar(fontsize = 7, col = ifelse(val > 0.5, "white", "black")))
    },
    border          = TRUE,
    column_title    = "Confusion matrix: cell type vs mmc_macro_area",
    column_title_gp = gpar(fontsize = 13, fontface = "bold")
  ), merge_legend = TRUE)
  dev.off()
  
  cat("Top ROI per cell type:\n")
  top_area <- apply(conf_norm, 1, function(x) names(which.max(x)))
  print(data.frame(CellType = names(top_area), TopROI = top_area))
  # Atteso:
  #   Ex_L5PT_* -> Motor o Visual
  #   Ex_L23_*  -> Visual, Somatosensory, Prefrontal
  #   Ex_L6CT_* -> Motor, Prefrontal
  #   Inh_*     -> distribuzione più ampia
  #   Astro/Micro/Oligo -> Unassigned
}

# ── CHECK 6: concordanza mmc_subclass vs mmc_macro_area ──────────────────────
cat("\n=== CHECK 6: concordanza subclass MapMyCells vs mmc_macro_area ===\n")

if ("mmc_subclass" %in% colnames(CTX@meta.data)) {
  subclass_roi <- CTX@meta.data %>%
    filter(!is.na(mmc_subclass), mmc_subclass != "other",
           !is.na(mmc_macro_area), mmc_macro_area != "Unassigned") %>%
    count(mmc_subclass, mmc_macro_area) %>%
    group_by(mmc_subclass) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    arrange(desc(n)) %>%
    select(mmc_subclass, top_macro_area = mmc_macro_area, pct, n_cells = n)
  
  cat("Top macro-area per subclass (prime 20):\n")
  print(head(subclass_roi, 20))
  cat("\nSubclass con distribuzione areale non dominante (top area < 40%):\n")
  print(subclass_roi %>% filter(pct < 40))
}

# ── SUMMARY ───────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════\n")
cat("RIEPILOGO VALIDATION CHECKS ROI (04D)\n")
cat("════════════════════════════════════════\n")
cat("Cellule totali:             ", ncol(CTX), "\n")
cat("Cellule con ROI assegnata:  ",
    sum(!is.na(CTX$mmc_macro_area) & CTX$mmc_macro_area != "Unassigned"),
    sprintf("(%.1f%%)\n",
            100 * mean(!is.na(CTX$mmc_macro_area) & CTX$mmc_macro_area != "Unassigned")))
if ("mmc_roi_purity" %in% colnames(CTX@meta.data))
  cat("Purezza media mapping:      ",
      round(mean(CTX$mmc_roi_purity, na.rm = TRUE), 1), "%\n")
cat("Distribuzione macro-aree:\n")
print(table(CTX$mmc_macro_area, useNA = "always"))

# ==============================================================================
# SALVATAGGIO FINALE
# ==============================================================================

save(CTX, file = file.path(outdir, "CTX_celltype_roi.RData"))

write.csv(
  CTX@meta.data %>%
    select(any_of(c("cell_type", "mmc_class", "mmc_subclass", "mmc_cluster",
                    "mmc_cluster_alias", "mmc_roi", "mmc_roi_purity",
                    "mmc_macro_area", "mmc_cluster_entropy"))),
  file.path(outdir, "04CD_roi_assignments_final.csv"),
  quote = FALSE
)

message("=== 04CD completato ===")
message("Output: CTX_celltype_roi.RData, 04CD_roi_assignments_final.csv")
