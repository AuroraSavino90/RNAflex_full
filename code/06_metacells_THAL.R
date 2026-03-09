################################################################################
##  06_metacells_THAL.R                                                        ##
##                                                                             ##
##  Input:  THAL/THAL_celltype.RData  (da 05_cell_type_annotation_THAL)      ##
##  Output: THAL/Metacells_THAL.RData  (SC_THAL, SC_THAL_seurat)             ##
##          THAL/metacells/*.png                                               ##
##                                                                             ##
##  Visualizzazioni (output identici a 06_metacells.R, area THAL):           ##
##    01. Purity — istogramma + violin per CT                                 ##
##    02. Dimensione metacelle — distribuzione per CT                         ##
##    03. UMAP metacelle (PCA → UMAP sul Seurat delle metacelle)              ##
##    04. Overlay metacelle sull'UMAP originale                               ##
##    05. Composizione CT × campione                                           ##
##    06. Heatmap espressione geni chiave talamici (z-score per CT)           ##
##    07. VlnPlot Fos (IEG principale) per metacelle split per trattamento   ##
##    08. Boxplot Ctrl vs LSD per geni chiave talamici                        ##
##    09. Correlazione singola cellula vs metacella                           ##
##    10. QC metacelle (nFeature, nCount, percent.mt)                        ##
################################################################################

source("code/00_config.R")

outdir   <- "results/20260224"
thal_dir <- file.path(outdir, "THAL")
dir.create(thal_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(thal_dir, "THAL_celltype.RData"))   # THAL

# ==============================================================================
# PARAMETRI
# ==============================================================================

GAMMA   <- 10   # ~1 metacella ogni GAMMA cellule

mc_dir <- file.path(thal_dir, "metacells")
dir.create(mc_dir, recursive = TRUE, showWarnings = FALSE)

savemc <- function(p, filename, width = 3000, height = 2500, res = 300) {
  png(file.path(mc_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

get_ct_colors <- function(ct_levels) {
  n <- length(ct_levels)
  setNames(colorRampPalette(pal_ct)(n), sort(ct_levels))
}

# ==============================================================================
# 06.1 — AGGREGAZIONE SUPERCELL
# ==============================================================================

message("=== 06.1 Metacells THAL ===")

SC_THAL <- SCimplify_from_embedding(
  X     = Embeddings(THAL, "harmony"),
  k.knn = 20,
  gamma = GAMMA
)

message(sprintf("  THAL: %d cellule → %d metacelle (gamma=%d)",
                ncol(THAL), SC_THAL$n.supercells, GAMMA))

sct_mat <- as.matrix(LayerData(THAL, assay = "SCT", layer = "data"))
SC.GE   <- supercell_GE(sct_mat, SC_THAL$membership)
colnames(SC.GE) <- paste0("MC_", seq_len(ncol(SC.GE)))

# Metadati per metacella
meta_cols <- c("cell_type","orig.ident","treatment","time",
               "mmc_class","mmc_subclass","mmc_roi","mmc_macro_area")
for (col in meta_cols) {
  if (col %in% colnames(THAL@meta.data))
    SC_THAL[[col]] <- supercell_assign(THAL@meta.data[[col]],
                                       SC_THAL$membership, method = "jaccard")
}

purity <- supercell_purity(
  clusters             = THAL$cell_type,
  supercell_membership = SC_THAL$membership,
  method               = "max_proportion"
)
SC_THAL$purity  <- purity
SC_THAL$mc_size <- as.numeric(table(SC_THAL$membership))

message(sprintf("  Purity THAL — mediana: %.2f | <0.7: %d metacelle (%.1f%%)",
                median(purity),
                sum(purity < 0.7), 100 * mean(purity < 0.7)))

SC_THAL_seurat <- CreateSeuratObject(
  counts       = SC.GE,
  project      = "Metacells_THAL",
  min.cells    = 0,
  min.features = 0
)
for (col in c(meta_cols, "purity", "mc_size")) {
  if (!is.null(SC_THAL[[col]]))
    SC_THAL_seurat@meta.data[[col]] <- SC_THAL[[col]]
}

SC_THAL_seurat <- NormalizeData(SC_THAL_seurat, verbose = FALSE)

# ==============================================================================
# 06.2 — VISUALIZZAZIONI QUALITÀ — THAL
# ==============================================================================

message("\n=== 06.2 Visualizzazioni metacelle THAL ===")

meta       <- SC_THAL_seurat@meta.data
purity     <- SC_THAL_seurat$purity
mc_size    <- SC_THAL_seurat$mc_size
ct_levels  <- sort(unique(meta$cell_type))
pal_ct_thal <- get_ct_colors(ct_levels)

# ------------------------------------------------------------------ #
#  1. PURITY                                                           #
# ------------------------------------------------------------------ #
message("  1. Purity...")

p_hist_pur <- ggplot(data.frame(purity = purity), aes(x = purity)) +
  geom_histogram(bins = 40, fill = "#27AE60", color = "white", alpha = 0.85) +
  geom_vline(xintercept = median(purity), color = "#E74C3C",
             linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0.7, color = "grey50",
             linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = median(purity) + 0.02, y = Inf,
           label = sprintf("mediana=%.2f", median(purity)),
           vjust = 2, hjust = 0, color = "#E74C3C", size = 3) +
  theme_bw(base_size = 10) +
  labs(title    = "Purity metacelle \u2014 THAL",
       subtitle = sprintf("%d metacelle | \u03b3=%d | <0.7: %.1f%%",
                          length(purity), GAMMA, 100 * mean(purity < 0.7)),
       x = "Purity (max proportion)", y = "N metacelle")

p_vln_pur <- ggplot(
  data.frame(purity = purity, cell_type = meta$cell_type),
  aes(x = reorder(cell_type, purity, median), y = purity, fill = cell_type)
) +
  geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = pal_ct_thal, guide = "none") +
  theme_bw(base_size = 9) +
  labs(title = "Purity per cell type \u2014 THAL",
       subtitle = "Linea tratt. = soglia 0.7", x = "", y = "Purity")

savemc(p_hist_pur | p_vln_pur, "01_purity_THAL.png",
       width  = 5500,
       height = max(2200, length(ct_levels) * 100 + 700))

# ------------------------------------------------------------------ #
#  2. DIMENSIONE METACELLE                                             #
# ------------------------------------------------------------------ #
message("  2. Dimensione metacelle...")

size_df <- data.frame(mc_size = mc_size, cell_type = meta$cell_type,
                      treatment = meta$treatment)

p_size_hist <- ggplot(size_df, aes(x = mc_size)) +
  geom_histogram(bins = 30, fill = "#8E44AD", color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean(mc_size), color = "#E74C3C", linetype = "dashed") +
  theme_bw(base_size = 10) +
  labs(title = "Dimensione metacelle \u2014 THAL",
       x = "N cellule per metacella", y = "N metacelle")

p_size_ct <- ggplot(size_df,
                    aes(x = reorder(cell_type, mc_size, median),
                        y = mc_size, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  scale_fill_manual(values = pal_ct_thal, guide = "none") +
  theme_bw(base_size = 9) +
  labs(title = "Dimensione per cell type \u2014 THAL", x = "", y = "N cellule per metacella")

savemc(p_size_hist | p_size_ct, "02_size_THAL.png",
       width  = 5500,
       height = max(2200, length(ct_levels) * 100 + 700))

# ------------------------------------------------------------------ #
#  3. UMAP METACELLE                                                   #
# ------------------------------------------------------------------ #
message("  3. UMAP metacelle...")

SC_THAL_seurat <- FindVariableFeatures(SC_THAL_seurat, nfeatures = 2000, verbose = FALSE)
SC_THAL_seurat <- ScaleData(SC_THAL_seurat, verbose = FALSE)
SC_THAL_seurat <- RunPCA(SC_THAL_seurat,
                         npcs    = min(30, ncol(SC_THAL_seurat) - 1),
                         verbose = FALSE)
SC_THAL_seurat <- RunUMAP(SC_THAL_seurat,
                          dims     = 1:min(20, ncol(SC_THAL_seurat) - 1),
                          verbose  = FALSE,
                          min.dist = 0.3)

p_umap_ct <- DimPlot(SC_THAL_seurat, group.by = "cell_type",
                     pt.size    = pmax(0.5, 3 / sqrt(ncol(SC_THAL_seurat))),
                     label      = TRUE, label.size = 3, repel = TRUE,
                     cols       = pal_ct_thal) +
  ggtitle("UMAP metacelle \u2014 cell type | THAL") +
  theme(legend.text = element_text(size = 8))

p_umap_trt <- DimPlot(SC_THAL_seurat, group.by = "treatment",
                      pt.size = pmax(0.5, 3 / sqrt(ncol(SC_THAL_seurat))),
                      cols    = pal_trt) +
  ggtitle("UMAP metacelle \u2014 trattamento | THAL")

umap_df <- as.data.frame(Embeddings(SC_THAL_seurat, "umap"))
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df <- cbind(umap_df, meta[, c("cell_type","treatment","purity","mc_size")])

p_umap_pur <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=purity, size=mc_size)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Purity") +
  scale_size_continuous(range = c(0.5, 4), name = "N cellule") +
  theme_bw(base_size = 9) +
  labs(title = "UMAP metacelle \u2014 purity e dimensione | THAL")

savemc(
  (p_umap_ct | p_umap_trt) / p_umap_pur +
    plot_annotation(title = "Metacelle UMAP \u2014 THAL"),
  "03_UMAP_panel_THAL.png",
  width = 5500, height = 4500
)

# ------------------------------------------------------------------ #
#  4. OVERLAY METACELLE SU UMAP ORIGINALE                             #
# ------------------------------------------------------------------ #
message("  4. Overlay UMAP...")

if ("umap" %in% names(THAL@reductions)) {

  sc_umap <- as.data.frame(Embeddings(THAL, "umap"))
  colnames(sc_umap) <- c("UMAP1","UMAP2")
  sc_umap$cell_type <- THAL$cell_type

  mc_centroids <- do.call(rbind, lapply(
    seq_len(ncol(SC_THAL_seurat)), function(i) {
      cells_i <- intersect(names(SC_THAL$membership)[SC_THAL$membership == i],
                           rownames(sc_umap))
      if (length(cells_i) == 0) return(NULL)
      colMeans(sc_umap[cells_i, c("UMAP1","UMAP2"), drop = FALSE])
    }
  ))
  mc_df <- as.data.frame(mc_centroids)
  colnames(mc_df) <- c("UMAP1","UMAP2")
  mc_df$cell_type <- meta$cell_type[seq_len(nrow(mc_df))]
  mc_df$purity    <- meta$purity[seq_len(nrow(mc_df))]
  mc_df$mc_size   <- meta$mc_size[seq_len(nrow(mc_df))]

  p_overlay <- ggplot() +
    geom_point(data = sc_umap,
               aes(x=UMAP1, y=UMAP2, color=cell_type),
               size = 0.15, alpha = 0.2) +
    geom_point(data = mc_df,
               aes(x=UMAP1, y=UMAP2, fill=cell_type, size=mc_size, alpha=purity),
               shape = 21, color = "white", stroke = 0.3) +
    scale_color_manual(values = pal_ct_thal, guide = "none") +
    scale_fill_manual(values  = pal_ct_thal, name = "Cell type") +
    scale_size_continuous(range = c(1, 6), name = "N cellule/MC") +
    scale_alpha_continuous(range = c(0.4, 1), name = "Purity") +
    theme_bw(base_size = 9) +
    guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title    = "Metacelle sovrapposte UMAP originale \u2014 THAL",
         subtitle = "Punti grandi = metacelle (dim \u221d N cell, opacity \u221d purity)",
         x = "UMAP 1", y = "UMAP 2")

  savemc(p_overlay, "04_UMAP_overlay_THAL.png", width = 3500, height = 3000)
}

# ------------------------------------------------------------------ #
#  5. COMPOSIZIONE CT × CAMPIONE                                       #
# ------------------------------------------------------------------ #
message("  5. Composizione...")

comp_df <- meta %>%
  dplyr::count(cell_type, orig.ident, treatment) %>%
  group_by(orig.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_comp_bar <- ggplot(comp_df, aes(x=orig.ident, y=frac, fill=cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ treatment, scales = "free_x") +
  scale_fill_manual(values = pal_ct_thal) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7)) +
  labs(title = "Composizione metacelle per campione \u2014 THAL",
       x = "", y = "Proporzione", fill = "Cell type")

savemc(p_comp_bar, "05_composition_bar_THAL.png", width = 3500, height = 2500)

ct_sample_tab <- table(meta$cell_type, meta$orig.ident)
if (nrow(ct_sample_tab) > 1) {
  png(file.path(mc_dir, "05_composition_heatmap_THAL.png"),
      width  = 3000,
      height = max(2000, nrow(ct_sample_tab) * 80 + 600),
      res    = 300)
  pheatmap(as.matrix(ct_sample_tab),
           display_numbers = TRUE, number_format = "%d",
           color           = colorRampPalette(c("white","#2E86C1"))(50),
           cluster_cols    = FALSE, fontsize = 9,
           main            = "N metacelle CT \u00d7 campione \u2014 THAL")
  dev.off()
}

# ------------------------------------------------------------------ #
#  6. HEATMAP ESPRESSIONE GENI CHIAVE TALAMICI (z-score per CT)       #
# ------------------------------------------------------------------ #
message("  6. Espressione geni chiave talamici...")

genes_key_thal <- filter_genes(
  c("Slc17a6", "Syt2", "Ntng1",            # pan-TC glutamatergico
    "Tnnt1",   "Prkcd",  "Scnn1a", "Grik4", # TC_FO first-order
    "Necab1",  "Cbln1",  "Cbln2",  "Calb1", # TC_HO higher-order
    "Cplx3",   "Htr2c",  "Penk",            # TC_MTX matrix/intralaminar
    "Gad1",    "Gad2",   "Pvalb",  "Etv1",  # TRN
    "Fos",     "Arc",    "Egr1",   "Nr4a1", # IEG
    "Aqp4",    "Mog",    "Cx3cr1", "Pdgfra","Cldn5"), # non-neuronali
  SC_THAL_seurat
)

if (length(genes_key_thal) >= 5) {
  avg_expr <- AverageExpression(
    SC_THAL_seurat, features = genes_key_thal,
    group.by = "cell_type", assays = "RNA",
    layer = "data", verbose = FALSE
  )$RNA

  mat_z <- t(scale(t(as.matrix(avg_expr))))
  mat_z[is.nan(mat_z)] <- 0
  lim     <- max(abs(mat_z), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))

  ht <- Heatmap(
    mat_z, name = "z-score", col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_names_gp     = gpar(fontsize = 8),
    column_names_gp  = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title     = "Espressione media metacelle (z-score) \u2014 THAL",
    heatmap_legend_param = list(title = "z-score")
  )
  png(file.path(mc_dir, "06_expr_heatmap_THAL.png"),
      width  = 3500,
      height = max(2000, length(genes_key_thal) * 55 + 600),
      res    = 300)
  draw(ht)
  dev.off()
}

# ------------------------------------------------------------------ #
#  7. VlnPlot Fos split per trattamento                               #
# ------------------------------------------------------------------ #
if ("Fos" %in% rownames(SC_THAL_seurat)) {
  message("  7. VlnPlot Fos...")
  p_fos <- VlnPlot(SC_THAL_seurat, features = "Fos",
                   group.by = "cell_type", split.by = "treatment",
                   pt.size = 0.3, cols = pal_trt) +
    ggtitle("Fos nelle metacelle \u2014 THAL") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  savemc(p_fos, "07_Fos_violin_THAL.png",
         width  = max(3000, length(ct_levels) * 280 + 800),
         height = 2500)
}

# ------------------------------------------------------------------ #
#  8. BOXPLOT Ctrl vs LSD — geni chiave talamici                      #
# ------------------------------------------------------------------ #
message("  8. Ctrl vs LSD nelle metacelle talamiche...")

genes_lsd_thal <- filter_genes(
  c("Fos",    "Arc",    "Egr1",   "Nr4a1",  "Fosb",   "Junb",
    "Htr2a",  "Htr2c",  "Htr2b",                       # recettori 5-HT (target LSD)
    "Tnnt1",  "Prkcd",  "Necab1", "Cplx3",             # identità TC per CT
    "Gad1",   "Pvalb",                                  # TRN
    "Bdnf",   "Ntrk2",  "Camk2a"),
  SC_THAL_seurat
)

if (length(genes_lsd_thal) > 0 && "treatment" %in% colnames(meta)) {
  expr_long <- as.data.frame(
    t(as.matrix(GetAssayData(SC_THAL_seurat, assay = "RNA",
                             layer = "data")[genes_lsd_thal, , drop = FALSE]))
  ) %>%
    mutate(treatment = meta$treatment,
           cell_type = meta$cell_type,
           mc_size   = meta$mc_size) %>%
    pivot_longer(cols = all_of(genes_lsd_thal),
                 names_to  = "gene",
                 values_to = "expr")

  p_box_lsd <- ggplot(expr_long, aes(x=treatment, y=expr, fill=treatment)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.4, alpha = 0.8) +
    geom_jitter(aes(size = mc_size), width = 0.15, alpha = 0.3, color = "grey40") +
    facet_grid(gene ~ cell_type, scales = "free_y") +
    scale_fill_manual(values = pal_trt) +
    scale_size_continuous(range = c(0.3, 2), guide = "none") +
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text  = element_text(size = 6),
          legend.position = "none") +
    labs(title = "Espressione geni chiave \u2014 Ctrl vs LSD | THAL",
         x = "", y = "Espressione (log-norm)")

  savemc(p_box_lsd, "08_expr_ctrl_lsd_THAL.png",
         width  = max(3000, length(ct_levels) * 300 + 800),
         height = max(3000, length(genes_lsd_thal) * 250 + 800))
}

# ------------------------------------------------------------------ #
#  9. CORRELAZIONE SINGOLA CELLULA vs METACELLA                       #
# ------------------------------------------------------------------ #
message("  9. Correlazione SC vs metacella THAL...")

genes_cor_thal <- filter_genes(
  c("Slc17a6","Tnnt1","Prkcd","Scnn1a","Grik4",   # TC_FO
    "Necab1","Cbln1","Calb1",                       # TC_HO
    "Cplx3","Htr2c",                                # TC_MTX
    "Gad1","Pvalb",                                 # TRN
    "Fos","Arc","Aqp4","Mog"),                      # IEG + gliali
  SC_THAL_seurat
)
genes_cor_thal <- intersect(genes_cor_thal, rownames(THAL))

if (length(genes_cor_thal) >= 5) {

  mc_expr <- as.matrix(
    GetAssayData(SC_THAL_seurat, assay = "RNA", layer = "data")[genes_cor_thal, , drop=FALSE]
  )
  sc_expr <- as.matrix(
    LayerData(THAL, assay = "SCT", layer = "data")[genes_cor_thal, , drop=FALSE]
  )

  mc_ids <- sort(unique(SC_THAL$membership))
  sc_avg <- do.call(cbind, lapply(mc_ids, function(i) {
    cells_i <- intersect(names(SC_THAL$membership)[SC_THAL$membership == i],
                         colnames(sc_expr))
    if (length(cells_i) == 0) return(rep(0, nrow(sc_expr)))
    rowMeans(sc_expr[, cells_i, drop = FALSE])
  }))
  colnames(sc_avg) <- paste0("MC_", mc_ids)

  shared_mc <- intersect(colnames(mc_expr), colnames(sc_avg))
  if (length(shared_mc) >= 10) {
    cor_vals <- sapply(genes_cor_thal, function(g) {
      tryCatch(cor(sc_avg[g, shared_mc], mc_expr[g, shared_mc], method = "pearson"),
               error = function(e) NA_real_)
    })

    cor_df <- data.frame(gene = genes_cor_thal, pearson_r = cor_vals) %>%
      arrange(desc(pearson_r))
    write.csv(cor_df, file.path(mc_dir, "09_SC_MC_correlation_THAL.csv"), row.names = FALSE)

    p_cor <- ggplot(cor_df, aes(x=reorder(gene, pearson_r), y=pearson_r, fill=pearson_r)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
      coord_flip() +
      scale_fill_gradient2(low="#2471A3", mid="grey90", high="#CB4335",
                           midpoint=0.5, name="r") +
      theme_bw(base_size = 9) +
      labs(title    = "Correlazione SC vs metacella \u2014 THAL",
           subtitle = "Linea tratt. = r=0.9",
           x = "", y = "Pearson r")

    savemc(p_cor, "09_SC_MC_correlation_THAL.png",
           width  = 2500,
           height = max(1800, length(genes_cor_thal) * 70 + 500))
  }
}

# ------------------------------------------------------------------ #
# 10. QC METACELLE                                                    #
# ------------------------------------------------------------------ #
message("  10. QC metacelle THAL...")

SC_THAL_seurat[["percent.mt"]] <- PercentageFeatureSet(SC_THAL_seurat, pattern = "^mt-")
p_qc <- VlnPlot(
  SC_THAL_seurat,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  group.by = "cell_type",
  pt.size  = 0.2, ncol = 3,
  cols     = pal_ct_thal
) + plot_annotation(
  title    = "QC metacelle \u2014 THAL",
  subtitle = sprintf("%d metacelle totali", ncol(SC_THAL_seurat))
)
savemc(p_qc, "10_QC_THAL.png",
       width  = max(3000, length(ct_levels) * 280 + 800),
       height = 2800)

message(sprintf("  Completato THAL: %d metacelle", ncol(SC_THAL_seurat)))

# ==============================================================================
# 06.3 — SALVATAGGIO
# ==============================================================================

message("=== 06.3 Salvataggio ===")

save(SC_THAL, SC_THAL_seurat,
     file = file.path(thal_dir, "Metacells_THAL.RData"))

message(sprintf("Salvato: THAL/Metacells_THAL.RData | %d metacelle", ncol(SC_THAL_seurat)))
