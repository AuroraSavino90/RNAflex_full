################################################################################
##  06_metacells.R                                                             ##
##                                                                             ##
##  Input:  CTX.RData, THAL_clean.RData  (da 05_cell_type_annotation)  ##
##  Output: Metacells_CTX.RData  (SC_CTX, SC_CTX_seurat)                     ##
##          Metacells_THAL.RData (SC_THAL, SC_THAL_seurat)                   ##
##                                                                             ##
##  Flusso:                                                                    ##
##    06.1  Aggregazione SuperCell (CTX + THAL)                               ##
##    06.2  Visualizzazioni qualità metacelle (per area):                     ##
##          1. Purity — istogramma + violin per CT                            ##
##          2. Dimensione metacelle — distribuzione per CT                    ##
##          3. UMAP metacelle (PCA → UMAP sul Seurat delle metacelle)         ##
##          4. Overlay metacelle sull'UMAP originale                          ##
##          5. Composizione — barplot + heatmap CT × campione                 ##
##          6. Heatmap espressione geni chiave (z-score per CT)               ##
##          7. VlnPlot Htr2a per metacelle split per trattamento              ##
##          8. Boxplot Ctrl vs LSD per geni chiave                            ##
##          9. Correlazione singola cellula vs metacella                      ##
##         10. QC metacelle (nFeature, nCount, percent.mt)                   ##
##    06.3  Salvataggio                                                        ##
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_celltype.RData"))    # CTX
#load(file.path(outdir, "THAL_clean.RData"))   # THAL_clean

# ==============================================================================
# PARAMETRI
# ==============================================================================

GAMMA <- 10   # fattore aggregazione: ~1 metacella ogni GAMMA cellule
# es. 10.000 cellule / GAMMA=20 → ~500 metacelle

# Sottodirectory per le visualizzazioni
mc_dir <- file.path(outdir, "metacells")
dir.create(mc_dir, recursive = TRUE, showWarnings = FALSE)

savemc <- function(p, filename, width = 3000, height = 2500, res = 300) {
  png(file.path(mc_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# Palette locale per trattamento (ridondante con 00_config ma esplicita)
get_ct_colors <- function(ct_levels) {
  n <- length(ct_levels)
  setNames(colorRampPalette(pal_ct)(n), sort(ct_levels))
}

# ==============================================================================
# 06.1 — AGGREGAZIONE SUPERCELL (CTX + THAL)
# ==============================================================================


  area   <- "CTX"       # "CTX" o "THAL"
  sc_nm  <- paste0("SC_", area)                # "SC_CTX"
  ser_nm <- paste0("SC_", area, "_seurat")     # "SC_CTX_seurat"
  obj    <- get("CTX")
  
  message(sprintf("=== 06.1 Metacells %s ===", area))
  
  # ---- Aggregazione --------------------------------------------------------
  SC <- SCimplify_from_embedding(
    X     = Embeddings(obj, "harmony"),
    k.knn = 20,
    gamma = GAMMA
  )
  
  message(sprintf("  %s: %d cellule → %d metacelle (gamma=%d)",
                  area, ncol(obj), SC$n.supercells, GAMMA))
  
  # ---- Matrice espressione metacelle (SCT data layer) ----------------------
  sct_mat <- as.matrix(LayerData(obj, assay = "SCT", layer = "data"))
  SC.GE   <- supercell_GE(sct_mat, SC$membership)
  colnames(SC.GE) <- paste0("MC_", seq_len(ncol(SC.GE)))
  
  # ---- Metadati per metacella (jaccard majority vote) ----------------------
  meta_cols <- c("cell_type","orig.ident","treatment","time","area",
                 "mmc_class","mmc_subclass","broad_anat","fine_anat")
  for (col in meta_cols) {
    if (col %in% colnames(obj@meta.data))
      SC[[col]] <- supercell_assign(obj@meta.data[[col]],
                                    SC$membership, method = "jaccard")
  }
  
  # ---- Purity --------------------------------------------------------------
  purity <- supercell_purity(
    clusters             = obj$cell_type,
    supercell_membership = SC$membership,
    method               = "max_proportion"
  )
  SC$purity   <- purity
  SC$mc_size  <- as.numeric(table(SC$membership))
  
  message(sprintf("  Purity %s — mediana: %.2f | <0.7: %d metacelle (%.1f%%)",
                  area, median(purity),
                  sum(purity < 0.7), 100 * mean(purity < 0.7)))
  
  # ---- Oggetto Seurat delle metacelle --------------------------------------
  SC_seurat <- CreateSeuratObject(
    counts       = SC.GE,
    project      = paste0("Metacells_", area),
    min.cells    = 0,
    min.features = 0
  )
  for (col in c(meta_cols, "purity", "mc_size")) {
    if (!is.null(SC[[col]]))
      SC_seurat@meta.data[[col]] <- SC[[col]]
  }
  
  SC_seurat <- NormalizeData(SC_seurat, verbose = FALSE)
  SC_CTX_seurat<-SC_seurat
 
# ==============================================================================
# 06.2 — VISUALIZZAZIONI QUALITÀ — CTX
# ==============================================================================

message("\n=== 06.2 Visualizzazioni metacelle — CTX ===")

meta        <- SC_CTX_seurat@meta.data
purity      <- SC_CTX_seurat$purity
mc_size     <- SC_CTX_seurat$mc_size
ct_levels   <- sort(unique(meta$cell_type))
pal_ct_ctx  <- get_ct_colors(ct_levels)

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
  labs(title    = "Purity metacelle — CTX",
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
  scale_fill_manual(values = pal_ct_ctx, guide = "none") +
  theme_bw(base_size = 9) +
  labs(title    = "Purity per cell type — CTX",
       subtitle = "Linea tratt. = soglia 0.7",
       x = "", y = "Purity")

savemc(p_hist_pur | p_vln_pur, "01_purity_CTX.png",
       width  = 5500,
       height = max(2200, length(ct_levels) * 100 + 700))

# ------------------------------------------------------------------ #
#  2. DIMENSIONE METACELLE                                             #
# ------------------------------------------------------------------ #
message("  2. Dimensione metacelle...")

size_df <- data.frame(mc_size   = mc_size,
                      cell_type = meta$cell_type,
                      treatment = meta$treatment)

p_size_hist <- ggplot(size_df, aes(x = mc_size)) +
  geom_histogram(bins = 30, fill = "#8E44AD", color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean(mc_size), color = "#E74C3C",
             linetype = "dashed") +
  theme_bw(base_size = 10) +
  labs(title = "Dimensione metacelle — CTX",
       x = "N cellule per metacella", y = "N metacelle")

p_size_ct <- ggplot(size_df,
                    aes(x = reorder(cell_type, mc_size, median),
                        y = mc_size, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  scale_fill_manual(values = pal_ct_ctx, guide = "none") +
  theme_bw(base_size = 9) +
  labs(title = "Dimensione per cell type — CTX",
       x = "", y = "N cellule per metacella")

savemc(p_size_hist | p_size_ct, "02_size_CTX.png",
       width  = 5500,
       height = max(2200, length(ct_levels) * 100 + 700))

# ------------------------------------------------------------------ #
#  3. UMAP METACELLE (spazio metacelle)                               #
# ------------------------------------------------------------------ #
message("  3. UMAP metacelle...")

SC_CTX_seurat <- FindVariableFeatures(SC_CTX_seurat, nfeatures = 2000, verbose = FALSE)
SC_CTX_seurat <- ScaleData(SC_CTX_seurat, verbose = FALSE)
SC_CTX_seurat <- RunPCA(SC_CTX_seurat,
                        npcs    = min(30, ncol(SC_CTX_seurat) - 1),
                        verbose = FALSE)
SC_CTX_seurat <- RunUMAP(SC_CTX_seurat,
                         dims     = 1:min(20, ncol(SC_CTX_seurat) - 1),
                         verbose  = FALSE,
                         min.dist = 0.3)

p_umap_ct <- DimPlot(SC_CTX_seurat, group.by = "cell_type",
                     pt.size    = pmax(0.5, 3 / sqrt(ncol(SC_CTX_seurat))),
                     label      = TRUE, label.size = 3, repel = TRUE,
                     cols       = pal_ct_ctx) +
  ggtitle("UMAP metacelle — cell type | CTX") +
  theme(legend.text = element_text(size = 8))

p_umap_trt <- DimPlot(SC_CTX_seurat, group.by = "treatment",
                      pt.size = pmax(0.5, 3 / sqrt(ncol(SC_CTX_seurat))),
                      cols    = pal_trt) +
  ggtitle("UMAP metacelle — trattamento | CTX")

umap_df <- as.data.frame(Embeddings(SC_CTX_seurat, "umap"))
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df <- cbind(umap_df,
                 meta[, c("cell_type","treatment","purity","mc_size")])

p_umap_pur <- ggplot(umap_df,
                     aes(x = UMAP1, y = UMAP2,
                         color = purity, size = mc_size)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Purity") +
  scale_size_continuous(range = c(0.5, 4), name = "N cellule") +
  theme_bw(base_size = 9) +
  labs(title = "UMAP metacelle — purity e dimensione | CTX")

savemc(
  (p_umap_ct | p_umap_trt) / p_umap_pur +
    plot_annotation(title = "Metacelle UMAP — CTX"),
  "03_UMAP_panel_CTX.png",
  width = 5500, height = 4500
)

# ------------------------------------------------------------------ #
#  4. OVERLAY METACELLE SU UMAP ORIGINALE                             #
# ------------------------------------------------------------------ #
message("  4. Overlay UMAP...")

if ("umap" %in% names(CTX@reductions)) {
  
  sc_umap <- as.data.frame(Embeddings(CTX, "umap"))
  colnames(sc_umap) <- c("UMAP1","UMAP2")
  sc_umap$cell_type <- CTX$cell_type
  
  mc_centroids <- do.call(rbind, lapply(
    seq_len(ncol(SC_CTX_seurat)), function(i) {
      cells_i <- intersect(names(SC$membership)[SC$membership == i],
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
               aes(x = UMAP1, y = UMAP2, color = cell_type),
               size = 0.15, alpha = 0.2) +
    geom_point(data = mc_df,
               aes(x = UMAP1, y = UMAP2,
                   fill  = cell_type,
                   size  = mc_size,
                   alpha = purity),
               shape = 21, color = "white", stroke = 0.3) +
    scale_color_manual(values = pal_ct_ctx, guide = "none") +
    scale_fill_manual(values  = pal_ct_ctx, name = "Cell type") +
    scale_size_continuous(range = c(1, 6), name = "N cellule/MC") +
    scale_alpha_continuous(range = c(0.4, 1), name = "Purity") +
    theme_bw(base_size = 9) +
    guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title    = "Metacelle sovrapposte UMAP originale — CTX",
         subtitle = "Punti grandi = metacelle (dim \u221d N cell, opacity \u221d purity)",
         x = "UMAP 1", y = "UMAP 2")
  
  savemc(p_overlay, "04_UMAP_overlay_CTX.png",
         width = 3500, height = 3000)
}

# ------------------------------------------------------------------ #
#  5. COMPOSIZIONE CT \u00d7 CAMPIONE                                 #
# ------------------------------------------------------------------ #
message("  5. Composizione...")

comp_df <- meta %>%
  count(cell_type, orig.ident, treatment) %>%
  group_by(orig.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_comp_bar <- ggplot(comp_df,
                     aes(x = orig.ident, y = frac, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ treatment, scales = "free_x") +
  scale_fill_manual(values = pal_ct_ctx) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7)) +
  labs(title = "Composizione metacelle per campione — CTX",
       x = "", y = "Proporzione", fill = "Cell type")

savemc(p_comp_bar, "05_composition_bar_CTX.png",
       width = 3500, height = 2500)

ct_sample_tab <- table(meta$cell_type, meta$orig.ident)
if (nrow(ct_sample_tab) > 1) {
  png(file.path(mc_dir, "05_composition_heatmap_CTX.png"),
      width  = 3000,
      height = max(2000, nrow(ct_sample_tab) * 80 + 600),
      res    = 300)
  pheatmap(as.matrix(ct_sample_tab),
           display_numbers = TRUE, number_format = "%d",
           color           = colorRampPalette(c("white","#2E86C1"))(50),
           cluster_cols    = FALSE, fontsize = 9,
           main            = "N metacelle CT \u00d7 campione — CTX")
  dev.off()
}

# ------------------------------------------------------------------ #
#  6. HEATMAP ESPRESSIONE GENI CHIAVE (z-score per CT)                #
# ------------------------------------------------------------------ #
message("  6. Espressione geni chiave...")

genes_key <- filter_genes(
  c("Htr2a","Gnaq","Plcb1","Prkca","Camk2a",
    "Fos","Arc","Egr1","Nr4a1","Fosb","Junb",
    "Slc17a7","Slc17a6","Gad1","Gad2",
    "Rorb","Fezf2","Foxp2","Pvalb","Sst",
    "Aqp4","Mog","Cx3cr1","Pdgfra"),
  SC_CTX_seurat
)

if (length(genes_key) >= 5) {
  avg_expr <- AverageExpression(
    SC_CTX_seurat,
    features = genes_key,
    group.by = "cell_type",
    assays   = "RNA",
    layer    = "data",
    verbose  = FALSE
  )$RNA
  
  mat_z <- t(scale(t(as.matrix(avg_expr))))
  mat_z[is.nan(mat_z)] <- 0
  lim     <- max(abs(mat_z), na.rm = TRUE)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("#2471A3","white","#CB4335"))
  
  ht <- Heatmap(
    mat_z,
    name             = "z-score",
    col              = col_fun,
    cluster_rows     = TRUE, cluster_columns = TRUE,
    row_names_gp     = gpar(fontsize = 8),
    column_names_gp  = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title     = "Espressione media metacelle (z-score) — CTX",
    heatmap_legend_param = list(title = "z-score")
  )
  png(file.path(mc_dir, "06_expr_heatmap_CTX.png"),
      width  = 3500,
      height = max(2000, length(genes_key) * 55 + 600),
      res    = 300)
  draw(ht)
  dev.off()
}

# ------------------------------------------------------------------ #
#  7. VlnPlot Htr2a split per trattamento                             #
# ------------------------------------------------------------------ #
if ("Htr2a" %in% rownames(SC_CTX_seurat)) {
  message("  7. VlnPlot Htr2a...")
  p_htr2a <- VlnPlot(SC_CTX_seurat, features = "Htr2a",
                     group.by = "cell_type", split.by = "treatment",
                     pt.size = 0.3, cols = pal_trt) +
    ggtitle("Htr2a nelle metacelle — CTX") +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  savemc(p_htr2a, "07_Htr2a_violin_CTX.png",
         width  = max(3000, length(ct_levels) * 280 + 800),
         height = 2500)
}

# ------------------------------------------------------------------ #
#  8. BOXPLOT Ctrl vs LSD — geni chiave per CT                        #
# ------------------------------------------------------------------ #
message("  8. Ctrl vs LSD nelle metacelle...")

genes_lsd <- filter_genes(
  c("Htr2a","Fos","Arc","Egr1","Nr4a1","Fosb","Junb",
    "Bdnf","Ntrk2","Camk2a","Gnaq","Plcb1"),
  SC_CTX_seurat
)

if (length(genes_lsd) > 0 && "treatment" %in% colnames(meta)) {
  expr_long <- as.data.frame(
    t(as.matrix(GetAssayData(SC_CTX_seurat, assay = "RNA",
                             layer = "data")[genes_lsd, , drop = FALSE]))
  ) %>%
    mutate(treatment = meta$treatment,
           cell_type = meta$cell_type,
           mc_size   = meta$mc_size) %>%
    pivot_longer(cols      = all_of(genes_lsd),
                 names_to  = "gene",
                 values_to = "expr")
  
  p_box_lsd <- ggplot(expr_long,
                      aes(x = treatment, y = expr, fill = treatment)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.4, alpha = 0.8) +
    geom_jitter(aes(size = mc_size), width = 0.15, alpha = 0.3,
                color = "grey40") +
    facet_grid(gene ~ cell_type, scales = "free_y") +
    scale_fill_manual(values = pal_trt) +
    scale_size_continuous(range = c(0.3, 2), guide = "none") +
    theme_bw(base_size = 7) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1),
          strip.text      = element_text(size = 6),
          legend.position = "none") +
    labs(title = "Espressione geni chiave — Ctrl vs LSD | CTX",
         x = "", y = "Espressione (log-norm)")
  
  savemc(p_box_lsd, "08_expr_ctrl_lsd_CTX.png",
         width  = max(3000, length(ct_levels) * 300 + 800),
         height = max(3000, length(genes_lsd) * 250 + 800))
}

# ------------------------------------------------------------------ #
#  9. CORRELAZIONE SINGOLA CELLULA vs METACELLA                       #
# ------------------------------------------------------------------ #
message("  9. Correlazione SC vs metacella...")

genes_cor <- filter_genes(
  c("Htr2a","Fos","Egr1","Arc","Gad1","Slc17a7","Aqp4","Mog",
    "Pvalb","Sst","Rorb","Fezf2","Foxp2"),
  SC_CTX_seurat
)
genes_cor <- intersect(genes_cor, rownames(CTX))

if (length(genes_cor) >= 5) {
  
  mc_expr <- as.matrix(
    GetAssayData(SC_CTX_seurat, assay = "RNA",
                 layer = "data")[genes_cor, , drop = FALSE]
  )
  sc_expr <- as.matrix(
    LayerData(CTX, assay = "SCT", layer = "data")[genes_cor, , drop = FALSE]
  )
  
  mc_ids <- sort(unique(SC$membership))
  sc_avg <- do.call(cbind, lapply(mc_ids, function(i) {
    cells_i <- intersect(names(SC$membership)[SC$membership == i],
                         colnames(sc_expr))
    if (length(cells_i) == 0) return(rep(0, nrow(sc_expr)))
    rowMeans(sc_expr[, cells_i, drop = FALSE])
  }))
  colnames(sc_avg) <- paste0("MC_", mc_ids)
  
  shared_mc <- intersect(colnames(mc_expr), colnames(sc_avg))
  if (length(shared_mc) >= 10) {
    cor_vals <- sapply(genes_cor, function(g) {
      tryCatch(cor(sc_avg[g, shared_mc], mc_expr[g, shared_mc],
                   method = "pearson"),
               error = function(e) NA_real_)
    })
    
    cor_df <- data.frame(gene = genes_cor, pearson_r = cor_vals) %>%
      arrange(desc(pearson_r))
    write.csv(cor_df,
              file.path(mc_dir, "09_SC_MC_correlation_CTX.csv"),
              row.names = FALSE)
    
    p_cor <- ggplot(cor_df,
                    aes(x = reorder(gene, pearson_r), y = pearson_r,
                        fill = pearson_r)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
      coord_flip() +
      scale_fill_gradient2(low = "#2471A3", mid = "grey90", high = "#CB4335",
                           midpoint = 0.5, name = "r") +
      theme_bw(base_size = 9) +
      labs(title    = "Correlazione SC vs metacella — CTX",
           subtitle = "Linea tratt. = r=0.9",
           x = "", y = "Pearson r")
    
    savemc(p_cor, "09_SC_MC_correlation_CTX.png",
           width  = 2500,
           height = max(1800, length(genes_cor) * 70 + 500))
  }
}

# ------------------------------------------------------------------ #
# 10. QC METACELLE                                                    #
# ------------------------------------------------------------------ #
message("  10. QC metacelle...")

SC_CTX_seurat[["percent.mt"]] <- PercentageFeatureSet(SC_CTX_seurat,
                                                      pattern = "^mt-")
p_qc <- VlnPlot(
  SC_CTX_seurat,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  group.by = "cell_type",
  pt.size  = 0.2, ncol = 3,
  cols     = pal_ct_ctx
) + plot_annotation(
  title    = "QC metacelle — CTX",
  subtitle = sprintf("%d metacelle totali", ncol(SC_CTX_seurat))
)
savemc(p_qc, "10_QC_CTX.png",
       width  = max(3000, length(ct_levels) * 280 + 800),
       height = 2800)

message(sprintf("  Completato CTX: %d metacelle", ncol(SC_CTX_seurat)))


# 06.3 — SALVATAGGIO
# ==============================================================================

message("=== 06.3 Salvataggio ===")

save(SC,  SC_CTX_seurat,
     file = file.path(outdir, "Metacells_CTX.RData"))

message("Salvati: Metacells_CTX.RData, Metacells_THAL.RData")
message(paste("  CTX:  %d metacelle",
                ncol(SC_CTX_seurat)))
