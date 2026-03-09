################################################################################
##  ST_02_norm_cluster.R
##
##  Input:  ST_01_experiment_merged.RData
##          ST_01_slices_filtered.RData
##  Output: ST_02_experiment_merged_clustered.RData
##          df_toplot.RData  — data.frame di plotting con coordinate spaziali
##
##  Flusso:
##    1. SCTransform + PCA
##    2. Integrazione Harmony per batch correction
##    3. Clustering Seurat (su riduzione Harmony)
##    4. UMAP
##    5. Costruzione df di plotting con coordinate spaziali e mirror
##    6. Plot diagnostici
##    7. Salvataggio
################################################################################

source("code/ST_pipeline/ST_00_config.R")

load(file.path(outdir, "ST_01_experiment_merged.RData"))
load(file.path(outdir, "ST_01_slices_filtered.RData"))

# ==============================================================================
# 1. SCTransform + PCA
# ==============================================================================
options(future.globals.maxSize = 8000 * 1024^2)

message("=== ST_02A. SCTransform + PCA ===")

experiment.merged <- SCTransform(
  experiment.merged,
  assay   = "Spatial",
  verbose = FALSE
) %>%
  RunPCA(verbose = FALSE)

# QC PCA: varianza spiegata
pdf(file.path(outdir, "ST_02_PCA_elbow.pdf"), 5, 4)
ElbowPlot(experiment.merged, ndims = 50)
dev.off()

# Check batch effect prima di Harmony
p_pca_pre <- DimPlot(experiment.merged,
                     reduction = "pca",
                     group.by  = "capture.area",
                     pt.size   = 0.3) +
  ggtitle("PCA pre-Harmony") +
  theme_classic(base_size = 10)
savepng(p_pca_pre, "ST_02_PCA_pre_harmony.png", width = 3000, height = 2500)


# ==============================================================================
# 2. Harmony batch correction
# ==============================================================================

message("=== ST_02B. Harmony integration ===")

experiment.merged.h <- experiment.merged %>%
  RunHarmony("capture.area", plot_convergence = FALSE, verbose = FALSE)

p_pca_post <- DimPlot(experiment.merged.h,
                      reduction = "harmony",
                      group.by  = "capture.area",
                      pt.size   = 0.3) +
  ggtitle("PCA post-Harmony") +
  theme_classic(base_size = 10)
savepng(p_pca_post, "ST_02_PCA_post_harmony.png", width = 3000, height = 2500)


# ==============================================================================
# 3. Clustering
# ==============================================================================

message("=== ST_02C. Clustering ===")

experiment.merged.h <- FindNeighbors(
  experiment.merged.h, reduction = "harmony", dims = HARMONY_DIMS
) %>%
  FindClusters(verbose = FALSE)

# Propaga cluster all'oggetto non-harmonized per compatibilità downstream
experiment.merged$cluster <- Idents(experiment.merged.h)
experiment.merged$cluster.condition <- paste(
  Idents(experiment.merged.h),
  experiment.merged.h$capture.area,
  sep = "_"
)
Idents(experiment.merged) <- "cluster"


# ==============================================================================
# 4. UMAP
# ==============================================================================

message("=== ST_02D. UMAP ===")

experiment.merged.h <- RunUMAP(
  experiment.merged.h,
  reduction = "harmony",
  dims      = UMAP_DIMS,
  verbose   = FALSE
)

p_umap <- DimPlot(experiment.merged.h,
                  group.by = c("ident", "capture.area"),
                  pt.size  = 0.3) &
  theme_classic(base_size = 9)
savepng(p_umap, "ST_02_UMAP_clusters.png", width = 6000, height = 2500)

# SpatialDimPlot (solo primo campione per brevità)
png(file.path(outdir, "ST_02_SpatialDimPlot.png"),
    width = 6000, height = 2000, res = 300)
SpatialDimPlot(experiment.merged.h)
dev.off()


# ==============================================================================
# 5. Costruzione df di plotting
# ==============================================================================
# df contiene: coordinate spaziali (imagecol, imagerow) già orientate/mirror,
# metadati spot (area, cluster, area.spot), e successivamente score di transfer.
# Le operazioni di mirror sono specifiche per questo esperimento — verificale
# visivamente con il plot qui sotto e modifica se necessario.
# ==============================================================================

message("=== ST_02E. Costruzione df di plotting ===")

# Estrai coordinate direttamente dal metadata di experiment.merged
# (Seurat le conserva in meta.data dopo il merge come imagerow/imagecol)
coords_all <- do.call(rbind, lapply(SAMPLE_NAMES, function(sname) {
  obj <- get(sname)
  coords <- GetTissueCoordinates(obj)
  # Forza numeric esplicitamente
  coords[] <- lapply(coords, as.numeric)
  coords
}))

# Verifica che i barcode coincidano
shared_cells <- intersect(Cells(experiment.merged), rownames(coords_all))
cat(sprintf(
  "Spot con coordinate: %d / %d\n",
  length(shared_cells), ncol(experiment.merged)
))

# ==============================================================================
# Costruzione df di plotting
# ==============================================================================

df <- data.frame(
  area      = experiment.merged$capture.area[match(shared_cells,
                                                   Cells(experiment.merged))],
  x         = as.numeric(coords_all[shared_cells, "x"]),
  y         = as.numeric(coords_all[shared_cells, "y"]),
  row.names = shared_cells,
  stringsAsFactors = FALSE
)
df$area <- factor(df$area, levels = SAMPLE_ORDER)

# Centra ogni sezione intorno allo zero
for (sname in SAMPLE_NAMES) {
  mask <- df$area == sname
  df$x[mask] <- df$x[mask] - mean(range(df$x[mask]))
  df$y[mask] <- df$y[mask] - mean(range(df$y[mask]))
}

# Mirror per allineare le sezioni
df$x[df$area == "CTRL90min_1"] <- -df$x[df$area == "CTRL90min_1"]
df$x[df$area == "LSD90min_2"]  <- -df$x[df$area == "LSD90min_2"]
df$y[df$area == "LSD90min_2"]  <- -df$y[df$area == "LSD90min_2"]
df$x[df$area == "CTRL90min_2"] <- -df$x[df$area == "CTRL90min_2"]
df$y[df$area == "CTRL90min_2"] <- -df$y[df$area == "CTRL90min_2"]
df$x[df$area == "CTRL24h_1"]   <- -df$x[df$area == "CTRL24h_1"]
df$x[df$area == "LSD24h_1"]    <- -df$x[df$area == "LSD24h_1"]
df$x[df$area == "CTRL24h_2"]   <- -df$x[df$area == "CTRL24h_2"]
df$y[df$area == "CTRL24h_2"]   <- -df$y[df$area == "CTRL24h_2"]


# Aggiungi cluster dopo FindClusters
df$cluster <- experiment.merged$cluster[match(rownames(df),
                                              Cells(experiment.merged))]

# ==============================================================================
# Plot verifica orientamento
# ==============================================================================

p_orient <- ggplot(df, aes(x = x, y = -y, colour = factor(cluster))) +
  geom_point(size = 0.6) +
  facet_wrap(~ area, ncol = 2) +
  theme_classic(base_size = 9) +
  labs(title = "Cluster ST nello spazio", colour = "cluster")
savepng(p_orient, "ST_02_clusters_space.png", width = 4000, height = 6000)




# Ripeti il plot con i cluster
p_orient_clustered <- ggplot(df, aes(x = x, y = -y, colour = factor(cluster))) +
  geom_point(size = 0.6) +
  facet_wrap(~ area, ncol = 2) +
  theme_classic(base_size = 9) +
  labs(title = "Cluster ST nello spazio", colour = "cluster")
savepng(p_orient_clustered, "ST_02_clusters_space_clustered.png", width = 4000, height = 5000)


# ==============================================================================
# 6. Plot diagnostici aggiuntivi
# ==============================================================================

# Numero cluster per sezione
cat("\nCluster per cattura:\n")
print(table(df$area, df$cluster))

# Heatmap cluster x sezione (proporzioni)
cluster_section_mat <- table(df$cluster, df$area)
cluster_section_prop <- prop.table(cluster_section_mat, margin = 2)

png(file.path(outdir, "ST_02_heatmap_cluster_section.png"),
    width = 3000, height = 2500, res = 300)
pheatmap(
  cluster_section_prop,
  color    = colorRampPalette(c("white", "#457B9D"))(50),
  fontsize = 8,
  main     = "Proporzione cluster per sezione"
)
dev.off()


# ==============================================================================
# 7. Salvataggio
# ==============================================================================

save(experiment.merged,
     experiment.merged.h,
     file = file.path(outdir, "ST_02_experiment_merged_clustered.RData"))

save(df, file = file.path(outdir, "ST_02_df_toplot.RData"))

message("=== ST_02 completato ===")
message("  Cluster trovati: ", nlevels(experiment.merged$cluster))
message("  Salvati: ST_02_experiment_merged_clustered.RData, ST_02_df_toplot.RData")
