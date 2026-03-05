################################################################################
##                                                                            ##
##   PIPELINE scRNA-seq — Seurat v5                                          ##
##   LSD vs Ctrl | 90min & 24h | Cortex + Thalamus                          ##
##                                                                            ##
##   Sezioni:                                                                 ##
##     01. Caricamento dati, metadati, QC metrics e plot pre-filtraggio      ##
##     02. Doublet detection (scDblFinder) + filtraggio QC                   ##
##     03. SCTransform, IntegrateLayers (Harmony), UMAP, clustering          ##
##         → oggetti separati CTX e THAL                                     ##
##                                                                            ##
##   Logica Seurat v5 rispetto a v4:                                         ##
##     - merge() crea layer separati per campione (counts.C1, counts.C2…)   ##
##     - SCTransform() processa ogni layer indipendentemente                 ##
##     - IntegrateLayers(HarmonyIntegration) sostituisce RunHarmony()        ##
##     - JoinLayers() compatta i layer prima di FindMarkers / GetAssayData   ##
##     - SplitObject() non serve più per normalizzare: si usa split = col    ##
##                                                                            ##
################################################################################

# ==============================================================================
# LIBRERIE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)            # >= 5.0
  library(SeuratObject)      # >= 5.0
  library(harmony)           # backend per IntegrateLayers
  library(scDblFinder)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

# ==============================================================================
# IMPOSTAZIONI GLOBALI
# ==============================================================================

set.seed(42)
options(Seurat.object.assay.version = "v5")  # assay v5 esplicito

date_str <- format(Sys.Date(), "%Y%m%d")
outdir   <- file.path("results", date_str)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

savepng <- function(p, filename, width = 2500, height = 2000, res = 300) {
  png(file.path(outdir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# Palette condivise
pal_trt  <- c("Ctrl" = "#3498DB", "LSD"       = "#E74C3C")
pal_time <- c("90min" = "#F39C12", "24h"      = "#8E44AD")
pal_area <- c("Cortex" = "#27AE60", "Thalamus" = "#2E86C1")

# ==============================================================================
# SEZIONE 01 — CARICAMENTO DATI E QC
# ==============================================================================

message("\n=== 01. Caricamento dati e QC ===\n")

# Campioni:  C = Cortex | T = Thalamus
# 1, 2 = 90min;  3, 4 = 24h
# 1, 3 = Ctrl;   2, 4 = LSD
dirs <- c(
  "../../3_Data/PoolC1/", "../../3_Data/PoolC2/",
  "../../3_Data/PoolC3/", "../../3_Data/PoolC4/",
  "../../3_Data/PoolT1/", "../../3_Data/PoolT2/",
  "../../3_Data/PoolT3/", "../../3_Data/PoolT4/"
)
proj <- c("C1","C2","C3","C4","T1","T2","T3","T4")

# --- Caricamento: un Seurat per campione ------------------------------------
raw_list    <- lapply(dirs, Read10X)
seurat_list <- mapply(
  function(r, p) CreateSeuratObject(counts = r, project = p,
                                    min.cells = 3, min.features = 200),
  raw_list, proj, SIMPLIFY = FALSE
)
names(seurat_list) <- proj
rm(raw_list); gc()

# --- Metadati biologici (assegnati per campione, prima del merge) -----------
treatment_map <- c(C1="Ctrl",C2="LSD", C3="Ctrl",C4="LSD",
                   T1="Ctrl",T2="LSD", T3="Ctrl",T4="LSD")
time_map      <- c(C1="90min",C2="90min",C3="24h",C4="24h",
                   T1="90min",T2="90min",T3="24h",T4="24h")
area_map      <- c(C1="Cortex",C2="Cortex",C3="Cortex",C4="Cortex",
                   T1="Thalamus",T2="Thalamus",T3="Thalamus",T4="Thalamus")

seurat_list <- lapply(names(seurat_list), function(nm) {
  obj           <- seurat_list[[nm]]
  obj$treatment <- treatment_map[[nm]]
  obj$time      <- time_map[[nm]]
  obj$area      <- area_map[[nm]]
  message(sprintf("  %s (%s | %s | %s): %d cellule, %d geni",
                  nm, area_map[nm], treatment_map[nm], time_map[nm],
                  ncol(obj), nrow(obj)))
  obj
})
names(seurat_list) <- proj

# --- QC metrics (per campione, prima del merge) -----------------------------
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]]       <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.rb"]]       <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["log10GenesPerUMI"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  obj
})

# --- Merge temporaneo SOLO per i plot QC ------------------------------------
brain_preqc <- merge(seurat_list[[1]], y = seurat_list[-1],
                     add.cell.ids = proj, project = "LSD_brain")
# JoinLayers necessario in v5 per accedere ai dati in modo unificato
brain_preqc <- JoinLayers(brain_preqc, assay = "RNA")

p_vln_pre <- VlnPlot(
  brain_preqc,
  features = c("nFeature_RNA","nCount_RNA","percent.mt","log10GenesPerUMI"),
  group.by = "orig.ident", pt.size = 0, ncol = 4
) + plot_annotation(title = "QC pre-filtraggio — tutti i campioni")
savepng(p_vln_pre, "01_QC_violin_prefilter.png", width = 7000, height = 2000)

p_scatter_ctx <- FeatureScatter(
  brain_preqc[, brain_preqc$area == "Cortex"],
  feature1 = "nFeature_RNA", feature2 = "percent.mt",
  group.by = "orig.ident"
) +
  geom_hline(yintercept = c(5, 10, 15),
             linetype = "dashed", color = c("green","orange","red")) +
  ggtitle("Cortex — nFeature vs %mt")

p_scatter_thal <- FeatureScatter(
  brain_preqc[, brain_preqc$area == "Thalamus"],
  feature1 = "nFeature_RNA", feature2 = "percent.mt",
  group.by = "orig.ident"
) +
  geom_hline(yintercept = c(5, 10, 15),
             linetype = "dashed", color = c("green","orange","red")) +
  ggtitle("Thalamus — nFeature vs %mt")

p_scatter_count <- FeatureScatter(
  brain_preqc,
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
  group.by = "area"
) +
  scale_color_manual(values = pal_area) +
  geom_hline(yintercept = c(200, 7500),
             linetype = "dashed", color = "red") +
  ggtitle("Count vs Feature — tutte le aree")

savepng(
  (p_scatter_ctx | p_scatter_thal) / p_scatter_count,
  "01_QC_scatter.png", width = 5000, height = 4000
)

rm(brain_preqc); gc()

# ==============================================================================
# SEZIONE 02 — DOUBLET DETECTION (scDblFinder) + FILTRAGGIO QC
# ==============================================================================

message("\n=== 02. Doublet detection ===\n")

# scDblFinder va eseguito PER CAMPIONE prima del merge:
# il tasso atteso di doublet dipende dal numero di cellule caricate per corsia
seurat_list <- lapply(names(seurat_list), function(nm) {
  obj <- seurat_list[[nm]]
  sce <- scDblFinder(as.SingleCellExperiment(obj), verbose = FALSE)
  obj$doublet_class <- sce$scDblFinder.class   # "singlet" / "doublet"
  obj$doublet_score <- sce$scDblFinder.score
  n_dbl <- sum(obj$doublet_class == "doublet")
  message(sprintf("  %s: %d doublet / %d cellule (%.1f%%)",
                  nm, n_dbl, ncol(obj), 100 * n_dbl / ncol(obj)))
  obj
})
names(seurat_list) <- proj

# --- Tabella riassuntiva ----------------------------------------------------
doublet_summary <- do.call(rbind, lapply(names(seurat_list), function(nm) {
  obj <- seurat_list[[nm]]
  data.frame(
    sample      = nm,
    area        = obj$area[1],
    treatment   = obj$treatment[1],
    time        = obj$time[1],
    n_total     = ncol(obj),
    n_doublet   = sum(obj$doublet_class == "doublet"),
    n_singlet   = sum(obj$doublet_class == "singlet"),
    pct_doublet = 100 * mean(obj$doublet_class == "doublet"),
    stringsAsFactors = FALSE
  )
}))
print(doublet_summary %>% mutate(pct_doublet = round(pct_doublet, 2)))

# --- Plot doublet -----------------------------------------------------------
brain_dbl <- merge(seurat_list[[1]], y = seurat_list[-1],
                   add.cell.ids = proj, project = "LSD_brain")
brain_dbl <- JoinLayers(brain_dbl, assay = "RNA")

doublet_df <- brain_dbl@meta.data %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, orig.ident, area, treatment, time,
         doublet_class, doublet_score)

p_pct <- ggplot(doublet_summary,
                aes(x = sample, y = pct_doublet, fill = area)) +
  geom_bar(stat = "identity", color = "white", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_doublet, n_doublet)),
            vjust = -0.3, size = 3) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed",
             color = c("orange","red"), linewidth = 0.5) +
  scale_fill_manual(values = pal_area)  +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(title    = "% Doublet per campione (scDblFinder)",
       subtitle = "Linee: 5% arancio | 10% rosso",
       x = "", y = "% doublet", fill = "Area")

p_abs <- doublet_df %>%
  mutate(class = ifelse(doublet_class == "singlet","Singlet","Doublet")) %>%
  ggplot(aes(x = orig.ident, fill = class)) +
  geom_bar(color = "white", width = 0.7) +
  geom_text(data = doublet_summary,
            aes(x = sample, y = n_total,
                label = formatC(n_total, format = "d", big.mark = ",")),
            inherit.aes = FALSE, vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Singlet" = "#82E0AA","Doublet" = "#E74C3C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_grid(~ area + time, scales = "free_x", space = "free") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(title = "Singlet vs Doublet — N assoluti",
       x = "", y = "N cellule", fill = "")

p_score <- ggplot(doublet_df,
                  aes(x = orig.ident, y = doublet_score, fill = doublet_class)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c("singlet" = "#82E0AA","doublet" = "#E74C3C")) +
  facet_grid(~ area + time, scales = "free_x", space = "free") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribuzione doublet score",
       x = "", y = "Score", fill = "Classe")

savepng(
  (p_pct / p_abs / p_score) +
    plot_annotation(
      title    = "QC Doublet Detection — scDblFinder",
      subtitle = sprintf(
        "Totale: %s cellule | Doublet: %s (%.1f%%)",
        formatC(sum(doublet_summary$n_total),   format = "d", big.mark = ","),
        formatC(sum(doublet_summary$n_doublet), format = "d", big.mark = ","),
        100 * sum(doublet_summary$n_doublet) / sum(doublet_summary$n_total)
      )
    ),
  "02_QC_doublets.png", width = 4000, height = 6000
)

rm(brain_dbl, doublet_df); gc()

# --- Filtraggio QC ----------------------------------------------------------
# <<< SOGLIE: aggiornare dopo aver visto i plot QC >>>
QC <- list(
  nFeature_min = 200,
  nFeature_max = 7500,
  percent_mt   = 15,
  log10ratio   = 0.80
)

seurat_list <- lapply(seurat_list, function(obj) {
  n_before <- ncol(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA     >  QC$nFeature_min &
      nFeature_RNA     <  QC$nFeature_max &
      percent.mt       <  QC$percent_mt   &
      log10GenesPerUMI >  QC$log10ratio   &
      doublet_class    == "singlet"
  )
  message(sprintf("  %s: %d → %d cellule (rimosso %d, %.1f%%)",
                  obj$orig.ident[1], n_before, ncol(obj),
                  n_before - ncol(obj),
                  100 * (n_before - ncol(obj)) / n_before))
  obj
})
names(seurat_list) <- proj

# ==============================================================================
# SEZIONE 03 — NORMALIZZAZIONE, INTEGRAZIONE E CLUSTERING (Seurat v5)
#
#  Flusso per ogni area (CTX, THAL):
#    1. merge() dei campioni dell'area    → layer separati per campione
#    2. SCTransform()                     → normalizza ogni layer
#    3. RunPCA()                          → PCA su geni variabili SCT
#    4. IntegrateLayers(Harmony)          → corregge batch tra campioni
#    5. JoinLayers()                      → compatta i layer
#    6. RunUMAP / FindNeighbors / FindClusters (su riduzione harmony)
# ==============================================================================

message("\n=== 03. Normalizzazione, integrazione e clustering ===\n")

# ---- Parametri -------------------------------------------------------------
N_FEATURES <- 3000   # geni variabili per SCTransform
N_PCS      <- 50     # PC da calcolare
N_PCS_CTX  <- 40     # <<< aggiornare dopo ElbowPlot CTX
N_PCS_THAL <- 40     # <<< aggiornare dopo ElbowPlot THAL
RES_CTX    <- 0.4    # resoluzione Leiden per CTX
RES_THAL   <- 0.4    # resoluzione Leiden per THAL

# ---- Funzione pipeline v5 --------------------------------------------------
run_pipeline_v5 <- function(obj_list, label, dims_use, resolution) {
  
  message(sprintf("  [%s] Avvio con %d campioni: %s",
                  label, length(obj_list),
                  paste(names(obj_list), collapse = ", ")))
  
  # 1. Merge — in Seurat v5 mantiene layer separati per campione
  #    Non serve SplitObject + SCTransform loop come in v4
  obj <- merge(
    obj_list[[1]],
    y            = obj_list[-1],
    add.cell.ids = names(obj_list),
    project      = paste0(label, "_project")
  )
  message(sprintf("  [%s] Dopo merge: %d cellule, %d geni",
                  label, ncol(obj), nrow(obj)))
  # Layer presenti: obj[["RNA"]]@layers — uno per campione
  
  # 2. SCTransform — normalizza ogni layer separatamente in automatico
  #    Equivalente al vecchio loop SCTransform per campione in v4
  #    vars.to.regress rimuove la variazione da %mt
  obj <- SCTransform(
    obj,
    vars.to.regress     = "percent.mt",
    variable.features.n = N_FEATURES,
    verbose             = FALSE
  )
  # DefaultAssay(obj) == "SCT"
  
  # 3. PCA
  obj <- RunPCA(obj, npcs = N_PCS, verbose = FALSE)
  
  savepng(
    ElbowPlot(obj, ndims = N_PCS) + ggtitle(paste("Elbow —", label)),
    paste0("03_elbow_", label, ".png")
  )
  
  # 4. Integrazione batch con Harmony via IntegrateLayers
  #    - NON è necessario RunHarmony() separato come in v4
  #    - group.by.vars = "orig.ident": corregge il batch per campione
  #    - orig.reduction = "pca": parte dallo spazio PCA
  #    - new.reduction = "harmony": salva la riduzione integrata
  obj <- IntegrateLayers(
    object         = obj,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    group.by.vars  = "orig.ident",
    verbose        = FALSE
  )
  message(sprintf("  [%s] IntegrateLayers (Harmony) completato", label))
  
  # 5. JoinLayers — compatta counts.C1 + counts.C2 + … → "counts"
  #    Obbligatorio prima di FindMarkers, GetAssayData, DESeq2
  obj <- JoinLayers(obj, assay = "RNA")
  message(sprintf("  [%s] JoinLayers completato", label))
  
  # 6. UMAP, vicini e clustering nello spazio Harmony
  obj <- RunUMAP(obj,     reduction = "harmony", dims = dims_use, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, algorithm = 4, verbose = FALSE)
  # algorithm = 4: Leiden (più stabile e riproducibile di Louvain)
  
  n_cl <- length(levels(obj$seurat_clusters))
  message(sprintf("  [%s] %d cluster (res = %.2f) | %d cellule",
                  label, n_cl, resolution, ncol(obj)))
  obj
}

# ---- 3a. Cortex -------------------------------------------------------------
message("\n--- 03a. Cortex ---")

CTX <- run_pipeline_v5(
  obj_list   = seurat_list[c("C1","C2","C3","C4")],
  label      = "CTX",
  dims_use   = 1:N_PCS_CTX,
  resolution = RES_CTX
)
CTX$initial_clusters <- CTX$seurat_clusters

# ---- 3b. Thalamus -----------------------------------------------------------
message("\n--- 03b. Thalamus ---")

THAL <- run_pipeline_v5(
  obj_list   = seurat_list[c("T1","T2","T3","T4")],
  label      = "THAL",
  dims_use   = 1:N_PCS_THAL,
  resolution = RES_THAL
)
THAL$initial_clusters <- THAL$seurat_clusters

rm(seurat_list); gc()

# ---- 3c. QC post-integrazione ----------------------------------------------

plot_integration_qc <- function(obj, label) {
  
  # PCA vs Harmony: verifica che il batch sia rimosso
  p_pca  <- DimPlot(obj, reduction = "pca",    group.by = "orig.ident",
                    pt.size = 0.3) +
    ggtitle(paste(label, "— PCA (pre-Harmony)")) +
    theme(legend.text = element_text(size = 7))
  p_harm <- DimPlot(obj, reduction = "harmony", group.by = "orig.ident",
                    pt.size = 0.3) +
    ggtitle(paste(label, "— Harmony (post-integrazione)")) +
    theme(legend.text = element_text(size = 7))
  savepng((p_pca | p_harm) +
            plot_annotation(title = paste("Integrazione batch —", label)),
          paste0("03_harmony_", label, ".png"),
          width = 5500, height = 2500)
  
  # UMAP: cluster, trattamento, timepoint, campione + QC metrics
  p_cl   <- DimPlot(obj, reduction = "umap", label = TRUE,
                    label.size = 3, pt.size = 0.3) +
    ggtitle(paste(label, "— Cluster"))
  p_trt  <- DimPlot(obj, reduction = "umap", group.by = "treatment",
                    pt.size = 0.3, cols = pal_trt) +
    ggtitle(paste(label, "— Trattamento"))
  p_time <- DimPlot(obj, reduction = "umap", group.by = "time",
                    pt.size = 0.3, cols = pal_time) +
    ggtitle(paste(label, "— Timepoint"))
  p_samp <- DimPlot(obj, reduction = "umap", group.by = "orig.ident",
                    pt.size = 0.3) +
    ggtitle(paste(label, "— Campione"))
  p_mt   <- FeaturePlot(obj, features = "percent.mt",
                        reduction = "umap", pt.size = 0.3) +
    scale_color_viridis_c(option = "magma") +
    ggtitle(paste(label, "— %MT"))
  p_feat <- FeaturePlot(obj, features = "nFeature_RNA",
                        reduction = "umap", pt.size = 0.3) +
    scale_color_viridis_c(option = "viridis") +
    ggtitle(paste(label, "— nFeature"))
  savepng(
    (p_cl | p_trt | p_time) / (p_samp | p_mt | p_feat) +
      plot_annotation(title = paste("UMAP —", label)),
    paste0("03_UMAP_panel_", label, ".png"),
    width = 7000, height = 4500
  )
  
  # Heatmap proporzione cluster × campione
  # Se un cluster è presente solo in un campione = probabile artefatto di batch
  ct_prop <- prop.table(table(Idents(obj), obj$orig.ident), margin = 2)
  png(file.path(outdir, paste0("03_cluster_composition_", label, ".png")),
      width = 3000, height = max(2000, nrow(ct_prop) * 80 + 500), res = 300)
  pheatmap(ct_prop,
           color           = colorRampPalette(c("white","#2E86C1"))(50),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize = 8, cluster_cols = FALSE,
           main = paste("Proporzione cluster per campione —", label))
  dev.off()
  
  message(sprintf("  %s: %d cluster | %d cellule",
                  label, length(levels(Idents(obj))), ncol(obj)))
}

plot_integration_qc(CTX,  "CTX")
plot_integration_qc(THAL, "THAL")

# ---- 3d. Panel comparativo CTX vs THAL -------------------------------------
savepng(
  (DimPlot(CTX,  reduction = "harmony", group.by = "orig.ident", pt.size = 0.3) +
     ggtitle("CTX — Harmony") |
     DimPlot(THAL, reduction = "harmony", group.by = "orig.ident", pt.size = 0.3) +
     ggtitle("THAL — Harmony")) /
    (DimPlot(CTX,  reduction = "umap", label = TRUE,
             label.size = 3, pt.size = 0.3) +
       ggtitle("CTX — Clustering iniziale") |
       DimPlot(THAL, reduction = "umap", label = TRUE,
               label.size = 3, pt.size = 0.3) +
       ggtitle("THAL — Clustering iniziale")) +
    plot_annotation(title = "Integrazione e clustering — CTX vs THAL"),
  "03_harmony_by_area.png", width = 6000, height = 5000
)

# ==============================================================================
# SALVATAGGIO
# ==============================================================================

message("\n=== Salvataggio ===")
save(CTX,  file = file.path(outdir, "CTX_initial.RData"))
save(THAL, file = file.path(outdir, "THAL_initial.RData"))

message(sprintf(
  "Fatto.\n  CTX:  %d cellule | %d cluster\n  THAL: %d cellule | %d cluster",
  ncol(CTX),  length(levels(Idents(CTX))),
  ncol(THAL), length(levels(Idents(THAL)))
))

################################################################################
#  NOTE SEURAT v4 → v5 (documentazione)
#
#  v4                                        v5
#  ───────────────────────────────────────── ──────────────────────────────────
#  SplitObject() + loop SCTransform          merge() → SCTransform() su layers
#  merge() + unione manuale VariableFeatures automatico dentro SCTransform()
#  RunHarmony(group.by.vars = "orig.ident")  IntegrateLayers(HarmonyIntegration,
#                                              group.by.vars = "orig.ident")
#  slot = "data" / SCT@data                  layer = "data" dopo JoinLayers()
#  FindMarkers senza preparazione            richiede JoinLayers() prima
#  FindIntegrationAnchors + IntegrateData    IntegrateLayers(CCAIntegration)
#  SplitObject per area → Harmony separato  merge() per area → IntegrateLayers
################################################################################