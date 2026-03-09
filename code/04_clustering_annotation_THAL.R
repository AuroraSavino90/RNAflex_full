################################################################################
##  04_clustering_annotation_THAL.R
##
##  Approccio alternativo: clustering non supervisionato sull'intero THAL
##  (senza filtro a priori per ROI talamica), poi annotazione dei cluster
##  tramite MapMyCells (class, subclass, cluster_alias) e DEG tra cluster.
##
##  Flusso:
##    A. Carica THAL_initial.RData, SCTransform + PCA + UMAP + clustering
##    B. Assegna a ogni cluster Seurat la label MapMyCells per plurality vote
##       (class_name, subclass_name, cluster_name)
##    C. Plot diagnostici (UMAP, composizione, purezza cluster)
##    D. DEG tra cluster annotati (FindAllMarkers)
##    E. Salvataggio
##
##  Input:
##    - results/20260224/THAL_initial.RData
##    - MAPPING_CSV (stesso file CSV MapMyCells già usato in 04_mapmycells_THAL.R)
##
##  Output:
##    - THAL_clustered_annotated.RData
##    - 04CL_UMAP_seurat_clusters.png
##    - 04CL_UMAP_mmc_class.png / _subclass.png / _cluster.png
##    - 04CL_cluster_annotation_table.csv
##    - 04CL_DEG_allmarkers.csv
##    - 04CL_dotplot_top_markers.png
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

library(dplyr)
library(ggplot2)
library(Seurat)
library(openxlsx)

# --------------------------------------------------------------------------
# PARAMETRI
# --------------------------------------------------------------------------
MAPPING_CSV <- file.path(outdir,
  "THAL_from_mapmycells/THAL_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1772014736467.csv")

N_PCS         <- 30
RESOLUTION    <- 0.5   # aumenta per più cluster, diminuisci per meno
MIN_PCT       <- 0.25
LOGFC_THRESH  <- 0.5
TOP_N_MARKERS <- 5     # marker per cluster nel dotplot

# ==============================================================================
# A. SCTransform + PCA + UMAP + Clustering
# ==============================================================================

message("=== A. Caricamento e preprocessing THAL_initial ===")

load(file.path(outdir, "THAL_initial.RData"))   # → THAL

# JoinLayers se necessario (Seurat v5 con oggetti merged)
if (inherits(THAL[["RNA"]], "Assay5")) {
  DefaultAssay(THAL) <- "RNA"
  THAL <- JoinLayers(THAL)
}
DefaultAssay(THAL) <- "RNA"

message(sprintf("  THAL_initial: %d cellule, %d geni",
                ncol(THAL), nrow(THAL)))

options(future.globals.maxSize = 4 * 1024^3)

THAL <- SCTransform(THAL, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:N_PCS, verbose = FALSE) %>%
  FindNeighbors(dims = 1:N_PCS, verbose = FALSE) %>%
  FindClusters(resolution = RESOLUTION, verbose = FALSE)

message(sprintf("  Cluster trovati: %d (resolution = %.2f)",
                nlevels(THAL$seurat_clusters), RESOLUTION))

# ==============================================================================
# B. Assegna label MapMyCells per plurality vote
# ==============================================================================

message("=== B. Annotazione cluster con MapMyCells ===")

mapping <- read.csv(MAPPING_CSV, comment.char = "#")
rownames(mapping) <- mapping$cell_id

# Trasferisci label MapMyCells alle cellule di THAL
THAL$mmc_class    <- mapping[colnames(THAL), "class_name"]
THAL$mmc_subclass <- mapping[colnames(THAL), "subclass_name"]
THAL$mmc_cluster  <- mapping[colnames(THAL), "cluster_name"]

if ("cluster_alias" %in% colnames(mapping)) {
  THAL$mmc_cluster_alias <- as.integer(
    mapping[colnames(THAL), "cluster_alias"]
  )
}

# Plurality vote per cluster Seurat: label più frequente per ogni livello
plurality_vote <- function(labels, cluster_ids) {
  clusters <- levels(factor(cluster_ids))
  result <- sapply(clusters, function(cl) {
    cells_cl <- cluster_ids == cl
    tab <- table(labels[cells_cl], useNA = "no")
    if (length(tab) == 0) return(NA_character_)
    names(which.max(tab))
  })
  setNames(result, clusters)
}

vote_class    <- plurality_vote(THAL$mmc_class,    THAL$seurat_clusters)
vote_subclass <- plurality_vote(THAL$mmc_subclass, THAL$seurat_clusters)
vote_cluster  <- plurality_vote(THAL$mmc_cluster,  THAL$seurat_clusters)

# Assegna label per cluster
# Sostituisci le tre righe di assegnazione con:
cl_ids <- as.character(THAL$seurat_clusters)

THAL$cluster_mmc_class    <- unname(vote_class[cl_ids])
THAL$cluster_mmc_subclass <- unname(vote_subclass[cl_ids])
THAL$cluster_mmc_cluster  <- unname(vote_cluster[cl_ids])
# Label compatta per plot: "N — Subclass"
THAL$cluster_label <- paste0(
  as.character(THAL$seurat_clusters), " — ",
  THAL$cluster_mmc_subclass
)

# Purezza: % cellule nel cluster con la label maggioritaria
purity_df <- THAL@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells           = n(),
    mmc_class_vote    = vote_class[as.character(seurat_clusters[1])],
    mmc_subclass_vote = vote_subclass[as.character(seurat_clusters[1])],
    mmc_cluster_vote  = vote_cluster[as.character(seurat_clusters[1])],
    purity_class      = mean(mmc_class    == vote_class[as.character(seurat_clusters[1])],    na.rm = TRUE),
    purity_subclass   = mean(mmc_subclass == vote_subclass[as.character(seurat_clusters[1])], na.rm = TRUE),
    purity_cluster    = mean(mmc_cluster  == vote_cluster[as.character(seurat_clusters[1])],  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(as.character(seurat_clusters)))
message("\nAnnotazione cluster (plurality vote):")
print(purity_df %>%
        dplyr::select(seurat_clusters, n_cells, mmc_subclass_vote,
                      purity_subclass, mmc_class_vote, purity_class) %>%
        mutate(across(starts_with("purity"), ~ round(.x * 100, 1))))
write.csv(purity_df,
          file.path(outdir, "04CL_cluster_annotation_table.csv"),
          row.names = FALSE)

# ==============================================================================
# C. Plot diagnostici
# ==============================================================================

message("=== C. Plot diagnostici ===")

# Palette subclass
subclass_levels <- sort(unique(na.omit(THAL$cluster_mmc_subclass)))
pal_subclass    <- setNames(
  scales::hue_pal()(length(subclass_levels)),
  subclass_levels
)

# Palette class
class_levels <- sort(unique(na.omit(THAL$cluster_mmc_class)))
pal_class    <- setNames(
  scales::hue_pal()(length(class_levels)),
  class_levels
)

# 1. UMAP cluster Seurat numerici
p_clust <- DimPlot(THAL, group.by = "seurat_clusters",
                   label = TRUE, label.size = 3, pt.size = 0.3,
                   repel = TRUE) +
  ggtitle(sprintf("Cluster Seurat (res=%.2f) — THAL_initial completo",
                  RESOLUTION)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
savepng(p_clust, "04CL_UMAP_seurat_clusters.png", width = 4000, height = 3500)

# 2. UMAP colorato per MMC class
p_class <- DimPlot(THAL, group.by = "cluster_mmc_class",
                   cols = pal_class, label = TRUE, label.size = 2.5,
                   pt.size = 0.3, repel = TRUE) +
  ggtitle("Cluster — MMC class (plurality vote)") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1))
savepng(p_class, "04CL_UMAP_mmc_class.png", width = 4500, height = 3500)

# 3. UMAP colorato per MMC subclass
p_sub <- DimPlot(THAL, group.by = "cluster_mmc_subclass",
                 cols = pal_subclass, label = TRUE, label.size = 2.2,
                 pt.size = 0.3, repel = TRUE) +
  ggtitle("Cluster — MMC subclass (plurality vote)") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1))
savepng(p_sub, "04CL_UMAP_mmc_subclass.png", width = 5000, height = 3500)

# 4. UMAP con label compatta (N — Subclass)
umap_centroids <- THAL@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  mutate(cl  = THAL$seurat_clusters,
         lbl = THAL$cluster_label) %>%
  group_by(cl) %>%
  summarise(umap_1 = median(umap_1),
            umap_2 = median(umap_2),
            lbl    = lbl[1],
            .groups = "drop")

p_label <- DimPlot(THAL, group.by = "seurat_clusters",
                   label = FALSE, pt.size = 0.3) +
  geom_text(data    = umap_centroids,
            mapping = aes(x = umap_1, y = umap_2, label = lbl),
            size = 2.5, fontface = "bold", inherit.aes = FALSE) +
  ggtitle("Cluster — label: N — MMC subclass") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
savepng(p_label, "04CL_UMAP_cluster_label.png", width = 5000, height = 3500)

# 5. Barplot purezza subclass per cluster
p_purity <- ggplot(purity_df,
                   aes(x = reorder(seurat_clusters, -purity_subclass),
                       y = purity_subclass * 100,
                       fill = purity_subclass)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%\n%s\n(n=%d)",
                                purity_subclass * 100,
                                mmc_subclass_vote, n_cells)),
            vjust = -0.2, size = 2.2, lineheight = 0.85) +
  scale_fill_gradient(low = "#F4A261", high = "#2D6A4F",
                      name = "Purezza") +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(title    = "Purezza MMC subclass per cluster Seurat",
       subtitle = "% cellule con label maggioritaria | etichetta = subclass (n cellule)",
       x = "Cluster Seurat", y = "% purezza subclass") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
savepng(p_purity, "04CL_barplot_cluster_purity.png", width = 4500, height = 2500)

# 6. Composizione class per cluster (stacked barplot)
comp_df <- THAL@meta.data %>%
  dplyr::count(seurat_clusters, mmc_class) %>%
  group_by(seurat_clusters) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df,
                 aes(x = seurat_clusters, y = pct, fill = mmc_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_class, na.value = "grey70",
                    name = "MMC class") +
  labs(title = "Composizione MMC class per cluster Seurat",
       x = "Cluster", y = "% cellule") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
savepng(p_comp, "04CL_barplot_class_composition.png",
        width = 5000, height = 2500)

# ==============================================================================
# D. DEG tra cluster annotati
# ==============================================================================

message("=== D. FindAllMarkers per cluster Seurat ===")

DefaultAssay(THAL) <- "SCT"
Idents(THAL)       <- "seurat_clusters"

# PrepSCTFindMarkers necessario con SCTransform v2 multi-sample
THAL <- PrepSCTFindMarkers(THAL)

markers_all <- FindAllMarkers(
  THAL,
  only.pos        = TRUE,
  min.pct         = MIN_PCT,
  logfc.threshold = LOGFC_THRESH,
  test.use        = "wilcox",
  verbose         = TRUE
)


# Aggiungi label subclass al cluster
markers_all <- markers_all %>%
  mutate(
    subclass_label = vote_subclass[as.character(cluster)],
    cluster_label  = paste0(cluster, " — ", subclass_label)
  ) %>%
  arrange(cluster, desc(avg_log2FC))

message(sprintf("  DEG totali: %d", nrow(markers_all)))
message("  DEG per cluster:")
print(table(markers_all$cluster))

write.csv(markers_all,
          file.path(outdir, "04CL_DEG_allmarkers.csv"),
          row.names = FALSE)

# Top N marker per cluster
top_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = TOP_N_MARKERS, with_ties = FALSE) %>%
  ungroup()

write.csv(top_markers,
          file.path(outdir, "04CL_DEG_top_markers.csv"),
          row.names = FALSE)

# Dotplot top marker
top_genes_ordered <- top_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

# Usa cluster_label come identity per il dotplot
THAL$cluster_label_factor <- factor(
  THAL$cluster_label,
  levels = unique(paste0(sort(as.integer(levels(THAL$seurat_clusters))),
                         " — ",
                         vote_subclass[as.character(sort(as.integer(
                           levels(THAL$seurat_clusters))))]))
)
Idents(THAL) <- "cluster_label_factor"

p_dot <- DotPlot(THAL,
                 features  = top_genes_ordered,
                 cols      = c("grey90", "#1D3557"),
                 dot.scale = 5) +
  coord_flip() +
  labs(title    = sprintf("Top %d markers per cluster (logFC ordinato)",
                          TOP_N_MARKERS),
       subtitle = "Cluster etichettati: N — MMC subclass",
       x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size  = 7),
    plot.title   = element_text(face  = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

savepng(p_dot, "04CL_dotplot_top_markers.png",
        width  = max(4000, length(levels(Idents(THAL))) * 300 + 2000),
        height = max(3000, length(top_genes_ordered) * 120 + 1000))

message("Salvato: 04CL_dotplot_top_markers.png")

# ==============================================================================
# E. Salvataggio
# ==============================================================================

message("=== E. Salvataggio ===")

# Ripristina Idents a seurat_clusters per coerenza downstream
Idents(THAL) <- "seurat_clusters"

save(THAL, markers_all, top_markers, purity_df, vote_subclass,
     file = file.path(outdir, "THAL_clustered_annotated.RData"))

message("=== 04_clustering_annotation_THAL completato ===")
message(sprintf("  Output: %s", outdir))
message("  THAL_clustered_annotated.RData")
message("  04CL_cluster_annotation_table.csv")
message("  04CL_DEG_allmarkers.csv")
message("  04CL_DEG_top_markers.csv")
message("  04CL_*.png")


################################################################################
##  SEZIONE DESeq2 — da attaccare a 04_clustering_annotation_THAL.R
##
##  Calcola DEG LSD vs Ctrl per ogni cluster Seurat × timepoint,
##  usando la stessa run_deseq2() di 07_deseq2_THAL.R.
##  La colonna identity è "seurat_clusters"; la label del cluster
##  (N — MMC subclass) viene aggiunta ai risultati per leggibilità.
################################################################################

library(DESeq2)

# --------------------------------------------------------------------------
# PARAMETRI (coerenti con 07_deseq2_THAL.R)
# --------------------------------------------------------------------------
FDR_THR             <- 0.05
LFC_THR             <- 0.25
MIN_CELLS           <- 20
MIN_COUNTS          <- 10
MIN_CELLS_EXPRESSED <- 5

lfc_type <- if (requireNamespace("apeglm", quietly = TRUE)) {
  "apeglm"
} else if (requireNamespace("ashr", quietly = TRUE)) {
  message("INFO: apeglm non disponibile, uso ashr"); "ashr"
} else {
  message("INFO: uso normal per lfcShrink"); "normal"
}
message(sprintf("LFC shrinkage: %s", lfc_type))

# --------------------------------------------------------------------------
# FUNZIONE run_deseq2 (identica a 07_deseq2_THAL.R,
# generalizzata per qualsiasi colonna identity)
# --------------------------------------------------------------------------
run_deseq2 <- function(obj, area_label, time_val,
                       cell_type_col = "seurat_clusters",
                       condition_col = "treatment",
                       assay_use     = "RNA",
                       min_cells     = MIN_CELLS) {
  
  cells_t <- colnames(obj)[obj$time == time_val]
  if (length(cells_t) == 0) {
    message(sprintf("  [%s | %s] nessuna cellula", area_label, time_val))
    return(list())
  }
  obj_t <- obj[, cells_t]
  
  cell_types <- sort(unique(as.character(unlist(obj_t[[cell_type_col]]))))
  degs_list  <- list()
  
  for (ct in cell_types) {
    cells_ct <- colnames(obj_t)[as.character(unlist(obj_t[[cell_type_col]])) == ct]
    if (length(cells_ct) < min_cells) {
      message(sprintf("    skip %s/%s/cluster%s: solo %d cellule (min=%d)",
                      area_label, time_val, ct, length(cells_ct), min_cells))
      next
    }
    
    obj_ct    <- obj_t[, cells_ct]
    DefaultAssay(obj_ct) <- assay_use
    count_mat <- LayerData(obj_ct, assay = assay_use, layer = "counts")
    count_mat <- as(count_mat, "dgCMatrix")
    count_mat <- count_mat[Matrix::rowSums(count_mat > 0) >= MIN_CELLS_EXPRESSED, ]
    
    if (nrow(count_mat) < 100) {
      message(sprintf("    skip %s/%s/cluster%s: solo %d geni dopo filtro",
                      area_label, time_val, ct, nrow(count_mat)))
      next
    }
    
    cond     <- factor(unlist(obj_ct[[condition_col]]), levels = c("Ctrl", "LSD"))
    cond_tab <- table(cond)
    
    if (any(cond_tab == 0) || length(unique(cond)) < 2) {
      message(sprintf("    skip %s/%s/cluster%s: condizione mancante (%s)",
                      area_label, time_val, ct,
                      paste(names(cond_tab[cond_tab == 0]), collapse = ",")))
      next
    }
    
    MIN_MC_PER_COND <- 3
    if (min(cond_tab) < MIN_MC_PER_COND) {
      message(sprintf(
        "    skip %s/%s/cluster%s: troppo poche unità (Ctrl=%d, LSD=%d)",
        area_label, time_val, ct, cond_tab["Ctrl"], cond_tab["LSD"]))
      next
    }
    
    col_data <- data.frame(condition = cond, row.names = colnames(count_mat))
    
    tryCatch({
      
      dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(count_mat)),
        colData   = col_data,
        design    = ~ condition
      )
      dds <- estimateSizeFactors(dds, type = "poscounts")
      dds <- dds[rowSums(counts(dds)) >= MIN_COUNTS, ]
      
      fit_ok <- FALSE
      
      dds <- tryCatch({
        d <- estimateDispersions(dds, fitType = "local", quiet = TRUE)
        fit_ok <<- TRUE; d
      }, error = function(e) {
        message(sprintf("    [cluster%s] local fallito, provo mean", ct)); dds
      })
      
      if (!fit_ok) {
        dds <- tryCatch({
          d <- estimateDispersions(dds, fitType = "mean", quiet = TRUE)
          fit_ok <<- TRUE; d
        }, error = function(e) {
          message(sprintf("    [cluster%s] mean fallito, provo gene-wise", ct)); dds
        })
      }
      
      if (!fit_ok) {
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        message(sprintf("    [cluster%s] dispersioni gene-wise (fallback)", ct))
      }
      
      if (fit_ok) {
        disp_fitted <- mcols(dds)$dispFit
        disp_gene   <- mcols(dds)$dispGeneEst
        if (!is.null(disp_fitted) && !is.null(disp_gene)) {
          ratio  <- disp_fitted / disp_gene
          ratio  <- ratio[is.finite(ratio) & ratio > 0]
          spread <- if (length(ratio) > 10) diff(range(log10(ratio))) else Inf
          if (spread < 2) {
            dispersions(dds) <- mcols(dds)$dispGeneEst
            message(sprintf("    [cluster%s] spread <2 OOM: gene-wise", ct))
          }
        }
      }
      
      dds <- nbinomWaldTest(dds, quiet = TRUE)
      
      res <- tryCatch({
        if (lfc_type %in% c("apeglm", "ashr")) {
          lfcShrink(dds, coef = "condition_LSD_vs_Ctrl",
                    type = lfc_type, quiet = TRUE)
        } else {
          lfcShrink(dds, contrast = c("condition", "LSD", "Ctrl"),
                    type = "normal", quiet = TRUE)
        }
      }, error = function(e) {
        message(sprintf("    [cluster%s] lfcShrink fallito, uso results() raw", ct))
        results(dds, contrast = c("condition", "LSD", "Ctrl"))
      })
      
      df              <- as.data.frame(res)
      df$gene         <- rownames(df)
      df$cluster      <- ct
      df$cluster_label <- paste0(ct, " — ", vote_subclass[ct])
      degs_list[[ct]] <- df
      
      n_sig <- sum(!is.na(df$padj) & df$padj < FDR_THR &
                     abs(df$log2FoldChange) > LFC_THR)
      message(sprintf("    %s | %s | cluster%-3s (%-30s) n=%4d | DEG=%d",
                      area_label, time_val, ct,
                      vote_subclass[ct],
                      length(cells_ct), n_sig))
      
    }, error = function(e) {
      message(sprintf("    WARN cluster%s/%s/%s: %s",
                      ct, area_label, time_val,
                      gsub("\n", " ", conditionMessage(e))))
    })
  }
  degs_list
}

# --------------------------------------------------------------------------
# Assicura layer RNA disponibile (JoinLayers se Assay5)
# --------------------------------------------------------------------------
DefaultAssay(THAL) <- "RNA"
if (inherits(THAL[["RNA"]], "Assay5")) {
  THAL <- JoinLayers(THAL)
}

# --------------------------------------------------------------------------
# RUN DESeq2 per cluster × timepoint
# --------------------------------------------------------------------------
message("=== DESeq2 cluster × timepoint — 90min ===")
DEGs_CL_90 <- run_deseq2(THAL, "THAL_clusters", "90min",
                         cell_type_col = "seurat_clusters")

message("=== DESeq2 cluster × timepoint — 24h ===")
DEGs_CL_24 <- run_deseq2(THAL, "THAL_clusters", "24h",
                         cell_type_col = "seurat_clusters")

# --------------------------------------------------------------------------
# Salva risultati grezzi
# --------------------------------------------------------------------------
save(DEGs_CL_90, DEGs_CL_24, THAL, markers_all, top_markers, purity_df,
     vote_subclass,
     file = file.path(outdir, "THAL_clustered_annotated.RData"))
message("Salvato: THAL_clustered_annotated.RData")

# --------------------------------------------------------------------------
# CSV flat con tutti i DEG
# --------------------------------------------------------------------------
deg_flat <- bind_rows(
  bind_rows(DEGs_CL_90) %>% mutate(timepoint = "90min"),
  bind_rows(DEGs_CL_24) %>% mutate(timepoint = "24h")
) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(cluster, timepoint, padj)

write.csv(deg_flat,
          file.path(outdir, "04CL_DEG_LSDvsCtrl_allclusters.csv"),
          row.names = FALSE)
message("Salvato: 04CL_DEG_LSDvsCtrl_allclusters.csv")

# --------------------------------------------------------------------------
# Barplot DEG up/down per cluster × timepoint
# --------------------------------------------------------------------------
plot_deg_barplot <- function(degs_list, title, fdr = FDR_THR, lfc = LFC_THR) {
  df <- do.call(rbind, lapply(names(degs_list), function(ct) {
    d <- degs_list[[ct]]
    d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange), ]
    rbind(
      data.frame(cluster       = ct,
                 cluster_label = paste0(ct, " — ", vote_subclass[ct]),
                 direction     = "Up",
                 n = sum(d$padj < fdr & d$log2FoldChange >  lfc)),
      data.frame(cluster       = ct,
                 cluster_label = paste0(ct, " — ", vote_subclass[ct]),
                 direction     = "Down",
                 n = sum(d$padj < fdr & d$log2FoldChange < -lfc))
    )
  }))
  if (nrow(df) == 0 || all(df$n == 0)) return(NULL)
  
  df$n_signed     <- ifelse(df$direction == "Down", -df$n, df$n)
  total           <- tapply(df$n, df$cluster_label, sum)
  df$cluster_label <- factor(df$cluster_label,
                             levels = names(sort(total, decreasing = TRUE)))
  
  ggplot(df, aes(x = cluster_label, y = n_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#CB4335", "Down" = "#2471A3")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top") +
    labs(title    = title,
         subtitle = sprintf("FDR < %.2f | |LFC| > %.2f  (LSD vs Ctrl)",
                            fdr, lfc),
         x = "", y = "N DEG (+ = up in LSD, - = down in LSD)")
}

p_bar_90 <- plot_deg_barplot(DEGs_CL_90, "DEG per cluster — 90min")
p_bar_24 <- plot_deg_barplot(DEGs_CL_24, "DEG per cluster — 24h")

if (!is.null(p_bar_90))
  savepng(p_bar_90, "04CL_DEG_barplot_90min.png",
          width = 3500,
          height = max(1800, length(DEGs_CL_90) * 100 + 600))

if (!is.null(p_bar_24))
  savepng(p_bar_24, "04CL_DEG_barplot_24h.png",
          width = 3500,
          height = max(1800, length(DEGs_CL_24) * 100 + 600))

# --------------------------------------------------------------------------
# Volcano per cluster con più DEG (90min e 24h)
# --------------------------------------------------------------------------
plot_volcano <- function(df, title, fdr = FDR_THR, lfc = LFC_THR, top_n = 15) {
  df <- df[!is.na(df$padj), ]
  df$color <- ifelse(df$padj < fdr & df$log2FoldChange >  lfc, "up",
                     ifelse(df$padj < fdr & df$log2FoldChange < -lfc, "down", "ns"))
  top <- df[df$color != "ns", ]
  top <- top[order(top$padj), ][seq_len(min(top_n, nrow(top))), ]
  
  ggplot(df, aes(x = log2FoldChange,
                 y = -log10(padj + 1e-300),
                 color = color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc, lfc),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(fdr),
               linetype = "dashed", color = "gray50") +
    ggrepel::geom_text_repel(data = top, aes(label = gene),
                             size = 2.8, max.overlaps = 20,
                             box.padding = 0.4) +
    scale_color_manual(
      values = c("up" = "#E74C3C", "down" = "#3498DB", "ns" = "gray70")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(title = title,
         x = paste0("log2FC (", lfc_type, ")"),
         y = "-log10(FDR)")
}

for (tp in c("90min", "24h")) {
  degs_l <- if (tp == "90min") DEGs_CL_90 else DEGs_CL_24
  if (length(degs_l) == 0) next
  n_sig  <- sapply(degs_l, function(d)
    sum(!is.na(d$padj) & d$padj < FDR_THR & abs(d$log2FoldChange) > LFC_THR))
  top_cl <- names(which.max(n_sig))
  if (length(top_cl) == 0) next
  p_v <- plot_volcano(
    degs_l[[top_cl]],
    sprintf("Volcano %s | cluster %s — %s", tp, top_cl, vote_subclass[top_cl])
  )
  savepng(p_v,
          sprintf("04CL_Volcano_%s_cluster%s.png", tp, top_cl),
          width = 3000, height = 2500)
}

message("=== DESeq2 cluster completato ===")

