################################################################################
##  10_subclustering.R                                                         ##
##                                                                             ##
##  Input:  CTX_clean.RData, THAL_clean.RData  (da sezione 05)               ##
##  Output: CTX_neurons.RData, THAL_neurons.RData                             ##
##          CTX_oligo.RData, THAL_oligo.RData                                 ##
##                                                                             ##
##  Seurat v5: IntegrateLayers + JoinLayers nel subclustering                 ##
##  Nessun SplitObject, nessun RunHarmony diretto                             ##
################################################################################

source("00_config.R")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_clean.RData"))   # CTX_clean
load(file.path(outdir, "THAL_clean.RData"))  # THAL_clean

# ==============================================================================
# PARAMETRI
# ==============================================================================

N_PCS_NEU_CTX  <- 30   # <<< aggiornare dopo elbow
N_PCS_NEU_THAL <- 20   # <<< aggiornare dopo elbow
N_PCS_OLIGO    <- 15   # <<< aggiornare dopo elbow

RES_NEURONS <- 0.5
RES_OLIGO   <- 0.4

# ==============================================================================
# FUNZIONE SUBCLUSTERING (Seurat v5)
# ==============================================================================

run_subclustering_v5 <- function(obj, cell_types_keep, label,
                                  n_pcs      = 30,
                                  resolution = 0.5,
                                  n_features = 2000) {

  # Seleziona cellule dei tipi di interesse
  cells_keep <- colnames(obj)[obj$cell_type %in% cell_types_keep]
  if (length(cells_keep) < 50) {
    message(sprintf("  %s: troppo poche cellule (%d) — skip",
                    label, length(cells_keep)))
    return(NULL)
  }

  sub_obj <- obj[, cells_keep]
  message(sprintf("  %s: %d cellule, %d tipi selezionati",
                  label, ncol(sub_obj),
                  length(unique(sub_obj$cell_type))))

  # Seurat v5: dopo il subset i layer sono ancora separati per campione
  # (il subset eredita la struttura layered dal genitore)
  # → SCTransform + IntegrateLayers + JoinLayers

  sub_obj <- SCTransform(
    sub_obj,
    vars.to.regress     = "percent.mt",
    variable.features.n = n_features,
    verbose             = FALSE
  )

  sub_obj <- RunPCA(sub_obj, npcs = 50, verbose = FALSE)

  savepng(
    ElbowPlot(sub_obj, ndims = 50) + ggtitle(paste("Elbow —", label)),
    paste0("10_elbow_", label, ".png")
  )

  sub_obj <- IntegrateLayers(
    object         = sub_obj,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    group.by.vars  = "orig.ident",
    verbose        = FALSE
  )

  # JoinLayers sull'assay RNA (non SCT)
  sub_obj <- JoinLayers(sub_obj, assay = "RNA")

  sub_obj <- RunUMAP(sub_obj, reduction = "harmony",
                     dims = 1:n_pcs, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony",
                            dims = 1:n_pcs, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolution,
                           algorithm = 4, verbose = FALSE)

  n_cl <- length(levels(sub_obj$seurat_clusters))
  message(sprintf("  %s: %d cluster (res=%.2f)", label, n_cl, resolution))

  # --- Silhouette per valutare la qualità del clustering ---------------------
  emb <- Embeddings(sub_obj, "harmony")[, 1:n_pcs]
  sil <- tryCatch(
    silhouette(as.integer(sub_obj$seurat_clusters), dist(emb)),
    error = function(e) NULL
  )
  if (!is.null(sil)) {
    sub_obj$silhouette <- sil[, 3]
    p_sil <- ggplot(
      data.frame(cluster   = sub_obj$seurat_clusters,
                 sil_width = sub_obj$silhouette),
      aes(x = cluster, y = sil_width, fill = cluster)
    ) +
      geom_boxplot(outlier.size = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      theme_bw(base_size = 9) + theme(legend.position = "none") +
      labs(title = paste("Silhouette —", label),
           x = "Cluster", y = "Silhouette width")
    savepng(p_sil, paste0("10_silhouette_", label, ".png"),
            width = 2500, height = 2000)
  }

  sub_obj
}

# ==============================================================================
# 10.1 — NEURONI CTX
# ==============================================================================

message("=== 10.1 Subclustering neuroni CTX ===")

# Identifica tipi cellulari neurali nel CTX (adattare ai nomi usati in 05)
neuron_ct_ctx <- grep(
  "Ex|Inh|L[2-6]|PV|SST|VIP|Lamp|Sncg|IT|PT|CT|NP",
  unique(CTX_clean$cell_type),
  value = TRUE, ignore.case = TRUE
)
message(sprintf("  CT neurali CTX selezionati: %s",
                paste(neuron_ct_ctx, collapse = ", ")))

CTX_neurons <- run_subclustering_v5(
  CTX_clean, neuron_ct_ctx, "CTX_neurons",
  n_pcs = N_PCS_NEU_CTX, resolution = RES_NEURONS
)

if (!is.null(CTX_neurons)) {

  # DotPlot marcatori
  genes_neu_ctx <- filter_genes(
    unique(c(markers_exc_ctx, markers_inh_ctx, markers_htr2a_signature)),
    CTX_neurons
  )
  p_dot <- DotPlot(CTX_neurons, features = genes_neu_ctx,
                   cluster.idents = TRUE) +
    RotatedAxis() +
    ggtitle("Marcatori neuroni — CTX")
  savepng(p_dot, "10_DotPlot_CTX_neurons.png",
          width = 7000, height = 3000)

  # UMAP: cluster + trattamento + Htr2a
  p_cl  <- DimPlot(CTX_neurons, label = TRUE, label.size = 3, pt.size = 0.3) +
    ggtitle("CTX Neurons — cluster")
  p_trt <- DimPlot(CTX_neurons, group.by = "treatment", pt.size = 0.3,
                   cols = pal_trt) +
    ggtitle("CTX Neurons — trattamento")

  htr2a_present <- "Htr2a" %in% rownames(CTX_neurons)
  p_htr2a <- if (htr2a_present) {
    FeaturePlot(CTX_neurons, features = "Htr2a",
                reduction = "umap", pt.size = 0.3) +
      scale_color_viridis_c(option = "inferno") +
      ggtitle("Htr2a — CTX neurons")
  } else {
    ggplot() + theme_void() + labs(title = "Htr2a non presente")
  }
  savepng((p_cl | p_trt | p_htr2a),
          "10_UMAP_CTX_neurons.png", width = 7500, height = 2500)

  # FeaturePlot Htr2a split per trattamento
  if (htr2a_present) {
    p_split <- FeaturePlot(CTX_neurons, features = "Htr2a",
                           split.by = "treatment", reduction = "umap",
                           pt.size = 0.3) +
      plot_annotation(title = "Htr2a — CTX neurons | Ctrl vs LSD")
    savepng(p_split, "10_Htr2a_split_CTX_neurons.png",
            width = 5000, height = 2500)
  }

  save(CTX_neurons, file = file.path(outdir, "CTX_neurons.RData"))
  message(sprintf("  Salvato: CTX_neurons.RData (%d cellule, %d cluster)",
                  ncol(CTX_neurons),
                  length(levels(CTX_neurons$seurat_clusters))))
}

# ==============================================================================
# 10.2 — NEURONI THAL
# ==============================================================================

message("=== 10.2 Subclustering neuroni THAL ===")

neuron_ct_thal <- grep(
  "TC|TRN|relay|neuron|Excit|Inhib",
  unique(THAL_clean$cell_type),
  value = TRUE, ignore.case = TRUE
)
message(sprintf("  CT neurali THAL selezionati: %s",
                paste(neuron_ct_thal, collapse = ", ")))

THAL_neurons <- run_subclustering_v5(
  THAL_clean, neuron_ct_thal, "THAL_neurons",
  n_pcs = N_PCS_NEU_THAL, resolution = RES_NEURONS
)

if (!is.null(THAL_neurons)) {

  genes_neu_thal <- filter_genes(
    unique(c(markers_thal, markers_htr2a_signature)), THAL_neurons
  )
  p_dot_thal <- DotPlot(THAL_neurons, features = genes_neu_thal,
                        cluster.idents = TRUE) +
    RotatedAxis() + ggtitle("Marcatori neuroni — THAL")
  savepng(p_dot_thal, "10_DotPlot_THAL_neurons.png",
          width = 6000, height = 3000)

  p_cl_t  <- DimPlot(THAL_neurons, label = TRUE, label.size = 3, pt.size = 0.3) +
    ggtitle("THAL Neurons — cluster")
  p_trt_t <- DimPlot(THAL_neurons, group.by = "treatment", pt.size = 0.3,
                     cols = pal_trt) +
    ggtitle("THAL Neurons — trattamento")
  savepng(p_cl_t | p_trt_t, "10_UMAP_THAL_neurons.png",
          width = 5000, height = 2500)

  save(THAL_neurons, file = file.path(outdir, "THAL_neurons.RData"))
  message(sprintf("  Salvato: THAL_neurons.RData (%d cellule, %d cluster)",
                  ncol(THAL_neurons),
                  length(levels(THAL_neurons$seurat_clusters))))
}

# ==============================================================================
# 10.3 — OLIGODENDROCITI (entrambe le aree)
# ==============================================================================

message("=== 10.3 Subclustering oligodendrociti ===")

oligo_ct_ctx  <- grep("Oligo|OPC|NF", unique(CTX_clean$cell_type),
                      value = TRUE, ignore.case = TRUE)
oligo_ct_thal <- grep("Oligo|OPC|NF", unique(THAL_clean$cell_type),
                      value = TRUE, ignore.case = TRUE)

CTX_oligo  <- run_subclustering_v5(CTX_clean,  oligo_ct_ctx,  "CTX_oligo",
                                    n_pcs = N_PCS_OLIGO, resolution = RES_OLIGO)
THAL_oligo <- run_subclustering_v5(THAL_clean, oligo_ct_thal, "THAL_oligo",
                                    n_pcs = N_PCS_OLIGO, resolution = RES_OLIGO)

for (obj_nm in c("CTX_oligo","THAL_oligo")) {
  obj <- get(obj_nm)
  if (is.null(obj)) next

  genes_ol <- filter_genes(markers_oligo_fine, obj)
  p_dot_ol <- DotPlot(obj, features = genes_ol,
                      cluster.idents = TRUE) +
    RotatedAxis() + ggtitle(paste("Marcatori oligo —", obj_nm))
  savepng(p_dot_ol, paste0("10_DotPlot_", obj_nm, ".png"),
          width = 5000, height = 2500)

  # Proporzione cluster per campione × trattamento
  frac_df <- as.data.frame(
    prop.table(table(obj$orig.ident, obj$seurat_clusters), margin = 1)
  )
  colnames(frac_df) <- c("sample","cluster","fraction")
  frac_df$treatment <- ifelse(frac_df$sample %in% c("C1","C3","T1","T3"),
                              "Ctrl","LSD")

  p_prop <- ggplot(frac_df, aes(x = cluster, y = fraction, fill = treatment)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = pal_trt) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Proporzioni cluster oligo —", obj_nm),
         x = "", y = "Frazione per campione")
  savepng(p_prop, paste0("10_oligo_proportions_", obj_nm, ".png"))

  save(list = obj_nm,
       file = file.path(outdir, paste0(obj_nm, ".RData")))
  message(sprintf("  Salvato: %s.RData (%d cellule, %d cluster)",
                  obj_nm, ncol(obj),
                  length(levels(obj$seurat_clusters))))
}

message("=== 10 completato ===")
