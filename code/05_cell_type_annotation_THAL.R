################################################################################
##  05_cell_type_annotation_THAL.R                                             ##
##                                                                             ##
##  Input:  THAL_annotated_MMC_fin.RData  (da 04_mapmycells_THAL)            ##
##          Oggetto THAL filtrato: TH_high + non-neuroni (broad_anat NA)      ##
##                                                                             ##
##  Output: THAL/THAL_celltype.RData                                          ##
##                                                                             ##
##  Flusso:                                                                    ##
##    PARTE A — ispezione visiva (eseguire sempre)                            ##
##      1. DotPlot marcatori per cluster con strip cell type                  ##
##      2. DotPlot QC per cluster (nCount, nFeature, percent.mt)             ##
##    PARTE B — compilare mappa rinomina, poi eseguire                        ##
##      3. Rimozione cluster bassa qualità / contaminanti                     ##
##      4. RenameIdents → colonna cell_type                                   ##
##    PARTE C — plot post-annotazione                                          ##
##      5. AddModuleScore per cell type                                        ##
##      6. FeaturePlot geni chiave con label CT                               ##
##      7. Heatmap cluster × MapMyCells subclass                              ##
##      8. UMAP annotato + composizione CT × campione                         ##
################################################################################

source("code/00_config.R")

outdir   <- "results/20260224"
thal_dir <- file.path(outdir, "THAL")
dir.create(thal_dir, recursive = TRUE, showWarnings = FALSE)

# Helper savepng nella sottodirectory THAL
savepng_thal <- function(p, filename, width = 3000, height = 2500, res = 300) {
  png(file.path(thal_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "THAL_annotated_MMC_fin.RData"))   # THAL
obj_nm <- "THAL"

message("\nLivelli cluster THAL: ", paste(levels(Idents(THAL)), collapse = " | "))

# ==============================================================================
# PARTE A — ISPEZIONE VISIVA DEI CLUSTER
# ==============================================================================

message("\n=== 05A. Ispezione visiva cluster THAL ===")

# ------------------------------------------------------------------------------
# Mappa gene → cell type per il talamo
# Ordine: gliali → endotelio → TC relay → TC higher-order → TC identità → TRN
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Mappa gene → cell type per il talamo
# Ordine: gliali → endotelio → TC relay → TC higher-order → TC identità → TRN
# ------------------------------------------------------------------------------
gene_ct_map_thal <- c(
  # -----------------------------------------------------------------------
  # CELLULE NON NEURONALI
  # -----------------------------------------------------------------------
  
  # Astrociti — marcatori stabili e specifici
  Aqp4    = "Astro",    Aldoc   = "Astro",    Mlc1    = "Astro",
  Acsbg1  = "Astro",    Gjb6    = "Astro",    S100b   = "Astro",
  
  # Oligodendrociti maturi
  Mog     = "Oligo",    Mbp     = "Oligo",    Mobp    = "Oligo",
  Plp1    = "Oligo",    Ermn    = "Oligo",    Aspa    = "Oligo",
  
  # NFOL (Newly Formed Oligodendrocyte)
  Bcas1   = "NFOL",     Enpp6   = "NFOL",     Tmem2   = "NFOL",
  
  # OPC
  Pdgfra  = "OPC",      Cspg4   = "OPC",      Ptprz1  = "OPC",
  
  # Microglia
  P2ry12  = "Micro",    Cx3cr1  = "Micro",    Tmem119 = "Micro",
  C1qa    = "Micro",    Hexb    = "Micro",
  
  # Endotelio + periciti
  Cldn5   = "Endo",     Flt1    = "Endo",     Pecam1  = "Endo",
  Rgs5    = "Peri",     Pdgfrb  = "Peri",
  
  # -----------------------------------------------------------------------
  # NEURONI TALAMICI GLUTAMATERGICI
  # -----------------------------------------------------------------------
  #
  # Framework da Phillips et al. 2019 (Nat Neurosci) — tre profili molecolari
  # lungo un continuum, validati da retrogradely-labeled scRNA-seq:
  #
  #   TC_FO  = First-order relay (VPM, VPL, dLGN, MGv, VB complex)
  #            → profilo "primario": Tnnt1+, Prkcd+, Scnn1a+
  #            Tnnt1 (troponin T slow): gene più selettivo per relay sensoriali
  #            nel topo adulto (Nat Neurosci 2019; Sci Reports 2025 spatial);
  #            Scnn1a: VPM/VPL/dLGN-core (non è L4 cortex nel talamo);
  #            Prkcd: distingue FO da HO in quasi tutti i sistemi
  #
  #   TC_HO  = Higher-order (LP, PO, MD, VA/VL motori, Pf/CM intralaminar)
  #            → profilo "secondario": Necab1+, Cbln1+, Cbln2+
  #            Necab1 (N-terminal EF-hand calcium binding 1): altamente
  #            selettivo per HO in LP, PO, MD (Phillips 2019; dLGN eLife 2021
  #            — shell vs core); Cbln1/2 forti in LP/PO/MD adulto
  #
  #   TC_MTX = Matrix/intralaminar (CM, Pf, CL, Re, Rh — diffuse projection)
  #            → profilo "terziario": Cplx3+, Htr2c+, Cd44+
  #            Cplx3 (complexin 3): espresso selettivamente nelle cellule
  #            matrix del talamo (Deschênes, Bhatt, Jones revisioni);
  #            Htr2c: recettore 5-HT2C, altamente arricchito nelle cellule
  #            intralaminar/matrix (rilevante anche per LSD)
  #
  # Nota: Slc17a6 (VGluT2) = pan-marker TC glutamatergici → utile come
  # controllo positivo per tutti i sottotipi TC ma non li discrimina.
  # Rorb/Hcn1/Cacna1g: ampiamente distribuiti, scarso potere discriminante
  # tra FO e HO nel topo adulto → rimossi come marcatori principali.
  
  # Pan-TC glutamatergico (controllo positivo, striscia separata)
  Slc17a6 = "TC_pan",   Ntng1   = "TC_pan",   Syt2    = "TC_pan",
  
  # TC First-Order relay sensoriali (VPM/VPL/dLGN/MGv)
  Tnnt1   = "TC_FO",    Prkcd   = "TC_FO",    Scnn1a  = "TC_FO",
  Kcnab2  = "TC_FO",    Grik4   = "TC_FO",    Slc6a1  = "TC_FO",
  
  # TC Higher-Order (LP/PO/MD/VA/VL motori)
  Necab1  = "TC_HO",    Cbln1   = "TC_HO",    Cbln2   = "TC_HO",
  Calb1   = "TC_HO",    Pou4f1  = "TC_HO",    Nnat    = "TC_HO",
  
  # TC Matrix / intralaminar (CM/Pf/CL/Re/Rh)
  Cplx3   = "TC_MTX",   Htr2c   = "TC_MTX",   Cd44    = "TC_MTX",
  Rxfp3   = "TC_MTX",   Penk    = "TC_MTX",
  
  # -----------------------------------------------------------------------
  # TRN (Thalamic Reticular Nucleus) — GABAergici
  # -----------------------------------------------------------------------
  # Pvalb è il marcatore più affidabile e specifico per TRN nel topo adulto
  # (Pinault 2004; Halassa 2011; confermato da ABCA WMB-10Xv3).
  # Gad1/Gad2: pan-GABAergici (separano TRN da tutti i TC glutamatergici).
  # Etv1 (ETS variant 1): espresso selettivamente nel TRN (Allen ISH).
  # Sst: usato con cautela — presente in TRN ma anche in alcuni TC motori.
  Gad1    = "TRN",      Gad2    = "TRN",
  Pvalb   = "TRN",      Etv1    = "TRN",      Cacna1i = "TRN"
)

gene_thal_2<-c(# Aggiunte a gene_ct_map_thal:
  
  # Identità MTX costitutiva
  Nrn1    = "TC_MTX",   Synpr   = "TC_MTX",
  Reln    = "TC_MTX",   Calb2   = "TC_MTX",
  
  # Pf-specifici
  Tacr1   = "TC_Pf",    Rxfp3   = "TC_Pf",    Nr4a2   = "TC_Pf",
  
  # CL-specifici  
  # (Calb1 e Pou4f1 già presenti in TC_HO — valuta se duplicarli qui
  #  dopo aver visto il DotPlot; non duplicare a priori)
  
  # Attività-dipendenti (late IEG / secondary response)
  # Late IEG / secondary response — picco 60-120min, AP-1 dipendenti
  Penk    = "Activity",   # già hai — sposta qui se l'ipotesi è confermata
  Bdnf    = "Activity",   # BDNF: indotto da attività neuronale, picco ~90min
  Vgf     = "Activity",   # VGF nerve growth factor inducible — fortemente AP-1 responsivo
  Nptx2   = "Activity",   # neuronal pentraxin 2 — indotto da attività sinaptica
  Nrg1    = "Activity",   # neuregulin 1 — modulato da attività in talamostriatal
  
  # Primary IEG (picco ~30-60min) — per confronto temporale
  Fos     = "IEG",        # già hai
  Arc     = "IEG",        # già hai
  Egr1    = "IEG",
  Nr4a1   = "IEG",        # Nur77 — IEG canonico, utile per timeline
  Fosb    = "IEG")        # più stabile di Fos, picco leggermente più tardivo)

# --- 05A.1 — DotPlot annotato con strip cell type ----------------------------

make_annotated_dotplot <- function(obj, gene_ct_map, area_label) {

  genes_ordered <- names(gene_ct_map)[names(gene_ct_map) %in% rownames(obj)]
  ct_per_gene   <- gene_ct_map[genes_ordered]

  if (length(genes_ordered) == 0) {
    message(sprintf("  %s: nessun gene della mappa trovato", area_label))
    return(NULL)
  }

  dp_tmp        <- DotPlot(obj, features = genes_ordered, cluster.idents = TRUE)
  cluster_order <- levels(dp_tmp$data$id)

  p_dot <- DotPlot(
    obj,
    features       = genes_ordered,
    cols           = c("lightgrey", "#C0392B"),
    dot.scale      = 5.5,
    cluster.idents = TRUE
  ) +
    scale_y_discrete(limits = rev(cluster_order)) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y     = element_text(size = 8),
      axis.title      = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 7),
      plot.margin     = margin(t = 0, r = 5, b = 5, l = 5)
    )

  ct_unique  <- unique(ct_per_gene)
  n_ct       <- length(ct_unique)
  strip_pal  <- setNames(colorRampPalette(pal_ct)(n_ct), ct_unique)

  strip_df <- data.frame(
    gene     = factor(genes_ordered, levels = genes_ordered),
    ct_group = factor(ct_per_gene, levels = ct_unique),
    x        = seq_along(genes_ordered),
    stringsAsFactors = FALSE
  )

  block_df <- strip_df %>%
    group_by(ct_group) %>%
    summarise(
      xmin = min(x) - 0.5,
      xmax = max(x) + 0.5,
      xmid = mean(x),
      n    = n(),
      .groups = "drop"
    )

  p_strip <- ggplot(block_df) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = ct_group),
      color = "white", linewidth = 0.5
    ) +
    geom_text(
      data = dplyr::filter(block_df, n >= 2),
      aes(x = xmid, y = 0.5, label = ct_group),
      size = 2.0, fontface = "bold", color = "white", vjust = 0.5
    ) +
    scale_fill_manual(values = strip_pal, guide = "none") +
    scale_x_continuous(
      limits = c(0.5, length(genes_ordered) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste(area_label, "\u2014 Marcatori per cluster (pre-annotazione)")) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 10, face = "bold",
                                     hjust = 0,
                                     margin = margin(b = 3)),
      plot.margin     = margin(t = 6, r = 5, b = 0, l = 5)
    )

  p_strip / p_dot + plot_layout(heights = c(0.07, 1))
}

p_dot_thal <- make_annotated_dotplot(THAL, gene_ct_map_thal, "THAL")
savepng_thal(p_dot_thal, "05_DotPlot_prefilt_THAL.png",
             width = 4500, height = 1800)
message("  Salvato: 05_DotPlot_prefilt_THAL.png")

# --- 05A.2 — DotPlot QC per cluster ------------------------------------------

meta_qc  <- THAL@meta.data
meta_qc$cluster <- as.character(Idents(THAL))

qc_summary <- meta_qc %>%
  group_by(cluster) %>%
  summarise(
    mean_nCount   = mean(nCount_RNA,   na.rm = TRUE),
    mean_nFeature = mean(nFeature_RNA, na.rm = TRUE),
    mean_mt       = mean(percent.mt,   na.rm = TRUE),
    n_cells       = n(),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(cluster,
                          levels = as.character(sort(as.numeric(unique(cluster))))))

qc_long <- qc_summary %>%
  tidyr::pivot_longer(
    cols      = c(mean_nCount, mean_nFeature, mean_mt),
    names_to  = "metric",
    values_to = "mean_val"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("mean_nCount", "mean_nFeature", "mean_mt"),
                    labels = c("nCount_RNA",  "nFeature_RNA",  "percent.mt"))
  ) %>%
  group_by(metric) %>%
  mutate(z = scale(mean_val)[, 1]) %>%
  ungroup()

thresholds <- meta_qc %>%
  summarise(
    nCount_RNA   = median(nCount_RNA,   na.rm = TRUE) + 2 * mad(nCount_RNA,   na.rm = TRUE),
    nFeature_RNA = median(nFeature_RNA, na.rm = TRUE) + 2 * mad(nFeature_RNA, na.rm = TRUE),
    percent.mt   = 5
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "threshold") %>%
  mutate(metric = factor(metric, levels = c("nCount_RNA","nFeature_RNA","percent.mt")))

p_qc <- ggplot(qc_long, aes(x = mean_val, y = cluster)) +
  geom_segment(
    aes(x = 0, xend = mean_val, y = cluster, yend = cluster),
    color = "grey88", linewidth = 0.4
  ) +
  geom_point(aes(color = z, size = n_cells), alpha = 0.92) +
  geom_vline(
    data     = thresholds,
    aes(xintercept = threshold),
    linetype = "dashed", color = "firebrick", linewidth = 0.65,
    na.rm    = TRUE
  ) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_color_gradient2(
    low = "#3498DB", mid = "grey88", high = "#E74C3C",
    midpoint = 0, name = "z-score"
  ) +
  scale_size_continuous(name = "N cellule", range = c(2, 9),
                        breaks = scales::pretty_breaks(4)) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  theme(
    strip.text         = element_text(face = "bold", size = 10),
    strip.background   = element_rect(fill = "grey92", color = NA),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 8, angle = 30, hjust = 1),
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.key.size    = unit(0.45, "cm"),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 8)
  ) +
  labs(title = "THAL \u2014 Metriche QC per cluster")

n_cl <- length(levels(qc_long$cluster))
savepng_thal(p_qc, "05_QC_dotplot_THAL.png",
             width  = 3500,
             height = max(1200, n_cl * 75 + 400))
message("  Salvato: 05_QC_dotplot_THAL.png")

message("
>>> PARTE A completata.
    Apri 05_DotPlot_prefilt_THAL.png e 05_QC_dotplot_THAL.png,
    poi compila la mappa di rinomina in PARTE B e ri-esegui.
")

# ==============================================================================
# PARTE B — RIMOZIONE CLUSTER + RINOMINA → cell_type
# ==============================================================================

message("=== 05B. Rinomina cluster THAL → cell_type ===")

# ------------------------------------------------------------------------------
# <<< ISTRUZIONI PER COMPILARE LA MAPPA >>>
#
# Apri 05_DotPlot_prefilt_THAL.png e identifica ogni cluster:
#
# Esempi tipici THAL:
#   "0"  = "TC_relay"    (Slc17a6+, Rorb+, Grik4+, Hcn1+)
#   "1"  = "TC_relay_2"  (Slc17a6+, Rorb+, Prkcg+)
#   "2"  = "TC_HO"       (Slc17a6+, Calb1+, Pou4f1+)
#   "3"  = "TC_HO_2"     (Slc17a6+, Calb2+, Nrgn+)
#   "4"  = "TRN"         (Gad1+, Gad2+, Pvalb+, Sst+)
#   "5"  = "Astro"       (Aqp4+, Aldoc+, Mlc1+)
#   "6"  = "Oligo"       (Mog+, Mbp+, Plp1+)
#   "7"  = "NFOL"        (Bcas1+, Enpp6+)
#   "8"  = "OPC"         (Pdgfra+, Cspg4+)
#   "9"  = "Micro"       (Cx3cr1+, P2ry12+, C1qa+)
#   "10" = "Endo"        (Cldn5+, Pecam1+, Flt1+)
#   "11" = "Peri"        (Rgs5+, Pdgfrb+)
#
# Cluster da scartare: bassa qualità (nFeature basso, mt alto), doublet evidenti
# ------------------------------------------------------------------------------

CLUSTERS_REMOVE_THAL <- c()   # <<< es. c("12", "13")

# <<< DA COMPILARE dopo aver visto i plot >>>
thal_rename <- NULL

# Esempio da decommentare e adattare:
thal_rename <- c(
  "1"  = "NA",
  "2"  = "TC_FO",
  "3"  = "TC_MTX",
  "4"  = "Oligo",
  "5"  = "NA2",
  "6"  = "Astro",
  "7"  = "TRN2",
  "8"  = "TC_FO2",
  "9"  = "Endo",
  "10" = "TC_FO3",
  "11" = "TRN3",
  "12" = "Micro",
  "13" = "NA3",
  "14" = "OPC",
  "15" = "Astro2",
  "16" = "TRN3",
  "17" = "TC_MTX",
  "18" = "Astro3",
  "19" = "TRN"
)

# ---- 1. Rimozione cluster espliciti -----------------------------------------
n_before <- ncol(THAL)
if (length(CLUSTERS_REMOVE_THAL) > 0) {
  bad        <- as.character(CLUSTERS_REMOVE_THAL)
  cells_keep <- !(as.character(Idents(THAL)) %in% bad)
  THAL       <- THAL[, cells_keep]
  message(sprintf("  THAL: rimossi cluster %s → -%d cellule",
                  paste(bad, collapse = ","),
                  n_before - ncol(THAL)))
}

# ---- 2. Rinomina cluster → label cell_type ----------------------------------
if (!is.null(thal_rename)) {
  missing_keys <- setdiff(names(thal_rename), as.character(levels(Idents(THAL))))
  if (length(missing_keys) > 0)
    warning(sprintf("  THAL: chiavi non trovate nella mappa: %s",
                    paste(missing_keys, collapse = ", ")))

  THAL <- RenameIdents(THAL, thal_rename)

  if ("Low_quality" %in% levels(Idents(THAL))) {
    n_lq <- sum(Idents(THAL) == "Low_quality")
    THAL <- THAL[, Idents(THAL) != "Low_quality"]
    message(sprintf("  THAL: rimossi %d cluster 'Low_quality'", n_lq))
  }

  message(sprintf("  THAL: cluster rinominati in %d tipi cellulari",
                  length(unique(as.character(Idents(THAL))))))
} else {
  message("  THAL: nessuna mappa fornita — cell_type = cluster numerico (provvisorio)")
}

# ---- 3. Salva cell_type come colonna ----------------------------------------
THAL$cell_type <- as.character(Idents(THAL))

# ==============================================================================
# PARTE C — PLOT POST-ANNOTAZIONE
# ==============================================================================

message("\n=== 05C. Plot post-annotazione THAL ===")

# --- 05C.0 — AddModuleScore per cell type ------------------------------------

ct_markers_thal <- split(names(gene_ct_map_thal), unname(gene_ct_map_thal))

ct_freq  <- sort(table(THAL$cell_type), decreasing = TRUE)
ct_order <- intersect(names(ct_freq), names(ct_markers_thal))
ct_order <- ct_order[vapply(ct_order, function(ct) {
  genes <- filter_genes(ct_markers_thal[[ct]], THAL)
  length(genes) >= 2
}, FUN.VALUE = logical(1))]

if (length(ct_order) > 0) {
  score_cols <- c()
  for (ct in ct_order) {
    genes    <- filter_genes(ct_markers_thal[[ct]], THAL)
    score_nm <- paste0("score_", gsub("[^A-Za-z0-9]", "_", ct))
    score_cols <- c(score_cols, score_nm)

    THAL <- AddModuleScore(THAL, features = list(genes),
                           name = score_nm, seed = 42)
    colnames(THAL@meta.data)[
      colnames(THAL@meta.data) == paste0(score_nm, "1")
    ] <- score_nm
  }

  n_cols    <- 4
  plot_list <- lapply(seq_along(ct_order), function(i) {
    ct       <- ct_order[i]
    score_nm <- score_cols[i]
    n_cells  <- sum(THAL$cell_type == ct, na.rm = TRUE)

    FeaturePlot(THAL, features = score_nm, reduction = "umap",
                pt.size = 0.2, combine = FALSE)[[1]] +
      scale_color_gradient2(low = "blue", mid = "white", high = "#C0392B",
                            midpoint = 0, name = "score") +
      ggtitle(sprintf("%s  (n=%d)", ct, n_cells)) +
      theme_void(base_size = 8) +
      theme(
        plot.title      = element_text(size = 7.5, face = "bold", hjust = 0.5),
        legend.key.size = unit(0.3, "cm"),
        legend.text     = element_text(size = 5),
        legend.title    = element_text(size = 6),
        plot.margin     = margin(3, 3, 3, 3)
      )
  })

  n_rows <- ceiling(length(plot_list) / n_cols)
  savepng_thal(
    wrap_plots(plot_list, ncol = n_cols) +
      plot_annotation(
        title    = "THAL \u2014 Score marcatori per cell type",
        subtitle = "Colore = AddModuleScore(marcatori specifici) | Scala divergente: grigio=0, rosso=alto",
        theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                         plot.subtitle = element_text(size = 8, color = "grey40"))
      ),
    "05_ModuleScore_per_CT_THAL.png",
    width  = n_cols * 1800,
    height = n_rows * 1800 + 400
  )
  message(sprintf("  Salvato: 05_ModuleScore_per_CT_THAL.png (%d pannelli)",
                  length(plot_list)))
}

# --- 05C.1 — FeaturePlot geni chiave talamici --------------------------------
# Pan-neuronali + identità TC + TRN + gliali

key_genes_thal <- filter_genes(
  c("Slc17a6", "Gad1", "Gad2", "Pvalb", "Sst",
    "Rorb", "Hcn1", "Grik4", "Calb1", "Calb2",
    "Lhx2", "Gbx2", "Nr4a2",
    "Aqp4", "Mog", "Pdgfra", "Cx3cr1", "Cldn5",
    "Fos", "Arc"),
  THAL
)

umap_coords <- as.data.frame(Embeddings(THAL, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords$cell_type <- THAL$cell_type

ct_centroids <- umap_coords %>%
  group_by(cell_type) %>%
  summarise(cx = mean(UMAP_1), cy = mean(UMAP_2), .groups = "drop")

avg_expr <- AverageExpression(
  THAL, features = key_genes_thal, assays = "SCT",
  layer = "data", group.by = "cell_type", verbose = FALSE
)$SCT

gene_labels <- lapply(key_genes_thal, function(gene) {
  if (!gene %in% rownames(avg_expr)) return(NULL)
  top_ct   <- names(which.max(avg_expr[gene, ]))
  centroid <- ct_centroids[ct_centroids$cell_type == top_ct, ]
  if (nrow(centroid) == 0) return(NULL)
  data.frame(gene = gene, ct_label = top_ct,
             cx = centroid$cx, cy = centroid$cy)
})
gene_labels <- do.call(rbind, Filter(Negate(is.null), gene_labels))

fp_list <- lapply(key_genes_thal, function(gene) {
  p <- FeaturePlot(THAL, features = gene, reduction = "umap",
                   pt.size = 0.2, combine = FALSE,
                   cols = c("grey90", viridis::inferno(100)))[[1]] +
    theme_bw(base_size = 8) +
    theme(
      plot.title      = element_text(size = 8, face = "bold"),
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 6),
      plot.margin     = margin(2, 2, 2, 2)
    )
  
  lbl_row <- gene_labels[!is.na(gene_labels$gene) & gene_labels$gene == gene, ]
  if (!is.null(lbl_row) && nrow(lbl_row) > 0) {
    p <- p +
      annotate("label",
               x = lbl_row$cx, y = lbl_row$cy,
               label = lbl_row$ct_label,
               size = 2.2, fill = "white", color = "black",
               alpha = 0.85, label.size = 0.3,
               label.padding = unit(0.15, "lines"))
  }
  p
})

n_cols_fp <- 4
savepng_thal(wrap_plots(fp_list, ncol = n_cols_fp),
             "05_FeaturePlot_key_genes_THAL.png",
             width  = 6000,
             height = ceiling(length(fp_list) / n_cols_fp) * 1600)
message("  Salvato: 05_FeaturePlot_key_genes_THAL.png")

# --- 05C.2 — Heatmap cluster × MapMyCells subclass ---------------------------

if ("mmc_subclass" %in% colnames(THAL@meta.data)) {
  ct_tab  <- table(THAL$cell_type, THAL$mmc_subclass)
  ct_prop <- prop.table(ct_tab, margin = 1)

  row_ct    <- rownames(ct_prop)
  row_order <- order(row_ct, rownames(ct_prop))
  ct_prop_ordered <- ct_prop[row_order, ]

  col_keep     <- apply(ct_prop_ordered, 2, max) >= 0.01
  ct_prop_plot <- ct_prop_ordered[, col_keep, drop = FALSE]

  png(file.path(thal_dir, "05_cluster_vs_MMC_THAL.png"),
      width  = max(3500, ncol(ct_prop_plot) * 65 + 700),
      height = max(2000, nrow(ct_prop_plot) * 60 + 400),
      res    = 300)
  pheatmap(
    ct_prop_plot,
    color           = colorRampPalette(c("white", "#1A5276", "#2ECC71"))(80),
    display_numbers = FALSE,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    fontsize        = 7, fontsize_row = 7, fontsize_col = 7,
    border_color    = NA,
    main            = "Cluster vs MapMyCells subclass (cell_type) \u2014 THAL"
  )
  dev.off()
  message("  Salvato: 05_cluster_vs_MMC_THAL.png")
}

# --- 05C.3 — UMAP annotato + trattamento + composizione ----------------------

ct_levels <- sort(unique(THAL$cell_type))
n_ct      <- length(ct_levels)
ct_colors <- setNames(colorRampPalette(pal_ct)(n_ct), ct_levels)

p_ct <- DimPlot(
  THAL, group.by = "cell_type", reduction = "umap",
  label = TRUE, repel = TRUE, label.size = 3,
  pt.size = 0.3, cols = ct_colors
) +
  ggtitle("THAL \u2014 Tipo cellulare") +
  theme(legend.text = element_text(size = 8))

p_trt <- DimPlot(
  THAL, group.by = "treatment", reduction = "umap",
  pt.size = 0.3, cols = pal_trt
) + ggtitle("THAL \u2014 Trattamento")

p_time <- DimPlot(
  THAL, group.by = "time", reduction = "umap",
  pt.size = 0.3, cols = pal_time
) + ggtitle("THAL \u2014 Timepoint")

p_samp <- DimPlot(
  THAL, group.by = "orig.ident", reduction = "umap",
  pt.size = 0.3
) + ggtitle("THAL \u2014 Campione")

savepng_thal(
  (p_ct | p_trt) / (p_time | p_samp) +
    plot_annotation(title = "Annotazione finale \u2014 THAL"),
  "05_UMAP_annotated_THAL.png",
  width = 6000, height = 5000
)

# --- 05C.4 — DotPlot marcatori per cell_type annotato ------------------------

markers_thal_use <- filter_genes(
  unique(names(gene_ct_map_thal)),
  THAL
)
markers_thal_use <- markers_thal_use[seq_len(min(60, length(markers_thal_use)))]

p_dot_ann <- DotPlot(
  THAL, features = markers_thal_use,
  group.by = "cell_type", cluster.idents = FALSE
) +
  RotatedAxis() +
  ggtitle("THAL \u2014 Marcatori per tipo cellulare")
savepng_thal(p_dot_ann, "05_DotPlot_annotated_THAL.png",
             width = 7000, height = max(2500, n_ct * 160 + 600))

# --- 05C.5 — FeaturePlot IEG + Fos ------------------------------------------

feat_thal <- filter_genes(c("Fos", "Arc", "Egr1", "Nr4a1", "Fosb", "Junb"), THAL)
if (length(feat_thal) > 0) {
  p_ieg <- FeaturePlot(
    THAL, features = feat_thal,
    reduction = "umap", ncol = 3, pt.size = 0.2
  ) & scale_color_viridis_c(option = "inferno") &
    theme(legend.position = "right", plot.title = element_text(size = 9))
  savepng_thal(p_ieg, "05_FeaturePlot_IEG_THAL.png",
               width  = 6000,
               height = ceiling(length(feat_thal) / 3) * 2000)
}

# --- 05C.6 — Heatmap proporzione cell_type × campione -----------------------

ct_prop_samp <- prop.table(table(THAL$cell_type, THAL$orig.ident), margin = 2)
png(file.path(thal_dir, "05_CellType_composition_THAL.png"),
    width  = 3000,
    height = max(2000, nrow(ct_prop_samp) * 80 + 600),
    res    = 300)
pheatmap(
  ct_prop_samp,
  color           = colorRampPalette(c("white", "#27AE60"))(50),
  display_numbers = TRUE, number_format = "%.2f",
  fontsize        = 8, cluster_cols = FALSE,
  main            = "Proporzione CT per campione \u2014 THAL"
)
dev.off()
message("  Salvato: 05_CellType_composition_THAL.png")

# ==============================================================================
# SALVATAGGIO
# ==============================================================================

message("\n=== 05. Salvataggio ===")
save(THAL, file = file.path(thal_dir, "THAL_celltype.RData"))
message("Salvato: THAL/THAL_celltype.RData")
