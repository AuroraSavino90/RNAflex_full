################################################################################
##  05_cell_type_annotation.R                                                  ##
##                                                                             ##
##  Input:  CTX_annotated.RData, THAL_annotated.RData  (da 04_mapmycells)    ##
##          Se MapMyCells non è stato eseguito:                               ##
##          CTX_initial.RData,   THAL_initial.RData   (da 0_QC)              ##
##                                                                             ##
##  Output: CTX_clean.RData, THAL_clean.RData                                 ##
##          (attesi da 06_metacells, 07_deseq2, 10_subclustering, 11_htr2a)  ##
##                                                                             ##
##  Flusso:                                                                    ##
##    PARTE A — ispezione visiva (eseguire sempre)                            ##
##      1. DotPlot marcatori per cluster                                      ##
##      2. FeaturePlot QC + geni chiave                                       ##
##      3. Heatmap cluster × MapMyCells subclass (se disponibile)            ##
##    PARTE B — compilare le mappe di rinomina, poi eseguire                  ##
##      4. Rimozione cluster di bassa qualità / contaminanti                  ##
##      5. RenameIdents → colonna cell_type                                   ##
##      6. Plot post-annotazione                                               ##
##      7. Salvataggio CTX_clean / THAL_clean                                 ##
################################################################################

source("code/00_config.R")
outdir<-"results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

# Prova prima la versione annotata (MapMyCells), poi quella iniziale

  f_ann  <- file.path(outdir, "CTX_annotated.RData")
  load(f_ann, envir = .GlobalEnv)
 
# ==============================================================================
# PARTE A — ISPEZIONE VISIVA DEI CLUSTER
# ==============================================================================

message("\n=== 05A. Ispezione visiva cluster ===")

# Stampa i livelli di cluster correnti per compilare le mappe sotto
message("\nLivelli cluster CTX:  ", paste(levels(Idents(CTX)),  collapse = " | "))

  # --- 05A.1 — DotPlot annotato: marcatori per cluster con strip cell type ----
  #
  # I geni sono ordinati per cell type di appartenenza (blocchi contigui),
  # così il plot mostra blocchi diagonali puliti quando i cluster
  # corrispondono alle identità attese.
  # Sopra le colonne viene disegnata una strip colorata che indica
  # a quale cell type appartiene ogni gruppo di marcatori.
  
  # Mappa gene → cell type — CTX
  # Ordine: gliali → endotelio → Ex (strato) → Inh (sottotipo)
  gene_ct_map_ctx <- c(
    # Astrociti
    Mlc1    = "Astro",     Gjb6     = "Astro",    Aqp4     = "Astro",
    Acsbg1  = "Astro",     Aldoc    = "Astro",    S100b    = "Astro",
    # Oligo maturi
    Ermn    = "Oligo",     Opalin   = "Oligo",    Mog      = "Oligo",
    Aspa    = "Oligo",     Mobp     = "Oligo",    Mbp      = "Oligo",
    Plp1    = "Oligo",     Ptgds    = "Oligo",
    # NFOL (Newly Formed OL)
    Bcas1   = "NFOL",      Enpp6    = "NFOL",     Neu4     = "NFOL",     Bmp4    = "NFOL",
    # OPC
    Pdgfra  = "OPC",       Cspg4    = "OPC",      Ptprz1   = "OPC",
    Cacng4  = "OPC",       Tmem100  = "OPC",      Matn4    = "OPC",
    # Microglia
    C1qa    = "Micro",     C1qb     = "Micro",    C1qc     = "Micro",
    Ctss    = "Micro",     P2ry12   = "Micro",    Cx3cr1   = "Micro",    Tmem119 = "Micro",
    # Endotelio
    Flt1    = "Endo",      Cldn5    = "Endo",     Tm4sf1   = "Endo",
    Pecam1  = "Endo",      Rgs5     = "Endo",
    # Eccitatori pan
    Snap25  = "Ex",        Neurod6  = "Ex",       Slc17a7  = "Ex",
    # Ex L2/3
    Rasgrf2 = "Ex_L23",    Cux1     = "Ex_L23",   Cux2     = "Ex_L23",
    Otof    = "Ex_L23",    Satb2    = "Ex_L23",
    # Ex L4
    Rorb    = "Ex_L4",     Rspo1    = "Ex_L4",    Il1rapl2 = "Ex_L4",
    # Ex L5 IT
    Galnt14 = "Ex_L5IT",   Sulf1    = "Ex_L5IT",  Bcl6     = "Ex_L5IT",
    # Ex L5 PT
    Fezf2   = "Ex_L5PT",   Bcl11b   = "Ex_L5PT",  Ldb2     = "Ex_L5PT",
    Slc17a8 = "Ex_L5PT",   Chrna6   = "Ex_L5PT",
    # Ex L6 CT
    Foxp2   = "Ex_L6CT",   Syt6     = "Ex_L6CT",  Hs3st4   = "Ex_L6CT",
    Tshz2   = "Ex_L6CT",   Nr4a2    = "Ex_L6CT",  Ctgf     = "Ex_L6CT",
    # Ex L6b
    Ostn    = "Ex_L6b",    Hpcal1   = "Ex_L6b",
    # Inibitori pan
    Gad1    = "Inh",       Gad2     = "Inh",      Slc32a1  = "Inh",
    # Inh PV
    Pvalb   = "Inh_PV",    Kcnc1    = "Inh_PV",   Syt2     = "Inh_PV",
    # Inh SST
    Sst     = "Inh_SST",   Chodl    = "Inh_SST",  Npy      = "Inh_SST",  Tac2  = "Inh_SST",
    # Inh VIP
    Vip     = "Inh_VIP",   Calb2    = "Inh_VIP",  Cck      = "Inh_VIP",  Tac1  = "Inh_VIP",
    # Inh Lamp5
    Lamp5   = "Inh_Lamp5", Ndnf     = "Inh_Lamp5", Lhx6   = "Inh_Lamp5",
    # Inh Sncg / Reln
    Sncg    = "Inh_Sncg",  Cnr1     = "Inh_Sncg",  Reln   = "Inh_Sncg"
  )
  
  # Mappa gene → cell type — THAL
  # Ordine: gliali → endotelio → TC (sottotipo) → TRN
  gene_ct_map_thal <- c(
    # Astrociti
    Mlc1    = "Astro",     Gjb6    = "Astro",    Aqp4    = "Astro",
    Acsbg1  = "Astro",     Aldoc   = "Astro",    S100b   = "Astro",
    # Oligo
    Ermn    = "Oligo",     Mog     = "Oligo",    Mobp    = "Oligo",
    Plp1    = "Oligo",     Mbp     = "Oligo",    Ptgds   = "Oligo",
    # NFOL
    Bcas1   = "NFOL",      Enpp6   = "NFOL",     Neu4    = "NFOL",     Bmp4   = "NFOL",
    # OPC
    Pdgfra  = "OPC",       Ptprz1  = "OPC",      Cspg4   = "OPC",
    # Microglia
    C1qa    = "Micro",     C1qb    = "Micro",    C1qc    = "Micro",
    Ctss    = "Micro",     P2ry12  = "Micro",    Cx3cr1  = "Micro",   Tmem119 = "Micro",
    # Endotelio
    Flt1    = "Endo",      Cldn5   = "Endo",     Pecam1  = "Endo",    Rgs5    = "Endo",
    # TC relay sensoriali (VPM/VPL/dLGN/MGv)
    Slc17a6 = "TC_relay",  Rorb    = "TC_relay", Grik4   = "TC_relay",
    Hcn1    = "TC_relay",  Cacna1g = "TC_relay", Kcnab2  = "TC_relay",
    # TC higher-order (PO/LP/MD/CM)
    Tbr1    = "TC_HO",     Pou4f1  = "TC_HO",   Calb1   = "TC_HO",
    Calb2   = "TC_HO",     Nrgn    = "TC_HO",   Prkcg   = "TC_HO",   Grid2   = "TC_HO",
    # Identità TH generale
    Nr4a2   = "TH_id",     Lhx2    = "TH_id",   Gbx2    = "TH_id",
    # TRN
    Gad1    = "TRN",       Gad2    = "TRN",     Pvalb   = "TRN",     Sst     = "TRN"
  )
  
  # Funzione principale: costruisce DotPlot + strip annotazione
  make_annotated_dotplot <- function(obj, gene_ct_map, area_label) {
    
    # 1. Filtra geni presenti nel dataset, mantieni ordine della mappa
    genes_ordered <- names(gene_ct_map)[names(gene_ct_map) %in% rownames(obj)]
    ct_per_gene   <- gene_ct_map[genes_ordered]
    
    if (length(genes_ordered) == 0) {
      message(sprintf("  %s: nessun gene della mappa trovato", area_label))
      return(NULL)
    }
    
    # 2. Estrai l'ordine dei cluster dal dendrogramma (cluster.idents)
    dp_tmp        <- DotPlot(obj, features = genes_ordered, cluster.idents = TRUE)
    cluster_order <- levels(dp_tmp$data$id)
    
    # 3. DotPlot finale — senza titolo (il titolo va nella strip)
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
    
    # 4. Strip colorata: un rettangolo per gruppo CT, larghezza = n geni del gruppo
    ct_unique <- unique(ct_per_gene)   # ordine di comparsa nella mappa
    n_ct      <- length(ct_unique)
    
    # Palette: interpola pal_ct se serve (gestisce n_ct > lunghezza di pal_ct)
    strip_pal <- setNames(colorRampPalette(pal_ct)(n_ct), ct_unique)
    
    strip_df <- data.frame(
      gene     = factor(genes_ordered, levels = genes_ordered),
      ct_group = factor(ct_per_gene, levels = ct_unique),
      x        = seq_along(genes_ordered),
      stringsAsFactors = FALSE
    )
    
    block_df <- strip_df %>%
      group_by(ct_group) %>%
      summarise(
        xmin  = min(x) - 0.5,
        xmax  = max(x) + 0.5,
        xmid  = mean(x),
        n     = n(),         # n geni nel blocco (serve per nascondere label blocchi piccoli)
        .groups = "drop"
      )
    
    p_strip <- ggplot(block_df) +
      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = ct_group),
        color = "white", linewidth = 0.5
      ) +
      # Etichetta solo sui blocchi >= 2 geni (evita sovrapposizioni)
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
    
    # 5. Combina strip (sopra, sottile) + DotPlot con patchwork
    p_combined <- p_strip / p_dot +
      plot_layout(heights = c(0.07, 1))
    
    p_combined
  }
  
  # Costruisci e salva
  p_dot_ctx  <- make_annotated_dotplot(CTX,  gene_ct_map_ctx,  "CTX")
  
    savepng(p_dot_ctx,  "05_DotPlot_prefilt_CTX.png",
            width = 4000, height = 1500)
   message("  Salvati: 05_DotPlot_prefilt_*.png")
  
    # --- 05A.2 — DotPlot QC per cluster: nCount, nFeature, percent.mt ----------
    #
    # Ogni punto = un cluster Seurat.
    # Colore   = valore medio della metrica (z-score per confronto cross-metrica)
    # Blu = basso, Rosso = alto rispetto alla media del dataset
    # Dimensione = numero di cellule nel cluster
    # I cluster con nCount/nFeature molto bassi o percent.mt alta
    # sono candidati alla rimozione in PARTE B.
    #
    # Un unico file PNG con tre pannelli affiancati (CTX e THAL separati).
    
    
      meta <- CTX@meta.data
      meta$cluster <- as.character(Idents(CTX))
      
      # Sommario per cluster
      qc_summary <- meta %>%
        group_by(cluster) %>%
        summarise(
          mean_nCount   = mean(nCount_RNA,   na.rm = TRUE),
          mean_nFeature = mean(nFeature_RNA, na.rm = TRUE),
          mean_mt       = mean(percent.mt,   na.rm = TRUE),
          n_cells       = n(),
          .groups = "drop"
        ) %>%
        # Ordina numericamente per un asse y pulito
        mutate(cluster = factor(cluster,
                                levels = as.character(
                                  sort(as.numeric(unique(cluster))))))
      
      # Formato long + z-score per metrica (colore relativo al range del dataset)
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
      
      # Soglie globali mediana ± 2 MAD (per metrica) → linee verticali di allerta
      thresholds <- meta %>%
        summarise(
          nCount_RNA   = median(nCount_RNA,   na.rm = TRUE) + 2 * mad(nCount_RNA,   na.rm = TRUE),
          nFeature_RNA = median(nFeature_RNA, na.rm = TRUE) + 2 * mad(nFeature_RNA, na.rm = TRUE),
          percent.mt   = 5   # soglia fissa percent.mt — modificare se necessario
        ) %>%
        tidyr::pivot_longer(everything(),
                            names_to = "metric", values_to = "threshold") %>%
        mutate(metric = factor(metric,
                               levels = c("nCount_RNA","nFeature_RNA","percent.mt")))
      
      p_qc <- ggplot(qc_long, aes(x = mean_val, y = cluster)) +
        # Linea di riferimento per ogni cluster (guida l'occhio)
        geom_segment(
          aes(x = 0, xend = mean_val, y = cluster, yend = cluster),
          color = "grey88", linewidth = 0.4
        ) +
        geom_point(aes(color = z, size = n_cells), alpha = 0.92) +
        # Linea di allerta verticale per metrica
        geom_vline(
          data     = thresholds,
          aes(xintercept = threshold),
          linetype = "dashed", color = "firebrick", linewidth = 0.65,
          na.rm    = TRUE
        ) +
        facet_wrap(~ metric, scales = "free_x", nrow = 1) +
        scale_color_gradient2(
          low      = "#3498DB",
          mid      = "grey88",
          high     = "#E74C3C",
          midpoint = 0,
          name     = "z-score"
        ) +
        scale_size_continuous(
          name   = "N cellule",
          range  = c(2, 9),
          breaks = scales::pretty_breaks(4)
        ) +
        scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
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
        labs(title = paste(obj_nm, "\u2014 Metriche QC per cluster"))
      
      n_cl <- length(levels(qc_long$cluster))
      savepng(p_qc,
              paste0("05_QC_dotplot_", obj_nm, ".png"),
              width  = 3500,
              height = max(1200, n_cl * 75 + 400))
    
    message("  Salvati: 05_QC_dotplot_CTX.png, 05_QC_dotplot_THAL.png")


message("
>>> PARTI A completata.
    Apri i file 05_DotPlot_prefilt_*.png e 05_cluster_vs_MMC_*.png,
    poi compila le mappe di rinomina in PARTE B e ri-esegui lo script.
")

# ==============================================================================
# PARTE B — RIMOZIONE CLUSTER + RINOMINA → cell_type
# ==============================================================================

message("=== 05B. Rinomina cluster → cell_type ===")

# ------------------------------------------------------------------------------
# <<< ISTRUZIONI PER COMPILARE LE MAPPE >>>
#
# 1. Apri 05_DotPlot_prefilt_CTX.png e guarda quali geni si accendono
#    in ogni cluster numerico.
# 2. Per ogni cluster numerico inserisci il nome del tipo cellulare
#    come valore del vettore (la chiave è il numero del cluster).
# 3. I cluster da scartare (bassa qualità, doublet, contaminanti non
#    rimossi in 04) vanno inseriti in CLUSTERS_REMOVE_*.
#    Oppure assegnali a "Low_quality" e verranno rimossi automaticamente.
#
# Esempi tipici CTX:
#   "0"  = "Ex_L23_IT"   (Slc17a7+, Rasgrf2+, Cux1+)
#   "1"  = "Ex_L4_IT"    (Slc17a7+, Rorb+)
#   "2"  = "Ex_L5_IT"    (Slc17a7+, Bcl6+)
#   "3"  = "Ex_L5_PT"    (Slc17a7+, Fezf2+, Bcl11b+)
#   "4"  = "Ex_L6_CT"    (Slc17a7+, Foxp2+, Syt6+)
#   "5"  = "Inh_PV"      (Gad1+, Pvalb+)
#   "6"  = "Inh_SST"     (Gad1+, Sst+)
#   "7"  = "Inh_VIP"     (Gad1+, Vip+)
#   "8"  = "Inh_Lamp5"   (Gad1+, Lamp5+)
#   "9"  = "Astro"       (Aqp4+, Aldoc+)
#   "10" = "Oligo"       (Mog+, Mbp+)
#   "11" = "OPC"         (Pdgfra+, Cspg4+)
#   "12" = "Micro"       (Cx3cr1+, P2ry12+)
#   "13" = "Endo"        (Cldn5+, Pecam1+)
#
# Esempi tipici THAL:
#   "0" = "TC_relay"     (Slc17a6+, Rorb+, Grik4+)
#   "1" = "TC_HO"        (Slc17a6+, Calb1+)
#   "2" = "TRN"          (Gad1+, Pvalb+, Sst+)
#   "3" = "Astro"        (Aqp4+)
#   "4" = "Oligo"        (Mog+)
#   "5" = "OPC"          (Pdgfra+)
#   "6" = "Micro"        (Cx3cr1+)
#   "7" = "Endo"         (Cldn5+)
# ------------------------------------------------------------------------------

# --- Cluster da rimuovere prima della rinomina --------------------------------
# (numeri come stringhe, vettore vuoto = nessuna rimozione)
CLUSTERS_REMOVE_CTX  <- c()   # <<< es. c("14","15")
CLUSTERS_REMOVE_THAL <- c()   # <<< es. c("8")

# --- Mappa rinomina CTX -------------------------------------------------------
# <<< DA COMPILARE: chiave = numero cluster (stringa), valore = nome CT >>>
ctx_rename <- c(
  
  "1"  = "Ex1",    
  "2"  = "Ex_L23_1",
  "3"  = "Ex_L23_2",
  "4"  = "Ex2",
  "5"  = "Astro",
  "6"  = "Oligo",
  "7"  = "Ex3",
  "8"  = "Inh_SST",
  "9"  = "Ex_L6CT_1",
  "10"  = "Ex_L5PT_1",
  "11"  = "Inh_PV",
  "12"  = "Inh_VIP",
  "13"  = "Endo",
  "14"  = "Micro",
  "15"  = "Ex_L5PT_2",
  "16"  = "OPC",
  "17"  = "Ex_L6CT_2",
  "18"  = "Ex_L6CT_3"
  
)

# Esempio da decommentare e adattare:
# ctx_rename <- c(
#   "0"  = "Ex_L23_IT",
#   "1"  = "Ex_L4_IT",
#   "2"  = "Ex_L5_IT",
#   "3"  = "Ex_L5_PT",
#   "4"  = "Ex_L6_CT",
#   "5"  = "Inh_PV",
#   "6"  = "Inh_SST",
#   "7"  = "Inh_VIP",
#   "8"  = "Inh_Lamp5",
#   "9"  = "Astro",
#   "10" = "Oligo",
#   "11" = "OPC",
#   "12" = "Micro",
#   "13" = "Endo"
# )

# --- Mappa rinomina THAL ------------------------------------------------------
# <<< DA COMPILARE >>>
thal_rename <- NULL

# Esempio da decommentare e adattare:
# thal_rename <- c(
#   "0" = "TC_relay",
#   "1" = "TC_HO",
#   "2" = "TRN",
#   "3" = "Astro",
#   "4" = "Oligo",
#   "5" = "OPC",
#   "6" = "Micro",
#   "7" = "Endo"
# )

# ==============================================================================
# ESECUZIONE RINOMINA + RIMOZIONE
# ==============================================================================


  clusters_remove <- if (obj_nm == "CTX") CLUSTERS_REMOVE_CTX else CLUSTERS_REMOVE_THAL
  rename_map      <- if (obj_nm == "CTX") ctx_rename          else thal_rename

  n_before <- ncol(CTX)

  # ---- 1. Rimozione cluster espliciti ----------------------------------------
  if (length(clusters_remove) > 0) {
    bad <- as.character(clusters_remove)
    cells_keep <- !(as.character(Idents(CTX)) %in% bad)
    CTX <- CTX[, cells_keep]
    message(sprintf("  %s: rimossi cluster %s → -%d cellule",
                    obj_nm,
                    paste(bad, collapse = ","),
                    n_before - ncol(CTX)))
  }

  # ---- 2. Rinomina cluster → label cell_type ---------------------------------
  if (!is.null(rename_map)) {

    # Verifica che tutte le chiavi corrispondano a cluster esistenti
    missing_keys <- setdiff(names(rename_map), as.character(levels(Idents(CTX))))
    if (length(missing_keys) > 0)
      warning(sprintf("  %s: chiavi non trovate nella mappa: %s",
                      obj_nm, paste(missing_keys, collapse = ", ")))

    CTX <- RenameIdents(CTX, rename_map)
    message(sprintf("  %s: cluster rinominati in %d tipi cellulari",
                    obj_nm, length(unique(as.character(Idents(CTX))))))

    # Rimuovi eventuali celle con label "Low_quality" (se assegnata in rename_map)
    if ("Low_quality" %in% levels(Idents(CTX))) {
      n_lq   <- sum(Idents(CTX) == "Low_quality")
      CTX    <- CTX[, Idents(CTX) != "Low_quality"]
      message(sprintf("  %s: rimossi %d cluster 'Low_quality'", obj_nm, n_lq))
    }

  } else {
    # Nessuna mappa: usa i numeri di cluster come cell_type provvisori
    message(sprintf(
      "  %s: nessuna mappa fornita — cell_type = cluster numerico (provvisorio)",
      obj_nm
    ))
  }

  # ---- 3. Salva cell_type come colonna in @meta.data -------------------------
  CTX$cell_type <- as.character(Idents(CTX))


# ==============================================================================
# PARTE C — PLOT POST-ANNOTAZIONE
# ==============================================================================
 
  message("\n=== 05C. Plot post-annotazione ===")
  
  # --- 05C.0 — UMAP score per cell type: un pannello per tipo cellulare ------
  #
  # Per ogni cell_type presente nel dataset viene calcolato uno score combinato
  # (AddModuleScore) dei suoi marcatori specifici, ricavati invertendo gene_ct_map_ctx.
  # Ogni pannello = un cell type; il colore mostra lo score sulla UMAP.
  # I pannelli sono ordinati per frequenza (cell type più abbondante prima).
  # Scopo: verificare che ogni tipo cellulare abbia espressione specifica
  #        nella regione UMAP attesa, dopo la rinomina.
  
  # ---- Mappa inversa cell_type → geni marcatori (da gene_ct_map_ctx) ---------
  ct_markers_ctx <- split(names(gene_ct_map_ctx),
                          unname(gene_ct_map_ctx))
  # Unifica gruppi simili (Ex/Ex_pan → tutti i marcatori eccitatori pan)
  # già separati nella mappa; nessuna modifica necessaria
  
  # Tieni solo i cell_type effettivamente presenti e con >= 2 marcatori nel dataset
  ct_present <- sort(unique(CTX$cell_type))
  
    ct_map_use  <- if (obj_nm == "CTX") gene_ct_map_ctx else gene_ct_map_thal
    ct_present  <- sort(unique(CTX$cell_type))
    
    # Mappa inversa per quest'area
    ct_markers <- split(names(ct_map_use), unname(ct_map_use))
    
    # Ordina cell_type per abbondanza (più cellule = pannello più importante)
    ct_freq <- sort(table(CTX$cell_type), decreasing = TRUE)
    ct_order <- intersect(names(ct_freq), names(ct_markers))
    
    # Filtra: tieni solo CT con >= 2 marcatori presenti nel dataset
    ct_order <- ct_order[sapply(ct_order, function(ct) {
      genes <- filter_genes(ct_markers[[ct]], CTX)
      length(genes) >= 2
    })]
    
    if (length(ct_order) == 0) {
      message(sprintf("  %s: nessun cell_type con marcatori sufficienti", obj_nm))
      next
    }
    
    # ---- Calcola AddModuleScore per ogni cell_type ----------------------------
    # Ogni score viene salvato come colonna "score_<CT>" in @meta.data
    # (AddModuleScore usa geni random come controllo e centra su 0)
    score_cols <- c()
    for (ct in ct_order) {
      genes     <- filter_genes(ct_markers[[ct]], CTX)
      score_nm  <- paste0("score_", gsub("[^A-Za-z0-9]", "_", ct))
      score_cols <- c(score_cols, score_nm)
      
      CTX <- AddModuleScore(
        CTX,
        features = list(genes),
        name     = score_nm,
        seed     = 42
      )
      # AddModuleScore aggiunge "1" alla fine del nome — rinomina
      colnames(CTX@meta.data)[
        colnames(CTX@meta.data) == paste0(score_nm, "1")
      ] <- score_nm
    }
    
    # ---- Un FeaturePlot per cell_type ----------------------------------------
    n_cols    <- 4
    plot_list <- lapply(seq_along(ct_order), function(i) {
      ct       <- ct_order[i]
      score_nm <- score_cols[i]
      n_cells  <- sum(CTX$cell_type == ct)
      
      FeaturePlot(
        CTX,
        features  = score_nm,
        reduction = "umap",
        pt.size   = 0.2,
        combine   = FALSE
      )[[1]] +
        scale_color_gradient2(
          low      = "blue",
          mid      = "white",
          high     = "#C0392B",
          midpoint = 0,
          name     = "score"
        ) +
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
    savepng(
      wrap_plots(plot_list, ncol = n_cols) +
        plot_annotation(
          title    = paste(obj_nm, "— Score marcatori per cell type"),
          subtitle = "Colore = AddModuleScore(marcatori specifici) | Scala divergente: grigio=0, arancio→rosso=alto",
          theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                           plot.subtitle = element_text(size = 8, color = "grey40"))
        ),
      paste0("05_ModuleScore_per_CT_", obj_nm, ".png"),
      width  = n_cols * 1800,
      height = n_rows * 1800 + 400
    )
    message(sprintf("  %s: salvato 05_ModuleScore_per_CT_%s.png (%d pannelli)",
                    obj_nm, obj_nm, length(plot_list)))
    

  
  # --- 05C.0 — FeaturePlot geni chiave con etichetta cell_type ---------------
  #
  # Ora che cell_type è assegnato, ogni pannello mostra il nome del tipo
  # cellulare dominante nel cluster con espressione più alta per quel gene.
  
  key_genes <- filter_genes(
    c("Slc17a7","Slc17a6","Gad1","Gad2","Pvalb","Sst","Vip","Aqp4",
      "Mog","Pdgfra","Cx3cr1","Cldn5","Htr2a","Fos","Arc"),
    CTX
  )
  
     # Coordinate UMAP + cell_type per ogni cellula
    umap_coords <- as.data.frame(Embeddings(CTX, "umap"))
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    umap_coords$cell_type <- CTX$cell_type
    
    # Centroide di ogni cell_type nello spazio UMAP
    ct_centroids <- umap_coords %>%
      group_by(cell_type) %>%
      summarise(cx = mean(UMAP_1), cy = mean(UMAP_2), .groups = "drop")
    
    # Espressione media per gene x cell_type
    assay_use <- "SCT" 
    avg_expr  <- AverageExpression(
      CTX,
      features = key_genes,
      assays   = assay_use,
      layer    = "data",
      group.by = "cell_type",
      verbose  = FALSE
    )[[assay_use]]
    
    # Per ogni gene: cell_type con espressione media massima → etichetta
    gene_labels <- lapply(key_genes, function(gene) {
      if (!gene %in% rownames(avg_expr)) return(NULL)
      top_ct   <- names(which.max(avg_expr[gene, ]))
      centroid <- ct_centroids[ct_centroids$cell_type == top_ct, ]
      if (nrow(centroid) == 0) return(NULL)
      data.frame(gene = gene, ct_label = top_ct,
                 cx = centroid$cx, cy = centroid$cy)
    })
    gene_labels <- do.call(rbind, Filter(Negate(is.null), gene_labels))
    
    plot_list <- lapply(key_genes, function(gene) {
      p <- FeaturePlot(
        CTX, features = gene,
        reduction = "umap", pt.size = 0.2, combine = FALSE
      )[[1]] +
        scale_color_viridis_c(option = "inferno") +
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
      
      lbl_row <- gene_labels[gene_labels$gene == gene, ]
      if (nrow(lbl_row) > 0) {
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
    
    n_cols <- 4
    savepng(wrap_plots(plot_list, ncol = n_cols),
            "05_FeaturePlot_key_genes_CTX.png",
            width  = 6000,
            height = ceiling(length(plot_list) / n_cols) * 1600)
  
  
  # --- 05C.1 — Heatmap cluster × MapMyCells subclass con cell_type reale ------
  #
  # Le righe (cluster originali) sono annotate con il cell_type finale
  # e ordinate per tipo — mostra quanto ogni cluster e' "puro" secondo MMC.
  
    # Usa cluster originale (Idents prima della rinomina) se ancora disponibile,
    # altrimenti raggruppa per cell_type direttamente
    group_var <- "cell_type"
    
    ct_tab  <- table(CTX@meta.data[[group_var]], CTX$mmc_subclass)
    ct_prop <- prop.table(ct_tab, margin = 1)
    
    # Annotazione righe con cell_type
      row_ct <- rownames(ct_prop)
    
    # Ordina righe per cell_type (blocchi contigui)
    row_order       <- order(row_ct, rownames(ct_prop))
    ct_prop_ordered <- ct_prop[row_order, ]
    row_ct_ordered  <- row_ct[row_order]
    
    # Palette annotazione
    ct_unique  <- unique(row_ct_ordered)
     
    row_ann <- data.frame(CellType = row_ct_ordered,
                          row.names = rownames(ct_prop_ordered))
    
    # Filtra subclass con < 1% in tutti i cluster (riduce rumore)
    col_keep     <- apply(ct_prop_ordered, 2, max) >= 0.01
    ct_prop_plot <- ct_prop_ordered[, col_keep, drop = FALSE]
    
    png(file.path(outdir, paste0("05_cluster_vs_MMC_", obj_nm, ".png")),
        width  = max(3500, ncol(ct_prop_plot) * 65 + 700),
        height = max(2000, nrow(ct_prop_plot) * 60 + 400),
        res    = 300)
    pheatmap(
      ct_prop_plot,
      color             = colorRampPalette(c("white", "#1A5276", "#2ECC71"))(80),
      display_numbers   = FALSE,
      cluster_rows      = T,
      cluster_cols      = TRUE,,
      fontsize          = 7, fontsize_row = 7, fontsize_col = 7,border_color = NA,
      main = paste("Cluster vs MapMyCells subclass (cell_type) —", obj_nm)
    )
    dev.off()
    message(sprintf("  Salvato: 05_cluster_vs_MMC_%s.png", obj_nm))
  
  
message("\n=== 05C. Plot post-annotazione ===")

  ct_levels <- sort(unique(CTX$cell_type))
  n_ct      <- length(ct_levels)
  ct_colors <- setNames(pal_ct[seq_len(n_ct)], ct_levels)

  # --- UMAP: cell_type + trattamento + timepoint + campione -----------------
  p_ct <- DimPlot(
    CTX, group.by = "cell_type", reduction = "umap",
    label = TRUE, repel = TRUE, label.size = 3,
    pt.size = 0.3, cols = ct_colors
  ) +
    ggtitle(paste(obj_nm, "— Tipo cellulare")) +
    theme(legend.text = element_text(size = 8))

  p_trt <- DimPlot(
    CTX, group.by = "treatment", reduction = "umap",
    pt.size = 0.3, cols = pal_trt
  ) + ggtitle(paste(obj_nm, "— Trattamento"))

  p_time <- DimPlot(
    CTX, group.by = "time", reduction = "umap",
    pt.size = 0.3, cols = pal_time
  ) + ggtitle(paste(obj_nm, "— Timepoint"))

  p_samp <- DimPlot(
    CTX, group.by = "orig.ident", reduction = "umap",
    pt.size = 0.3
  ) + ggtitle(paste(obj_nm, "— Campione"))

  savepng(
    (p_ct | p_trt) / (p_time | p_samp) +
      plot_annotation(title = paste("Annotazione finale —", obj_nm)),
    paste0("05_UMAP_annotated_", obj_nm, ".png"),
    width = 6000, height = 5000
  )

  # --- DotPlot marcatori per cell_type (cluster ordinati) -------------------
  markers_use <- if (obj_nm == "CTX") {
    filter_genes(unique(c(markers_broad, markers_exc_ctx, markers_inh_ctx)), CTX)
  } else {
    filter_genes(unique(c(markers_broad, markers_thal)), CTX)
  }
  # Massimo 60 geni per leggibilità
  markers_use <- markers_use[seq_len(min(60, length(markers_use)))]

  p_dot_ann <- DotPlot(
    CTX, features = markers_use,
    group.by = "cell_type", cluster.idents = FALSE
  ) +
    RotatedAxis() +
    ggtitle(paste(obj_nm, "— Marcatori per tipo cellulare"))
  savepng(p_dot_ann, paste0("05_DotPlot_annotated_", obj_nm, ".png"),
          width = 7000, height = max(2500, n_ct * 160 + 600))

  # --- FeaturePlot Htr2a + IEG ---------------------------------------------
  feat_show <- filter_genes(c("Htr2a","Fos","Arc","Egr1","Nr4a1"), CTX)
  if (length(feat_show) > 0) {
    p_feat <- FeaturePlot(
      CTX, features = feat_show,
      reduction = "umap", ncol = 3, pt.size = 0.2
    ) & scale_color_viridis_c(option = "inferno") &
      theme(legend.position = "right", plot.title = element_text(size = 9))
    savepng(p_feat, paste0("05_FeaturePlot_Htr2a_IEG_", obj_nm, ".png"),
            width = 6000,
            height = ceiling(length(feat_show) / 3) * 2000)
  }

  # --- Heatmap proporzione cell_type × campione ----------------------------
  ct_prop <- prop.table(table(CTX$cell_type, CTX$orig.ident), margin = 2)
  png(file.path(outdir, paste0("05_CellType_composition_", obj_nm, ".png")),
      width  = 3000,
      height = max(2000, nrow(ct_prop) * 80 + 600),
      res    = 300)
  pheatmap(
    ct_prop,
    color           = colorRampPalette(c("white","#27AE60"))(50),
    display_numbers = TRUE, number_format = "%.2f",
    fontsize         = 8,  cluster_cols = FALSE,
    main             = paste("Proporzione CT per campione —", obj_nm)
  )
  dev.off()


# ==============================================================================
# SALVATAGGIO
# ==============================================================================

message("\n=== 05. Salvataggio ===")

save(CTX,  file = file.path(outdir, "CTX_celltype.RData"))

