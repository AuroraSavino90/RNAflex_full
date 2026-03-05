################################################################################
##  04_mapmycells.R                                                            ##
##                                                                             ##
##  Input:  CTX_initial.RData, THAL_initial.RData  (da sezione 03)            ##
##  Output: brain_combined.h5ad  (da caricare su MapMyCells)                  ##
##          CTX_annotated.RData, THAL_annotated.RData  (post-MapMyCells)      ##
##                                                                             ##
##  Flusso:                                                                    ##
##    PARTE A — esegui fino all'export h5ad, poi vai su MapMyCells online     ##
##    PARTE B — decommentare dopo aver scaricato il CSV risultati             ##
################################################################################

source("code/00_config.R")
outdir<-"results/20260224"
# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_initial.RData"))   # CTX

# ==============================================================================
# PARTE A — EXPORT H5AD PER MAPMYCELLS
# ==============================================================================

message("=== 04A. Export h5ad per MapMyCells ===")

# Seurat v5: LayerData() al posto di @assays$RNA@counts
# I conteggi raw sono nel layer "counts" dell'assay RNA
# (dopo JoinLayers in sezione 03 è un layer unico)

export_h5ad <- function(obj, path, label) {

  # Estrai matrice counts RNA (layer "counts" dopo JoinLayers)
  counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")

  # AnnData vuole cellule come righe, geni come colonne (transposta)
  counts_t <- Matrix::t(counts_mat)

  ad <- AnnData(
    X   = counts_t,
    obs = data.frame(
      sample    = obj$orig.ident,
      treatment = obj$treatment,
      time      = obj$time,
      row.names = rownames(counts_t)
    ),
    var = data.frame(
      gene_name = colnames(counts_t),
      row.names = colnames(counts_t)
    )
  )

  write_h5ad(ad, path, compression = "gzip")
  message(sprintf("  %s: salvato %s (%.1f MB, %d cellule, %d geni)",
                  label, path, file.size(path) / 2^20,
                  nrow(counts_t), ncol(counts_t)))
}

export_h5ad(CTX,  file.path(outdir, "CTX_for_mapmycells.h5ad"),  "CTX")


message("
>>> PROSSIMI PASSI (PARTE A completata):
  1. Vai su https://portal.brain-map.org/atlases-and-data/bkp/mapmycells
  2. Carica CTX_for_mapmycells.h5ad 
  3. Seleziona: Whole Mouse Brain (CCN20230722), Hierarchical mapping
  4. Scarica il CSV dei risultati (es. 'MapMyCells_results.csv')
  5. Scarica il metadata Excel 'cl.df_CCN20230722.xlsx' dal sito AllenBrain
  6. Copia i file in: ", outdir, "
  7. Passa alla parte B
")

# ==============================================================================
# PARTE B — PROCESSAMENTO RISULTATI MAPMYCELLS
# (decommentare dopo aver scaricato il CSV)
# ==============================================================================

# --- Parametri ---------------------------------------------------------------
MAPPING_CSV <- file.path(outdir, "CTX_from_mapmycells/CTX_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1772014289213.csv")
META_XLSX   <- file.path("data/cl.df_CCN202307220.xlsx")

# Soglie: cellule con meno di N assegnazioni alla stessa classe → "other"
MIN_CELLS_CLASS    <- 100
MIN_CELLS_SUBCLASS <- 20

# --- Esegui solo se il file esiste -------------------------------------------

  message("=== 04B. Processamento MapMyCells ===")

  mapping <- read.csv(MAPPING_CSV, comment.char = "#")
  meta_df <- read.xlsx(META_XLSX)

  # Verifica colonne attese
  expected_cols <- c("cell_id","class_name","subclass_name","cluster_name")
  if (!all(expected_cols %in% colnames(mapping))) {
    stop(sprintf(
      "Colonne mancanti nel CSV MapMyCells. Attese: %s. Trovate: %s",
      paste(expected_cols, collapse = ", "),
      paste(colnames(mapping), collapse = ", ")
    ))
  }
  rownames(mapping) <- mapping$cell_id

  # --- 04B.1 — Classifica "other" per classi rare ---------------------------
  class_freq    <- table(mapping$class_name)
  sub_freq      <- table(mapping$subclass_name)

  keep_class <- names(class_freq[class_freq >= MIN_CELLS_CLASS])
  keep_sub   <- names(sub_freq[sub_freq   >= MIN_CELLS_SUBCLASS])

  mapping$class_clean    <- ifelse(mapping$class_name    %in% keep_class,
                                   mapping$class_name,    "other")
  mapping$subclass_clean <- ifelse(mapping$subclass_name %in% keep_sub,
                                   mapping$subclass_name, "other")

  # --- 04B.2 — Aggiungi annotazioni agli oggetti Seurat ---------------------
  add_mapmycells <- function(obj, mapping_df, label) {

    # Solo le cellule dell'oggetto (in caso di export combinato)
    shared_cells <- intersect(colnames(obj), rownames(mapping_df))
    if (length(shared_cells) == 0) {
      message(sprintf("  WARN: nessuna cellula condivisa tra %s e mapping", label))
      return(obj)
    }
    missing <- setdiff(colnames(obj), rownames(mapping_df))
    if (length(missing) > 0)
      message(sprintf("  WARN %s: %d cellule senza mapping", label, length(missing)))

    m <- mapping_df[shared_cells, ]

    obj$mmc_class    <- NA_character_
    obj$mmc_subclass <- NA_character_
    obj$mmc_cluster  <- NA_character_

    obj$mmc_class[   match(shared_cells, colnames(obj))] <- m$class_clean
    obj$mmc_subclass[match(shared_cells, colnames(obj))] <- m$subclass_clean
    obj$mmc_cluster[ match(shared_cells, colnames(obj))] <- m$cluster_name

    message(sprintf("  %s: %d/%d cellule annotate | %d classi | %d subclassi",
                    label, length(shared_cells), ncol(obj),
                    length(unique(m$class_clean[m$class_clean != "other"])),
                    length(unique(m$subclass_clean[m$subclass_clean != "other"]))))
    obj
  }

  CTX  <- add_mapmycells(CTX,  mapping, "CTX")
  obj_nm<-"CTX"#nome dell'area nei file di results
  
  # --- 04B.3 — Plot: cluster Seurat vs MapMyCells subclass ------------------
    if (!"mmc_subclass" %in% colnames(CTX@meta.data)) next

    # Heatmap: cluster numerico (Seurat) × subclass (MapMyCells)
    ct_tab <- table(Idents(CTX), CTX$mmc_subclass)
    ct_prop <- prop.table(ct_tab, margin = 1)  # proporzione per cluster Seurat

    n_sub <- ncol(ct_prop)
    n_cl  <- nrow(ct_prop)

    png(file.path(outdir, paste0("04_cluster_vs_MMC_", obj_nm, ".png")),
        width  = max(3000, n_sub * 80 + 1000),
        height = max(2000, n_cl  * 80 + 600),
        res    = 300)
    pheatmap(
      ct_prop,
      color           = colorRampPalette(c("white","#27AE60"))(50),
      display_numbers = FALSE,
      fontsize        = 7,
      main            = paste("Cluster Seurat vs MapMyCells subclass —", obj_nm)
    )
    dev.off()

    # UMAP colorato per subclass MapMyCells
    # "other" in grigio, disegnato per primo (rimane sotto le classi reali)
    sub_levels   <- sort(unique(CTX$mmc_subclass[!is.na(CTX$mmc_subclass)]))
    sub_real     <- sub_levels[sub_levels != "other"]
    pal_sub      <- setNames(
      c(scales::hue_pal()(length(sub_real)), "grey70"),
      c(sub_real, "other")
    )
    
    p_mmc <- DimPlot(CTX, group.by = "mmc_subclass",
                     label = TRUE, repel = TRUE, label.size = 2.5,
                     pt.size = 0.3,
                     cols  = pal_sub,
                     order = c(sub_real, "other")  # "other" disegnato per primo → sotto
    ) +
      ggtitle(paste("MapMyCells subclass —", obj_nm)) +
      theme(legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"))
    savepng(p_mmc, paste0("04_UMAP_MMC_subclass_", obj_nm, ".png"),
            width = 6500, height = 3500)
    
    # UMAP colorato per class MapMyCells — "other" in grigio
    cl_levels <- sort(unique(CTX$mmc_class[!is.na(CTX$mmc_class)]))
    cl_real   <- cl_levels[cl_levels != "other"]
    pal_cl    <- setNames(
      c(scales::hue_pal()(length(cl_real)), "grey70"),
      c(cl_real, "other")
    )
    
    p_mmc_cl <- DimPlot(CTX, group.by = "mmc_class",
                        label = TRUE, repel = TRUE, label.size = 3,
                        pt.size = 0.3,
                        cols  = pal_cl,
                        order = c(cl_real, "other")
    ) +
      ggtitle(paste("MapMyCells class —", obj_nm))
    savepng(p_mmc_cl, paste0("04_UMAP_MMC_class_", obj_nm, ".png"),
            width = 4000, height = 3000)
    
    
    # --- 04B.4 — Annotazione anatomica broad da meta_df --------------------------
    # broad_anat viene estratto dalla colonna CCF_broad.freq del metadata Excel
    # (formato "Isocortex:..."); prendiamo solo la parte prima dei ":" per avere
    # la macroarea CCF (es. "Isocortex", "Thalamus", "Striatum", ecc.)
    # fine_anat è l'annotazione anatomica fine-grained (colonna anatomical_annotation)
    
    fine_anat  <- meta_df$anatomical_annotation
    broad_anat <- sub("\\:.*", "", meta_df$CCF_broad.freq)
    
    # Lookup table cluster_id → broad_anat / fine_anat
    # Usiamo cluster_id_label come chiave (stessa colonna usata per match in mapping)
    anat_lookup <- data.frame(
      cluster_id = meta_df$cluster_id_label,
      broad_anat = broad_anat,
      fine_anat  = fine_anat,
      stringsAsFactors = FALSE
    )
    
    # Aggiungi broad_anat e fine_anat a ogni oggetto tramite mmc_cluster

      
      idx <- match(CTX$mmc_cluster, anat_lookup$cluster_id)
      
      CTX$broad_anat <- anat_lookup$broad_anat[idx]   # NA se cluster non trovato
      CTX$fine_anat  <- anat_lookup$fine_anat[idx]
      
      n_mapped <- sum(!is.na(CTX$broad_anat))
      message(sprintf("  %s: broad_anat assegnato a %d/%d cellule",
                      obj_nm, n_mapped, ncol(obj)))
    
    # Stampa le macroaree distinte trovate (utile per compilare le soglie sotto)
    cat("\n=== Macroaree CCF presenti (broad_anat) ===\n")
    cat("CTX: "); cat(sort(unique(CTX$broad_anat)), sep=", "); cat("\n")
    
    # --- 04B.5 — Analisi contaminazione anatomica --------------------------------
    # Testa quante cellule di CTX non sono corticali e quante di THAL non sono
    # talamiche secondo la macroarea CCF assegnata da MapMyCells (broad_anat).
    #
    # <<< Modificare le liste se le etichette nel tuo dataset differiscono >>>
    # (le etichette esatte compaiono nel printout qui sopra)
    CORTEX_CLASSES  <- c(
      "CTXsp",
      "Isocortex",
      "NA",
      "OLF"
    )
    THALAMUS_CLASSES <- c(
      "Thalamus",
      "Epithalamus"
    )
    
    CONTAM_THRESHOLD <- 0.05   # soglia accettabile per il test binomiale (5%)
    
   
      
      if (!"broad_anat" %in% colnames(CTX@meta.data)) next
      
      expected_classes <- if (obj_nm == "CTX") CORTEX_CLASSES else THALAMUS_CLASSES
      area_label       <- if (obj_nm == "CTX") "Cortex"       else "Thalamus"
      
      meta <- CTX@meta.data
      # Cellule senza mapping → "unmapped"
      meta$broad_anat_clean <- ifelse(is.na(meta$broad_anat), "unmapped",
                                      meta$broad_anat)
      
      # ---- Proporzioni per macroarea ----------------------------------------
      class_tab <- as.data.frame(table(meta$broad_anat_clean)) %>%
        rename(broad_anat = Var1, n_cells = Freq) %>%
        mutate(
          pct         = round(100 * n_cells / sum(n_cells), 2),
          is_expected = broad_anat %in% expected_classes
        ) %>%
        arrange(desc(n_cells))
      
      cat(sprintf("\n=== Proporzione macroaree CCF (broad_anat) — %s ===\n", obj_nm))
      print(class_tab)
      
      n_expected    <- sum(class_tab$n_cells[ class_tab$is_expected])
      n_contaminant <- sum(class_tab$n_cells[!class_tab$is_expected])
      pct_expected  <- round(100 * n_expected    / nrow(meta), 2)
      pct_contam    <- round(100 * n_contaminant / nrow(meta), 2)
      
      cat(sprintf(
        "\n  Attese (%s): %d (%.2f%%)\n  Contaminanti: %d (%.2f%%)\n",
        paste(expected_classes, collapse = " / "),
        n_expected, pct_expected,
        n_contaminant, pct_contam
      ))
      
      # ---- Test chi-quadro: proporzione contaminanti omogenea tra campioni? --
      meta$is_contaminant <- !(meta$broad_anat_clean %in% expected_classes)
      cont_tab <- table(meta$orig.ident, meta$is_contaminant)
      
      # Assicura che ci siano entrambe le colonne (possibile se 0 contaminanti)
      if (ncol(cont_tab) == 2) {
        colnames(cont_tab) <- c("expected", "contaminant")
        
        cont_df <- as.data.frame(cont_tab) %>%
          rename(sample = Var1, status = Var2, n = Freq) %>%
          group_by(sample) %>%
          mutate(pct = round(100 * n / sum(n), 2)) %>%
          ungroup() %>%
          filter(status == "contaminant") %>%
          select(sample, n, pct) %>%
          arrange(desc(pct))
        
        cat(sprintf("\n--- Contaminanti per campione — %s ---\n", obj_nm))
        print(cont_df)
        
        if (nrow(cont_tab) > 1 && all(cont_tab > 0)) {
          chi_res <- chisq.test(cont_tab)
          cat(sprintf(
            "\n  Chi-quadro tra campioni: X²=%.2f, df=%d, p=%.4f\n",
            chi_res$statistic, chi_res$parameter, chi_res$p.value
          ))
          if (chi_res$p.value < 0.05)
            cat("  *** Proporzione contaminanti significativamente diversa tra campioni\n")
          else
            cat("  Proporzione contaminanti omogenea tra campioni (p >= 0.05)\n")
        }
      } else {
        # Nessun contaminante trovato
        cont_df <- data.frame(sample = rownames(cont_tab), n = 0, pct = 0)
        cat(sprintf("\n  Nessun contaminante trovato in %s\n", obj_nm))
      }
      

      # ---- Plot 1: barplot proporzioni macroaree broad_anat ------------------
      # Classi < 1% e non-attese vengono raggruppate in "altre (<1%)"
      class_plot <- class_tab %>%
        mutate(
          label = ifelse(pct < 1 & !is_expected,
                         "altre (<1%)", as.character(broad_anat)),
          fill  = case_when(
            is_expected             ~ "#27AE60",
            label == "altre (<1%)" ~ "grey70",
            label == "unmapped"    ~ "grey40",
            TRUE                   ~ "#E74C3C"
          )
        ) %>%
        group_by(label, fill) %>%
        summarise(n_cells = sum(n_cells), .groups = "drop") %>%
        mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
        arrange(desc(n_cells))
      
      fill_pal <- setNames(class_plot$fill, class_plot$label)
      
      p_bar_anat <- ggplot(class_plot,
                           aes(x = reorder(label, n_cells),
                               y = pct, fill = label)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = sprintf("%.1f%%", pct)),
                  hjust = -0.1, size = 3) +
        coord_flip() +
        scale_fill_manual(values = fill_pal, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_bw(base_size = 10) +
        labs(
          title    = paste("Macroaree CCF (broad_anat) —", obj_nm),
          subtitle = sprintf(
            "Verde = atteso (%s) | Rosso = contaminante | n=%d cellule",
            area_label, nrow(meta)),
          x = "", y = "% cellule"
        )
      savepng(p_bar_anat,
              paste0("04_anatomy_barplot_", obj_nm, ".png"),
              width = 3500, height = max(2000, nrow(class_plot) * 120 + 700))
      
     
      # ---- Plot 3: UMAP con contaminanti evidenziati -------------------------
      # Cellule attese in grigio chiaro (sotto), contaminanti in colori forti
      CTX$contam_label <- ifelse(
        meta$broad_anat_clean %in% expected_classes, "Atteso",
        ifelse(meta$broad_anat_clean == "unmapped",  "Non mappato",
               meta$broad_anat_clean)
      )
      
      contam_classes <- setdiff(sort(unique(CTX$contam_label)),
                                c("Atteso", "Non mappato"))
      n_contam <- length(contam_classes)
      
      # Palette scalabile: combina Set1 (9) + Dark2 (8) + Paired (12) poi
      # interpola con colorRampPalette se servono ancora più colori
      base_cols <- unique(c(
        RColorBrewer::brewer.pal(9,  "Set1"),
        RColorBrewer::brewer.pal(8,  "Dark2"),
        RColorBrewer::brewer.pal(12, "Paired")
      ))
      if (n_contam > length(base_cols)) {
        base_cols <- colorRampPalette(base_cols)(n_contam)
      }
      pal_contam <- c(
        "Atteso"      = "grey85",
        "Non mappato" = "grey50",
        setNames(base_cols[seq_len(n_contam)], contam_classes)
      )
      order_contam <- c("Atteso", "Non mappato", contam_classes)
      
      p_umap_contam <- DimPlot(
        CTX, group.by = "contam_label",
        pt.size = 0.3,
        cols    = pal_contam,
        order   = order_contam
      ) +
        ggtitle(paste("Contaminanti anatomici (broad_anat) —", obj_nm)) +
        theme(legend.text = element_text(size = 8))
      
      savepng(p_umap_contam,
              paste0("04_UMAP_contamination_", obj_nm, ".png"),
              width = 4000, height = 3500)
      
      # ---- Plot 4: UMAP con fine_anat (dettaglio anatomico) ------------------
      # Utile per capire da quale struttura specifica provengono i contaminanti
  
        # Mostra solo le regioni con >= 10 cellule; le altre → "altre"
        fine_freq  <- table(CTX$fine_anat)
        fine_keep  <- names(fine_freq[fine_freq >= 10])
        CTX$fine_anat_plot <- ifelse(
          CTX$fine_anat %in% fine_keep,
          CTX$fine_anat, "altre (<10 cellule)"
        )
        CTX$fine_anat_plot[is.na(CTX$fine_anat_plot)] <- "unmapped"
        
        fine_levels <- sort(unique(CTX$fine_anat_plot))
        fine_real   <- setdiff(fine_levels, c("altre (<10 cellule)", "unmapped"))
        pal_fine    <- setNames(
          c(scales::hue_pal()(length(fine_real)), "grey70", "grey40"),
          c(fine_real, "altre (<10 cellule)", "unmapped")
        )
        
        p_umap_fine <- DimPlot(
          CTX, group.by = "fine_anat_plot",
          pt.size = 0.2,
          label   = FALSE,
          cols    = pal_fine,
          order   = c(fine_real, "altre (<10 cellule)", "unmapped")
        ) +
          ggtitle(paste("Annotazione anatomica fine (fine_anat) —", obj_nm)) +
          theme(legend.text     = element_text(size = 7),
                legend.key.size = unit(0.35, "cm"))
        savepng(p_umap_fine,
                paste0("04_UMAP_fine_anat_", obj_nm, ".png"),
                width = 6000, height = 3500)
      
      
      # Salva tabelle
      write.csv(class_tab,
                file.path(outdir, paste0("04_anatomy_proportions_", obj_nm, ".csv")),
                row.names = FALSE)
      write.csv(cont_df,
                file.path(outdir, paste0("04_contamination_by_sample_", obj_nm, ".csv")),
                row.names = FALSE)
      

      
      # --- 04B.6 — Filtraggio: tieni solo cellule nelle macroaree attese ----------
      # Rimuove i contaminanti anatomici identificati nell'analisi broad_anat.
      # Le cellule con broad_anat NA (unmapped) vengono rimosse per default:
      # impostare KEEP_UNMAPPED = TRUE per tenerle se si vuole essere conservativi.
      
      KEEP_UNMAPPED <- FALSE   # <<< TRUE per tenere le cellule non mappate
      
       expected_classes <- if (obj_nm == "CTX") CORTEX_CLASSES else THALAMUS_CLASSES
        n_before <- ncol(CTX)
        
        # Maschera: cellule da tenere
        is_expected <- CTX$broad_anat %in% expected_classes
        is_unmapped <- is.na(CTX$broad_anat)
        
        keep_mask <- if (KEEP_UNMAPPED) is_expected | is_unmapped else is_expected
        
        obj_filt  <- CTX[, keep_mask]
        n_after   <- ncol(obj_filt)
        n_removed <- n_before - n_after
        
        message(sprintf(
          "  %s: %d → %d cellule (rimossi %d contaminanti, %.1f%%)",
          obj_nm, n_before, n_after, n_removed,
          100 * n_removed / n_before
        ))
        
        # Riepilogo per campione dopo il filtraggio
        cat(sprintf("\n  Cellule per campione dopo filtraggio — %s:\n", obj_nm))
        print(table(obj_filt$orig.ident))
        
        CTX<-obj_filt
      
  # --- 04B.4 — Salvataggio --------------------------------------------------
  save(CTX,  file = file.path(outdir, "CTX_annotated.RData"))
  message("Salvati: CTX_annotated.RData")

