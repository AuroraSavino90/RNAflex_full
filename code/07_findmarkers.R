################################################################################
##  07_findmarkers.R                                                           ##
##                                                                             ##
##  Input:  CTX_celltype.RData   (da sezione 05)                              ##
##          Metacells_CTX.RData  (da sezione 06)                               ##
##                                                                             ##
##  Output: DEGs_SingleCell.RData                                              ##
##          DEGs_Metacell.RData                                                ##
##          DEG_concordance.RData                                              ##
##                                                                             ##
##  Per ogni area x timepoint x tipo cellulare:                               ##
##    - FindMarkers (Wilcoxon) su singole cellule                              ##
##    - FindMarkers (Wilcoxon) su metacelle                                    ##
##    - Concordanza SC intersect MC                                            ##
##    - Volcano e barplot DEG                                                  ##
##                                                                             ##
##  Test: Wilcoxon rank-sum (default Seurat)                                  ##
##    Pro: non parametrico, robusto a sparsita' e distribuzioni non normali   ##
##    Contro: non modella struttura gerarchica (animale come effetto random)  ##
##    Alternativa per pubblicazione: pseudo-bulk + DESeq2 (07_deseq2.R)      ##
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_celltype.RData"))   # -> CTX
load(file.path(outdir, "Metacells_CTX.RData"))  # -> SC_CTX, SC_CTX_seurat

# ==============================================================================
# PARAMETRI
# ==============================================================================

FDR_THR  <- 0.05   # soglia FDR (p_val_adj in FindMarkers)
LFC_THR  <- 0.25   # soglia |log2FC| (logfc.threshold in FindMarkers)
MIN_CELLS <- 20    # minimo cellule per gruppo per eseguire il test
MIN_PCT   <- 0.10  # gene espresso in almeno questa % di cellule (min.pct)

# ==============================================================================
# 07.1 — FUNZIONE run_findmarkers
# ==============================================================================
# Wrappa FindMarkers per iterare su ogni cell type x timepoint.
# Produce un data.frame con le stesse colonne usate nei plot/concordance a valle:
#   gene, log2FoldChange, padj, cell_type
# (FindMarkers chiama i valori avg_log2FC e p_val_adj — li rinominiamo)

run_findmarkers <- function(obj,
                            area_label,
                            time_val,
                            cell_type_col = "cell_type",
                            condition_col = "treatment",
                            assay_use     = "RNA",
                            min_cells     = MIN_CELLS) {
  
  # Subset per timepoint
  cells_t <- colnames(obj)[obj$time == time_val]
  if (length(cells_t) == 0) {
    message(sprintf("  [%s | %s] nessuna cellula", area_label, time_val))
    return(list())
  }
  obj_t <- obj[, cells_t]
  
  # Verifica che entrambe le condizioni siano presenti nel timepoint
  cond_all <- unlist(obj_t[[condition_col]])
  names(cond_all) <- colnames(obj_t)
  if (!all(c("Ctrl","LSD") %in% unique(cond_all))) {
    message(sprintf("  [%s | %s] condizioni mancanti", area_label, time_val))
    return(list())
  }
  
  # Metadati allineati ai barcode: usa match esplicito per evitare
  # disallineamenti tra colnames e slot dei metadati
  meta_t     <- obj_t@meta.data
  ct_vec     <- meta_t[[cell_type_col]]   # vettore allineato a rownames(meta_t)
  cond_vec   <- meta_t[[condition_col]]   # idem
  names(ct_vec)   <- rownames(meta_t)
  names(cond_vec) <- rownames(meta_t)
  
  cell_types <- sort(unique(ct_vec))
  degs_list  <- list()
  
  for (ct in cell_types) {
    
    # Barcode delle celle di questo CT (allineati tramite meta_t)
    cells_ct <- rownames(meta_t)[ct_vec == ct]
    # Verifica che tutti i barcode esistano nell'oggetto
    cells_ct <- intersect(cells_ct, colnames(obj_t))
    if (length(cells_ct) == 0) next
    
    # Condizione per ciascuna cella del CT
    cond_ct <- cond_vec[cells_ct]
    n_ctrl  <- sum(cond_ct == "Ctrl", na.rm = TRUE)
    n_lsd   <- sum(cond_ct == "LSD",  na.rm = TRUE)
    
    if (n_ctrl < min_cells || n_lsd < min_cells) {
      message(sprintf("    skip %s/%s/%s: Ctrl=%d LSD=%d (min=%d)",
                      area_label, time_val, ct,
                      n_ctrl, n_lsd, min_cells))
      next
    }
    
    tryCatch({
      
      # Subset al CT e imposta Idents = condizione sull'oggetto ridotto
      # (evita subscript out of bounds che si ha facendo subset su obj_t intero
      #  con barcode che potrebbero non corrispondere all'ordine interno)
      obj_ct <- obj_t[, cells_ct]
      Idents(obj_ct) <- condition_col
      
      res <- FindMarkers(
        object          = obj_ct,
        ident.1         = "LSD",       # positivo = up in LSD vs Ctrl
        ident.2         = "Ctrl",
        assay           = assay_use,
        slot            = "counts",    # counts raw per Wilcoxon
        test.use        = "wilcox",
        logfc.threshold = LFC_THR,
        min.pct         = MIN_PCT,
        min.cells.group = min_cells,
        only.pos        = FALSE,
        verbose         = FALSE
      )
      
      if (nrow(res) == 0) {
        message(sprintf("    %s | %s | %-20s Ctrl=%4d LSD=%4d | 0 geni testati",
                        area_label, time_val, ct, n_ctrl, n_lsd))
        next
      }
      
      # Rinomina colonne per compatibilita' con plot/concordance a valle
      # FindMarkers restituisce: p_val, avg_log2FC, pct.1, pct.2, p_val_adj
      res$gene           <- rownames(res)
      res$log2FoldChange <- res$avg_log2FC
      res$padj           <- res$p_val_adj
      res$cell_type      <- ct
      
      degs_list[[ct]] <- res
      
      n_sig <- sum(res$padj < FDR_THR & abs(res$log2FoldChange) > LFC_THR,
                   na.rm = TRUE)
      message(sprintf("    %s | %s | %-20s Ctrl=%4d LSD=%4d | FDR<%.2f & |LFC|>%.2f: %d geni",
                      area_label, time_val, ct,
                      n_ctrl, n_lsd, FDR_THR, LFC_THR, n_sig))
      
    }, error = function(e) {
      message(sprintf("    WARN %s/%s/%s: %s",
                      area_label, time_val, ct, conditionMessage(e)))
    })
  }
  
  degs_list
}

# ==============================================================================
# 07.2 — FindMarkers SINGOLA CELLULA
# ==============================================================================

message("=== 07.2 FindMarkers — singola cellula ===")

message("--- CTX 90min ---")
DEGs_SC_CTX_90 <- run_findmarkers(CTX, "CTX", "90min")

message("--- CTX 24h ---")
DEGs_SC_CTX_24 <- run_findmarkers(CTX, "CTX", "24h")

save(DEGs_SC_CTX_90, DEGs_SC_CTX_24,
     file = file.path(outdir, "DEGs_SingleCell_markers.RData"))
message("Salvato: DEGs_SingleCell_markers.RData")

# ==============================================================================
# 07.3 — FindMarkers METACELLE
# ==============================================================================

message("=== 07.3 FindMarkers — metacelle ===")

message("--- MC CTX 90min ---")
DEGs_MC_CTX_90 <- run_findmarkers(SC_CTX_seurat, "MC_CTX", "90min",
                                  min_cells = 3)

message("--- MC CTX 24h ---")
DEGs_MC_CTX_24 <- run_findmarkers(SC_CTX_seurat, "MC_CTX", "24h",
                                  min_cells = 3)

save(DEGs_MC_CTX_90, DEGs_MC_CTX_24,
     file = file.path(outdir, "DEGs_Metacell_markers.RData"))
message("Salvato: DEGs_Metacell.RData")

# ==============================================================================
# 07.4 — CONCORDANZA SC intersect MC
# ==============================================================================

message("=== 07.4 Concordanza SC vs MC ===")

concordance <- function(sc_list, mc_list, label,
                        fdr = FDR_THR, lfc = LFC_THR) {
  
  common_ct <- intersect(names(sc_list), names(mc_list))
  if (length(common_ct) == 0) {
    message(sprintf("  %s: nessun CT in comune tra SC e MC", label))
    return(NULL)
  }
  
  by_ct <- lapply(setNames(common_ct, common_ct), function(ct) {
    sc_d   <- sc_list[[ct]]
    mc_d   <- mc_list[[ct]]
    sc_sig <- sc_d$gene[!is.na(sc_d$padj) & sc_d$padj < fdr &
                          abs(sc_d$log2FoldChange) > lfc]
    mc_sig <- mc_d$gene[!is.na(mc_d$padj) & mc_d$padj < fdr &
                          abs(mc_d$log2FoldChange) > lfc]
    list(
      concordant = intersect(sc_sig, mc_sig),
      sc_only    = setdiff(sc_sig, mc_sig),
      mc_only    = setdiff(mc_sig, sc_sig),
      n_sc       = length(sc_sig),
      n_mc       = length(mc_sig),
      n_shared   = length(intersect(sc_sig, mc_sig))
    )
  })
  
  summ <- do.call(rbind, lapply(names(by_ct), function(ct) {
    sc_d <- sc_list[[ct]]
    mc_d <- mc_list[[ct]]
    n_union <- length(union(
      sc_d$gene[!is.na(sc_d$padj) & sc_d$padj < fdr],
      mc_d$gene[!is.na(mc_d$padj) & mc_d$padj < fdr]
    ))
    data.frame(
      cell_type = ct,
      n_sc      = by_ct[[ct]]$n_sc,
      n_mc      = by_ct[[ct]]$n_mc,
      n_shared  = by_ct[[ct]]$n_shared,
      jaccard   = by_ct[[ct]]$n_shared / max(1, n_union)
    )
  }))
  
  message(sprintf("  %s — media Jaccard: %.2f", label, mean(summ$jaccard)))
  print(summ)
  list(by_celltype = by_ct, summary = summ)
}

concordance_CTX_90 <- concordance(DEGs_SC_CTX_90, DEGs_MC_CTX_90, "CTX_90")
concordance_CTX_24 <- concordance(DEGs_SC_CTX_24, DEGs_MC_CTX_24, "CTX_24")

save(concordance_CTX_90, concordance_CTX_24,
     file = file.path(outdir, "DEG_concordance_markers.RData"))
message("Salvato: DEG_concordance.RData")

# ==============================================================================
# 07.5 — BARPLOT DEG UP/DOWN per CT
# ==============================================================================

message("=== 07.5 Plot DEG barplot ===")

plot_deg_barplot <- function(degs_list, title,
                             fdr = FDR_THR, lfc = LFC_THR) {
  df <- do.call(rbind, lapply(names(degs_list), function(ct) {
    d <- degs_list[[ct]]
    d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange), ]
    rbind(
      data.frame(cell_type = ct, direction = "Up",
                 n = sum(d$padj < fdr & d$log2FoldChange >  lfc)),
      data.frame(cell_type = ct, direction = "Down",
                 n = sum(d$padj < fdr & d$log2FoldChange < -lfc))
    )
  }))
  if (nrow(df) == 0 || all(df$n == 0)) return(NULL)
  
  df$n_signed  <- ifelse(df$direction == "Down", -df$n, df$n)
  total        <- tapply(df$n, df$cell_type, sum)
  df$cell_type <- factor(df$cell_type,
                         levels = names(sort(total, decreasing = TRUE)))
  
  ggplot(df, aes(x = cell_type, y = n_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#CB4335", "Down" = "#2471A3")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top") +
    labs(title    = title,
         subtitle = sprintf("FDR < %.2f | |LFC| > %.2f (Wilcoxon)", fdr, lfc),
         x = "", y = "N DEG (positivo = up in LSD, negativo = down in LSD)")
}

deg_sets <- list(
  SC_CTX_90 = DEGs_SC_CTX_90,
  SC_CTX_24 = DEGs_SC_CTX_24,
  MC_CTX_90 = DEGs_MC_CTX_90,
  MC_CTX_24 = DEGs_MC_CTX_24
)

for (nm in names(deg_sets)) {
  p <- plot_deg_barplot(deg_sets[[nm]], paste("DEG —", nm))
  if (!is.null(p))
    savepng(p, paste0("07_DEG_barplot_", nm, ".png"),
            width  = 2500,
            height = max(1800, length(deg_sets[[nm]]) * 90 + 600))
}

# ==============================================================================
# 07.6 — VOLCANO PLOT
# ==============================================================================

message("=== 07.6 Volcano ===")

plot_volcano <- function(df, title, fdr = FDR_THR, lfc = LFC_THR,
                         top_n = 15) {
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
    scale_color_manual(values = c("up"   = "#E74C3C",
                                  "down" = "#3498DB",
                                  "ns"   = "gray70")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(title = title,
         x     = "log2FC (LSD vs Ctrl, Wilcoxon)",
         y     = "-log10(FDR)")
}

for (nm in c("SC_CTX_90","SC_CTX_24","MC_CTX_90","MC_CTX_24")) {
  degs_l <- deg_sets[[nm]]
  if (length(degs_l) == 0) next
  
  n_sig  <- sapply(degs_l, function(d)
    sum(!is.na(d$padj) & d$padj < FDR_THR & abs(d$log2FoldChange) > LFC_THR))
  top_ct <- names(which.max(n_sig))
  if (length(top_ct) == 0) next
  
  p_v <- plot_volcano(degs_l[[top_ct]],
                      paste("Volcano", nm, "|", top_ct))
  savepng(p_v, paste0("07_Volcano_", nm, "_", top_ct, ".png"))
}

message("=== 07 completato ===")
message(sprintf("Output in: %s", outdir))