################################################################################
##  08_ieg_summary.R                                                           ##
##                                                                             ##
##  Input:  DEGs_SingleCell.RData  (da sezione 07)                            ##
##  Output: tabelle IEG a video + PNG heatmap e dotplot                       ##
##                                                                             ##
##  Per ogni area × timepoint:                                                 ##
##    - Conta DEG up/down per tipo cellulare                                  ##
##    - Identifica IEG upregolati per tipo cellulare                          ##
##    - Heatmap LFC IEG × CT (asterischi per FDR < 0.05)                     ##
##    - Dotplot: n IEG up + LFC max + FDR min                                ##
################################################################################

source("code/00_config.R")
outdir<-"results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "DEGs_SingleCell.RData"))
# Oggetti: DEGs_SC_CTX_90, DEGs_SC_CTX_24, DEGs_SC_THAL_90, DEGs_SC_THAL_24

# ==============================================================================
# PARAMETRI
# ==============================================================================

FDR_THR <- 0.05
LFC_THR <- 0.0

IEG_GENES <- c(
  "Fos","Fosb","Fosl1","Fosl2",
  "Jun","Junb","Jund",
  "Arc","Egr1","Egr2","Egr3",
  "Nr4a1","Nr4a2","Nr4a3",
  "Npas4","Bdnf","Ntrk2",
  "Dusp1","Dusp5","Dusp6",
  "Rheb","Homer1","Narp"
)

# ==============================================================================
# FUNZIONE: analisi IEG su una lista DEG
# ==============================================================================

analyze_ieg <- function(degs_list, label) {

  message(sprintf("\n=== IEG — %s ===", label))

  # --- 1. Conteggio DEG per CT -------------------------------------------
  deg_counts <- do.call(rbind, lapply(names(degs_list), function(ct) {
    df <- degs_list[[ct]]
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    data.frame(
      cell_type = ct,
      n_tested  = nrow(df),
      n_up      = sum(df$padj < FDR_THR & df$log2FoldChange >  LFC_THR),
      n_down    = sum(df$padj < FDR_THR & df$log2FoldChange < -LFC_THR),
      n_sig     = sum(df$padj < FDR_THR & abs(df$log2FoldChange) > LFC_THR)
    )
  })) %>% arrange(desc(n_sig))

  cat(sprintf("\n--- DEG %s (FDR<%.2f, |LFC|>%.2f) ---\n",
              label, FDR_THR, LFC_THR))
  print(deg_counts)

  # --- 2. IEG upregolati per CT -------------------------------------------
  ieg_results <- do.call(rbind, lapply(names(degs_list), function(ct) {
    df <- degs_list[[ct]]
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    ieg_up <- df[df$gene %in% IEG_GENES &
                   df$padj < FDR_THR &
                   df$log2FoldChange > LFC_THR, ]
    ieg_up <- ieg_up[order(ieg_up$padj), ]
    data.frame(
      cell_type = ct,
      n_ieg_up  = nrow(ieg_up),
      ieg_genes = if (nrow(ieg_up) > 0) paste(ieg_up$gene, collapse=", ") else "",
      max_lfc   = if (nrow(ieg_up) > 0) round(max(ieg_up$log2FoldChange), 2) else NA_real_,
      min_padj  = if (nrow(ieg_up) > 0) signif(min(ieg_up$padj), 2)         else NA_real_
    )
  })) %>% arrange(desc(n_ieg_up))

  cat(sprintf("\n--- IEG upregolati %s ---\n", label))
  print(ieg_results)
  cat(sprintf("\n%d/%d CT con almeno 1 IEG upregolato\n",
              sum(ieg_results$n_ieg_up > 0), nrow(ieg_results)))

  # --- 3. Tabella dettagliata LFC × CT × gene -----------------------------
  ieg_detail <- do.call(rbind, lapply(names(degs_list), function(ct) {
    df <- degs_list[[ct]]
    df <- df[df$gene %in% IEG_GENES & !is.na(df$padj), ]
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(cell_type = ct,
             sig_up    = padj < FDR_THR & log2FoldChange > LFC_THR) %>%
      select(cell_type, gene, log2FoldChange, padj, sig_up)
  }))

  cat(sprintf("\n--- Dettaglio IEG sig up %s ---\n", label))
  if (!is.null(ieg_detail) && any(ieg_detail$sig_up, na.rm = TRUE)) {
    print(ieg_detail[ieg_detail$sig_up, ] %>%
            arrange(cell_type, desc(log2FoldChange)))
  } else {
    cat("  Nessun IEG significativamente upregolato.\n")
  }

  # --- 4. Barplot DEG up/down ----------------------------------------------
  deg_long <- deg_counts %>%
    mutate(n_down_signed = -n_down) %>%
    select(cell_type, n_up, n_down_signed) %>%
    pivot_longer(c(n_up, n_down_signed),
                 names_to  = "direction",
                 values_to = "n") %>%
    mutate(direction = recode(direction,
                              n_up          = "Up (LSD > Ctrl)",
                              n_down_signed = "Down (LSD < Ctrl)"))

  p_bar <- ggplot(deg_long,
                  aes(x = reorder(cell_type, abs(n)),
                      y = n, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Up (LSD > Ctrl)"   = "#CB4335",
                                 "Down (LSD < Ctrl)" = "#2471A3"),
                      name = "") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top") +
    labs(title    = paste("DEG per CT —", label),
         subtitle = sprintf("FDR < %.2f | |LFC| > %.2f", FDR_THR, LFC_THR),
         x = "", y = "N geni")
  savepng(p_bar,
          paste0("08_DEG_barplot_", label, ".png"),
          width = 2500, height = max(1800, nrow(deg_counts) * 80 + 600))

  # --- 5. Heatmap IEG LFC × CT -------------------------------------------
  if (!is.null(ieg_detail) && nrow(ieg_detail) > 0) {

    ieg_mat_df <- ieg_detail %>%
      select(cell_type, gene, log2FoldChange) %>%
      pivot_wider(names_from  = gene,
                  values_from = log2FoldChange,
                  values_fill = 0)

    ieg_mat <- as.matrix(ieg_mat_df[, -1])
    rownames(ieg_mat) <- ieg_mat_df$cell_type

    # Ordina CT per n IEG up
    ct_order <- ieg_results$cell_type[order(-ieg_results$n_ieg_up)]
    ct_order <- ct_order[ct_order %in% rownames(ieg_mat)]
    ieg_mat  <- ieg_mat[ct_order, , drop = FALSE]

    # Maschera asterischi significatività
    sig_mask <- matrix("", nrow = nrow(ieg_mat), ncol = ncol(ieg_mat),
                       dimnames = dimnames(ieg_mat))
    for (i in seq_len(nrow(ieg_mat))) {
      ct_i <- rownames(ieg_mat)[i]
      for (j in seq_len(ncol(ieg_mat))) {
        g <- colnames(ieg_mat)[j]
        row_ij <- ieg_detail[ieg_detail$cell_type == ct_i &
                                ieg_detail$gene == g, ]
        if (nrow(row_ij) > 0 && !is.na(row_ij$sig_up[1]) &&
            row_ij$sig_up[1]) sig_mask[i, j] <- "*"
      }
    }

    lim     <- max(abs(ieg_mat), na.rm = TRUE)
    col_fun <- colorRampPalette(c("#2471A3","white","#CB4335"))(100)
    breaks  <- seq(-lim, lim, length.out = 101)

    graphics.off()
    png(file.path(outdir, paste0("08_IEG_heatmap_", label, ".png")),
        width  = max(2800, ncol(ieg_mat) * 120 + 800),
        height = max(2000, nrow(ieg_mat) * 90  + 600),
        res    = 300)
    pheatmap::pheatmap(
      ieg_mat,
      color           = col_fun,
      breaks          = breaks,
      cluster_rows    = FALSE,
      cluster_cols    = TRUE,
      display_numbers = sig_mask,
      number_color    = "black",
      fontsize_number = 11,
      fontsize_row    = 9,
      fontsize_col    = 9,
      angle_col       = 45,
      main            = paste0("IEG log2FC (LSD vs Ctrl) — ", label,
                               "\n* FDR < ", FDR_THR, " & LFC > ", LFC_THR)
    )
    dev.off()
    message(sprintf("  Salvato: 08_IEG_heatmap_%s.png", label))

    # --- 6. Dotplot: n IEG up + max LFC + min FDR -------------------------
    ieg_pos <- ieg_results %>% filter(n_ieg_up > 0)
    if (nrow(ieg_pos) > 0) {
      p_dot <- ggplot(ieg_pos,
                      aes(x     = reorder(cell_type, n_ieg_up),
                          y     = n_ieg_up,
                          size  = max_lfc,
                          color = -log10(min_padj + 1e-300))) +
        geom_point() +
        coord_flip() +
        scale_color_gradient(low = "gray70", high = "#CB4335",
                             name = "-log10(FDR min)") +
        scale_size_continuous(range = c(3, 10), name = "LFC max") +
        theme_bw(base_size = 10) +
        labs(title    = paste("IEG upregolati —", label),
             subtitle = "Dimensione = LFC max | Colore = -log10(FDR min)",
             x = "", y = "N IEG upregolati")
      savepng(p_dot,
              paste0("08_IEG_dotplot_", label, ".png"),
              width = 2800, height = max(1800, nrow(ieg_pos) * 100 + 600))
      message(sprintf("  Salvato: 08_IEG_dotplot_%s.png", label))
    }
  }

  invisible(list(deg_counts = deg_counts,
                 ieg_results = ieg_results,
                 ieg_detail  = ieg_detail))
}

# ==============================================================================
# ESECUZIONE
# ==============================================================================

ieg_CTX_90  <- analyze_ieg(DEGs_SC_CTX_90,  "CTX_90min")
ieg_CTX_24  <- analyze_ieg(DEGs_SC_CTX_24,  "CTX_24h")

message(sprintf("\n=== 08 completato — output in: %s ===", outdir))
