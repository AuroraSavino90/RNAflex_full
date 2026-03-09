################################################################################
##  07_deseq2_THAL.R                                                           ##
##                                                                             ##
##  Input:  THAL/THAL_celltype.RData       (da 05_cell_type_annotation_THAL) ##
##          THAL/Metacells_THAL.RData      (da 06_metacells_THAL)             ##
##  Output: THAL/DEGs_SingleCell_THAL.RData                                   ##
##          THAL/DEGs_Metacell_THAL.RData                                     ##
##                                                                             ##
##  Per ogni timepoint × tipo cellulare talamica:                             ##
##    - DESeq2 su singole cellule (pseudo-bulk intra-CT)                      ##
##    - DESeq2 su metacelle                                                    ##
##    - Concordanza SC ∩ MC                                                    ##
##    - Volcano e barplot DEG                                                  ##
################################################################################

source("code/00_config.R")

outdir   <- "results/20260224"
thal_dir <- file.path(outdir, "THAL")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(thal_dir, "THAL_celltype.RData"))       # THAL
load(file.path(thal_dir, "Metacells_THAL.RData"))      # SC_THAL, SC_THAL_seurat

# ==============================================================================
# PARAMETRI DESeq2
# ==============================================================================

FDR_THR              <- 0.05
LFC_THR              <- 0.25
MIN_CELLS            <- 20
MIN_COUNTS           <- 10
MIN_CELLS_EXPRESSED  <- 5

# ==============================================================================
# 07.1 — FUNZIONE DESeq2
# ==============================================================================

lfc_type <- if (requireNamespace("apeglm", quietly = TRUE)) {
  "apeglm"
} else if (requireNamespace("ashr", quietly = TRUE)) {
  message("INFO: apeglm non disponibile, uso ashr")
  "ashr"
} else {
  message("INFO: uso normal per lfcShrink")
  "normal"
}
message(sprintf("LFC shrinkage: %s", lfc_type))

run_deseq2 <- function(obj, area_label, time_val,
                       cell_type_col = "cell_type",
                       condition_col = "treatment",
                       assay_use     = "RNA",
                       min_cells     = MIN_CELLS) {

  cells_t <- colnames(obj)[obj$time == time_val]
  if (length(cells_t) == 0) {
    message(sprintf("  [%s | %s] nessuna cellula", area_label, time_val))
    return(list())
  }
  obj_t <- obj[, cells_t]

  cell_types <- sort(unique(unlist(obj_t[[cell_type_col]])))
  degs_list  <- list()

  for (ct in cell_types) {
    cells_ct <- colnames(obj_t)[unlist(obj_t[[cell_type_col]]) == ct]
    if (length(cells_ct) < min_cells) {
      message(sprintf("    skip %s/%s/%s: solo %d cellule (min=%d)",
                      area_label, time_val, ct, length(cells_ct), min_cells))
      next
    }

    obj_ct    <- obj_t[, cells_ct]
    count_mat <- LayerData(obj_ct, assay = assay_use, layer = "counts")
    count_mat <- as(count_mat, "dgCMatrix")
    count_mat <- count_mat[Matrix::rowSums(count_mat > 0) >= MIN_CELLS_EXPRESSED, ]

    if (nrow(count_mat) < 100) {
      message(sprintf("    skip %s/%s/%s: solo %d geni dopo filtro",
                      area_label, time_val, ct, nrow(count_mat)))
      next
    }

    cond     <- factor(unlist(obj_ct[[condition_col]]), levels = c("Ctrl","LSD"))
    cond_tab <- table(cond)

    if (any(cond_tab == 0) || length(unique(cond)) < 2) {
      message(sprintf("    skip %s/%s/%s: condizione mancante (%s)",
                      area_label, time_val, ct,
                      paste(names(cond_tab[cond_tab == 0]), collapse = ",")))
      next
    }

    MIN_MC_PER_COND <- 3
    if (min(cond_tab) < MIN_MC_PER_COND) {
      message(sprintf(
        "    skip %s/%s/%s: troppo poche unità per condizione (min=%d, Ctrl=%d, LSD=%d)",
        area_label, time_val, ct, MIN_MC_PER_COND,
        cond_tab["Ctrl"], cond_tab["LSD"]))
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
        message(sprintf("    [%s/%s/%s] local fallito, provo mean", area_label, time_val, ct))
        dds
      })

      if (!fit_ok) {
        dds <- tryCatch({
          d <- estimateDispersions(dds, fitType = "mean", quiet = TRUE)
          fit_ok <<- TRUE; d
        }, error = function(e) {
          message(sprintf("    [%s/%s/%s] mean fallito, provo gene-wise", area_label, time_val, ct))
          dds
        })
      }

      if (!fit_ok) {
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        message(sprintf("    [%s/%s/%s] dispersioni gene-wise (fallback finale)", area_label, time_val, ct))
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
            message(sprintf("    [%s/%s/%s] spread <2 OOM: gene-wise", area_label, time_val, ct))
          }
        }
      }

      dds <- nbinomWaldTest(dds, quiet = TRUE)

      res <- tryCatch({
        if (lfc_type %in% c("apeglm","ashr")) {
          lfcShrink(dds, coef = "condition_LSD_vs_Ctrl", type = lfc_type, quiet = TRUE)
        } else {
          lfcShrink(dds, contrast = c("condition","LSD","Ctrl"), type = "normal", quiet = TRUE)
        }
      }, error = function(e) {
        message(sprintf("    [%s/%s/%s] lfcShrink fallito, uso results() raw",
                        area_label, time_val, ct))
        results(dds, contrast = c("condition","LSD","Ctrl"))
      })

      df           <- as.data.frame(res)
      df$gene      <- rownames(df)
      df$cell_type <- ct
      degs_list[[ct]] <- df

      n_sig <- sum(!is.na(df$padj) & df$padj < FDR_THR & abs(df$log2FoldChange) > LFC_THR)
      message(sprintf("    %s | %s | %-20s n=%4d | FDR<%.2f & |LFC|>%.2f: %d DEG",
                      area_label, time_val, ct, length(cells_ct),
                      FDR_THR, LFC_THR, n_sig))

    }, error = function(e) {
      message(sprintf("    WARN %s/%s/%s: %s",
                      area_label, time_val, ct,
                      gsub("\n", " ", conditionMessage(e))))
    })
  }
  degs_list
}

# ==============================================================================
# 07.2 — DESeq2 SINGOLA CELLULA — THAL
# ==============================================================================

message("=== 07.2 DESeq2 — singola cellula THAL ===")

message("--- THAL 90min ---")
DEGs_SC_THAL_90 <- run_deseq2(THAL, "THAL", "90min")

message("--- THAL 24h ---")
DEGs_SC_THAL_24 <- run_deseq2(THAL, "THAL", "24h")

save(DEGs_SC_THAL_90, DEGs_SC_THAL_24,
     file = file.path(thal_dir, "DEGs_SingleCell_THAL.RData"))
message("Salvato: THAL/DEGs_SingleCell_THAL.RData")

# ==============================================================================
# 07.3 — DESeq2 METACELLE — THAL
# ==============================================================================

message("=== 07.3 DESeq2 — metacelle THAL ===")
MIN_CELLS <- 3

message("--- MC THAL 90min ---")
DEGs_MC_THAL_90 <- run_deseq2(SC_THAL_seurat, "MC_THAL", "90min")

message("--- MC THAL 24h ---")
DEGs_MC_THAL_24 <- run_deseq2(SC_THAL_seurat, "MC_THAL", "24h")

save(DEGs_MC_THAL_90, DEGs_MC_THAL_24,
     file = file.path(thal_dir, "DEGs_Metacell_THAL.RData"))
message("Salvato: THAL/DEGs_Metacell_THAL.RData")

# ==============================================================================
# 07.4 — CONCORDANZA SC ∩ MC — THAL
# ==============================================================================

message("=== 07.4 Concordanza SC vs MC THAL ===")

concordance <- function(sc_list, mc_list, label, fdr = FDR_THR, lfc = LFC_THR) {
  common_ct <- intersect(names(sc_list), names(mc_list))
  if (length(common_ct) == 0) {
    message(sprintf("  %s: nessun CT in comune tra SC e MC", label))
    return(NULL)
  }

  by_ct <- lapply(setNames(common_ct, common_ct), function(ct) {
    sc_sig <- sc_list[[ct]] %>%
      filter(!is.na(padj), padj < fdr, abs(log2FoldChange) > lfc) %>%
      pull(gene)
    mc_sig <- mc_list[[ct]] %>%
      filter(!is.na(padj), padj < fdr, abs(log2FoldChange) > lfc) %>%
      pull(gene)
    list(
      concordant = intersect(sc_sig, mc_sig),
      sc_only    = setdiff(sc_sig, mc_sig),
      mc_only    = setdiff(mc_sig, sc_sig),
      n_sc       = length(sc_sig),
      n_mc       = length(mc_sig),
      n_shared   = length(intersect(sc_sig, mc_sig))
    )
  })

  summ <- do.call(rbind, lapply(names(by_ct), function(ct)
    data.frame(cell_type = ct,
               n_sc      = by_ct[[ct]]$n_sc,
               n_mc      = by_ct[[ct]]$n_mc,
               n_shared  = by_ct[[ct]]$n_shared,
               jaccard   = by_ct[[ct]]$n_shared /
                 max(1, length(union(
                   sc_list[[ct]]$gene[!is.na(sc_list[[ct]]$padj) & sc_list[[ct]]$padj < fdr],
                   mc_list[[ct]]$gene[!is.na(mc_list[[ct]]$padj) & mc_list[[ct]]$padj < fdr]
                 ))))))

  message(sprintf("  %s — media Jaccard: %.2f", label, mean(summ$jaccard)))
  print(summ)
  list(by_celltype = by_ct, summary = summ)
}

concordance_THAL_90 <- concordance(DEGs_SC_THAL_90, DEGs_MC_THAL_90, "THAL_90")
concordance_THAL_24 <- concordance(DEGs_SC_THAL_24, DEGs_MC_THAL_24, "THAL_24")

# ==============================================================================
# 07.5 — PLOT BARPLOT DEG UP/DOWN per CT
# ==============================================================================

message("=== 07.5 Plot DEG barplot THAL ===")

plot_deg_barplot <- function(degs_list, title, fdr = FDR_THR, lfc = LFC_THR) {
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
  df$cell_type <- factor(df$cell_type, levels = names(sort(total, decreasing = TRUE)))

  ggplot(df, aes(x = cell_type, y = n_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#CB4335", "Down" = "#2471A3")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top") +
    labs(title    = title,
         subtitle = sprintf("FDR < %.2f | |LFC| > %.2f", fdr, lfc),
         x = "", y = "N DEG (positivo = up, negativo = down)")
}

deg_sets_thal <- list(
  SC_THAL_90 = DEGs_SC_THAL_90,
  SC_THAL_24 = DEGs_SC_THAL_24
)

for (nm in names(deg_sets_thal)) {
  p <- plot_deg_barplot(deg_sets_thal[[nm]], paste("DEG —", nm))
  if (!is.null(p))
    savepng(p, file.path(thal_dir, paste0("07_DEG_barplot_", nm, ".png")),
            width = 2500,
            height = max(1800, length(deg_sets_thal[[nm]]) * 90 + 600))
}

# ==============================================================================
# 07.6 — VOLCANO PLOT
# ==============================================================================

message("=== 07.6 Volcano THAL ===")

plot_volcano <- function(df, title, fdr = FDR_THR, lfc = LFC_THR, top_n = 15) {
  df <- df[!is.na(df$padj), ]
  df$color <- ifelse(df$padj < fdr & df$log2FoldChange >  lfc, "up",
              ifelse(df$padj < fdr & df$log2FoldChange < -lfc, "down", "ns"))
  top <- df[df$color != "ns", ]
  top <- top[order(top$padj), ][seq_len(min(top_n, nrow(top))), ]

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj + 1e-300), color = color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(fdr),  linetype = "dashed", color = "gray50") +
    ggrepel::geom_text_repel(data = top, aes(label = gene),
                              size = 2.8, max.overlaps = 20, box.padding = 0.4) +
    scale_color_manual(values = c("up" = "#E74C3C", "down" = "#3498DB", "ns" = "gray70")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(title = title,
         x = paste0("log2FC (", lfc_type, ")"),
         y = "-log10(FDR)")
}

for (nm in c("SC_THAL_90","SC_THAL_24")) {
  degs_l <- deg_sets_thal[[nm]]
  if (length(degs_l) == 0) next

  n_sig  <- sapply(degs_l, function(d)
    sum(!is.na(d$padj) & d$padj < FDR_THR & abs(d$log2FoldChange) > LFC_THR))
  top_ct <- names(which.max(n_sig))
  if (length(top_ct) == 0) next

  p_v <- plot_volcano(degs_l[[top_ct]], paste("Volcano", nm, "|", top_ct))
  savepng(p_v, file.path(thal_dir, paste0("07_Volcano_", nm, "_", top_ct, ".png")))
}

message("=== 07 THAL completato ===")
message(sprintf("Output in: %s", thal_dir))
