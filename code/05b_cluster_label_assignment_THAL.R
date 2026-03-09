################################################################################
##  05b_cluster_label_assignment_THAL.R                                       ##
##                                                                             ##
##  Input:  THAL/THAL_celltype.RData  (da 05_cell_type_annotation_THAL)      ##
##          MAPPING_CSV_THAL  (CSV gerarchico da MapMyCells per THAL)         ##
##                                                                             ##
##  Output: THAL/THAL_cluster_labels.csv                                      ##
##    THAL/05b_1_crosstab_THAL.png   heatmap cluster x mmc_roi  (proporzioni)##
##    THAL/05b_2_crosstab_THAL.png   heatmap cluster x mmc_class (proporzioni)##
##    THAL/05b_3_reliability_THAL.png  lollipop reliability per cluster       ##
##    THAL/05b_4_scatter_THAL.png      scatter purity_roi vs purity_class     ##
##                                                                             ##
##  Reliability score (0-1):                                                  ##
##    purity_roi   = % cellule con mmc_roi dominante                          ##
##    purity_class = % cellule con mmc_class dominante                        ##
##    concordance  = % cellule con ENTRAMBI i label dominanti                 ##
##    conf_median  = bootstrapping probability mediana dal CSV MMC            ##
##    score        = W_ROI*purity_roi + W_CLASS*purity_class +               ##
##                   W_CONC*concordance                                       ##
##    grade        = HIGH (>=80%) / MEDIUM / LOW (<50%) su concordance       ##
##                                                                             ##
##  Nota: usa mmc_roi (ROI talamica da MapMyCells) al posto di broad_anat,   ##
##  coerentemente con la pipeline THAL che non usa broad_anat.               ##
################################################################################

source("code/00_config.R")

outdir   <- "results/20260224"
thal_dir <- file.path(outdir, "THAL")
dir.create(thal_dir, recursive = TRUE, showWarnings = FALSE)

savepng_thal <- function(p, filename, width = 3000, height = 2500, res = 300) {
  png(file.path(thal_dir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# ==============================================================================
# PARAMETRI
# ==============================================================================

MAPPING_CSV_THAL <- file.path(outdir,
  "THAL_from_mapmycells/THAL_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC.csv")

# Classi attese per il talamo nell'atlas ABCA
EXPECTED_CLASS_THAL <- c("267 Glutamatergic", "267 GABAergic",
                          "Astrocyte", "Oligodendrocyte", "OPC",
                          "Microglia", "Vascular")

# ROI talamiche attese (label ABCA)
EXPECTED_ROI_THAL <- c("TH", "EPI")   # Thalamus + Epithalamus (habenula)

THR_PURITY_HIGH  <- 0.80
THR_PURITY_LOW   <- 0.50
THR_CONCORD_WARN <- 0.60
THR_CONF_LOW     <- 0.70

W_ROI   <- 0.30
W_CLASS <- 0.30
W_CONC  <- 0.40

grade_pal <- c(HIGH = "#27AE60", MEDIUM = "#F39C12", LOW = "#E74C3C")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(thal_dir, "THAL_celltype.RData"))   # THAL

if (file.exists(MAPPING_CSV_THAL)) {
  mapping_raw <- read.csv(MAPPING_CSV_THAL, comment.char = "#")
  rownames(mapping_raw) <- mapping_raw$cell_id
  has_mapping <- TRUE
  message(sprintf("  CSV MMC caricato: %d cellule", nrow(mapping_raw)))
} else {
  has_mapping <- FALSE
  message("  WARN: MAPPING_CSV_THAL non trovato — conf_median sarà NA")
}

# ==============================================================================
# METADATI
# ==============================================================================

meta         <- THAL@meta.data
meta$cluster <- as.character(Idents(THAL))

if (has_mapping) {
  conf_col <- grep("class.*bootstrapping_probability",
                   colnames(mapping_raw), value = TRUE,
                   ignore.case = TRUE)[1]
  shared <- intersect(rownames(meta), rownames(mapping_raw))
  meta$conf_class <- NA_real_
  meta[shared, "conf_class"] <- as.numeric(mapping_raw[shared, conf_col])
  message(sprintf("  Confidence: '%s' | mediana globale = %.2f",
                  conf_col, median(meta$conf_class, na.rm = TRUE)))
}

cluster_order <- as.character(sort(unique(meta$cell_type)))

# ==============================================================================
# FUNZIONI HELPER
# ==============================================================================

shannon_entropy <- function(tab) {
  tab <- tab[tab > 0]
  if (length(tab) <= 1) return(0)
  p <- tab / sum(tab)
  -sum(p * log2(p))
}
top1      <- function(v) { tab <- sort(table(v[!is.na(v)]), decreasing=TRUE); if (!length(tab)) return(NA_character_); names(tab)[1] }
top1_pct  <- function(v) { tab <- table(v[!is.na(v)]); if (!length(tab)) return(NA_real_); round(max(tab)/sum(tab),4) }
top2      <- function(v) { tab <- sort(table(v[!is.na(v)]), decreasing=TRUE); if (length(tab)<2) return(NA_character_); names(tab)[2] }
top2_pct  <- function(v) { tab <- sort(table(v[!is.na(v)]), decreasing=TRUE); if (length(tab)<2) return(0); round(as.numeric(tab[2])/sum(tab),4) }

# ==============================================================================
# CALCOLO METRICHE PER CLUSTER
# ==============================================================================

message("=== 05b. Label assignment THAL ===")

rows <- lapply(cluster_order, function(cl) {

  m  <- meta[meta$cell_type == cl, , drop = FALSE]
  nc <- nrow(m)

  # Usa mmc_roi al posto di broad_anat
  dom_roi   <- top1(m$mmc_roi)
  dom_class <- top1(m$mmc_class)
  pur_roi   <- top1_pct(m$mmc_roi)
  pur_class <- top1_pct(m$mmc_class)

  ent_roi   <- round(shannon_entropy(table(na.omit(m$mmc_roi))),   3)
  ent_class <- round(shannon_entropy(table(na.omit(m$mmc_class))), 3)

  rup_roi      <- top2(m$mmc_roi)
  rup_roi_pct  <- top2_pct(m$mmc_roi)
  rup_class    <- top2(m$mmc_class)
  rup_class_pct <- top2_pct(m$mmc_class)

  n_roi   <- length(unique(na.omit(m$mmc_roi)))
  n_class <- length(unique(na.omit(m$mmc_class)))

  # Concordance: % cellule con ENTRAMBI mmc_roi=dom E mmc_class=dom
  valid   <- !is.na(m$mmc_roi) & !is.na(m$mmc_class)
  n_valid <- sum(valid)
  concordance <- if (n_valid > 0 && !is.na(dom_roi) && !is.na(dom_class))
    round(sum(m$mmc_roi[valid]   == dom_roi &
                m$mmc_class[valid] == dom_class) / n_valid, 4)
  else NA_real_

  # Mismatch biologico: ROI attesa TH/EPI ma classe inattesa
  roi_ok   <- !is.na(dom_roi)   && dom_roi   %in% EXPECTED_ROI_THAL
  class_ok <- !is.na(dom_class) && dom_class %in% EXPECTED_CLASS_THAL
  mismatch <- roi_ok && !class_ok

  # Confidence bootstrapping
  cv           <- m$conf_class[!is.na(m$conf_class)]
  conf_median  <- if (length(cv)) round(median(cv), 3)              else NA_real_
  conf_pct_low <- if (length(cv)) round(mean(cv < THR_CONF_LOW), 3) else NA_real_

  # Score composito
  score <- round(
    W_ROI   * ifelse(is.na(pur_roi),     0, pur_roi)   +
    W_CLASS * ifelse(is.na(pur_class),   0, pur_class)  +
    W_CONC  * ifelse(is.na(concordance), 0, concordance),
    4)

  base_g <- if (!is.na(concordance)) concordance else
    mean(c(pur_roi, pur_class), na.rm = TRUE)
  grade <- if (is.na(base_g))          NA_character_ else
    if (base_g >= THR_PURITY_HIGH) "HIGH"    else
    if (base_g <  THR_PURITY_LOW)  "LOW"     else
    "MEDIUM"

  flags <- character(0)
  if (!is.na(pur_roi)     && pur_roi     < THR_PURITY_LOW)   flags <- c(flags, "MIXED_ROI")
  if (!is.na(pur_class)   && pur_class   < THR_PURITY_LOW)   flags <- c(flags, "MIXED_CLASS")
  if (!is.na(concordance) && concordance < THR_CONCORD_WARN) flags <- c(flags, "LOW_CONCORDANCE")
  if (mismatch)                                               flags <- c(flags, "ROI_CLASS_MISMATCH")
  if (!is.na(conf_median) && conf_median < THR_CONF_LOW)     flags <- c(flags, "LOW_CONF")

  list(
    cell_type             = cl,
    n_cells               = nc,
    mmc_roi_label         = dom_roi,
    mmc_roi_purity        = pur_roi,
    mmc_roi_entropy       = ent_roi,
    mmc_roi_n             = n_roi,
    mmc_roi_runner_up     = rup_roi,
    mmc_roi_runner_up_pct = rup_roi_pct,
    mmc_class_label       = dom_class,
    mmc_class_purity      = pur_class,
    mmc_class_entropy     = ent_class,
    mmc_class_n           = n_class,
    mmc_class_runner_up   = rup_class,
    mmc_class_runner_up_pct = rup_class_pct,
    concordance           = concordance,
    conf_median           = conf_median,
    conf_pct_low          = conf_pct_low,
    reliability_score     = score,
    reliability_grade     = grade,
    warning_flags         = paste(flags, collapse = "|")
  )
})

label_df <- data.frame(
  cell_type               = sapply(rows, `[[`, "cell_type"),
  n_cells                 = as.integer(sapply(rows, `[[`, "n_cells")),
  mmc_roi_label           = sapply(rows, `[[`, "mmc_roi_label"),
  mmc_roi_purity          = as.numeric(sapply(rows, `[[`, "mmc_roi_purity")),
  mmc_roi_entropy         = as.numeric(sapply(rows, `[[`, "mmc_roi_entropy")),
  mmc_roi_n               = as.integer(sapply(rows, `[[`, "mmc_roi_n")),
  mmc_roi_runner_up       = sapply(rows, `[[`, "mmc_roi_runner_up"),
  mmc_roi_runner_up_pct   = as.numeric(sapply(rows, `[[`, "mmc_roi_runner_up_pct")),
  mmc_class_label         = sapply(rows, `[[`, "mmc_class_label"),
  mmc_class_purity        = as.numeric(sapply(rows, `[[`, "mmc_class_purity")),
  mmc_class_entropy       = as.numeric(sapply(rows, `[[`, "mmc_class_entropy")),
  mmc_class_n             = as.integer(sapply(rows, `[[`, "mmc_class_n")),
  mmc_class_runner_up     = sapply(rows, `[[`, "mmc_class_runner_up"),
  mmc_class_runner_up_pct = as.numeric(sapply(rows, `[[`, "mmc_class_runner_up_pct")),
  concordance             = as.numeric(sapply(rows, `[[`, "concordance")),
  conf_median             = as.numeric(sapply(rows, `[[`, "conf_median")),
  conf_pct_low            = as.numeric(sapply(rows, `[[`, "conf_pct_low")),
  reliability_score       = as.numeric(sapply(rows, `[[`, "reliability_score")),
  reliability_grade       = sapply(rows, `[[`, "reliability_grade"),
  warning_flags           = sapply(rows, `[[`, "warning_flags"),
  stringsAsFactors        = FALSE
)

message(sprintf("  Cluster: %d | concordance mediana: %.2f",
                nrow(label_df), median(label_df$concordance, na.rm = TRUE)))
for (g in c("HIGH", "MEDIUM", "LOW")) {
  n <- sum(label_df$reliability_grade == g, na.rm = TRUE)
  message(sprintf("    %6s: %d (%.0f%%)", g, n, 100 * n / nrow(label_df)))
}

write.csv(label_df,
          file.path(thal_dir, "THAL_cluster_labels.csv"),
          row.names = FALSE)
message("  Salvato: THAL/THAL_cluster_labels.csv")

# ==============================================================================
# PLOT 1 & 2 — HEATMAP: cluster x mmc_roi / cluster x mmc_class
# ==============================================================================

message("=== 05b Plot 1-2: heatmap proporzionali THAL ===")

make_crosstab_heatmap <- function(col, filename, plot_title) {

  raw_tab  <- table(meta$cell_type, meta[[col]])
  ok_cols  <- colnames(raw_tab)[colnames(raw_tab) != "NA" & !is.na(colnames(raw_tab))]
  raw_tab  <- raw_tab[, ok_cols, drop = FALSE]
  prop_mat <- sweep(as.matrix(raw_tab), 1, rowSums(as.matrix(raw_tab)), "/")

  row_ord  <- cluster_order[cluster_order %in% rownames(prop_mat)]
  col_ord  <- names(sort(colSums(prop_mat), decreasing = TRUE))
  prop_mat <- prop_mat[row_ord, col_ord, drop = FALSE]

  display_mat <- matrix(
    ifelse(prop_mat >= 0.05, as.character(round(prop_mat * 100)), ""),
    nrow = nrow(prop_mat), ncol = ncol(prop_mat),
    dimnames = dimnames(prop_mat)
  )

  idx     <- match(row_ord, label_df$cell_type)
  row_ann <- data.frame(
    Grade   = factor(label_df$reliability_grade[idx], levels = c("HIGH","MEDIUM","LOW")),
    Warning = ifelse(is.na(label_df$warning_flags[idx]) |
                       label_df$warning_flags[idx] == "", "no", "si"),
    row.names = row_ord, stringsAsFactors = FALSE
  )

  ann_colors <- list(
    Grade   = c(HIGH = "#27AE60", MEDIUM = "#F39C12", LOW = "#E74C3C"),
    Warning = c(no = "grey90",    si     = "#E74C3C")
  )

  graphics.off()
  png(file.path(thal_dir, filename),
      width  = max(3000, ncol(prop_mat) * 220 + 1200),
      height = max(2000, nrow(prop_mat) * 80  + 800),
      res    = 300)
  print(pheatmap(
    prop_mat,
    color             = colorRampPalette(c("grey97","#2980B9","#1A5276"))(80),
    breaks            = c(seq(0, 0.99, length.out = 80), 1),
    display_numbers   = display_mat,
    number_color      = "white",
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_row    = row_ann,
    annotation_colors = ann_colors,
    fontsize          = 8, fontsize_row = 7, fontsize_col = 8,
    fontsize_number   = 7, angle_col = "45", border_color = "white",
    main              = plot_title
  ))
  dev.off()
  message(sprintf("  Salvato: %s", filename))
}

make_crosstab_heatmap(
  col        = "mmc_roi",
  filename   = "05b_1_crosstab_THAL.png",
  plot_title = "Proporzione mmc_roi per cluster - THAL\n(colore = % celle nel cluster)"
)

make_crosstab_heatmap(
  col        = "mmc_class",
  filename   = "05b_2_crosstab_THAL.png",
  plot_title = "Proporzione mmc_class per cluster - THAL\n(colore = % celle nel cluster)"
)

# ==============================================================================
# PLOT 3 — RELIABILITY LOLLIPOP
# ==============================================================================

message("=== 05b Plot 3: reliability lollipop THAL ===")

cl_by_score <- label_df$cell_type[order(label_df$reliability_score, decreasing = FALSE)]

comp_long <- rbind(
  data.frame(
    cell_type         = label_df$cell_type, n_cells = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    roi_label         = label_df$mmc_roi_label,
    class_label       = label_df$mmc_class_label,
    componente        = "Purity mmc_roi",
    valore            = label_df$mmc_roi_purity,
    stringsAsFactors  = FALSE
  ),
  data.frame(
    cell_type         = label_df$cell_type, n_cells = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    roi_label         = label_df$mmc_roi_label,
    class_label       = label_df$mmc_class_label,
    componente        = "Purity mmc_class",
    valore            = label_df$mmc_class_purity,
    stringsAsFactors  = FALSE
  ),
  data.frame(
    cell_type         = label_df$cell_type, n_cells = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    roi_label         = label_df$mmc_roi_label,
    class_label       = label_df$mmc_class_label,
    componente        = "Concordance (roi+class)",
    valore            = label_df$concordance,
    stringsAsFactors  = FALSE
  )
)

comp_long$cell_type    <- factor(comp_long$cell_type,    levels = cl_by_score)
comp_long$componente   <- factor(comp_long$componente,
                                 levels = c("Purity mmc_roi", "Purity mmc_class",
                                            "Concordance (roi+class)"))
comp_long$reliability_grade <- factor(comp_long$reliability_grade,
                                      levels = c("HIGH","MEDIUM","LOW"))
comp_long$has_warning <- comp_long$warning_flags != ""

ann_base       <- label_df[!duplicated(label_df$cell_type), ]
ann_base$cell_type  <- factor(ann_base$cell_type, levels = cl_by_score)
ann_base$label_text <- paste0(ann_base$mmc_roi_label, "  /  ", ann_base$mmc_class_label)

bg_df   <- comp_long[comp_long$componente == "Purity mmc_roi", ]
warn_df <- comp_long[comp_long$has_warning &
                       comp_long$componente == "Concordance (roi+class)", ]
has_conf <- any(!is.na(label_df$conf_median))
conf_df  <- ann_base

comp_colors <- c(
  "Purity mmc_roi"         = "#2980B9",
  "Purity mmc_class"       = "#8E44AD",
  "Concordance (roi+class)" = "#27AE60"
)

p_rel <- ggplot(comp_long, aes(x = valore, y = cell_type, color = componente)) +
  geom_segment(data = bg_df,
               aes(x = 0, xend = 1, y = cell_type, yend = cell_type),
               color = "grey92", linewidth = 4, inherit.aes = FALSE) +
  geom_segment(data = bg_df,
               aes(x = 0, xend = reliability_score,
                   y = cell_type, yend = cell_type, color = reliability_grade),
               linewidth = 4, alpha = 0.18, inherit.aes = FALSE) +
  geom_segment(aes(x = 0, xend = valore, yend = cell_type),
               linewidth = 1.4, lineend = "round",
               position  = position_dodge(width = 0.7)) +
  geom_point(aes(size = n_cells),
             position = position_dodge(width = 0.7), alpha = 0.95) +
  geom_vline(xintercept = THR_PURITY_LOW,
             linetype = "dashed",   color = "#E74C3C", linewidth = 0.45, alpha = 0.8) +
  geom_vline(xintercept = THR_CONCORD_WARN,
             linetype = "longdash", color = "#F39C12", linewidth = 0.4,  alpha = 0.7) +
  geom_vline(xintercept = THR_PURITY_HIGH,
             linetype = "dotted",   color = "#27AE60", linewidth = 0.45, alpha = 0.8) +
  geom_point(data = warn_df, aes(x = 1.06, y = cell_type),
             shape = 17, color = "#E74C3C", size = 2.8, inherit.aes = FALSE) +
  geom_text(data = warn_df,
            aes(x = 1.08, y = cell_type, label = gsub("\\|", "\n", warning_flags)),
            hjust = 0, vjust = 0.4, size = 1.7, color = "#E74C3C",
            lineheight = 0.85, inherit.aes = FALSE) +
  geom_text(data = ann_base,
            aes(x = -0.01, y = cell_type, label = label_text),
            hjust = 1, size = 2.1, color = "grey30", inherit.aes = FALSE) +
  { if (has_conf)
    geom_point(data = conf_df,
               aes(x = 1.16, y = cell_type, fill = conf_median),
               shape = 22, size = 3.5, color = "white", stroke = 0.3,
               inherit.aes = FALSE)
  } +
  scale_color_manual(
    values = c(comp_colors, HIGH = "#27AE60", MEDIUM = "#F39C12", LOW = "#E74C3C"),
    breaks = names(comp_colors), name = "Componente") +
  { if (has_conf)
    scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
                         midpoint = THR_CONF_LOW, limits = c(0,1), name = "Conf.\nMMC")
  } +
  scale_size_continuous(range = c(1.5, 5.5), name = "N cellule",
                        guide = guide_legend(override.aes = list(color="grey60", shape=16))) +
  scale_x_continuous(
    limits = c(-0.40, if (has_conf) 1.30 else 1.12),
    breaks = c(0, THR_PURITY_LOW, THR_CONCORD_WARN, THR_PURITY_HIGH, 1),
    labels = c("0",
               sprintf("%.0f%%", THR_PURITY_LOW   * 100),
               sprintf("%.0f%%", THR_CONCORD_WARN * 100),
               sprintf("%.0f%%", THR_PURITY_HIGH  * 100),
               "100%"),
    expand = c(0, 0)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_blank(), legend.position = "right",
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7)) +
  labs(title    = "Reliability per cluster - THAL",
       subtitle = sprintf(
         "Score = %.1fx pur_roi + %.1fx pur_class + %.1fx concordance | ordinati per score crescente",
         W_ROI, W_CLASS, W_CONC),
       x = "Proporzione celle")

savepng_thal(p_rel, "05b_3_reliability_THAL.png",
             width  = 5200,
             height = max(2200, length(cluster_order) * 110 + 900))

# ==============================================================================
# PLOT 4 — SCATTER: purity_roi vs purity_class
# ==============================================================================

message("=== 05b Plot 4: scatter purity_roi vs purity_class THAL ===")

scat_df <- label_df
scat_df$reliability_grade <- factor(scat_df$reliability_grade, levels = c("HIGH","MEDIUM","LOW"))
scat_df$has_warning <- scat_df$warning_flags != ""
warn_scat <- scat_df[scat_df$has_warning, ]

p_scat <- ggplot(scat_df,
                 aes(x = mmc_roi_purity, y = mmc_class_purity,
                     color = concordance, size = n_cells,
                     shape = reliability_grade)) +
  annotate("rect", xmin = THR_PURITY_HIGH, xmax = 1,
           ymin = THR_PURITY_HIGH, ymax = 1, fill = "#27AE60", alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = THR_PURITY_LOW,
           ymin = 0, ymax = THR_PURITY_LOW, fill = "#E74C3C", alpha = 0.06) +
  geom_hline(yintercept = c(THR_PURITY_LOW, THR_PURITY_HIGH),
             linetype = c("dashed","dotted"), color = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = c(THR_PURITY_LOW, THR_PURITY_HIGH),
             linetype = c("dashed","dotted"), color = "grey60", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "grey75", linewidth = 0.4) +
  geom_point(alpha = 0.85, stroke = 0.8) +
  geom_point(data = warn_scat,
             aes(x = mmc_roi_purity, y = mmc_class_purity),
             shape = 1, color = "#E74C3C", size = 5, stroke = 1.2,
             inherit.aes = FALSE) +
  ggrepel::geom_text_repel(aes(label = cell_type),
                            size = 2.5, color = "grey20",
                            max.overlaps = 20,
                            segment.size = 0.3, segment.color = "grey70") +
  scale_color_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
                        midpoint = 0.65, limits = c(0,1), name = "Concordance") +
  scale_size_continuous(range = c(2,8), name = "N cellule") +
  scale_shape_manual(values = c(HIGH = 16, MEDIUM = 17, LOW = 15), name = "Grade") +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  annotate("text", x=0.99, y=0.99, label="ZONA HIGH", color="#27AE60",
           size=2.4, fontface="bold", hjust=1, vjust=1) +
  annotate("text", x=0.01, y=0.01, label="ZONA LOW",  color="#E74C3C",
           size=2.4, fontface="bold", hjust=0, vjust=0) +
  theme_bw(base_size = 10) +
  theme(legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 7)) +
  labs(title    = "Purity mmc_roi vs mmc_class per cluster - THAL",
       subtitle = "Colore = concordance | cerchio rosso = cluster con warning",
       x = "Purity mmc_roi", y = "Purity mmc_class")

savepng_thal(p_scat, "05b_4_scatter_THAL.png", width = 3800, height = 3500)

message("  Salvati: 05b_1..4 THAL")
message("=== 05b THAL completato ===")
