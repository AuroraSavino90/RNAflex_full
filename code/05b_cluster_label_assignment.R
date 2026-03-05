################################################################################
##  05b_cluster_label_assignment.R                                            ##
##                                                                            ##
##  Input:  CTX_celltype.RData  (da 05_cell_type_annotation)                 ##
##          MAPPING_CSV         (CSV gerarchico da MapMyCells)                ##
##                                                                            ##
##  Output: CTX_cluster_labels.csv                                            ##
##    05b_1_crosstab_CTX.png   heatmap cluster x broad_anat (proporzioni)    ##
##    05b_2_crosstab_CTX.png   heatmap cluster x mmc_class  (proporzioni)    ##
##    05b_3_reliability_CTX.png  lollipop reliability per cluster             ##
##    05b_4_scatter_CTX.png      scatter purity_anat vs purity_class          ##
##                                                                            ##
##  Reliability score (0-1):                                                  ##
##    purity_anat  = % celle con broad_anat dominante                        ##
##    purity_class = % celle con mmc_class dominante                         ##
##    concordance  = % celle con ENTRAMBI i label dominanti                  ##
##    conf_median  = bootstrapping probability mediana dal CSV MMC           ##
##    score        = W_ANAT*purity_anat + W_CLASS*purity_class +             ##
##                   W_CONC*concordance                                       ##
##    grade        = HIGH (>=80%) / MEDIUM / LOW (<50%) su concordance       ##
##                                                                            ##
##  Warning flags:                                                            ##
##    MIXED_ANAT          purity_anat  < THR_PURITY_LOW                     ##
##    MIXED_CLASS         purity_class < THR_PURITY_LOW                     ##
##    LOW_CONCORDANCE     concordance  < THR_CONCORD_WARN                   ##
##    ANAT_CLASS_MISMATCH anat attesa in CTX ma class inattesa               ##
##    LOW_CONF            bootstrapping mediana < THR_CONF_LOW               ##
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

# ==============================================================================
# PARAMETRI
# ==============================================================================

MAPPING_CSV <- file.path(outdir,
                         "CTX_from_mapmycells/CTX_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1772014289213.csv")

EXPECTED_CLASS_CTX <- c("Glutamatergic", "GABAergic",
                        "Astrocyte", "Oligodendrocyte", "OPC",
                        "Microglia", "Vascular")
EXPECTED_ANAT_CTX  <- c("Isocortex", "Claustrum", "Cortical subplate")

THR_PURITY_HIGH  <- 0.80
THR_PURITY_LOW   <- 0.50
THR_CONCORD_WARN <- 0.60
THR_CONF_LOW     <- 0.70

W_ANAT  <- 0.30
W_CLASS <- 0.30
W_CONC  <- 0.40

grade_pal <- c(HIGH = "#27AE60", MEDIUM = "#F39C12", LOW = "#E74C3C")

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "CTX_celltype.RData"))   # -> CTX

if (file.exists(MAPPING_CSV)) {
  mapping_raw <- read.csv(MAPPING_CSV, comment.char = "#")
  rownames(mapping_raw) <- mapping_raw$cell_id
  has_mapping <- TRUE
  message(sprintf("  CSV MMC caricato: %d cellule | colonne: %s",
                  nrow(mapping_raw),
                  paste(colnames(mapping_raw), collapse = ", ")))
} else {
  has_mapping <- FALSE
  message("  WARN: MAPPING_CSV non trovato - conf_median sara' NA")
}

# ==============================================================================
# METADATI
# ==============================================================================

meta         <- CTX@meta.data
meta$cluster <- as.character(Idents(CTX))

# Confidence score per classe dal CSV MMC

  conf_col <- grep("class.*bootstrapping_probability",
                   colnames(mapping_raw), value = TRUE,
                   ignore.case = TRUE)[1]
  
    shared <- intersect(rownames(meta), rownames(mapping_raw))
    meta$conf_class <- NA_real_
    meta[shared, "conf_class"] <- as.numeric(mapping_raw[shared, conf_col])
    message(sprintf("  Confidence: '%s' | mediana globale = %.2f",
                    conf_col, median(meta$conf_class, na.rm = TRUE)))
  


cluster_order <- as.character(sort((unique(meta$cell_type))))

# ==============================================================================
# FUNZIONI HELPER
# ==============================================================================

shannon_entropy <- function(tab) {
  tab <- tab[tab > 0]
  if (length(tab) <= 1) return(0)
  p <- tab / sum(tab)
  -sum(p * log2(p))
}

top1 <- function(v) {
  tab <- sort(table(v[!is.na(v)]), decreasing = TRUE)
  if (!length(tab)) return(NA_character_)
  names(tab)[1]
}

top1_pct <- function(v) {
  tab <- table(v[!is.na(v)])
  if (!length(tab)) return(NA_real_)
  round(max(tab) / sum(tab), 4)
}

top2 <- function(v) {
  tab <- sort(table(v[!is.na(v)]), decreasing = TRUE)
  if (length(tab) < 2) return(NA_character_)
  names(tab)[2]
}

top2_pct <- function(v) {
  tab <- sort(table(v[!is.na(v)]), decreasing = TRUE)
  if (length(tab) < 2) return(0)
  round(as.numeric(tab[2]) / sum(tab), 4)
}

# ==============================================================================
# CALCOLO METRICHE PER CLUSTER
# ==============================================================================

message("=== 05b. Label assignment - CTX ===")

rows <- lapply(cluster_order, function(cl) {
  
  m  <- meta[meta$cell_type == cl, , drop = FALSE]
  nc <- nrow(m)
  
  dom_anat  <- top1(m$broad_anat)
  dom_class <- top1(m$mmc_class)
  pur_anat  <- top1_pct(m$broad_anat)
  pur_class <- top1_pct(m$mmc_class)
  
  ent_anat  <- round(shannon_entropy(table(na.omit(m$broad_anat))), 3)
  ent_class <- round(shannon_entropy(table(na.omit(m$mmc_class))),  3)
  
  rup_anat      <- top2(m$broad_anat)
  rup_anat_pct  <- top2_pct(m$broad_anat)
  rup_class     <- top2(m$mmc_class)
  rup_class_pct <- top2_pct(m$mmc_class)
  
  n_anat  <- length(unique(na.omit(m$broad_anat)))
  n_class <- length(unique(na.omit(m$mmc_class)))
  
  # Concordance: % celle con ENTRAMBI broad_anat=dom E mmc_class=dom
  valid   <- !is.na(m$broad_anat) & !is.na(m$mmc_class)
  n_valid <- sum(valid)
  concordance <- if (n_valid > 0 && !is.na(dom_anat) && !is.na(dom_class))
    round(sum(m$broad_anat[valid] == dom_anat &
                m$mmc_class[valid]  == dom_class) / n_valid, 4)
  else NA_real_
  
  # Mismatch biologico
  anat_ok  <- !is.na(dom_anat)  && dom_anat  %in% EXPECTED_ANAT_CTX
  class_ok <- !is.na(dom_class) && dom_class %in% EXPECTED_CLASS_CTX
  mismatch <- anat_ok && !class_ok
  
  # Confidence bootstrapping
  cv           <- m$conf_class[!is.na(m$conf_class)]
  conf_median  <- if (length(cv)) round(median(cv), 3)               else NA_real_
  conf_pct_low <- if (length(cv)) round(mean(cv < THR_CONF_LOW), 3)  else NA_real_
  
  # Score composito
  score <- round(
    W_ANAT  * ifelse(is.na(pur_anat),    0, pur_anat)   +
      W_CLASS * ifelse(is.na(pur_class),   0, pur_class)  +
      W_CONC  * ifelse(is.na(concordance), 0, concordance),
    4)
  
  # Grade su concordance
  base_g <- if (!is.na(concordance)) concordance else
    mean(c(pur_anat, pur_class), na.rm = TRUE)
  grade <- if (is.na(base_g))          NA_character_ else
    if (base_g >= THR_PURITY_HIGH) "HIGH"    else
      if (base_g <  THR_PURITY_LOW)  "LOW"     else
        "MEDIUM"
  
  # Warning flags
  flags <- character(0)
  if (!is.na(pur_anat)    && pur_anat    < THR_PURITY_LOW)   flags <- c(flags, "MIXED_ANAT")
  if (!is.na(pur_class)   && pur_class   < THR_PURITY_LOW)   flags <- c(flags, "MIXED_CLASS")
  if (!is.na(concordance) && concordance < THR_CONCORD_WARN) flags <- c(flags, "LOW_CONCORDANCE")
  if (mismatch)                                               flags <- c(flags, "ANAT_CLASS_MISMATCH")
  if (!is.na(conf_median) && conf_median < THR_CONF_LOW)     flags <- c(flags, "LOW_CONF")
  
  list(
    cell_type                  = cl,
    n_cells                  = nc,
    broad_anat_label         = dom_anat,
    broad_anat_purity        = pur_anat,
    broad_anat_entropy       = ent_anat,
    broad_anat_n             = n_anat,
    broad_anat_runner_up     = rup_anat,
    broad_anat_runner_up_pct = rup_anat_pct,
    mmc_class_label          = dom_class,
    mmc_class_purity         = pur_class,
    mmc_class_entropy        = ent_class,
    mmc_class_n              = n_class,
    mmc_class_runner_up      = rup_class,
    mmc_class_runner_up_pct  = rup_class_pct,
    concordance              = concordance,
    conf_median              = conf_median,
    conf_pct_low             = conf_pct_low,
    reliability_score        = score,
    reliability_grade        = grade,
    warning_flags            = paste(flags, collapse = "|")
  )
})

# Costruzione label_df: ogni colonna esplicitamente estratta dalla lista
# (evita problemi di tipo e nomi che si hanno con bind_rows/as.data.frame)
label_df <- data.frame(
  cell_type                  = sapply(rows, `[[`, "cell_type"),
  n_cells                  = as.integer(sapply(rows, `[[`, "n_cells")),
  broad_anat_label         = sapply(rows, `[[`, "broad_anat_label"),
  broad_anat_purity        = as.numeric(sapply(rows, `[[`, "broad_anat_purity")),
  broad_anat_entropy       = as.numeric(sapply(rows, `[[`, "broad_anat_entropy")),
  broad_anat_n             = as.integer(sapply(rows, `[[`, "broad_anat_n")),
  broad_anat_runner_up     = sapply(rows, `[[`, "broad_anat_runner_up"),
  broad_anat_runner_up_pct = as.numeric(sapply(rows, `[[`, "broad_anat_runner_up_pct")),
  mmc_class_label          = sapply(rows, `[[`, "mmc_class_label"),
  mmc_class_purity         = as.numeric(sapply(rows, `[[`, "mmc_class_purity")),
  mmc_class_entropy        = as.numeric(sapply(rows, `[[`, "mmc_class_entropy")),
  mmc_class_n              = as.integer(sapply(rows, `[[`, "mmc_class_n")),
  mmc_class_runner_up      = sapply(rows, `[[`, "mmc_class_runner_up"),
  mmc_class_runner_up_pct  = as.numeric(sapply(rows, `[[`, "mmc_class_runner_up_pct")),
  concordance              = as.numeric(sapply(rows, `[[`, "concordance")),
  conf_median              = as.numeric(sapply(rows, `[[`, "conf_median")),
  conf_pct_low             = as.numeric(sapply(rows, `[[`, "conf_pct_low")),
  reliability_score        = as.numeric(sapply(rows, `[[`, "reliability_score")),
  reliability_grade        = sapply(rows, `[[`, "reliability_grade"),
  warning_flags            = sapply(rows, `[[`, "warning_flags"),
  stringsAsFactors         = FALSE
)

# Sommario console
message(sprintf("  Cluster: %d | concordance mediana: %.2f",
                nrow(label_df), median(label_df$concordance, na.rm = TRUE)))
for (g in c("HIGH", "MEDIUM", "LOW")) {
  n <- sum(label_df$reliability_grade == g, na.rm = TRUE)
  message(sprintf("    %6s: %d (%.0f%%)", g, n, 100 * n / nrow(label_df)))
}
message(sprintf("  Warning: %d cluster",
                sum(label_df$warning_flags != "", na.rm = TRUE)))

write.csv(label_df,
          file.path(outdir, "CTX_cluster_labels.csv"),
          row.names = FALSE)
message("  Salvato: CTX_cluster_labels.csv")

# ==============================================================================
# PLOT 1 & 2 — HEATMAP: cluster x broad_anat / cluster x mmc_class
# ==============================================================================
# Asse Y = cluster (ordine numerico)
# Asse X = categoria, ordinata per frequenza globale decrescente
# Colore = proporzione celle in quel cluster con quella categoria
# Annotazione riga = reliability grade + warning
# Celle con >= 5% mostrano il valore numerico
# ==============================================================================

message("=== 05b Plot 1-2: heatmap proporzionali ===")

make_crosstab_heatmap <- function(col, filename, plot_title) {
  
  # Tabella assoluta -> matrice numerica, normalizzata per riga
  raw_tab  <- table(meta$cell_type, meta[[col]])
  ok_cols  <- colnames(raw_tab)[colnames(raw_tab) != "NA" & !is.na(colnames(raw_tab))]
  raw_tab  <- raw_tab[, ok_cols, drop = FALSE]
  prop_mat <- sweep(as.matrix(raw_tab), 1, rowSums(as.matrix(raw_tab)), "/")
  
  # Ordina: righe per cluster numerico, colonne per frequenza globale
  row_ord  <- cluster_order[cluster_order %in% rownames(prop_mat)]
  col_ord  <- names(sort(colSums(prop_mat), decreasing = TRUE))
  prop_mat <- prop_mat[row_ord, col_ord, drop = FALSE]
  
  # Testo celle: % intera solo se >= 5%
  display_mat <- matrix(
    ifelse(prop_mat >= 0.05, as.character(round(prop_mat * 100)), ""),
    nrow = nrow(prop_mat), ncol = ncol(prop_mat),
    dimnames = dimnames(prop_mat)
  )
  
  # Annotazione righe allineata a row_ord
  idx     <- match(row_ord, label_df$cell_type)
  row_ann <- data.frame(
    Grade   = factor(label_df$reliability_grade[idx],
                     levels = c("HIGH","MEDIUM","LOW")),
    Warning = ifelse(is.na(label_df$warning_flags[idx]) |
                       label_df$warning_flags[idx] == "", "no", "si"),
    row.names        = row_ord,
    stringsAsFactors = FALSE
  )
  
  ann_colors <- list(
    Grade   = c(HIGH = "#27AE60", MEDIUM = "#F39C12", LOW = "#E74C3C"),
    Warning = c(no = "grey90", si = "#E74C3C")
  )
  
  graphics.off()
  png(file.path(outdir, filename),
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
    fontsize          = 8,
    fontsize_row      = 7,
    fontsize_col      = 8,
    fontsize_number   = 7,
    angle_col         = "45",
    border_color      = "white",
    main              = plot_title
  ))
  dev.off()
  message(sprintf("  Salvato: %s", filename))
}

make_crosstab_heatmap(
  col        = "broad_anat",
  filename   = "05b_1_crosstab_CTX.png",
  plot_title = "Proporzione broad_anat per cluster - CTX\n(colore = % celle nel cluster | annotaz. destra = grade + warning)"
)

make_crosstab_heatmap(
  col        = "mmc_class",
  filename   = "05b_2_crosstab_CTX.png",
  plot_title = "Proporzione mmc_class per cluster - CTX\n(colore = % celle nel cluster | annotaz. destra = grade + warning)"
)

# ==============================================================================
# PLOT 3 — RELIABILITY LOLLIPOP
# ==============================================================================
# Ogni riga = un cluster, ordinato per reliability_score crescente
# (peggiori in basso).
# Tre segmenti per riga in dodge verticale:
#   blu   = purity broad_anat
#   viola = purity mmc_class
#   verde = concordance (metrica chiave)
# Sfondo tenue = reliability_score composito
# Triangolo rosso + testo = warning flags
# Quadratino a destra = conf bootstrapping MMC (se disponibile)
# ==============================================================================

message("=== 05b Plot 3: reliability lollipop ===")

# Ordine cluster: dal peggiore (score piu' basso) al migliore (in alto)
cl_by_score <- label_df$cell_type[order(label_df$reliability_score,
                                      decreasing = FALSE)]

# Data frame long per i tre componenti (senza pivot_longer)
comp_long <- rbind(
  data.frame(
    cell_type           = label_df$cell_type,
    n_cells           = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    broad_anat_label  = label_df$broad_anat_label,
    mmc_class_label   = label_df$mmc_class_label,
    componente        = "Purity broad_anat",
    valore            = label_df$broad_anat_purity,
    stringsAsFactors  = FALSE
  ),
  data.frame(
    cell_type           = label_df$cell_type,
    n_cells           = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    broad_anat_label  = label_df$broad_anat_label,
    mmc_class_label   = label_df$mmc_class_label,
    componente        = "Purity mmc_class",
    valore            = label_df$mmc_class_purity,
    stringsAsFactors  = FALSE
  ),
  data.frame(
    cell_type           = label_df$cell_type,
    n_cells           = label_df$n_cells,
    reliability_score = label_df$reliability_score,
    reliability_grade = label_df$reliability_grade,
    warning_flags     = label_df$warning_flags,
    conf_median       = label_df$conf_median,
    broad_anat_label  = label_df$broad_anat_label,
    mmc_class_label   = label_df$mmc_class_label,
    componente        = "Concordance (anat+class)",
    valore            = label_df$concordance,
    stringsAsFactors  = FALSE
  )
)

# Fattori con ordinamento esplicito
comp_long$cell_type    <- factor(comp_long$cell_type,    levels = cl_by_score)
comp_long$componente <- factor(comp_long$componente,
                               levels = c("Purity broad_anat",
                                          "Purity mmc_class",
                                          "Concordance (anat+class)"))
comp_long$reliability_grade <- factor(comp_long$reliability_grade,
                                      levels = c("HIGH","MEDIUM","LOW"))
comp_long$has_warning <- comp_long$warning_flags != ""

# df una-riga-per-cluster per le annotazioni
ann_base <- label_df[!duplicated(label_df$cell_type), ]
ann_base$cell_type    <- factor(ann_base$cell_type, levels = cl_by_score)
ann_base$label_text <- paste0(ann_base$broad_anat_label,
                              "  /  ",
                              ann_base$mmc_class_label)

# df per sfondo reliability_score (una riga per cluster, non triplicata)
bg_df <- comp_long[comp_long$componente == "Purity broad_anat", ]

# df per warning (una riga per cluster con warning, sulla riga concordance)
warn_df <- comp_long[comp_long$has_warning &
                       comp_long$componente == "Concordance (anat+class)", ]

# df per confidence (se disponibile)
has_conf <- any(!is.na(label_df$conf_median))
conf_df  <- ann_base

comp_colors <- c(
  "Purity broad_anat"      = "#2980B9",
  "Purity mmc_class"       = "#8E44AD",
  "Concordance (anat+class)" = "#27AE60"
)

p_rel <- ggplot(comp_long,
                aes(x = valore, y = cell_type, color = componente)) +
  
  # Sfondo grigio 0-1 per ogni cluster
  geom_segment(
    data        = bg_df,
    aes(x = 0, xend = 1, y = cell_type, yend = cell_type),
    color       = "grey92",
    linewidth   = 4,
    inherit.aes = FALSE) +
  
  # Sfondo tenue = reliability score composito
  geom_segment(
    data        = bg_df,
    aes(x = 0, xend = reliability_score,
        y = cell_type, yend = cell_type,
        color = reliability_grade),
    linewidth   = 4,
    alpha       = 0.18,
    inherit.aes = FALSE) +
  
  # Segmenti componenti (dodge verticale)
  geom_segment(
    aes(x = 0, xend = valore, yend = cell_type),
    linewidth = 1.4,
    lineend   = "round",
    position  = position_dodge(width = 0.7)) +
  
  # Punti terminali (dimensione = n cellule)
  geom_point(
    aes(size = n_cells),
    position = position_dodge(width = 0.7),
    alpha    = 0.95) +
  
  # Soglie verticali
  geom_vline(xintercept = THR_PURITY_LOW,
             linetype = "dashed",  color = "#E74C3C",
             linewidth = 0.45, alpha = 0.8) +
  geom_vline(xintercept = THR_CONCORD_WARN,
             linetype = "longdash", color = "#F39C12",
             linewidth = 0.4,  alpha = 0.7) +
  geom_vline(xintercept = THR_PURITY_HIGH,
             linetype = "dotted", color = "#27AE60",
             linewidth = 0.45, alpha = 0.8) +
  
  # Triangolo warning a destra
  geom_point(
    data        = warn_df,
    aes(x = 1.06, y = cell_type),
    shape       = 17,
    color       = "#E74C3C",
    size        = 2.8,
    inherit.aes = FALSE) +
  
  # Testo warning flags
  geom_text(
    data        = warn_df,
    aes(x = 1.08, y = cell_type,
        label = gsub("\\|", "\n", warning_flags)),
    hjust       = 0,
    vjust       = 0.4,
    size        = 1.7,
    color       = "#E74C3C",
    lineheight  = 0.85,
    inherit.aes = FALSE) +
  
  # Etichetta label dominanti a sinistra
  geom_text(
    data        = ann_base,
    aes(x = -0.01, y = cell_type, label = label_text),
    hjust       = 1,
    size        = 2.1,
    color       = "grey30",
    inherit.aes = FALSE) +
  
  # Quadratino confidence MMC (se disponibile)
  { if (has_conf)
    geom_point(
      data        = conf_df,
      aes(x = 1.16, y = cell_type, fill = conf_median),
      shape       = 22,
      size        = 3.5,
      color       = "white",
      stroke      = 0.3,
      inherit.aes = FALSE)
  } +
  
  scale_color_manual(
    values = c(
      comp_colors,
      HIGH   = "#27AE60",
      MEDIUM = "#F39C12",
      LOW    = "#E74C3C"),
    breaks = names(comp_colors),
    name   = "Componente") +
  
  { if (has_conf)
    scale_fill_gradient2(
      low      = "#E74C3C",
      mid      = "#F39C12",
      high     = "#27AE60",
      midpoint = THR_CONF_LOW,
      limits   = c(0, 1),
      name     = "Conf.\nMMC")
  } +
  
  scale_size_continuous(
    range = c(1.5, 5.5),
    name  = "N cellule",
    guide = guide_legend(override.aes = list(color = "grey60", shape = 16))) +
  
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
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_blank(),
    legend.position    = "right",
    legend.key.size    = unit(0.4, "cm"),
    legend.text        = element_text(size = 7)
  ) +
  labs(
    title    = "Reliability per cluster - CTX",
    subtitle = sprintf(
      paste0("Score = %.1fx pur_anat + %.1fx pur_class + %.1fx concordance",
             " | ordinati per score crescente (peggiori in basso)\n",
             "Grade su concordance: HIGH>=%.0f%% | MEDIUM %.0f-%.0f%%",
             " | LOW<%.0f%% | triangolo = warning"),
      W_ANAT, W_CLASS, W_CONC,
      THR_PURITY_HIGH * 100,
      THR_PURITY_LOW  * 100,
      THR_PURITY_HIGH * 100,
      THR_PURITY_LOW  * 100),
    x = "Proporzione celle"
  )

savepng(p_rel, "05b_3_reliability_CTX.png",
        width  = 5200,
        height = max(2200, length(cluster_order) * 110 + 900))

# ==============================================================================
# PLOT 4 — SCATTER: purity_anat vs purity_class
# ==============================================================================
# Ogni punto = un cluster
# Colore = concordance | Dimensione = n cellule | Forma = grade
# Cerchio rosso aggiuntivo = cluster con warning
# Quadranti delimitati da soglie HIGH/LOW
# ==============================================================================

message("=== 05b Plot 4: scatter purity_anat vs purity_class ===")

scat_df <- label_df
scat_df$reliability_grade <- factor(scat_df$reliability_grade,
                                    levels = c("HIGH","MEDIUM","LOW"))
scat_df$has_warning <- scat_df$warning_flags != ""

warn_scat <- scat_df[scat_df$has_warning, ]

p_scat <- ggplot(scat_df,
                 aes(x     = broad_anat_purity,
                     y     = mmc_class_purity,
                     color = concordance,
                     size  = n_cells,
                     shape = reliability_grade)) +
  
  # Sfondo zone HIGH / LOW
  annotate("rect",
           xmin = THR_PURITY_HIGH, xmax = 1,
           ymin = THR_PURITY_HIGH, ymax = 1,
           fill = "#27AE60", alpha = 0.06) +
  annotate("rect",
           xmin = 0, xmax = THR_PURITY_LOW,
           ymin = 0, ymax = THR_PURITY_LOW,
           fill = "#E74C3C", alpha = 0.06) +
  
  # Soglie
  geom_hline(yintercept = THR_PURITY_LOW,
             linetype = "dashed",  color = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = THR_PURITY_HIGH,
             linetype = "dotted",  color = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = THR_PURITY_LOW,
             linetype = "dashed",  color = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = THR_PURITY_HIGH,
             linetype = "dotted",  color = "grey60", linewidth = 0.5) +
  
  # Diagonale parita'
  geom_abline(slope = 1, intercept = 0,
              color = "grey75", linewidth = 0.4) +
  
  # Punti
  geom_point(alpha = 0.85, stroke = 0.8) +
  
  # Cerchio warning
  geom_point(data        = warn_scat,
             aes(x = broad_anat_purity, y = mmc_class_purity),
             shape       = 1,
             color       = "#E74C3C",
             size        = 5,
             stroke      = 1.2,
             inherit.aes = FALSE) +
  
  # Etichette cluster
  ggrepel::geom_text_repel(
    aes(label = cell_type),
    size          = 2.5,
    color         = "grey20",
    max.overlaps  = 20,
    segment.size  = 0.3,
    segment.color = "grey70") +
  
  scale_color_gradient2(
    low      = "#E74C3C",
    mid      = "#F39C12",
    high     = "#27AE60",
    midpoint = 0.65,
    limits   = c(0, 1),
    name     = "Concordance") +
  scale_size_continuous(range = c(2, 8), name = "N cellule") +
  scale_shape_manual(
    values = c(HIGH = 16, MEDIUM = 17, LOW = 15),
    name   = "Grade") +
  
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  
  annotate("text", x = 0.99, y = 0.99,
           label = "ZONA HIGH", color = "#27AE60",
           size = 2.4, fontface = "bold", hjust = 1, vjust = 1) +
  annotate("text", x = 0.01, y = 0.01,
           label = "ZONA LOW",  color = "#E74C3C",
           size = 2.4, fontface = "bold", hjust = 0, vjust = 0) +
  
  theme_bw(base_size = 10) +
  theme(
    legend.key.size = unit(0.45, "cm"),
    legend.text     = element_text(size = 7)
  ) +
  labs(
    title    = "Purity broad_anat vs mmc_class per cluster - CTX",
    subtitle = paste0(
      "Colore = concordance (% celle con entrambi i label dominanti)",
      " | cerchio rosso = cluster con warning\n",
      "Diagonale = parity line | soglie: tratteggio=50%, punteggiato=80%"),
    x = "Purity broad_anat",
    y = "Purity mmc_class"
  )

savepng(p_scat, "05b_4_scatter_CTX.png",
        width  = 3800,
        height = 3500)

message("  Salvati: 05b_1_crosstab_CTX.png | 05b_2_crosstab_CTX.png | 05b_3_reliability_CTX.png | 05b_4_scatter_CTX.png")
message("=== 05b completato ===")