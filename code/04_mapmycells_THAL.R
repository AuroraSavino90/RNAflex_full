################################################################################
##  04_mapmycells_THAL.R                                                       ##
##                                                                             ##
##  Input:  THAL_initial.RData  (da sezione 03)                                ##
##  Output: THAL_for_mapmycells.h5ad  (da caricare su MapMyCells)             ##
##          THAL_annotated.RData  (post-MapMyCells)                            ##
##          THAL_annotated_MMC.RData  (con ROI talamiche)                      ##
##                                                                             ##
##  Flusso:                                                                    ##
##    PARTE A — esegui fino all'export h5ad, poi vai su MapMyCells online      ##
##    PARTE B — decommentare dopo aver scaricato il CSV risultati              ##
##    PARTE C — assegnazione ROI talamica da cluster MapMyCells               ##
################################################################################

source("code/00_config.R")
outdir <- "results/20260224"

# ==============================================================================
# CARICAMENTO
# ==============================================================================

load(file.path(outdir, "THAL_initial.RData"))   # THAL

# ==============================================================================
# PARTE A — EXPORT H5AD PER MAPMYCELLS
# ==============================================================================

message("=== 04A. Export h5ad per MapMyCells (THAL) ===")

# Seurat v5: LayerData() al posto di @assays$RNA@counts
# I conteggi raw sono nel layer "counts" (dopo JoinLayers in sezione 03)

export_h5ad <- function(obj, path, label) {
  counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")
  counts_t   <- Matrix::t(counts_mat)

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

export_h5ad(THAL, file.path(outdir, "THAL_for_mapmycells.h5ad"), "THAL")

message("
>>> PROSSIMI PASSI (PARTE A completata):
  1. Vai su https://portal.brain-map.org/atlases-and-data/bkp/mapmycells
  2. Carica THAL_for_mapmycells.h5ad
  3. Seleziona: Whole Mouse Brain (CCN20230722), Hierarchical mapping
  4. Scarica il CSV dei risultati
  5. Scarica il metadata Excel cl.df_CCN20230722.xlsx dal sito AllenBrain
  6. Copia i file in: ", outdir, "
  7. Passa alla parte B
")

# ==============================================================================
# PARTE B — PROCESSAMENTO RISULTATI MAPMYCELLS
# (decommentare dopo aver scaricato il CSV)
# ==============================================================================

# --- Parametri ---------------------------------------------------------------
MAPPING_CSV <- file.path(outdir, "THAL_from_mapmycells/THAL_for_mapmycells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1772014736467.csv")
META_XLSX   <- file.path("data/cl.df_CCN202307220.xlsx")

MIN_CELLS_CLASS    <- 100
MIN_CELLS_SUBCLASS <- 20

obj    <- THAL
obj_nm <- "THAL"

message("=== 04B. Processamento MapMyCells (THAL) ===")

mapping <- read.csv(MAPPING_CSV, comment.char = "#")
meta_df <- read.xlsx(META_XLSX)

# Verifica colonne attese
expected_cols <- c("cell_id", "class_name", "subclass_name", "cluster_name")
if (!all(expected_cols %in% colnames(mapping))) {
  stop(sprintf(
    "Colonne mancanti nel CSV MapMyCells. Attese: %s. Trovate: %s",
    paste(expected_cols, collapse = ", "),
    paste(colnames(mapping), collapse = ", ")
  ))
}
rownames(mapping) <- mapping$cell_id

clust_meta <- read.csv("data/cluster_to_cluster_annotation_membership.csv")
meta <- read.csv("data/cell_metadata.csv", row.names = 1)



# ==============================================================================
# Arricchimento TH per cluster MapMyCells
# ==============================================================================
# Per ogni cluster (cluster_alias) presente nel tuo THAL:
#   - conta quante cellule dell'atlas provengono da TH vs altre ROI
#   - normalizza per il numero totale di cellule campionate per ROI nell'atlas
#   - calcola % TH grezza e % TH normalizzata (arricchimento relativo)
# ==============================================================================

library(dplyr)

# --- 1. Cellule per ROI nell'atlas (denominatore normalizzazione) ------------
# Quante cellule totali sono state campionate per ogni ROI nel WMB atlas?
# Questo è il fattore di correzione: TH ha meno cellule di Isocortex nell'atlas.

roi_totals <- meta %>%
  filter(!is.na(region_of_interest_acronym)) %>%
  count(region_of_interest_acronym, name = "n_roi_total")

cat("Cellule totali per ROI nell'atlas (prime 15):\n")
print(head(arrange(roi_totals, desc(n_roi_total)), 15))


# --- 2. Composizione per cluster: quante cellule per ROI -------------------
cluster_roi_counts <- meta %>%
  filter(!is.na(region_of_interest_acronym), !is.na(cluster_alias)) %>%
  count(cluster_alias, region_of_interest_acronym, name = "n_cells")


# --- 3. Aggiungi denominatore e calcola proporzione normalizzata ------------
cluster_roi_norm <- cluster_roi_counts %>%
  left_join(roi_totals, by = "region_of_interest_acronym") %>%
  mutate(
    # proporzione grezza: % cellule di quel cluster che vengono da questa ROI
    n_cluster_total = sum(n_cells),         # tot cellule nel cluster (per ROI)
    .by = cluster_alias
  ) %>%
  mutate(
    pct_raw = 100 * n_cells / n_cluster_total,
    # score normalizzato: (n_cells / n_roi_total) pesato sul cluster
    # equivale a: quanto è arricchita questa ROI rispetto alla sua dimensione nell'atlas?
    norm_score = (n_cells / n_roi_total)
  ) %>%
  # normalizza norm_score in modo che la somma per cluster = 1
  mutate(
    norm_score_sum = sum(norm_score),
    pct_norm = 100 * norm_score / norm_score_sum,
    .by = cluster_alias
  ) %>%
  select(cluster_alias, region_of_interest_acronym,
         n_cells, n_cluster_total, n_roi_total,
         pct_raw, pct_norm)


# --- 4. Estrai statistiche TH per cluster -----------------------------------
cluster_TH_stats <- cluster_roi_norm %>%
  group_by(cluster_alias) %>%
  summarise(
    n_cells_total    = first(n_cluster_total),
    n_cells_TH       = sum(n_cells[region_of_interest_acronym == "TH"]),
    pct_TH_raw       = sum(pct_raw[region_of_interest_acronym == "TH"]),
    pct_TH_norm      = sum(pct_norm[region_of_interest_acronym == "TH"]),
    top_roi          = region_of_interest_acronym[which.max(n_cells)],
    top_roi_pct_raw  = max(pct_raw),
    top_roi_pct_norm = pct_norm[which.max(pct_norm)],  # usa pct_norm invece di norm_score
    .groups = "drop"
  ) %>%
  arrange(desc(pct_TH_norm))

cat("\nCluster con maggior arricchimento TH (normalizzato) — prime 20 righe:\n")
print(head(cluster_TH_stats, 20))


# --- 5. Associa a ogni cellula di THAL le statistiche del suo cluster -------
# Estrai cluster_alias dal CSV MapMyCells (stesso approccio di 04C)

if ("cluster_alias" %in% colnames(mapping)) {
  THAL$mmc_cluster_alias <- mapping$cluster_alias[
    match(colnames(THAL), mapping$cell_id)
  ]
} else {
  THAL$mmc_cluster_alias <- as.character(
    as.integer(sub("^(\\d+).*", "\\1", THAL$mmc_cluster))
  )
}

# Join: cellula -> cluster_alias -> statistiche TH
THAL$cluster_pct_TH_raw  <- cluster_TH_stats$pct_TH_raw[
  match(THAL$mmc_cluster_alias, cluster_TH_stats$cluster_alias)
]
THAL$cluster_pct_TH_norm <- cluster_TH_stats$pct_TH_norm[
  match(THAL$mmc_cluster_alias, cluster_TH_stats$cluster_alias)
]
THAL$cluster_top_roi <- cluster_TH_stats$top_roi[
  match(THAL$mmc_cluster_alias, cluster_TH_stats$cluster_alias)
]

cat("\nDistribuzione cluster_pct_TH_norm nelle cellule di THAL:\n")
print(summary(THAL$cluster_pct_TH_norm))

cat("\nCellule per bin di arricchimento TH (normalizzato):\n")
print(table(cut(THAL$cluster_pct_TH_norm,
                breaks = c(0, 10, 25, 50, 75, 90, 100),
                include.lowest = TRUE)))


# --- 6. Plot diagnostici ---------------------------------------------------

# Istogramma: distribuzione pct_TH_norm per cellula
p_hist <- ggplot(THAL@meta.data, aes(x = cluster_pct_TH_norm)) +
  geom_histogram(binwidth = 5, fill = "#2980B9", color = "white") +
  geom_vline(xintercept = c(25, 50), linetype = "dashed", color = c("orange", "red")) +
  labs(title = "Arricchimento TH del cluster MapMyCells (normalizzato)",
       subtitle = "Ogni barra = cellule di THAL; linee = soglie 25% e 50%",
       x = "% TH normalizzata del cluster", y = "n cellule") +
  theme_classic(base_size = 12)
savepng(p_hist, "04D_hist_TH_enrichment.png", width = 3000, height = 2000)

# UMAP colorato per pct_TH_norm
p_umap <- FeaturePlot(THAL, features = "cluster_pct_TH_norm",
                      pt.size = 0.3, order = TRUE) +
  scale_color_gradientn(
    colors = c("grey90", "#F4A261", "#E63946", "#9B0000"),
    name   = "% TH norm"
  ) +
  ggtitle("Arricchimento TH del cluster (normalizzato)") +
  theme_classic(base_size = 12)
savepng(p_umap, "04D_UMAP_TH_enrichment.png", width = 3500, height = 3000)

# UMAP con soglia binaria (es. >50% TH normalizzato = "talamica")
THAL$is_TH_cluster <- ifelse(
  !is.na(THAL$cluster_pct_TH_norm) & THAL$cluster_pct_TH_norm >= 50,
  "TH-enriched", "other"
)

p_binary <- DimPlot(THAL, group.by = "is_TH_cluster",
                    cols = c("TH-enriched" = "#2980B9", "other" = "grey80"),
                    pt.size = 0.3, order = "TH-enriched") +
  ggtitle("Cluster talamici (pct_TH_norm >= 50%)") +
  theme_classic(base_size = 12)
savepng(p_binary, "04D_UMAP_TH_binary.png", width = 3500, height = 3000)

# Palette: colori distinti per ROI, grigio per TH (le cellule "corrette")
roi_levels <- sort(unique(na.omit(THAL$cluster_top_roi)))

# Costruisci palette: TH in blu THAL, tutto il resto con colori distinti
roi_pal <- setNames(
  scales::hue_pal()(length(roi_levels)),
  roi_levels
)
roi_pal["TH"] <- "#2980B9"   # sovrascrivi TH con il blu della pipeline

p_roi <- DimPlot(
  THAL,
  group.by = "cluster_top_roi",
  cols     = roi_pal,
  pt.size  = 0.3,
  order    = roi_levels[roi_levels != "TH"]   # porta ROI non-TH in primo piano
) +
  ggtitle("ROI predominante nel cluster MapMyCells") +
  labs(colour = "Top ROI (pct_norm)") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               ncol = 2))

savepng(p_roi, "04D_UMAP_top_roi.png", width = 4000, height = 3000)


cat(sprintf(
  "\nCellule TH-enriched (>=50%% norm): %d / %d (%.1f%%)\n",
  sum(THAL$is_TH_cluster == "TH-enriched", na.rm = TRUE),
  ncol(THAL),
  100 * mean(THAL$is_TH_cluster == "TH-enriched", na.rm = TRUE)
))


# --- 7. Tabella riassuntiva per cluster presenti in THAL -------------------
clusters_in_THAL <- unique(na.omit(THAL$mmc_cluster_alias))

summary_table <- cluster_TH_stats %>%
  filter(cluster_alias %in% clusters_in_THAL) %>%
  left_join(
    clust_meta %>% select(cluster_alias, cluster_annotation_term_name),
    by = "cluster_alias"
  ) %>%
  arrange(desc(pct_TH_norm)) %>%
  select(cluster_alias, cluster_annotation_term_name,
         n_cells_total, n_cells_TH,
         pct_TH_raw, pct_TH_norm, top_roi)

cat("\nTabella cluster presenti in THAL con arricchimento TH:\n")
print(summary_table, n = 40)

write.csv(summary_table,
          file.path(outdir, "04D_cluster_TH_enrichment.csv"),
          row.names = FALSE)

message("=== 04D completato ===")


# ==============================================================================
# Analisi cellule con basso arricchimento TH
# ==============================================================================

# --- 1. Definisci soglia e crea gruppi ---------------------------------------
# Guarda prima la distribuzione per scegliere una soglia informata
quantile(THAL$cluster_pct_TH_norm, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

TH_SOGLIA <- 0   # cellule TH_low = cluster con pct_TH_norm == 0 (nessun origine talamica nell'atlas)

THAL$TH_group <- case_when(
  is.na(THAL$cluster_pct_TH_norm)          ~ "no_mapping",
  THAL$cluster_pct_TH_norm > TH_SOGLIA     ~ "TH_high",
  TRUE                                       ~ "TH_low"
)

cat("Cellule per gruppo:\n")
print(table(THAL$TH_group))


# --- 2. Per le cellule TH_low: quali cluster sono? --------------------------
low_cells <- colnames(THAL)[THAL$TH_group == "TH_low"]

low_cluster_summary <- cluster_TH_stats %>%
  filter(cluster_alias %in% THAL$mmc_cluster_alias[THAL$TH_group == "TH_low"]) %>%
  left_join(clust_meta %>% select(cluster_alias, cluster_annotation_term_name),
            by = "cluster_alias") %>%
  mutate(n_cells_in_THAL = table(THAL$mmc_cluster_alias[THAL$TH_group == "TH_low"])[
    as.character(cluster_alias)
  ]) %>%
  arrange(pct_TH_norm) %>%
  select(cluster_alias, cluster_annotation_term_name,
         n_cells_in_THAL, pct_TH_norm, pct_TH_raw, top_roi)

cat("\nCluster a basso arricchimento TH presenti in THAL:\n")
print(low_cluster_summary, n = 40)


# --- 3. Composizione ROI completa per i cluster TH_low ----------------------
# Per ogni cluster TH_low: distribuzione completa delle ROI nell'atlas
low_aliases <- unique(THAL$mmc_cluster_alias[THAL$TH_group == "TH_low"])

low_roi_detail <- cluster_roi_norm %>%
  filter(cluster_alias %in% low_aliases) %>%
  left_join(clust_meta %>% select(cluster_alias, cluster_annotation_term_name),
            by = "cluster_alias") %>%
  arrange(cluster_alias, desc(pct_norm)) %>%
  select(cluster_alias, cluster_annotation_term_name,
         region_of_interest_acronym, n_cells, pct_raw, pct_norm)

cat("\nComposizione ROI per cluster TH_low (top ROI per cluster):\n")
low_roi_detail %>%
  group_by(cluster_alias, cluster_annotation_term_name) %>%
  slice_max(pct_norm, n = 3) %>%
  print(n = 60)


# --- 4. Aggregato: quali ROI dominano le cellule TH_low? --------------------
# Pesa ogni ROI per il numero di cellule TH_low che cadono in quel cluster
low_roi_agg <- cluster_roi_norm %>%
  filter(cluster_alias %in% low_aliases) %>%
  left_join(
    data.frame(
      cluster_alias   = as.integer(names(table(THAL$mmc_cluster_alias[THAL$TH_group == "TH_low"]))),
      n_cells_in_THAL = as.integer(table(THAL$mmc_cluster_alias[THAL$TH_group == "TH_low"]))
    ),
    by = "cluster_alias"
  ) %>%
  group_by(region_of_interest_acronym) %>%
  summarise(
    weighted_pct_norm = sum(pct_norm * n_cells_in_THAL, na.rm = TRUE) /
      sum(n_cells_in_THAL, na.rm = TRUE),
    weighted_pct_raw  = sum(pct_raw  * n_cells_in_THAL, na.rm = TRUE) /
      sum(n_cells_in_THAL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(weighted_pct_norm))

cat("\nROI che definiscono le cellule TH_low (peso normalizzato):\n")
print(low_roi_agg, n = 30)


# --- 5. Plot -----------------------------------------------------------------

# Barplot ROI aggregate per TH_low
p_bar_low <- ggplot(low_roi_agg %>% filter(weighted_pct_norm > 0.5),
                    aes(x = reorder(region_of_interest_acronym, weighted_pct_norm),
                        y = weighted_pct_norm,
                        fill = region_of_interest_acronym == "TH")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2980B9", "FALSE" = "#E74C3C"),
                    guide = "none") +
  labs(title    = sprintf("ROI di origine per cellule TH_low (soglia %.0f%%)", TH_SOGLIA),
       subtitle = "% normalizzata pesata per n cellule nel tuo dataset",
       x = "", y = "% normalizzata pesata") +
  theme_classic(base_size = 11)
savepng(p_bar_low, "04E_barplot_TH_low_ROI.png", width = 3500, height = 2500)

# UMAP: evidenzia TH_low vs TH_high
p_umap_group <- DimPlot(
  THAL, group.by = "TH_group",
  cols   = c("TH_high" = "#2980B9", "TH_low" = "#E74C3C", "no_mapping" = "grey70"),
  pt.size = 0.3,
  order  = c("TH_low", "no_mapping", "TH_high")
) +
  ggtitle(sprintf("Cellule TH_low vs TH_high (soglia %.0f%% norm)", TH_SOGLIA)) +
  theme_classic(base_size = 12)
savepng(p_umap_group, "04E_UMAP_TH_low_high.png", width = 3500, height = 3000)

# Heatmap cluster TH_low x top ROI
top_rois <- low_roi_agg %>% slice_max(weighted_pct_norm, n = 15) %>%
  pull(region_of_interest_acronym)

heatmap_mat <- cluster_roi_norm %>%
  filter(cluster_alias %in% low_aliases,
         region_of_interest_acronym %in% top_rois) %>%
  select(cluster_alias, region_of_interest_acronym, pct_norm) %>%
  tidyr::pivot_wider(names_from  = region_of_interest_acronym,
                     values_from = pct_norm,
                     values_fill = 0) %>%
  tibble::column_to_rownames("cluster_alias") %>%
  as.matrix()

png(file.path(outdir, "04E_heatmap_TH_low_clusters.png"),
    width = max(3000, ncol(heatmap_mat) * 120 + 800),
    height = max(2000, nrow(heatmap_mat) * 80 + 600),
    res = 300)
pheatmap(heatmap_mat,
         color    = colorRampPalette(c("white", "#E74C3C"))(50),
         fontsize = 7,
         main     = "Composizione ROI (norm) — cluster TH_low")
dev.off()

write.csv(low_cluster_summary,
          file.path(outdir, "04E_TH_low_cluster_summary.csv"),
          row.names = FALSE)
write.csv(low_roi_agg,
          file.path(outdir, "04E_TH_low_ROI_aggregated.csv"),
          row.names = FALSE)

message("=== 04E completato ===")

# ==============================================================================
# FILTRO TH — tieni TH_high + cellule no_mapping (NA in cluster_pct_TH_norm)
# ==============================================================================
# TH_high  = pct_TH_norm >= TH_SOGLIA  → cellule talamiche vere
# no_mapping = cluster_pct_TH_norm è NA → oligodendrociti, astrociti ecc.
#              che MapMyCells non ha potuto classificare; vengono tenuti
#              perché la loro esclusione è basata sull'assenza di dati,
#              non su evidenza di contaminazione.
# TH_low   = pct_TH_norm < TH_SOGLIA  → escluse, analizzate in ST_05
# ==============================================================================

# Salva THAL completo (con TH_low) prima del filtro — serve per 04F
THAL_all_cells <- THAL

n_before  <- ncol(THAL)
keep_mask <- is.na(THAL$cluster_pct_TH_norm) |
             THAL$cluster_pct_TH_norm >= TH_SOGLIA

THAL_filt <- THAL[, keep_mask]
n_after   <- ncol(THAL_filt)

message(sprintf(
  "Filtro TH: %d -> %d cellule (rimossi %d TH_low, %.1f%%)",
  n_before, n_after, n_before - n_after,
  100 * (n_before - n_after) / n_before
))

cat("\nCellule per campione dopo filtro:\n")
print(table(THAL_filt$orig.ident))

cat("\nCellule per TH_group nel dataset filtrato:\n")
print(table(THAL_filt$TH_group, useNA = "always"))

cat("\nCellule per TH_group rimosse:\n")
print(table(THAL$TH_group[!keep_mask]))

THAL <- THAL_filt
rm(THAL_filt)

# Salva conteggio cellule TH_low per ROI PRIMA del filtro
# (serve per 04F scatter — dopo il filtro le cellule TH_low non esistono più in THAL)
low_cells_per_roi_counts <- table(THAL_all_cells$mmc_cluster_alias[
  !is.na(THAL_all_cells$TH_group) & THAL_all_cells$TH_group == "TH_low"
])


# ==============================================================================
# Lookup nomi + distanza fisica da TH per ogni ROI
# Fonte: Allen CCFv3 Structure API (ontologia ufficiale)
# ==============================================================================

library(httr)
library(dplyr)

# --- 1. Scarica ontologia CCFv3 via Allen API --------------------------------
# Restituisce tutte le strutture con acronimo, nome completo e centroide 3D (µm)

resp <- GET(
  "https://api.brain-map.org/api/v2/data/Structure/query.json",
  query = list(
    criteria    = "ontology[name$eq'Mouse Brain Atlas']",
    num_rows    = 2000,
    start_row   = 0,
    include     = "structure_id_path",
    # campi utili: acronym, name, id, centroid x/y/z
    only        = paste(
      "acronym", "name", "id",
      "neuro_name_structure_id_path",
      "safe_name",
      sep = ","
    )
  )
)

# Centroide non è restituito da questo endpoint — usiamo l'endpoint mesh/summary
# che ha coordinate 3D dei centroidi per ogni struttura

resp_centroid <- GET(
  "https://api.brain-map.org/api/v2/data/StructureMesh/query.json",
  query = list(
    num_rows = 2000,
    only     = "structure_id,storage_directory"
  )
)

# Endpoint più diretto: structure con average coordinate dalla reference space
ontology_resp <- GET(
  "http://api.brain-map.org/api/v2/structure_graph_download/1.json"
)

# Parsing: l'ontologia CCFv3 è un albero JSON ricorsivo
parse_structure_tree <- function(node) {
  row <- data.frame(
    id      = node$id      %||% NA_integer_,
    acronym = node$acronym %||% NA_character_,
    name    = if (!is.null(node$`safe-name`) && length(node$`safe-name`) > 0)
      node$`safe-name`
    else if (!is.null(node$name) && length(node$name) > 0)
      node$name
    else NA_character_,
    stringsAsFactors = FALSE
  )
  children <- node$children
  if (length(children) > 0) {
    child_rows <- bind_rows(lapply(children, parse_structure_tree))
    row <- bind_rows(row, child_rows)
  }
  row
}

# Operatore %||% se non hai rlang caricato
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

ont_data  <- content(ontology_resp, as = "parsed")
structure_df <- parse_structure_tree(ont_data$msg[[1]])

cat("Strutture caricate:", nrow(structure_df), "\n")

# --- 2. Lookup acronimi delle tue ROI ----------------------------------------
rois_to_lookup <- low_roi_agg$region_of_interest_acronym

# Alcune ROI ABCA sono combinazioni (es. "AUD-TEa-PERI-ECT") — splittale
# per il lookup e riagganciala al nome composto
roi_lookup <- data.frame(
  region_of_interest_acronym = rois_to_lookup,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    # prendi il primo acronimo del gruppo come chiave principale
    primary_acronym = strsplit(region_of_interest_acronym, "-")[[1]][1],
    full_name = {
      hit <- structure_df$name[structure_df$acronym == primary_acronym]
      if (length(hit) == 0) NA_character_ else hit[1]
    }
  ) %>%
  ungroup()

# Aggiunge nomi completi per i casi semplici; 
# per i ROI compositi (es. "AUD-TEa-PERI-ECT") costruisce nome concatenato
roi_lookup <- roi_lookup %>%
  rowwise() %>%
  mutate(
    full_name = if (!is.na(full_name)) full_name else {
      parts  <- strsplit(region_of_interest_acronym, "-")[[1]]
      hits   <- sapply(parts, function(a) {
        h <- structure_df$name[structure_df$acronym == a]
        if (length(h) == 0) a else h[1]
      })
      paste(hits, collapse = " / ")
    }
  ) %>%
  ungroup()

cat("\nTabella acronimo -> nome completo:\n")
print(roi_lookup %>% select(region_of_interest_acronym, full_name), n = 30)


# --- 3. Distanza fisica da TH ------------------------------------------------
# Usiamo le coordinate 3D dei centroidi da MERFISH (cell_metadata ha x, y)
# oppure calcoliamo i centroidi delle ROI direttamente dal cell_metadata dell'atlas
# che contiene le colonne x, y (coordinate CCFv3 2D per sezione)
# Per la distanza 3D usiamo la coppia (x, y) disponibile + assumiamo z medio

# Calcola centroide 2D per ogni ROI nel cell_metadata dell'atlas
roi_centroids <- meta %>%
  filter(!is.na(region_of_interest_acronym),
         !is.na(x), !is.na(y)) %>%
  group_by(region_of_interest_acronym) %>%
  summarise(
    centroid_x = median(x, na.rm = TRUE),
    centroid_y = median(y, na.rm = TRUE),
    n_cells    = n(),
    .groups    = "drop"
  )

cat("\nCentroidi ROI (prime 10):\n")
print(head(roi_centroids, 10))

# Centroide di TH
TH_centroid <- roi_centroids %>%
  filter(region_of_interest_acronym == "TH")

if (nrow(TH_centroid) == 0) stop("TH non trovato in roi_centroids — verifica il label")

cat(sprintf("\nCentroide TH: x=%.3f, y=%.3f\n",
            TH_centroid$centroid_x, TH_centroid$centroid_y))

# Distanza euclidea da TH per ogni ROI
roi_distances <- roi_centroids %>%
  mutate(
    dist_from_TH = sqrt(
      (centroid_x - TH_centroid$centroid_x)^2 +
        (centroid_y - TH_centroid$centroid_y)^2
    )
  ) %>%
  arrange(dist_from_TH)

cat("\nDistanza da TH per ROI (ordinate):\n")
print(roi_distances, n = 30)



# --- 4. Tabella finale: combina tutto ----------------------------------------
final_table <- low_roi_agg %>%
  left_join(roi_lookup   %>% select(region_of_interest_acronym, full_name),
            by = "region_of_interest_acronym") %>%
  left_join(roi_distances %>% select(region_of_interest_acronym,
                                     centroid_x, centroid_y, dist_from_TH),
            by = "region_of_interest_acronym") %>%
  arrange(desc(weighted_pct_norm)) %>%
  select(region_of_interest_acronym, full_name,
         weighted_pct_norm, weighted_pct_raw,
         dist_from_TH, centroid_x, centroid_y)

cat("\n=== Tabella finale: ROI, nome, arricchimento e distanza da TH ===\n")
print(final_table, n = 30)


# Aggiungi n_cells_THAL_in_low_clusters alla final_table
# Usa low_cells_per_roi_counts salvato prima del filtro TH_high
low_cells_per_roi <- cluster_roi_norm %>%
  filter(cluster_alias %in% low_aliases) %>%
  left_join(
    data.frame(
      cluster_alias   = as.integer(names(low_cells_per_roi_counts)),
      n_cells_in_THAL = as.integer(low_cells_per_roi_counts)
    ),
    by = "cluster_alias"
  ) %>%
  group_by(region_of_interest_acronym) %>%
  summarise(
    n_cells_THAL_weighted = sum(pct_norm / 100 * n_cells_in_THAL, na.rm = TRUE),
    .groups = "drop"
  )

final_table <- final_table %>%
  left_join(low_cells_per_roi, by = "region_of_interest_acronym")

write.csv(final_table,
          file.path(outdir, "04F_ROI_names_distance_TH.csv"),
          row.names = FALSE)


# Rimuovi TH stesso e righe con NA
cor_data <- final_table %>%
  filter(region_of_interest_acronym != "TH",
         !is.na(dist_from_TH),
         !is.na(n_cells_THAL_weighted),
         !is.na(weighted_pct_norm))

# Pearson (relazione lineare)
cor_pearson <- cor.test(cor_data$dist_from_TH, cor_data$n_cells_THAL_weighted,
                        method = "pearson")

# Spearman (più robusto a outlier e non-linearità)
cor_spearman <- cor.test(cor_data$dist_from_TH, cor_data$n_cells_THAL_weighted,
                         method = "spearman")

cat("=== Correlazione distanza da TH ~ n cellule THAL nei cluster TH_low ===\n\n")
cat(sprintf("Pearson:  r = %.3f,  p = %.4f\n",
            cor_pearson$estimate,  cor_pearson$p.value))
cat(sprintf("Spearman: rho = %.3f, p = %.4f\n",
            cor_spearman$estimate, cor_spearman$p.value))


# --- 5. Plot: arricchimento vs distanza da TH --------------------------------
p_scatter <- ggplot(
  final_table %>% filter(region_of_interest_acronym != "TH"),
  aes(x = dist_from_TH, y = n_cells_THAL_weighted,
      label = region_of_interest_acronym)
) +
  geom_point(aes(size = weighted_pct_norm),
             color = "#E74C3C", alpha = 0.7) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  scale_size_continuous(
    name  = "% ROI normalizzata\nnei cluster TH_low",
    range = c(2, 10)
  ) +
  labs(
    title    = "Cellule THAL a rischio vs distanza fisica da TH",
    subtitle = "Y = n cellule THAL attribuibili a quella ROI | Dimensione = % ROI normalizzata nei cluster TH_low",
    x        = "Distanza dal centroide TH (unità CCFv3)",
    y        = "N cellule THAL nei cluster TH_low"
  ) +
  theme_classic(base_size = 12)


cor_label <- sprintf("Pearson r = %.2f (p = %.3f)\nSpearman rho = %.2f (p = %.3f)",
                     cor_pearson$estimate,  cor_pearson$p.value,
                     cor_spearman$estimate, cor_spearman$p.value)

p_scatter <- p_scatter +
  annotate("text", x = Inf, y = Inf, label = cor_label,
           hjust = 1.05, vjust = 1.5, size = 3.2, color = "grey30")

savepng(p_scatter, "04F_scatter_enrichment_vs_distance.png",
        width = 3500, height = 3000)

savepng(p_scatter, "04F_scatter_enrichment_vs_distance.png",
        width = 3500, height = 3000)

message("=== 04F completato ===")


##############UMAP

roi_name_map <- setNames(structure_df$name, structure_df$acronym)

# ==============================================================================
# 2. Aggiungi colonna nome esteso al metadata
# ==============================================================================

roi_sigla <- THAL$cluster_top_roi

THAL$cluster_top_roi_name <- ifelse(
  roi_sigla %in% names(roi_name_map),
  roi_name_map[roi_sigla],
  roi_sigla   # fallback: sigla originale se non trovata nell'ontologia
)

missing_roi <- unique(roi_sigla[!roi_sigla %in% names(roi_name_map)])
if (length(missing_roi) > 0) {
  message("ROI non trovate in structure_df: ", paste(missing_roi, collapse = ", "))
}

cat("ROI presenti nel dataset:\n")
print(sort(table(THAL$cluster_top_roi_name), decreasing = TRUE))

# ==============================================================================
# 3. Palette colori
# ==============================================================================

name_levels <- sort(unique(na.omit(THAL$cluster_top_roi_name)))
n_roi        <- length(name_levels)

if (n_roi <= 8) {
  pal_cols <- RColorBrewer::brewer.pal(max(n_roi, 3), "Set2")
} else if (n_roi <= 12) {
  pal_cols <- RColorBrewer::brewer.pal(n_roi, "Paired")
} else {
  pal_cols <- scales::hue_pal()(n_roi)
}
roi_pal_named <- setNames(pal_cols[seq_along(name_levels)], name_levels)

# Thalamus sempre in blu della pipeline
thal_label <- roi_name_map["TH"]
if (!is.na(thal_label) && thal_label %in% name_levels) {
  roi_pal_named[thal_label] <- "#2980B9"
}

# ==============================================================================
# 4. UMAP
# ==============================================================================

front_labels <- name_levels[name_levels != thal_label]

p_roi_named <- DimPlot(
  THAL,
  group.by = "cluster_top_roi_name",
  cols     = roi_pal_named,
  pt.size  = 0.3,
  order    = front_labels
) +
  ggtitle("ROI predominante nel cluster MapMyCells") +
  labs(colour = "Area (CCFv3)") +
  theme_classic(base_size = 12) +
  theme(legend.text = element_text(size = 9)) +
  guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1))

savepng(p_roi_named, "04D_UMAP_top_roi_names.png",
        width = 4500, height = 3000)

save(THAL, cluster_TH_stats, cluster_roi_norm, low_aliases, clust_meta,
     file = file.path(outdir, "THAL_annotated_MMC.RData"))

###Doublecheck che i vari cluster mappino sulle regioni predette

















# --- 04B.1 — Classifica "other" per classi rare ----------------------------
class_freq <- table(mapping$class_name)
sub_freq   <- table(mapping$subclass_name)

keep_class <- names(class_freq[class_freq >= MIN_CELLS_CLASS])
keep_sub   <- names(sub_freq[sub_freq   >= MIN_CELLS_SUBCLASS])

mapping$class_clean    <- ifelse(mapping$class_name    %in% keep_class,
                                 mapping$class_name,    "other")
mapping$subclass_clean <- ifelse(mapping$subclass_name %in% keep_sub,
                                 mapping$subclass_name, "other")

# --- 04B.2 — Aggiungi annotazioni agli oggetti Seurat ----------------------
add_mapmycells <- function(obj, mapping_df, label) {
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

THAL <- add_mapmycells(THAL, mapping, "THAL")

# --- 04B.3 — Plot: cluster Seurat vs MapMyCells subclass -------------------
if ("mmc_subclass" %in% colnames(THAL@meta.data)) {

  ct_tab  <- table(Idents(THAL), THAL$mmc_subclass)
  ct_prop <- prop.table(ct_tab, margin = 1)
  n_sub   <- ncol(ct_prop)
  n_cl    <- nrow(ct_prop)

  png(file.path(outdir, paste0("04_cluster_vs_MMC_", obj_nm, ".png")),
      width  = max(3000, n_sub * 80 + 1000),
      height = max(2000, n_cl  * 80 + 600),
      res    = 300)
  pheatmap(
    ct_prop,
    color           = colorRampPalette(c("white", "#2980B9"))(50),
    display_numbers = FALSE,
    fontsize        = 7,
    main            = paste("Cluster Seurat vs MapMyCells subclass —", obj_nm)
  )
  dev.off()

  sub_levels <- sort(unique(THAL$mmc_subclass[!is.na(THAL$mmc_subclass)]))
  sub_real   <- sub_levels[sub_levels != "other"]
  pal_sub    <- setNames(
    c(scales::hue_pal()(length(sub_real)), "grey70"),
    c(sub_real, "other")
  )

  p_mmc <- DimPlot(THAL, group.by = "mmc_subclass",
                   label = TRUE, repel = TRUE, label.size = 2.5, pt.size = 0.3,
                   cols = pal_sub, order = c(sub_real, "other")) +
    ggtitle(paste("MapMyCells subclass —", obj_nm)) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"))
  savepng(p_mmc, paste0("04_UMAP_MMC_subclass_", obj_nm, ".png"),
          width = 6500, height = 3500)

  cl_levels <- sort(unique(THAL$mmc_class[!is.na(THAL$mmc_class)]))
  cl_real   <- cl_levels[cl_levels != "other"]
  pal_cl    <- setNames(
    c(scales::hue_pal()(length(cl_real)), "grey70"),
    c(cl_real, "other")
  )

  p_mmc_cl <- DimPlot(THAL, group.by = "mmc_class",
                      label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.3,
                      cols = pal_cl, order = c(cl_real, "other")) +
    ggtitle(paste("MapMyCells class —", obj_nm))
  savepng(p_mmc_cl, paste0("04_UMAP_MMC_class_", obj_nm, ".png"),
          width = 4000, height = 3000)
}

# --- 04B.4 — Distribuzione cellule per campione (diagnostica) ---------------
cat(sprintf("\n  Cellule per campione — %s:\n", obj_nm))
print(table(THAL$orig.ident))

# --- 04B.5 — UMAP contaminanti stimati da cluster_top_roi --------------------
# Colorazione basata sulla ROI di origine del cluster (pct_norm):
# TH = grigio (cellule talamiche attese), tutto il resto = colori distinti

if ("cluster_top_roi" %in% colnames(THAL@meta.data)) {

  roi_levels_all <- sort(unique(na.omit(THAL$cluster_top_roi)))
  non_TH_rois    <- roi_levels_all[roi_levels_all != "TH"]

  cat(sprintf("\nROI di origine nei cluster (pct_norm): %s\n",
              paste(roi_levels_all, collapse = ", ")))
  cat(sprintf("Cellule con top_roi == TH: %d (%.1f%%)\n",
              sum(THAL$cluster_top_roi == "TH", na.rm = TRUE),
              100 * mean(THAL$cluster_top_roi == "TH", na.rm = TRUE)))
  cat(sprintf("Cellule con top_roi != TH: %d (%.1f%%)\n",
              sum(THAL$cluster_top_roi != "TH", na.rm = TRUE),
              100 * mean(THAL$cluster_top_roi != "TH", na.rm = TRUE)))

  # Palette: TH in grigio chiaro, non-TH con colori distinti
  base_cols <- unique(c(
    RColorBrewer::brewer.pal(9,  "Set1"),
    RColorBrewer::brewer.pal(8,  "Dark2"),
    RColorBrewer::brewer.pal(12, "Paired")
  ))
  if (length(non_TH_rois) > length(base_cols))
    base_cols <- colorRampPalette(base_cols)(length(non_TH_rois))

  pal_roi_contam <- c(
    "TH" = "grey85",
    setNames(base_cols[seq_len(length(non_TH_rois))], non_TH_rois)
  )
  if (any(is.na(THAL$cluster_top_roi)))
    pal_roi_contam["NA"] <- "grey50"

  p_umap_roi_contam <- DimPlot(
    THAL,
    group.by = "cluster_top_roi",
    cols     = pal_roi_contam,
    pt.size  = 0.3,
    order    = c("TH", non_TH_rois)   # TH sotto, non-TH sopra
  ) +
    ggtitle(paste("ROI di origine cluster (pct_norm) —", obj_nm)) +
    labs(colour = "Top ROI") +
    theme_classic(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 2))

  savepng(p_umap_roi_contam,
          paste0("04B5_UMAP_roi_origin_", obj_nm, ".png"),
          width = 4500, height = 3500)

  # Barplot: % cellule per ROI di origine
  roi_tab <- as.data.frame(table(THAL$cluster_top_roi, useNA = "always")) %>%
    rename(roi = Var1, n_cells = Freq) %>%
    mutate(
      pct    = round(100 * n_cells / sum(n_cells), 2),
      is_TH  = roi == "TH",
      roi    = as.character(roi)
    ) %>%
    arrange(desc(n_cells))

  fill_pal_bar <- setNames(
    ifelse(roi_tab$is_TH, "grey85",
           ifelse(is.na(roi_tab$roi), "grey50", "#E74C3C")),
    roi_tab$roi
  )

  p_bar_roi <- ggplot(roi_tab,
                      aes(x = reorder(roi, n_cells), y = pct, fill = roi)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_manual(values = fill_pal_bar, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 10) +
    labs(
      title    = paste("ROI di origine cluster (pct_norm) —", obj_nm),
      subtitle = "Grigio = TH (atteso) | Rosso = altra ROI",
      x = "", y = "% cellule"
    )
  savepng(p_bar_roi, paste0("04B5_barplot_roi_origin_", obj_nm, ".png"),
          width = 3500, height = max(2000, nrow(roi_tab) * 120 + 700))

  write.csv(roi_tab,
            file.path(outdir, paste0("04B5_roi_origin_proportions_", obj_nm, ".csv")),
            row.names = FALSE)
}

# --- Salvataggio post-04B (prima del filtro TH) ------------------------------
save(THAL, file = file.path(outdir, "THAL_annotated.RData"))
message("Salvato: THAL_annotated.RData")


# ==============================================================================
# PARTE C — ASSEGNAZIONE ROI TALAMICA DA MAPMYCELLS CLUSTER
# ==============================================================================
# Riusa cluster_roi_norm già calcolata nella prima parte (coerenza con TH_high/low).
# La ROI assegnata è quella con pct_norm massima (normalizzata per dimensione ROI
# nell'atlas), non la ROI con più cellule in assoluto.
# ==============================================================================

message("=== 04C. Assegnazione ROI talamica da cluster MapMyCells ===")

# --- 1. Mapping cluster_alias -> top_roi (per pct_norm) ----------------------
# Riusa cluster_roi_norm dalla prima parte — stesso criterio di TH_high/TH_low
cluster_roi_map <- cluster_roi_norm %>%
  group_by(cluster_alias) %>%
  slice_max(pct_norm, n = 1, with_ties = FALSE) %>%
  select(cluster_alias,
         roi_majority    = region_of_interest_acronym,
         pct_majority    = pct_norm,       # normalizzata
         pct_majority_raw = pct_raw,       # grezza per confronto
         n_cells_in_roi  = n_cells) %>%
  ungroup()

cat("\nMapping cluster_alias -> ROI (top pct_norm):\n")
print(head(cluster_roi_map, 20))
cat("Cluster totali con ROI assegnata:", nrow(cluster_roi_map), "\n")

cat("\nPurezza del mapping (pct_norm nella ROI maggioritaria):\n")
print(summary(cluster_roi_map$pct_majority))
cat("Cluster con pct_norm >50:", sum(cluster_roi_map$pct_majority > 50), "\n")
cat("Cluster con pct_norm >25:", sum(cluster_roi_map$pct_majority > 25), "\n")

# --- 2. Assegna ROI a ogni cellula di THAL -----------------------------------
THAL$mmc_roi <- cluster_roi_map$roi_majority[
  match(THAL$mmc_cluster_alias, cluster_roi_map$cluster_alias)
]

# Purezza: pct_norm della ROI assegnata (coerente con soglia TH_high/low)
THAL$mmc_roi_purity <- cluster_roi_map$pct_majority[
  match(THAL$mmc_cluster_alias, cluster_roi_map$cluster_alias)
]

n_mapped   <- sum(!is.na(THAL$mmc_roi))
n_unmapped <- sum(is.na(THAL$mmc_roi))
cat(sprintf("\nCellule con ROI assegnata: %d (%.1f%%)\n",
            n_mapped, 100 * n_mapped / ncol(THAL)))
cat(sprintf("Cellule senza ROI:         %d (%.1f%%)\n",
            n_unmapped, 100 * n_unmapped / ncol(THAL)))

cat("\nDistribuzione ROI talamiche:\n")
print(table(THAL$mmc_roi, useNA = "always"))

# --- 3. Aggrega ROI -> macro-area talamica -----------------------------------
macro_areas <- list(
  Ventral     = c("VENT", "VPL", "VPM", "VPLpc", "VPMpc", "VM", "VAL", "VL"),
  Mediodorsal = c("MED", "MD", "CM", "PF", "CL", "PC", "IMD"),
  Anterior    = c("ANT", "AV", "AM", "AD", "LD"),
  Lateral     = c("LAT", "LP", "PO", "PoT"),
  Posterior   = c("POST", "PIL", "SPA", "SPFp", "SPFm", "SGN"),
  Geniculate  = c("GEN", "LGd", "LGv", "MG"),
  Midline     = c("MID", "RE", "RH", "Xi", "PT"),
  Epithalamus = c("EPITH", "MH", "LH", "SMT")
)

roi_to_macro <- stack(macro_areas)
colnames(roi_to_macro) <- c("roi", "macro_area")
roi_to_macro <- roi_to_macro[!duplicated(roi_to_macro$roi), ]

THAL$mmc_macro_area <- as.character(roi_to_macro$macro_area[
  match(THAL$mmc_roi, roi_to_macro$roi)
])
THAL$mmc_macro_area[is.na(THAL$mmc_macro_area)] <- "Unassigned"

cat("\nDistribuzione macro-aree talamiche:\n")
print(table(THAL$mmc_macro_area))

# --- 4. Plot diagnostici -----------------------------------------------------
area_colors <- c(
  Ventral     = "#E63946",
  Mediodorsal = "#457B9D",
  Anterior    = "#2D6A4F",
  Lateral     = "#F4A261",
  Posterior   = "#9B5DE5",
  Geniculate  = "#00B4D8",
  Midline     = "#F77F00",
  Epithalamus = "#588157",
  Unassigned  = "#CCCCCC"
)

p_roi <- DimPlot(THAL, group.by = "mmc_macro_area",
                 cols = area_colors, label = TRUE, repel = TRUE,
                 pt.size = 0.3) +
  ggtitle("Nucleo Talamico (da MapMyCells cluster -> ROI, pct_norm)") +
  theme_classic(base_size = 12)
savepng(p_roi, "04C_UMAP_mmc_macro_area_THAL.png", width = 4000, height = 3500)

roi_levels <- sort(unique(na.omit(THAL$mmc_roi)))
pal_roi    <- setNames(scales::hue_pal()(length(roi_levels)), roi_levels)

p_roi_fine <- DimPlot(THAL, group.by = "mmc_roi",
                      cols = pal_roi, label = TRUE, repel = TRUE,
                      pt.size = 0.3, label.size = 2.5) +
  ggtitle("ROI ABCA talamica (da MapMyCells cluster, pct_norm)") +
  theme_classic(base_size = 12)
savepng(p_roi_fine, "04C_UMAP_mmc_roi_fine_THAL.png", width = 5000, height = 3500)

if ("cell_type" %in% colnames(THAL@meta.data)) {
  comp_df <- THAL@meta.data %>%
    filter(!is.na(mmc_roi)) %>%
    dplyr::count(cell_type, mmc_roi) %>%
    group_by(cell_type) %>%
    mutate(pct = n / sum(n) * 100)

  p_bar <- ggplot(comp_df, aes(x = cell_type, y = pct, fill = mmc_roi)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Composizione ROI talamica per cell type (pct_norm)",
         x = "Cell Type", y = "% cellule", fill = "ROI") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  savepng(p_bar, "04C_barplot_roi_by_celltype_THAL.png", width = 4000, height = 2500)
}

cat("\nPurezza media del mapping ROI per cellula (pct_norm):",
    round(mean(THAL$mmc_roi_purity, na.rm = TRUE), 1), "\n")

# --- 5. Salvataggio ----------------------------------------------------------
save(THAL, file = file.path(outdir, "THAL_annotated_MMC_fin.RData"))

write.csv(
  THAL@meta.data %>%
    dplyr::select(any_of(c("mmc_class", "mmc_subclass", "mmc_cluster",
                    "mmc_cluster_alias", "mmc_roi", "mmc_roi_purity",
                    "mmc_macro_area"))),
  file.path(outdir, "04C_roi_assignments_THAL.csv"),
  quote = FALSE
)

message("=== 04C completato: ROI talamica assegnata a ", n_mapped,
        " cellule (", round(100 * n_mapped / ncol(THAL), 1), "%) ===")
