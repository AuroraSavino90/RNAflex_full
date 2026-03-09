################################################################################
##  ST_01_load_qc.R
##
##  Input:  dati Visium 10X (filtered_feature_bc_matrix.h5 per campione)
##  Output: experiment.merged.RData  — oggetto Seurat merged pre-processing
##          ST_01_QC_plots/          — plot QC per campione
##
##  Flusso:
##    1. Carica ogni sezione Visium
##    2. Fix coordinate integer
##    3. Merge in experiment.merged
##    4. QC plots (nUMI, nGene, %MT per spot)
##    5. Filtraggio spot e geni
##    6. Salvataggio
################################################################################

source("code/ST_pipeline/ST_00_config.R")

# ==============================================================================
# 1. Carica ogni sezione Visium
# ==============================================================================

message("=== ST_01. Caricamento sezioni Visium ===")
# In ST_01, aggiungi questo blocco PRIMA del lapply di caricamento

# In ST_01, prima del lapply di caricamento

# Load normale
experiment.slices <- lapply(SAMPLE_NAMES, function(sname) {
  dir_path <- SAMPLE_DIRS[[sname]]
  
  obj <- Load10X_Spatial(
    data.dir      = dir_path,
    filename      = "filtered_feature_bc_matrix.h5",
    filter.matrix = TRUE,
    slice         = sname
  )
  
  obj <- fix_coords(obj)
  obj <- SetIdent(obj, value = sname)
  obj <- RenameCells(
    obj,
    new.names = paste(
      sapply(strsplit(Cells(obj), split = "-"), "[[", 1),
      sname, sep = "-"
    )
  )
  obj
})
names(experiment.slices) <- SAMPLE_NAMES

# Assegna in environment globale per accesso individuale
for (sname in SAMPLE_NAMES) assign(sname, experiment.slices[[sname]])

message("  Sezioni caricate: ", paste(SAMPLE_NAMES, collapse = ", "))


# ==============================================================================
# 2. Merge in experiment.merged
# ==============================================================================

message("=== Merge in oggetto unico ===")

experiment.merged <- merge(
  experiment.slices[[1]],
  experiment.slices[2:length(experiment.slices)]
)

# Aggiungi metadato capture.area
area_vec <- rep(NA_character_, ncol(experiment.merged))
for (sname in SAMPLE_NAMES) {
  area_vec[grep(sname, Cells(experiment.merged))] <- sname
}
experiment.merged <- AddMetaData(experiment.merged,
                                 metadata = area_vec,
                                 col.name = "capture.area")

experiment.merged$capture.area <- factor(
  experiment.merged$capture.area, levels = SAMPLE_ORDER
)

cat(sprintf("\nexperiment.merged: %d spot, %d geni\n",
            ncol(experiment.merged), nrow(experiment.merged)))
print(table(experiment.merged$capture.area))


# ==============================================================================
# 3. QC metrics
# ==============================================================================

message("=== Calcolo metriche QC ===")

experiment.merged[["pct_mt"]] <- PercentageFeatureSet(
  experiment.merged, pattern = "^mt-", assay = "Spatial"
)

# Riepilogo per campione
qc_summary <- experiment.merged@meta.data %>%
  group_by(capture.area) %>%
  summarise(
    n_spots       = n(),
    median_nUMI   = median(nCount_Spatial),
    median_nGene  = median(nFeature_Spatial),
    median_pct_mt = median(pct_mt),
    pct_low_UMI   = round(100 * mean(nCount_Spatial < MIN_UMI_SPOT), 1),
    pct_high_mt   = round(100 * mean(pct_mt > MAX_MT_PCT), 1),
    .groups = "drop"
  )

cat("\nRiepilogo QC per campione:\n")
print(qc_summary)
write.csv(qc_summary, file.path(outdir, "ST_01_QC_summary.csv"), row.names = FALSE)


# ==============================================================================
# 4. QC plots
# ==============================================================================

message("=== QC plots ===")

dir.create(file.path(outdir, "ST_01_QC_plots"), showWarnings = FALSE)

# Violin nUMI, nGene, %MT per campione
p_vln <- VlnPlot(
  experiment.merged,
  features  = c("nCount_Spatial", "nFeature_Spatial", "pct_mt"),
  group.by  = "capture.area",
  pt.size   = 0,
  ncol      = 3
) & theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

savepng(p_vln,
        file.path(outdir, "ST_01_QC_plots", "QC_violin.png"),
        width = 5000, height = 2500)

# Scatter nUMI vs nGene per campione
p_scatter_qc <- ggplot(
  experiment.merged@meta.data,
  aes(x = nCount_Spatial, y = nFeature_Spatial,
      colour = pct_mt, alpha = 0.5)
) +
  geom_point(size = 0.5) +
  facet_wrap(~ capture.area, ncol = 4) +
  scale_colour_viridis_c(name = "% MT") +
  geom_vline(xintercept = MIN_UMI_SPOT, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = MIN_GENES_SPOT, linetype = "dashed", colour = "red") +
  labs(title = "QC: nUMI vs nGene (rosso = soglie filtraggio)",
       x = "nUMI", y = "nGene") +
  theme_classic(base_size = 9) +
  guides(alpha = "none")

savepng(p_scatter_qc,
        file.path(outdir, "ST_01_QC_plots", "QC_scatter_nUMI_nGene.png"),
        width = 6000, height = 3000)

# Istogramma nUMI per campione
p_hist_umi <- ggplot(experiment.merged@meta.data,
                     aes(x = nCount_Spatial, fill = capture.area)) +
  geom_histogram(bins = 80, colour = "white", linewidth = 0.2) +
  facet_wrap(~ capture.area, ncol = 4, scales = "free_y") +
  geom_vline(xintercept = MIN_UMI_SPOT, linetype = "dashed", colour = "red") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Distribuzione nUMI per campione",
       x = "nUMI per spot", y = "N spot") +
  theme_classic(base_size = 9)

savepng(p_hist_umi,
        file.path(outdir, "ST_01_QC_plots", "QC_hist_nUMI.png"),
        width = 6000, height = 2500)


# ==============================================================================
# 5. Filtraggio spot e geni
# ==============================================================================

message("=== Filtraggio QC ===")

n_before <- ncol(experiment.merged)

selected_c <- colnames(experiment.merged)[
  experiment.merged$nCount_Spatial   >= MIN_UMI_SPOT  &
  experiment.merged$nFeature_Spatial >= MIN_GENES_SPOT &
  experiment.merged$pct_mt           <= MAX_MT_PCT
]

experiment.merged <- JoinLayers(experiment.merged)

selected_f <- rownames(experiment.merged)[
  rowSums(GetAssayData(experiment.merged, assay = "Spatial", layer = "counts") > 0) >= 3
]

# Rimuovi geni contaminanti
selected_f <- setdiff(selected_f, CONTAMINATED_GENES)

experiment.merged <- subset(experiment.merged,
                            features = selected_f,
                            cells    = selected_c)

n_after <- ncol(experiment.merged)
message(sprintf(
  "  Spot: %d -> %d (rimossi %d, %.1f%%)",
  n_before, n_after, n_before - n_after,
  100 * (n_before - n_after) / n_before
))
message(sprintf("  Geni: %d", nrow(experiment.merged)))

# Aggiorna oggetti individuali con i filtri del merged
for (sname in SAMPLE_NAMES) {
  obj <- get(sname)
  obj_filt <- subset(
    obj,
    features = selected_f,
    cells    = intersect(Cells(obj), selected_c)
  )
  assign(sname, obj_filt)
  message(sprintf("  %s: %d spot dopo filtraggio", sname, ncol(obj_filt)))
}


# ==============================================================================
# 6. Salvataggio
# ==============================================================================

save(experiment.merged,
     file = file.path(outdir, "ST_01_experiment_merged.RData"))

# Salva anche gli slice individuali filtrati (servono per GetTissueCoordinates)
for (sname in SAMPLE_NAMES) {
  assign(sname, get(sname))
}
save(list = SAMPLE_NAMES,
     file = file.path(outdir, "ST_01_slices_filtered.RData"))

message("=== ST_01 completato ===")
message("  Salvati: ST_01_experiment_merged.RData, ST_01_slices_filtered.RData")
