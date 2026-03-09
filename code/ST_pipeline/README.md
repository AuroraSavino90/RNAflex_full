# Pipeline Spatial Transcriptomics — LSD Thalamus

Pipeline R modulare per l'analisi di dati Visium 10X Spatial Transcriptomics.
Il focus finale è mappare i cluster TH_low (dalla pipeline scRNA MapMyCells) sui
dati spaziali e testare se la loro localizzazione ST corrisponde alla ROI
"scorretta" identificata dall'atlas ABCA.

---

## Struttura file

| Script | Input | Output principale |
|---|---|---|
| `ST_00_config.R` | — | Variabili globali, helper |
| `ST_01_load_qc.R` | Visium .h5 | `ST_01_experiment_merged.RData` |
| `ST_02_norm_cluster.R` | ST_01 | `ST_02_*clustered.RData`, `df_toplot` |
| `ST_03_label_transfer_ABA.R` | ST_02, ABA reference | `area.spot` per spot |
| `ST_04_sc_transfer_thalamus.R` | ST_03, THAL scRNA | Subcluster talamici ST |
| `ST_05_TH_low_ST_mapping.R` | ST_04, THAL_MMC | Test + plot TH_low mapping |
| `ST_RUN_ALL.R` | — | Esegue tutti gli step |

---

## Struttura directory attesa

```
progetto/
├── ST_00_config.R
├── ST_01_load_qc.R
├── ...
├── data/
│   ├── Exp1/
│   │   ├── A1/filtered_feature_bc_matrix.h5   (LSD90min_1)
│   │   ├── B1/                                 (CTRL90min_1)
│   │   ├── C1/                                 (LSD24h_1)
│   │   └── D1/                                 (CTRL24h_1)
│   ├── Exp2/
│   │   ├── A1/ ... D1/
│   └── Reference/
│       ├── expr_raw_counts_table.tsv
│       └── meta_table.tsv
└── results/
    └── YYYYMMDD/   ← creata automaticamente
```

---

## Parametri da modificare in ST_00_config.R

- `BASE_DIR`, `DATA_DIR` — percorsi progetto
- `SAMPLE_DIRS` — mappa campione -> cartella Visium
- `ABA_SECTIONS` — sezioni atlas vicine alle tue sezioni sperimentali
- `MIN_UMI_SPOT`, `MIN_GENES_SPOT`, `MAX_MT_PCT` — soglie QC
- `CONTAMINATED_GENES` — geni da escludere

## Parametri da modificare in ST_02

- Operazioni di mirror coordinate (sezione 5) — verifica con il plot di orientamento

## Parametri da modificare in ST_04

- `THALAMIC_CLUSTERS_ST` — cluster Seurat talamici (suggeriti automaticamente,
  ma verifica visiva raccomandata)
- `THAL_ANNOTATED_PATH` — percorso THAL_annotated.RData

## Parametri da modificare in ST_05

- `roi_to_areaspot` — mapping sigla ROI ABCA -> label area.spot nel tuo dataset
  (verifica dopo `unique(df$area.spot)`)
- `MMC_PATH` — percorso THAL_annotated_MMC.RData

---

## Dipendenze R

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "ggrepel",
                   "pheatmap", "Matrix", "scales", "RColorBrewer"))

# Bioconductor
BiocManager::install(c("harmony", "spacexr"))
```
