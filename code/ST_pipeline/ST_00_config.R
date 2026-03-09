################################################################################
##  ST_00_config.R
##
##  Configurazione centralizzata per la pipeline Spatial Transcriptomics.
##  Source questo file all'inizio di ogni script della pipeline ST.
################################################################################

# ==============================================================================
# LIBRERIE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(Matrix)
  library(harmony)
  library(spacexr)
  library(RColorBrewer)
  library(scales)
})

# ==============================================================================
# PERCORSI
# ==============================================================================

# <<< Modifica questi percorsi in base alla tua struttura directory >>>

BASE_DIR    <- "D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/Psychedelics/3-Spatial transcriptomics/3_Data_analysis/LSD_Spatial_Transcriptomics"                          # root del progetto ST
DATA_DIR    <- file.path(BASE_DIR, "data")  # dati raw e reference
RESULTS_DIR <- file.path(BASE_DIR, "results")

# Cartella risultati datata
date    <- "20260304"
outdir  <- file.path(RESULTS_DIR, date)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Percorsi dati processati 10X Visium
# Struttura attesa:
#   data/Exp1/A1/filtered_feature_bc_matrix.h5
#   data/Exp1/B1/...  ecc.
PROCESSED_DIR <- "D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/Psychedelics/3-Spatial transcriptomics/2_Processed_data"
EXP1_DIR <- file.path(PROCESSED_DIR, "Exp1")
EXP2_DIR <- file.path(PROCESSED_DIR, "Exp2")


# Reference atlas ABA per label transfer anatomico
ABA_REF_COUNTS <- file.path(DATA_DIR, "Reference", "expr_raw_counts_table.tsv")
ABA_REF_META   <- file.path(DATA_DIR, "Reference", "meta_table.tsv")

# Geni contaminanti (specifici per il tuo dataset)
CONTAMINATED_GENES <- c("Ttr", "Pmch", "Lars2", "Prkcd", "Enpp2", "Nrgn")

# Sezioni ABA da usare come reference (vicine alle tue sezioni sperimentali)
ABA_SECTIONS <- c(paste(20:29, "A", sep = ""), paste(20:29, "B", sep = ""))

# ==============================================================================
# CAMPIONI
# ==============================================================================

# Mappa: nome campione -> sottocartella Visium
SAMPLE_DIRS <- list(
  LSD90min_1  = file.path(EXP1_DIR, "A1"),
  CTRL90min_1 = file.path(EXP1_DIR, "B1"),
  LSD24h_1    = file.path(EXP1_DIR, "C1"),
  CTRL24h_1   = file.path(EXP1_DIR, "D1"),
  LSD90min_2  = file.path(EXP2_DIR, "A1"),
  CTRL90min_2 = file.path(EXP2_DIR, "B1"),
  LSD24h_2    = file.path(EXP2_DIR, "C1"),
  CTRL24h_2   = file.path(EXP2_DIR, "D1")
)

SAMPLE_NAMES <- names(SAMPLE_DIRS)

# Ordine fattore per i plot
SAMPLE_ORDER <- c("LSD90min_1", "CTRL90min_1", "LSD90min_2", "CTRL90min_2",
                  "LSD24h_1",   "CTRL24h_1",   "LSD24h_2",   "CTRL24h_2")

# ==============================================================================
# PARAMETRI QC
# ==============================================================================

MIN_UMI_SPOT    <- 500     # UMI minime per spot
MIN_GENES_SPOT  <- 300     # geni minimi per spot
MAX_MT_PCT      <- 25      # % mitocondri massima

MIN_UMI_REF     <- 1000    # per reference ABA
MIN_EXPR_REF    <- 100     # geni con expr totale > soglia nel reference

# ==============================================================================
# PARAMETRI ANALISI
# ==============================================================================

N_PCS          <- 30       # componenti PCA per UMAP/clustering/transfer
HARMONY_DIMS   <- 1:30     # dims Harmony
UMAP_DIMS      <- 1:30

# Soglia per "ambiguous" nel label transfer anatomico
LABEL_TRANSFER_MIN_SCORE <- 0.5

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Salva PNG con dimensioni standard
savepng <- function(p, filename, width = 3600, height = 3000, res = 300) {
  path <- if (grepl("/", filename)) filename
          else file.path(outdir, filename)
  png(path, width = width, height = height, res = res)
  print(p)
  dev.off()
  invisible(path)
}

# Operatore null-coalescing
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Forza integer le coordinate immagine di un oggetto Seurat Spatial
# In ST_00_config.R — sostituisci fix_coords con questa versione

fix_coords <- function(obj) {
  img <- obj@images[[1]]
  
  if (inherits(img, "VisiumV2")) {
    # Seurat v5: le coordinate sono in @boundaries$centroids@coords
    coords <- img@boundaries$centroids@coords
    coords[, 1] <- as.numeric(coords[, 1])
    coords[, 2] <- as.numeric(coords[, 2])
    obj@images[[1]]@boundaries$centroids@coords <- coords
    
  } else if (inherits(img, "VisiumV1") || inherits(img, "SlideSeq")) {
    # Seurat v4: le coordinate sono in @coordinates
    for (coord in c("tissue", "row", "col", "imagerow", "imagecol")) {
      obj@images[[1]]@coordinates[[coord]] <-
        as.integer(obj@images[[1]]@coordinates[[coord]])
    }
  }
  obj
}

message("ST_00_config.R caricato — outdir: ", outdir)
