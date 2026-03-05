################################################################################
##  00_config.R — Configurazione condivisa                                    ##
##  Source questo file all'inizio di ogni sezione:                            ##
##    source("00_config.R")                                                   ##
################################################################################

# ==============================================================================
# LIBRERIE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SuperCell)
  library(DESeq2)
  library(fgsea)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(anndata)
  library(openxlsx)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(Matrix)
  library(cluster)
  library(matrixStats)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

# ==============================================================================
# IMPOSTAZIONI GLOBALI
# ==============================================================================

set.seed(42)
options(Seurat.object.assay.version = "v5")

# Output directory — uguale per tutti gli script
date_str <- format(Sys.Date(), "%Y%m%d")
outdir   <- "results/20260224"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PALETTE CONDIVISE
# ==============================================================================

pal_trt  <- c("Ctrl" = "#3498DB", "LSD"       = "#E74C3C")
pal_time <- c("90min" = "#F39C12", "24h"      = "#8E44AD")
pal_area <- c("Cortex" = "#27AE60", "Thalamus" = "#2E86C1")

# Palette discreta per tipi cellulari (fino a 20 CT)
pal_ct <- c(
  RColorBrewer::brewer.pal(8,  "Set2"),
  RColorBrewer::brewer.pal(8,  "Set1"),
  RColorBrewer::brewer.pal(4,  "Dark2")
)

# ==============================================================================
# HELPER: salva PNG con impostazioni standard
# ==============================================================================

savepng <- function(p, filename, width = 2500, height = 2000, res = 300) {
  png(file.path(outdir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

# ==============================================================================
# HELPER: filtra geni presenti in un oggetto Seurat
# ==============================================================================

filter_genes <- function(genes, obj) intersect(genes, rownames(obj))

# ==============================================================================
# GENI MARCATORI — usati in sezioni 04, 05, 11, 14
# ==============================================================================

markers_broad <- c(
  "Mlc1","Gjb6","Aqp4","Acsbg1","Aldoc","S100b",
  "Ermn","Opalin","Mog","Aspa","Mobp","Ptgds","Plp1","Mbp",
  "Pdgfra","Cacng4","Tmem100","Matn4","Ptprz1","Cspg4",
  "C1qa","C1qb","C1qc","Ctss","P2ry12","Cx3cr1","Tmem119",
  "Snap25","Neurod6","Slc17a7","Slc17a6",
  "Gad1","Gad2","Slc32a1",
  "Flt1","Cldn5","Tm4sf1","Pecam1","Rgs5",
  "Neu4","Bmp4","Enpp6","Bcas1"
)

markers_exc_ctx <- c(
  "Rasgrf2","Calb1","Cux1","Cux2","Otof","Satb2",
  "Rorb","Rspo1","Il1rapl2",
  "Galnt14","Sulf1","Bcl6","Hs3st2",
  "Fezf2","Bcl11b","Ldb2","Chrna6","Slc17a8","Pld5","Coro6","Parm1",
  "Foxp2","Syt6","Hs3st4","Nxph4","Tshz2","Nr4a2","Ctgf",
  "Ostn","Hpcal1",
  "Car3","Sla2","Fam84b"
)

markers_inh_ctx <- c(
  "Pvalb","Kcnc1","Syt2",
  "Sst","Chodl","Npy","Hpse","Tac2",
  "Vip","Calb2","Cck","Tac1",
  "Lamp5","Ndnf","Lhx6","Nkx2-1",
  "Sncg","Cnr1","Reln",
  "Nr2f2","Prox1","Adarb2"
)

markers_thal <- c(
  "Slc17a6","Rorb","Grik4",
  "Tbr1","Pou4f1","Calb1","Calb2",
  "Nrgn","Kcnab2","Hcn1","Cacna1g",
  "Prkcg","Grid2",
  "Nr4a2","Lhx2","Gbx2",
  "Gad1","Gad2","Pvalb","Sst"
)

markers_oligo_fine <- c(
  "Pdgfra","Cspg4","Sox10","Olig2",
  "Bcas1","Bmp4","Neu4","Enpp6",
  "Mog","Mobp","Plp1","Mbp","Aspa",
  "Opalin","Ermn","Ptgds"
)

markers_htr2a_signature <- c(
  "Htr2a",
  "Fos","Arc","Egr1","Nr4a1","Npas4","Fosb","Junb","Dusp1","Dusp5",
  "Bdnf","Ntrk2",
  "Gnaq","Plcb1","Prkca",
  "Drd1","Drd2","Drd5",
  "Slc6a4","Tph2",
  "Camk2a","Camk4"
)

message("00_config.R caricato.")
