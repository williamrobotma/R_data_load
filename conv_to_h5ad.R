renv::init(bare=TRUE, bioconductor = TRUE)
if (!require("BiocManager", quietly = TRUE))
  renv::install("BiocManager")
renv::install("bioc::BiocVersion")
################################################################################
renv::install("gitcreds")
renv::install(c("Seurat", "mojaveazure/seurat-disk", "satijalab/seurat-data", "bioc::rhdf5", "jsonlite", "bioc::spatialLIBD", "rlist", "dmcable/spacexr"))
renv::snapshot()

library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Chijimatsu, Ryota. (2022). Integrated Data of Single cell RNA sequencing for
# Human Pancreatic Adenocarcinoma (Version 1) [Data set]. Zenodo.
# https://doi.org/10.5281/zenodo.6024273


# download
url <- "https://zenodo.org/records/6024273/files/pk_all.rds"
# fname <- "pk_all.rds"
fname <- basename(url)
if (!file.exists(basename(url))) {
    curl::curl_download(url, fname, quiet=FALSE)
}

# convert to h5ad
pk_all <- readRDS(fname)
SaveH5Seurat(pk_all, filename = "pk_all.h5Seurat")
Convert("pk_all.h5Seurat", dest = "h5ad")
unlink("pk_all.h5Seurat")
