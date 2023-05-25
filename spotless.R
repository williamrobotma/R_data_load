library(Seurat)
library(SeuratData)
library(SeuratDisk)
# library(rhdf5)
# library(feather)
library(rlist)

standards_dir <- "./spotless/standards/"

if (!dir.exists(standards_dir)){
  dir.create(standards_dir)
}

url <- "https://zenodo.org/record/7781323/files/standards.tar.gz"
if (!file.exists(basename(url)) && !dir.exists(file.path(standards_dir, "reference/"))) {
  curl::curl_download(url, basename(url))
}

if (!dir.exists(file.path(standards_dir, "reference/"))) {
  untar(basename(url), exdir=standards_dir)
}

conv_2_h5ad <- function(gs_dir) {
  gs_fnames <- list.files(gs_dir, pattern=".rds")
  
  gs_slides <- list()
  for (name in gs_fnames) {
    gs_slides[[tools::file_path_sans_ext(name)]] <- readRDS(file.path(gs_dir, name))
  }
  
  gs_slides
  gs_seurat <- list()
  for (name in names(gs_slides)) {
    print(name)
    metadata <- list()
    for (subname in names(gs_slides[[name]])) {
      print(subname)
      if (subname == "counts") {
        counts <- gs_slides[[name]][[subname]]
      }
      else if (subname != "dataset_properties") {
        i<-length(metadata) + 1
        metadata[[i]] <- as.data.frame(gs_slides[[name]][[subname]])
        colnames(metadata[[i]]) <- paste(subname, colnames(metadata[[i]]), sep=".")
      }
    }
    gs_seurat[[name]] <- CreateSeuratObject(counts, meta.data = list.cbind(metadata))
  }
  
  
  for (name in names(gs_seurat)) {
    SaveH5Seurat(gs_seurat[[name]], filename = file.path(gs_dir, paste(name, ".h5Seurat", sep="")), overwrite = TRUE)
    Convert(file.path(gs_dir, paste(name, ".h5Seurat", sep="")), dest = "h5ad", overwrite = TRUE)
    unlink(file.path(gs_dir, paste(name, ".h5Seurat", sep="")))
  }
}

gs_dirs <- grep("gold", list.dirs(standards_dir), value=TRUE)
for (gs_dir in gs_dirs) {
  if (basename(gs_dir) != "reference") {
    print("Converting")
    print(gs_dir)
    conv_2_h5ad(gs_dir)
  }
}


gs_ref_dir <- "./spotless/standards/reference/"
gs_refnames <- list.files(gs_ref_dir, pattern="^gold.*.rds$")
for (name in gs_refnames) {
  name_sans_ext <- tools::file_path_sans_ext(name)
  name_seurat <- paste(name_sans_ext, ".h5Seurat", sep="")
  gs_ref <- readRDS(file.path(gs_ref_dir, name))
  SaveH5Seurat(gs_ref, filename = file.path(gs_ref_dir, name_seurat), overwrite = TRUE)
  Convert(file.path(gs_ref_dir, name_seurat), dest = "h5ad", overwrite = TRUE)
  unlink(file.path(gs_ref_dir, name_seurat))
}

