library(Seurat)
library(SeuratData)
library(SeuratDisk)
# library(rhdf5)
# library(feather)
library(rlist)

gs_1_dir <- "./spotless/gold_standard_1"

gs_1_fnames <- list.files(gs_1_dir, pattern="*.rds")

gs_1_slides <- list()
for (name in gs_1_fnames) {
  gs_1_slides[[tools::file_path_sans_ext(name)]] <- readRDS(file.path(gs_1_dir, name))
}

gs_1_slides
gs_1_seurat <- list()
for (name in names(gs_1_slides)) {
  metadata <- list()
  for (subname in names(gs_1_slides[[name]])) {
    if (subname == "counts") {
      counts <- gs_1_slides[[name]][[subname]]
    }
    else if (subname != "dataset_properties") {
      i<-length(metadata) + 1
      metadata[[i]] <- as.data.frame(gs_1_slides[[name]][[subname]])
      # print(metadata[[i]])
      # print(colnames(metadata[[i]]))
      colnames(metadata[[i]]) <- paste(subname, colnames(metadata[[i]]), sep=".")
    }
  }
  gs_1_seurat[[name]] <- CreateSeuratObject(counts, meta.data = list.cbind(metadata))
}


for (name in names(gs_1_seurat)) {
  SaveH5Seurat(gs_1_seurat[[name]], filename = file.path(gs_1_dir, paste(name, ".h5seurat", sep="")), overwrite = TRUE)
  Convert(file.path(gs_1_dir, paste(name, ".h5seurat", sep="")), dest = "h5ad", overwrite = TRUE)
  unlink(file.path(gs_1_dir, paste(name, ".h5seurat", sep="")))
}

