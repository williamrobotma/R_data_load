library(Seurat)
library(SeuratData)
library(SeuratDisk)
# library(rhdf5)
# library(feather)
library(rlist)

data_dir <- "./data/"

if (!dir.exists(data_dir)) {
    dir.create(data_dir)
}

url <- "https://zenodo.org/records/7750930/files/data.zip"
dest_file <- file.path(data_dir, basename(url))
dest_dir <- file.path(data_dir, "data/")
if (!file.exists(dest_file) && !dir.exists(dest_dir)) {
    curl::curl_download(url, dest_file, quiet=FALSE)
}

if (!dir.exists(dest_dir)) {
    unzip(dest_file, exdir = data_dir)
    unlink(dest_file)
}

gs_fnames <- list.files(dest_dir, pattern = ".rds")

gs_slides <- list()
for (name in gs_fnames) {
    gs_slides[[tools::file_path_sans_ext(gsub('(.*)_\\w+', '\\1', name))]] <- readRDS(file.path(dest_dir, name))
}

gs_seurat <- list()
for (name in names(gs_slides)) {
    print(name)
    metadata <- list()
    # for (subname in names(gs_slides[[name]])) {
    #     print(subname)
    #     if (subname == "counts") {
    #         counts <- gs_slides[[name]][[subname]]
    #     } else if (subname != "dataset_properties") {
    #         i <- length(metadata) + 1
    #         metadata[[i]] <- as.data.frame(gs_slides[[name]][[subname]])
    #         colnames(metadata[[i]]) <- paste(subname, colnames(metadata[[i]]), sep = ".")
    #     }
    # }

    # gs_seurat[[name]] <- CreateSeuratObject(counts, meta.data = list.cbind(metadata))
    # assay <- "RNA"
    # if (!is.null(mainExpName(gs_slides[[name]]))) {
    #     assay <- mainExpName(gs_slides[[name]])
    # }
    counts <- "counts"
    data<-"logcounts"
    if (!("logcounts" %in% names(gs_slides[[name]]@assays))) {
        data<-NULL
    }
    if ("counts" %in% names(gs_slides[[name]]@assays)) {
        gs_seurat[[name]] <- as.Seurat(gs_slides[[name]], counts =counts, data=data, assay=NULL)
        print("done as seurat")
        # make it a SeuratV4 object so that SeuratDisk works
        assay <- gs_seurat[[name]]@active.assay
        gs_seurat[[name]][[assay]] <- CreateAssayObject(counts = gs_seurat[[name]][[assay]]$counts)
    }



}



for (name in names(gs_seurat)) {
    out_path <- file.path(dest_dir, paste(name, ".h5Seurat", sep = ""))
    print(out_path)
    SaveH5Seurat(gs_seurat[[name]], filename = out_path, overwrite = TRUE)
    Convert(out_path, dest = "h5ad", overwrite = TRUE)
    unlink(out_path)
}

