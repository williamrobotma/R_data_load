library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(spacexr)
library(rhdf5)

data_dir <- "preprocessed_data"
DSETS = c("dlpfc", "pdac", "spotless")

MAX_CORES = 16


for (dset in DSETS) {
    dset_dir = file.path(data_dir, dset, "all/raw_counts")

    print(paste("Running RCTD for", dset))

    # SC
    print("Loading SC counts")
    # Loading HDF5 in R already transposes because R uses FORTRAN indexing
    counts <- h5read(file.path(dset_dir, "sc_train.h5"), "X")
    rownames(counts) <- h5read(file.path(dset_dir, "sc_train.h5"), "genes")
    colnames(counts) <- h5read(file.path(dset_dir, "sc_train.h5"), "samples")
    counts <- data.frame(counts)

    print("Loading SC cell types")
    cell_types <- h5read(file.path(dset_dir, "sc_train.h5"), "cell_type")
    names(cell_types) <- colnames(counts)
    cell_types <- as.factor(cell_types)

    print("Generating reference")
    reference <- Reference(counts, cell_types)

    print(dim(reference@counts))
    table(reference@cell_types)
    saveRDS(reference, file.path(dset_dir, 'SCRef.rds'))


    # ST
    print("Loading ST counts")
    # Loading HDF5 in R already transposes because R uses FORTRAN ordering
    counts <- h5read(file.path(dset_dir, "unscaled/st_train.h5"),"X")
    rownames(counts) <- h5read(file.path(dset_dir, "unscaled/st_train.h5"), "genes")
    colnames(counts) <- h5read(file.path(dset_dir, "unscaled/st_train.h5"), "samples")
    counts <- data.frame(counts)

    print("Loading ST coords")
    # Transpose back from FORTRAN ordering
    coords <- t(h5read(file.path(dset_dir, "unscaled/st_train.h5"), "coords"))
    rownames(coords) <- colnames(counts)
    colnames(coords) <- c("x", "y")
    coords <- data.frame(coords)

    ### Create SpatialRNA object
    print("Generating puck")
    puck <- SpatialRNA(coords, counts)

    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts)) # observe Digital Gene Expression matrix

    print(head(puck@coords)) # start of coordinate data.frame
    barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names).

    myRCTD <- create.RCTD(puck, reference, max_cores=MAX_CORES)
    myRCTD <- run.RCTD(myRCTD, doublet_mode='full')
    results <- myRCTD@results

    norm_weights = normalize_weights(results$weights)

    # cell_type_names <- myRCTD@cell_type_info$info[[2]]
    h5createFile(file.path(dset_dir, "st_predictions.h5"))
    h5write(
        t(as.matrix(norm_weights)), # transposing to C-order
        file.path(dset_dir, "st_predictions.h5"),
        "predictions",
    )

    write.table(
        as.matrix(colnames(norm_weights)),
        file=file.path(dset_dir, "pred_columns.csv"),
        row.names=F,
        col.names=F,
    )
    write.table(
        as.matrix(rownames(norm_weights)),
        file=file.path(dset_dir, "pred_rows.csv"),
        row.names=F,
        col.names=F,
    )

}
