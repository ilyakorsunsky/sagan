#' @export
Read10X_Image_v2 <- function (
    image.dir, 
    image.name, 
    filter.matrix = TRUE, ...
) {
    res.type <- gsub('.*(hires|lowres)_image.png', '\\1', image.name)
    
    image <- readPNG(source = file.path(image.dir, image.name))
    scale.factors <- fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
    tissue.positions <- read.csv(file = file.path(image.dir, 
        "tissue_positions_list.csv"), col.names = c("barcodes", 
        "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, 
        as.is = TRUE, row.names = 1)
    if (filter.matrix) {
        tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
    }
    
    if (res.type == 'hires') {
        res.diff <- scale.factors$tissue_hires_scalef / scale.factors$tissue_lowres_scalef
        scale.factors$tissue_hires_scalef <- scale.factors$tissue_hires_scalef * res.diff 
        scale.factors$tissue_lowres_scalef <- scale.factors$tissue_lowres_scalef * res.diff
    } else if (res.type == 'lowres') {
      ## do nothing, this is the default Seurat expects   
    } else {
        stop('Could not resolve image resolution type')
    }
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
        scale.factors$tissue_lowres_scalef        
    
    spot.radius <- unnormalized.radius/max(dim(x = image))
    return(new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
        fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
        scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, 
        spot.radius = spot.radius))
}
environment(Read10X_Image_v2) <- environment(Seurat::Read10X_Image)


# ## This function generalizes Read10X and Read10X_h5, so you don't need to specify which one to use 
# ## TODO: decide whether we would ever want to use the unfiltered .h5 file?
# Read10X_generic <- function(data.dir, assay, to.upper=FALSE, ...) {    
#     ## Search for putative filenames
#     filename_h5 <- normalizePath(list.files(data.dir, pattern='filtered.*h5$', recursive=TRUE, full=TRUE))
#     filename_filtered_mtx <- normalizePath(list.files(data.dir, pattern='filtered.*mtx', recursive=TRUE, full=TRUE))
#     filename_mtx <- normalizePath(list.files(data.dir, pattern='mtx', recursive=TRUE, full=TRUE))
    
    
#     if (length(filename_h5) > 0) {
#         if (length(filename_h5) > 1) stop(glue('Multiple h5 files found in {data.dir}'))
#         data <- Read10X_h5(filename = filename_h5, ...)
#     } else if (length(filename_filtered_mtx)) {
#         if (length(filename_filtered_mtx) > 1) stop(glue('Multiple filtered mtx files found in {data.dir}'))
#         filtered_mtx_dir <- paste(head(strsplit(filename_filtered_mtx, '/')[[1]], -1), collapse='/')
#         data <- Read10X(data.dir = filtered_mtx_dir, ...)

#     } else if (length(filename_mtx)) {
#         if (length(filename_mtx) > 1) stop(glue('Multiple mtx files found in {data.dir}'))
#         mtx_dir <- paste(head(strsplit(filename_mtx, '/')[[1]], -1), collapse='/')
#         data <- Read10X(data.dir = mtx_dir, ...)
#     } else {
#         stop(glue('Could not find any h5 or mtx files in {data.dir}'))
#     }
#     if (to.upper) {
#         rownames(x = data) <- toupper(x = rownames(x = data))
#     }
#     object <- CreateSeuratObject(counts = data, assay = assay)
#     return(object)
# }
# environment(Read10X_generic) <- environment(Seurat::Read10X)


# Load10X_Spatial_v2 <- function(
#     data.dir, #filename = "filtered_feature_bc_matrix.h5", 
#     assay = "Spatial", slice = "slice1", filter.matrix = TRUE, 
#     to.upper = FALSE, image = NULL, img_tag = 'hires', ...
# ) {
#     ## (1) Check inputs
#     if (length(x = data.dir) > 1) {
#         warning("'Load10X_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
#         data.dir <- data.dir[1]
#     }

#     object <- Read10X_generic(data.dir, assay, to.upper=FALSE, ...)
    
#     ## (3) Read image 
#     if (is.null(x = image)) {
#         filename_img_full <- list.files(data.dir, pattern=as.character(glue('{img_tag}.*png')), recursive=TRUE, full.names=TRUE)        
#         if (length(filename_img_full) == 0) stop(glue('Could not find image PNG file with pattern {img_tag}'))
#         dirname_img <- paste(head(strsplit(filename_img_full, '/')[[1]], -1), collapse='/')
#         filename_img <- tail(strsplit(filename_img_full, '/')[[1]], 1)
#         image <- Read10X_Image_v2(image.dir = dirname_img, image.name = filename_img,  filter.matrix = filter.matrix)
#     } else {
#         if (!inherits(x = image, what = "VisiumV1")) 
#             stop("Image must be an object of class 'VisiumV1'.")
#     }
#     image <- image[Cells(x = object)]
#     DefaultAssay(object = image) <- assay
#     object[[slice]] <- image
#     return(object)    
# }

# environment(Load10X_Spatial_v2) <- environment(Seurat::Load10X_Spatial)


## TODO: handle case when feature names are different 
#' @export
Read10X_generic <- function(data.dirs, assay, to.upper=FALSE, libnames=NULL, strip.suffix=TRUE, add.prefix=FALSE, ...) {
    ## Map over multiple libraries 
    
    if (is.null(libnames)) {
        if (length(data.dirs) == 1) {
            libnames <- 'slice1'
        } else {
            libnames <- purrr::map_chr(strsplit(normalizePath(data.dirs), '/'), tail, 1)    
        }
    }
    libnames <- make.names(libnames)
    
    if (length(data.dirs) > 1) {
        add.prefix <- TRUE
    }
    
    data_list <- purrr::map2(data.dirs, libnames, function(data.dir, libname) {
        ## Search for putative filenames
#         filename_h5 <- normalizePath(list.files(data.dir, pattern='filtered.*h5$', recursive=TRUE, full=TRUE))
#         filename_filtered_mtx <- normalizePath(list.files(data.dir, pattern='filtered.*mtx', recursive=TRUE, full=TRUE))
#         filename_mtx <- normalizePath(list.files(data.dir, pattern='mtx', recursive=TRUE, full=TRUE))
        filename_h5 <- unique(normalizePath(list.files(data.dir, pattern = "h5$", recursive = TRUE, full = TRUE)))
        filename_h5 <- grep('filtered', filename_h5, value = TRUE)
        filename_mtx <- unique(normalizePath(list.files(data.dir, pattern = "mtx", recursive = TRUE, full = TRUE)))
        filename_filtered_mtx <- grep('filtered', filename_mtx, value = TRUE)
        filename_mtx <- setdiff(filename_mtx, filename_filtered_mtx)

        if (length(filename_h5) > 0) {
            if (length(filename_h5) > 1) stop(glue('Multiple h5 files found in {data.dir}'))
            data <- Read10X_h5(filename = filename_h5, ...)
        } else if (length(filename_filtered_mtx)) {
            if (length(filename_filtered_mtx) > 1) stop(glue('Multiple filtered mtx files found in {data.dir}'))
            filtered_mtx_dir <- paste(head(strsplit(filename_filtered_mtx, '/')[[1]], -1), collapse='/')
            data <- Read10X(data.dir = filtered_mtx_dir, ...)
        } else if (length(filename_mtx)) {
            if (length(filename_mtx) > 1) stop(glue('Multiple mtx files found in {data.dir}'))
            mtx_dir <- paste(head(strsplit(filename_mtx, '/')[[1]], -1), collapse='/')
            data <- Read10X(data.dir = mtx_dir, ...)
        } else {
            stop(glue('Could not find any h5 or mtx files in {data.dir}'))
        }
        if (to.upper) {
            rownames(x = data) <- toupper(x = rownames(x = data))
        }
        if (strip.suffix) {
            colnames(data) <- gsub('-\\d$', '', colnames(data))    
        }        
        if (add.prefix) {
            colnames(data) <- paste0(libname, '_', colnames(data))        
        }
        return(data)
    })
    
    ## If genes are different among libraries, deal with it somehow
    ## For now, just take common genes
    ## Eventually, include option to pad with zeros 
    if (length(data_list) == 1) {
        object <- CreateSeuratObject(counts = data_list[[1]], assay = assay)        
    } else {
        if (!Reduce(identical, purrr::map(data_list, rownames))) {
            genes_use <- Reduce(intersect, purrr::map(data_list, rownames))
            data_list <- purrr::map(data_list, function(X) X[genes_use, ])
        }    
        object <- CreateSeuratObject(counts = Reduce(Matrix::cbind2, data_list), assay = assay)
    }
    
    return(object)
}
environment(Read10X_generic) <- environment(Seurat::Read10X)

                                    
## NOTE: automatically looks like images in data.dir(s)
##       will include option to pass image(s)
## NOTE: removed "slice" parameter: will this hurt us? 
#' @export
Load10X_Spatial_v2 <- function(
    data.dirs, #filename = "filtered_feature_bc_matrix.h5", 
    assay = "Spatial", filter.matrix = TRUE, strip.suffix = TRUE, add.prefix=FALSE, 
    to.upper = FALSE, img_tag = 'lowres', libnames=NULL, ...
) {
    if (is.null(libnames)) {
        if (length(data.dirs) == 1) {
            libnames <- 'slice1'
        } else {
            libnames <- purrr::map_chr(strsplit(normalizePath(data.dirs), '/'), tail, 1)    
        }
    }

    if (length(data.dirs) > 1) {
        add.prefix <- TRUE
    }
    
#     ## (1) Check inputs
#     if (length(x = data.dir) > 1) {
#         warning("'Load10X_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
#         data.dir <- data.dir[1]
#     }
    object <- Read10X_generic(data.dirs, assay, to.upper=FALSE, libnames=libnames, strip.suffix=strip.suffix, add.prefix=add.prefix, ...)

    
    ## (3) Read images
    images <- purrr::map2(data.dirs, libnames, function(data.dir, libname) {
        ## Search for directory with spatial information
        filename_img_full <- list.files(data.dir, pattern=as.character(glue('{img_tag}.*png')), recursive=TRUE, full.names=TRUE)        
        if (length(filename_img_full) == 0) stop(glue('Could not find image PNG file with pattern {img_tag}'))
        dirname_img <- paste(head(strsplit(filename_img_full, '/')[[1]], -1), collapse='/')
        filename_img <- tail(strsplit(filename_img_full, '/')[[1]], 1)
        
        ## Read spatial information and attach to object
        image <- Read10X_Image_v2(image.dir = dirname_img, image.name = filename_img,  filter.matrix = filter.matrix)
        ## Remove suffix 
        if (strip.suffix) {
            rownames(image@coordinates) <- gsub('-\\d$', '', rownames(image@coordinates))            
        }
        if (add.prefix) {
            rownames(image@coordinates) <- paste0(libname, '_', rownames(image@coordinates))            
        }
        image <- image[intersect(Cells(object), Cells(image))]
        DefaultAssay(image) <- assay
        return(image)
    })
        
    names(images) <- libnames
    for (libname in libnames) {
        object[[libname]] <- images[[libname]]        
    }
    
    ## PREVIOUS CODE FOR ONE IMAGE 
#     if (is.null(x = image)) {
#         filename_img_full <- list.files(data.dir, pattern=as.character(glue('{img_tag}.*png')), recursive=TRUE, full.names=TRUE)        
#         if (length(filename_img_full) == 0) stop(glue('Could not find image PNG file with pattern {img_tag}'))
#         dirname_img <- paste(head(strsplit(filename_img_full, '/')[[1]], -1), collapse='/')
#         filename_img <- tail(strsplit(filename_img_full, '/')[[1]], 1)
#         image <- Read10X_Image_v2(image.dir = dirname_img, image.name = filename_img,  filter.matrix = filter.matrix)
#     } else {
#         if (!inherits(x = image, what = "VisiumV1")) 
#             stop("Image must be an object of class 'VisiumV1'.")
#     }
#     image <- image[Cells(x = object)]
#     DefaultAssay(object = image) <- assay
#     object[[slice]] <- image
    return(object)    
}

environment(Load10X_Spatial_v2) <- environment(Seurat::Load10X_Spatial)

                                    
