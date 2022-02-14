## Main function
#' @export
FindSpatiallyVariableFeatures.SpatialDE2 <- function(obj, features=NULL, sizefactors_colname='nCount_Spatial', images=NULL, verbose=FALSE) {    
    load_py_modules()
    if (verbose) {
        warning('To Use FindSpatiallyVariableFeatures.SpatialDE2, make sure you have SpatialDE2 from PMBio/SpatialDE, not pip.')        
    }
    if (is.null(images)) images <- names(obj@images) ## by default, use all images
    if (length(images) == 0) {
        stop('Could not find any images in this object')
    }
    res <- map(images, function(image) {
        assay <- obj@images[[image]]@assay
        if (verbose) message(glue('Start SVG identification in image {image} and assay {assay}'))
        if (is.null(features)) features <- rownames(obj@assays[[assay]])

        if (length(features) == 0) stop('0 features selected')
        if (verbose) message(glue('Using {length(features)} features'))
        
        if (length(features) < 150) {
            if (verbose) message('Doing SpatialDE2 SVG with reticulate (1 core)')
            svg_sde2_serial(obj, features, sizefactors_colname, image, anndata=anndata, SpatialDE=SpatialDE)
        } else {
            if (verbose) message('Doing SpatialDE2 SVG in parallel')
            svg_sde2_parallel(obj, features, sizefactors_colname, image, anndata=anndata)
        }        
    })
    if (length(res) == 1) {
        return(dplyr::arrange(res[[1]], padj))
    } else {
        res <- map(res, dplyr::arrange, padj)
        names(res) <- images
        return(res)
        
    }    
}

## Check that python packages are installed 
#' @export
load_py_modules <- function(env=parent.frame()) {
    with(env, {
        ## Check SpatialDE version 
        tryCatch({        
            SpatialDE <- reticulate::import('SpatialDE')
        }, error = function(e) {
            stop('Must install SpatialDE in python before running SpatialDE2 SVG.')
        })
        ## Check SpatialDE version 
        tryCatch({        
            anndata <- reticulate::import('anndata')
        }, error = function(e) {
            stop('Must install anndata in python before running SpatialDE2 SVG.')
        })
    })
    return(TRUE)        
}


## In one shot (no parallel)
#' @export
svg_sde2_serial <- function(obj, features=NULL, sizefactors_colname='nCount_Spatial', image=NULL, anndata, SpatialDE) {
    ## Set up Seurat object data access fields 
    if (is.null(image)) image <- names(obj@images)[1] ## by default, use first image 
    assay <- obj@images[[image]]@assay
    if (is.null(features)) features <- rownames(obj@assays[[assay]])
    cells_use <- rownames(obj@images[[image]]@coordinates)  
    adata <- anndata$AnnData(
        X = Matrix::t(obj@assays[[assay]]@counts[features, cells_use]),
        obsm = list('spatial' = as.matrix(dplyr::select(obj@images[[image]]@coordinates, imagerow, imagecol))),
        var=data.frame(gene=features, row.names=features), 
        obs=obj@meta.data[cells_use, ]
    )
    res <- SpatialDE$test(adata, omnibus=TRUE, sizefactors=adata$obs[[sizefactors_colname]])[[1]] 
    return(res)
}

#' @export
write_sde_program <- function(python_path) {
    svg_sde_program_str <- list(
        paste0('#!', python_path),
        'import SpatialDE',
        'import anndata',
        'import sys',
        'adata_fname=sys.argv[1]',
        'sizefactors_colname=sys.argv[2]', 
        'adata = anndata.read(adata_fname)',
        'svg_res, _ = SpatialDE.test(adata, omnibus=True, sizefactors=adata.obs[[sizefactors_colname]])',
        "svg_res.to_csv('{}_out'.format(adata_fname), index=False)"
    ) %>% 
        paste(collapse = '\n')
    
    fname_program <- tempfile(fileext = '.py')
    writeLines(svg_sde_program_str, fname_program)
    system(glue('chmod a+x {fname_program}'))
    return(fname_program)
}

## Function needs to be defined outside of future_map
## to avoid passing all data structures from parent environment 
## See: https://github.com/DavisVaughan/furrr/commit/5ed2b067be842dfd034105913791da676dcc8862
#' @export
do_svg_sde2_py_external <- function(fname_tmp, fname_program, sizefactors_colname) {
    cmd <- as.character(glue('python {fname_program} {fname_tmp} "{sizefactors_colname}"'))
    system(cmd)
    res <- read.csv(paste0(fname_tmp, '_out'))
    return(res)
}

#' @export
svg_sde2_parallel <- function(obj, features=NULL, sizefactors_colname='nCount_Spatial', image=NULL, ncores=NULL, python_path=NULL, anndata) {
    ## Write python program to run downstream
    if (is.null(python_path)) {
        python_path <- system('which python', intern=TRUE)
    }    
    fname_program <- write_sde_program(python_path)
    ## Set up number of cores
    if (is.null(ncores)) {
        ncores <- min(
            ceiling(length(features) / 100), ## no use sending fewer than 100 genes to one job
            round(parallel::detectCores() * .8) 
        )        
    }

    ## Set up Seurat object data access fields 
    if (is.null(image)) image <- names(obj@images)[1] ## by default, use first image 
    assay <- obj@images[[image]]@assay
    if (is.null(features)) features <- rownames(obj@assays[[assay]])

    ## First, initialize the data structures and then pass them to workers 
    ## This avoids passsing the whole Seurat object to each worker 
    cells_use <- rownames(obj@images[[image]]@coordinates)  
    adata_list <- split(features, rep(seq_len(ncores), length.out = length(features))) %>% 
        map(function(.features) {
            env$anndata$AnnData(
                X = Matrix::t(obj@assays[[assay]]@counts[.features, cells_use]),
                obsm = list('spatial' = as.matrix(dplyr::select(obj@images[[image]]@coordinates, imagerow, imagecol))),
                # var=obj@assays$Spatial@meta.features[.features, ], ## cannot be empty, so this doesn't always work
                var=data.frame(gene=.features, row.names=.features), 
                obs=obj@meta.data[cells_use, ] 
            )
        })

    ## NOTE: this part must be done serially, because it uses reticulate
    fnames <- map(adata_list, function(adata) {
        fname_tmp <- tempfile()
        adata$write(fname_tmp)
        return(fname_tmp)
    })

    ## This part is parallel 
    ## Calls python externally, not with reticulate
    # plan(multisession) ## CAUTION: fails with multicore plan 
    # plan(multicore, workers=ncores) ## CAUTION: fails with multicore plan 
    future::plan(future::multisession, workers=ncores) ## CAUTION: fails with multicore plan 
    system.time({
        sde_res <- furrr::future_map(
            fnames, ## list
            do_svg_sde2_py_external, ## function
            fname_program=fname_program, sizefactors_colname=sizefactors_colname ## params
        ) %>% 
            bind_rows()    
    })

    return(sde_res)
} 

