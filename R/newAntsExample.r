#' This is simply the script newAntsExample.sh as an R function
#'
#' This is simply the script newAntsExample.sh as an R function
#' @param reference path to moving image
#' @param target path to fixed image
#' @param setting presets options ar e fast and production
#' @param itkthreads integer: specify number of threads to be used
#' @export
newAntsExample <- function(reference, target, setting=c("fast","production"),itkthreads=4) {
    Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=itkthreads)

    if (setting[1] == "production") {
        percent <- 0.3
        its <- "10000x111110x11110"
        syn <- "[100x100x50,-0.01,5]"
    } else if (setting[1] == "fast") {
        percent <- 0.1
        its <- "10000x0x0"
        syn <- "[100x0x0,0,5]"
    }
    nm=paste0(target,"_fixed_",reference,"_moving_",setting[1])
    antsargs <- list(d=3,
                     m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="translation[0.1]", c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="6x4x2",l="1",## translation step
                     m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="rigid[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="3x2x1",##rigid alignment
                     m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="affine[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox",  f="3x2x1",
                     m=paste0("mattes[",target,",",reference,",0.5,32]"),m=paste0("cc[",target,",",reference,",0.5,4]"),t="SyN[0.2,3,0]", c=syn,  s="1x0.5x0vox", f="4x2x1", u="1",l="1",z="1",float="0",
                     o=paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]"))

    registration <- antsRegistration(antsargs)
    return(paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]"))
}



#' calls antsApplyTransformsToPoints
#'
#' calls antsApplyTransformsToPoints
#' @param mat a matrix or mesh3d
#' @param affine path to affine transformation
#' @param warp path to warp
#' @param antsdir path to directory containing antsApplyTransformsToPoints
#' @param IJK2RAS 4x4 matrix containing the transform between image and pointspace. For 2D pointclouds this must be the transformation for cbind(mat,0) - the z-coordinate supplemented with zeros.
#' @importFrom Morpho vert2points applyTransform
#' @importFrom Rvcg vcgUpdateNormals
#' @rdname antsTransformPoints
#' @export
antsTransformPoints <- function(mat,transformlist,IJK2RAS=diag(c(-1,-1,1,1)),...)UseMethod("antsTransformPoints")

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.matrix <- function(mat,transformlist, whichtoinvert=NA,IJK2RAS=diag(c(-1,-1,1,1)),...) {
    ptsdim <- ncol(mat)
    if (ptsdim == 3) {
        pts <- applyTransform(mat,IJK2RAS)
    } else if (ptsdim == 2) {
        pts <- mat
        pts <- cbind(pts,0)
        pts <- applyTransform(pts,IJK2RAS)
    } else
        stop("only 2D and 3D points allowed")
    #readit <- as.matrix(antsApplyTransformsToPoints(ptsdim,pts,transformlist,whichtoinvert))
    print(dim(pts))
    readit <- applyTransform(as.matrix(antsApplyTransformsToPoints(ptsdim,pts,transformlist,whichtoinvert)),trafo=IJK2RAS,inverse=T)
    #readit <- applyTransform(as.matrix(read.csv(output)[,1:3]),trafo=IJK2RAS,inverse=T)##convert back to RAS space
    readit <- readit[,1:ptsdim]
    
    return(readit)
}

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.mesh3d <- function(mat,transformlist, whichtoinvert=NA,IJK2RAS=diag(c(-1,-1,1,1)),...) {
    mesh <- mat
    mat <- vert2points(mat)
    readit <- antsTransformPoints(mat,transformlist,whichtoinvert,IJK2RAS,...)
    mesh$vb[1:3,] <- t(readit)
    mesh <- vcgUpdateNormals(mesh)
return(mesh)
}


#' interface to generate arguments for antsRegistration command
#'
#' interface to generate arguments for antsRegistration command
#'
#' @param reference path to moving image
#' @param target path to fixed image
#' @param setting name of the setting to be attached to output
#' @param percent parameter for ants
#' @param searchradius searchradius for CC metric
#' @param mattesbins binsize for Mattes MI metric
#' @param its parameters for non-syn registration
#' @param synargs character: convergence valzes MxNxO [MxNxO,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]
#' @param affine character vector of non-elastic transform to use
#' @param elastic character vector of elastic procedures (including parameters). E.g "SyN[0.2,3,0]"
#' @param mattesweight weight to attribute to mattes metric (used only in elastic registration)
#' @param ccweight weight to attribute to CC metric (used only in elastic registration)
#' @param dims dimensionality of data
#' @param elasticS smoothing sigmas (see command args of antsRegistration)
#' @param elasticF shrink factors (see command args of antsRegistration)
#' @param initTransform run initial coarse alignment. These features include using the geometric center of the
#' images (=0), the image intensities (=1), or the origin of the images (=2).
#' @param itkthreads integer: specify number of threads to be used
#' @return returns a list with
#' \item{antsargs}{arguments passed to antsRegistration}
#' \item{outname}{basename of outputname}
#' @export
createAntsArgs <- function(reference,target,setting="custom",percent=0.1,synpercent=0.1,searchradius=4,mattesbins=32,its="10000x10000x10000",synargs ="[100x10x1,0,5]", affine=c("trans","rigid","similarity","affine"),elastic="SyN[0.2,3,0]",mattesweight=0.5,ccweight=0.5,dims=3,initTransform=NULL,elasticS="3x1x0vox",elasticF="3x2x1",itkthreads=parallel::detectCores()) {
    if (nargs() == 0) {
        antsRegistration()
        return()
    }
    Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=itkthreads)
    nm=paste0(target,"_fixed_",reference,"_moving_",setting[1])
    weights <- c(mattesweight,ccweight)
    weights <- weights/sum(weights)
    mattesweight <- weights[1]
    ccweight <- weights[2]
    antsargs <- list(d=dims)
    affinedef <- c("translation","rigid","similarity","affine")
    if (!is.null(initTransform))
        antsargs <- append(antsargs,list(r=paste0("[",target,",",reference,",",initTransform,"]")))
    
    affineprefix <- list(m=paste0("mattes[",target,",",reference,",1,",mattesbins,",regular,",percent,"]"))
    affinesuffix <- list(c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="6x4x2",l="1")
    for (i in affine) {
        curaff <- match.arg(i,affinedef)
        curargs <- append(affineprefix,list(t=paste0(curaff,"[0.1]")))
        curargs <- append(curargs,affinesuffix)
        antsargs <- append(antsargs,curargs)
                          
    }
    if (length(elastic)) {
        elasticargs <- list()
        if (mattesweight > 0)
            elasticargs <- append(elasticargs,list(m=paste0("mattes[",target,",",reference,",0.5,", mattesbins,",regular,",synpercent,"]")))
        if (ccweight > 0 )
            elasticargs <- append(elasticargs,list(m=paste0("cc[",target,",",reference,",0.5,", searchradius, "]")))
        elasticsuffix <- list(c=synargs,s=elasticS,f=elasticF)
        for (i in elastic) {
            elmid <- list(t=i)
            antsargs <- append(antsargs,unlist(list(elasticargs,elmid,elasticsuffix)))
        }
    }    
    antsargs <- append(antsargs,list(u="1",l="1",z="1",float="1",
                                     o=paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]")))
    
    return(list(antsargs=antsargs,outname=basename(nm)))
}




cmdargs <- function(antsargs) {
    x <- antsargs$antsargs
    out <- paste0("-",names(x)," ",x)
    out <- gsub("-float","--float",out)
    cat(out)
    cat("\n")
}
