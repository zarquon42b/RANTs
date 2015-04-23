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
antsTransformPoints <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/",IJK2RAS=diag(c(-1,-1,1,1)),...)UseMethod("antsTransformPoints")

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.matrix <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/",IJK2RAS=diag(c(-1,-1,1,1)),...) {
    ptsdim <- ncol(mat)
    if (ptsdim == 3) {
        pts <- applyTransform(mat,IJK2RAS)
        pts <- cbind(pts,0)
    } else if (ptsdim == 2) {
          pts <- mat
          pts <- cbind(pts,0)
          pts <- applyTransform(pts,IJK2RAS)
          pts <- cbind(pts,0)
      } else
        stop("only 2D and 3D points allowed")
    colnames(pts) <- c("x","y","z","t")
    ptsname <- paste0(tempdir(),"/pts.csv")
    write.csv(pts,file=ptsname,row.names=F)
    output <- paste0(tempdir(),"/ptsDeformed2.csv")
    cmd <- paste0(antsdir,"antsApplyTransformsToPoints -d " ,ptsdim," -i ",ptsname," -o ",output)
    if (!missing(affine))
        cmd <- paste0(cmd," -t ",affine)

    if (!missing(warp))
        cmd <- paste0(cmd," -t ",warp)
    print(cmd)
    system(cmd)
    
    readit <- applyTransform(as.matrix(read.csv(output)[,1:3]),trafo=IJK2RAS,inverse=T)##convert back to RAS space
    readit <- readit[,1:ptsdim]
    
    return(readit)
}

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.mesh3d <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/",IJK2RAS=diag(c(-1,-1,1,1)),...) {
    mesh <- mat
    mat <- vert2points(mat)
    readit <- antsTransformPoints(mat,affine=affine,warp=warp,antsdir=antsdir,...)
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
#' @param synargs parameters for syn registration
#' @param itkthreads integer: specify number of threads to be used
#' @param trans run translation
#' @param rigid run rigid matching
#' @param affine run affine matching
#' @param syn run syn matching
#' @param BSplineDisp use B-Spline deformations
#' @param gauss use gaussian Smoothed deformation fields
#' @param expo use Exponential deformation type
#' @param mattesweight weight to attribute to mattes metric
#' @param ccweight weight to attribute to CC metric
#' @param dims dimensionality of data
#' @param initTransform run initial coarse alignment. These features include using the geometric center of the
#' images (=0), the image intensities (=1), or the origin of the images (=2).
#' @param PSE not supported yet
#' @return returns a list with
#' \item{antsargs}{arguments passed to antsRegistration}
#' \item{outname}{basename of outputname}
#' @export
createAntsArgs <- function(reference,target,setting="custom",percent=0.1,synpercent=0.1,searchradius=4,mattesbins=32,its="10000x10000x10000",synargs ="[100x10x1,0,5]",cc=FALSE,itkthreads=parallel::detectCores(),trans=T,rigid=T,similarity=T,affine=T,BSplineDisp=F,gauss=F,expo=F,syn=T,mattesweight=0.5,ccweight=0.5,dims=3,initTransform=NULL,PSE=NULL) {
    Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=itkthreads)
    nm=paste0(target,"_fixed_",reference,"_moving_",setting[1])
     antsargs <- list(d=dims)
    if (!is.null(initTransform))
        antsargs <- append(antsargs,list(r=paste0("[",target,",",reference,",",initTransform,"]")))
    if (trans)
        antsargs <- append(antsargs,list(m=paste0("mattes[",target,",",reference,",1,",mattesbins,",regular,",percent,"]"), t="translation[0.1]", c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="6x4x2",l="1"))
                                 ## translation step
    if (similarity)
        antsargs <- append(antsargs,list(
                     m=paste0("mattes[",target,",",reference,",1,",mattesbins,",regular,",percent,"]"), t="similarity[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="3x2x1"))
    ##rigid alignment
    if (rigid)
        antsargs <- append(antsargs,list(
                     m=paste0("mattes[",target,",",reference,",1,",mattesbins,",regular,",percent,"]"), t="rigid[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="3x2x1"))
    ##rigid alignment
    if (affine)
           antsargs <- append(antsargs,list(                 
               m=paste0("mattes[",target,",",reference,",1,",mattesbins,",regular,",percent,"]"), t="affine[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox",  f="3x2x1"))
    
    if (gauss)
        antsargs <- append(antsargs,list(         
                     m=paste0("mattes[",target,",",reference,",0.5,", mattesbins,",regular,",synpercent,"]"),m=paste0("cc[",target,",",reference,",0.5,", searchradius, "]"),t="GaussianDisplacementField[0.2,3,0]", c=synargs,  s="1x0.5x0vox", f="4x2x1"))

    if (syn)
        antsargs <- append(antsargs,list(         
            m=paste0("mattes[",target,",",reference,",",mattesweight,",",mattesbins,",regular,",synpercent,"]"),m=paste0("cc[",target,",",reference,",",ccweight ,",", searchradius, "]"),t="SyN[0.2,3,0]", c=synargs,  s="1x0.5x0vox", f="4x2x1"))

    if (expo)
        antsargs <- append(antsargs,list(         
            m=paste0("mattes[",target,",",reference,",0.5,32]"),m=paste0("cc[",target,",",reference,",0.5,4]"),t="Exponential[0.2,3,0.5,10]", c=synargs,  s="1x0.5x0vox", f="4x2x1"))
    if (BSplineDisp)
        antsargs <- append(antsargs,list(         
            m=paste0("mattes[",target,",",reference,",0.5,32]"),m=paste0("cc[",target,",",reference,",0.5,4]"),t="BSplineDisplacementField[0.2,3,0.5,3]", c=synargs,  s="1x0.5x0vox", f="4x2x1"))
    
    antsargs <- append(antsargs,list(u="1",l="1",z="1",float="0",
                                     o=paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]")))
    
        
     return(list(antsargs=antsargs,outname=basename(nm)))
}




