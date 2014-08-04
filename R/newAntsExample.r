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
#' @importFrom Morpho vert2points
#' @rdname antsTransformPoints
#' @export
antsTransformPoints <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/")UseMethod("antsTransformPoints")

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.matrix <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/") {
    pts <- mat%*% diag(c(-1,-1,1))
    write.csv(pts,file="pts.csv",row.names=F)
    cmd <- paste0(antsdir,"antsApplyTransformsToPoints -d 3 -i pts.csv -o ptsDeformed2.csv -t ",affine," -t ",warp)
    print(cmd)
    system(cmd)
    readit <- as.matrix(read.csv("ptsDeformed2.csv")[,1:3])%*%diag(c(-1,-1,1))##convert back to RAS space
    return(readit)
}

#' @rdname antsTransformPoints
#' @export
antsTransformPoints.mesh3d <- function(mat,affine,warp,antsdir="~/GIT/DEV/ANTSbuild/bin/") {
    mesh <- mat
    mat <- vert2points(mat)
    pts <- mat%*% diag(c(-1,-1,1))
    ptsname <- paste0(tempdir(),"/pts.csv")
    write.csv(pts,file=ptsname,row.names=F)
    output <- paste0(tempdir(),"/ptsDeformed2.csv")
    
    cmd <- paste0(antsdir,"antsApplyTransformsToPoints -d 3 -i ",ptsname," -o ",output," -t ",affine," -t ",warp)
    print(cmd)
    system(cmd)
    readit <- as.matrix(read.csv(output)[,1:3])%*%diag(c(-1,-1,1))##convert back to RAS space
    mesh$vb[1:3,] <- t(readit)
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
#' @param its parameters for non-syn registration
#' @param synargs parameters for syn registration
#' @param itkthreads integer: specify number of threads to be used
#' @param trans run translation
#' @param rigid run rigid matching
#' @param affine run affine matching
#' @param syn run syn matching
#' @export
createAntsArgs <- function(reference,target,setting="custom",percent=0.1,its="10000x10000x10000",synargs ="[100x10x1,0,5]",cc=FALSE,itkthreads=4,trans=T,rigid=T,affine=T,syn=T) {
    Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=itkthreads)
    nm=paste0(target,"_fixed_",reference,"_moving_",setting[1])
     antsargs <- list(d=3)
    if (trans)
              antsargs <- append(antsargs,list(m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="translation[0.1]", c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="6x4x2",l="1"))
                                 ## translation step
    if (rigid)
        antsargs <- append(antsargs,list(
                     m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="rigid[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox", f="3x2x1"))
    ##rigid alignment
    if (affine)
           antsargs <- append(antsargs,list(                 
               m=paste0("mattes[",target,",",reference,",1,32,regular,",percent,"]"), t="affine[0.1]",c=paste0("[",its,",1e-8,20]"), s="4x2x1vox",  f="3x2x1"))
    if (syn)
        antsargs <- append(antsargs,list(         
                     m=paste0("mattes[",target,",",reference,",0.5,32]"),m=paste0("cc[",target,",",reference,",0.5,4]"),t="SyN[0.2,3,0]", c=synargs,  s="1x0.5x0vox", f="4x2x1"))
    antsargs <- append(antsargs,list(u="1",l="1",z="1",float="0",
                                     o=paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]")))
     return(antsargs)
}



