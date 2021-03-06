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

