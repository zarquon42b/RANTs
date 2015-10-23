#' convert an affine transform file, generated by ANTsR to a 4x4 transformation matrix
#'
#' convert an affine transform file, generated by ANTsR to a 4x4 transformation matrix
#' @param file affine transform file, generated by ANTsR
#' @param IJK2RAS IJK to RAS transform
#' @return 4x4 matrix
#' @importFrom R.matlab readMat
#' @export
getAffineMat <- function(file,IJK2RAS=diag(c(-1,-1,1,1))) {
    m <- 3
    aff <- readMat(file)
    affine <- t(matrix(aff$AffineTransform.float.3.3[1:9],3,3,byrow = T))
    trans <- aff$AffineTransform.float.3.3[10:12]
    hgamm <- rbind(cbind(affine, 0), 0)
    hgamm[m + 1, m + 1] <- 1
    htrans <- diag(m + 1)
    htrans[1:m, m + 1] <- c(-aff$fixed)
    htrans2 <- diag(m + 1)
    htrans2[1:m, m + 1] <- c(trans)
    hall <- solve(IJK2RAS)%*%htrans2%*%solve(htrans)%*%t(hgamm) %*%htrans%*%IJK2RAS
    return(hall)
}

#' Create an itk transform from landmarks (rigid, similarity and affine)
#'
#' Create an itk transform from landmarks (rigid, similarity and affine) to transform an image
#' @param lmFix fix landmarks
#' @param lmMoving moving landmarks
#' @param type of transform (can be "rigid","similarity" or "affine")
#' @param IJK2RAS transform from point to image space.
#' @details
#' lmMoving are landmarks placed on the image to be transformed to the position of lmFix
#' @return
#' writes a transform to a tempfile and returns the path to that transform
#' @importFrom RcppOctave .CallOctave o_load
#' @importFrom Morpho rotonto applyTransform computeTransform
#' @export
landmarkTransform <- function(lmFix,lmMoving,type=c("rigid","similarity","affine"),IJK2RAS=diag(c(-1,-1,1,1)[seq_len(ncol(lmFix)+1)])) {
    m <- ncol(lmFix)
    file <- (tempfile(pattern = "transform"))
    file <- paste0(file,".mat")
    lmFix <- applyTransform(lmFix,IJK2RAS)
    lmMoving <- applyTransform(lmMoving,IJK2RAS)
    scale <- FALSE
    type = match.arg(type,c("rigid","affine","similarity"))
    #afftrans <- computeTransform(lmFix,lmMoving,type=type)
    afftrans0 <- computeTransform(lmMoving,lmFix,type=type)
    ## afftrans <- afftrans0[1:3,]
    ## AffineTransform_double_3_3 <- c(t(afftrans[1:m,1:m]))
    ## AffineTransform_double_3_3 <- c(AffineTransform_double_3_3,afftrans[,m+1])
    ## print(AffineTransform_double_3_3)
    ## fixed <- rep(0,m)
    ## o_load(list(AffineTransform_double_3_3=matrix(AffineTransform_double_3_3,12,1),fixed=as.matrix(fixed)))
    ## .CallOctave("save", "-v4", file, "AffineTransform_double_3_3", "fixed")
    transform2mat(afftrans0,file)
    return (file)
}

transform2mat <- function(affinemat,file,IJK2RAS=diag(c(1,1,1,1)[seq_len(nrow(affinemat))])) {
    affinemat <- IJK2RAS%*%affinemat%*%IJK2RAS
    m <- nrow(affinemat)-1
    affinemat <- affinemat[1:m,]
    AffineTransform_double_3_3 <- c(t(affinemat[1:m,1:m]))
    AffineTransform_double_3_3 <- c(AffineTransform_double_3_3,affinemat[,m+1])
    fixed <- rep(0,m)
    if (m == 2) {
       o_load(list(AffineTransform_double_2_2=matrix(AffineTransform_double_3_3,6,1),fixed=as.matrix(fixed)))
       .CallOctave("save", "-v4", file, "AffineTransform_double_2_2", "fixed") 
    } else {
        o_load(list(AffineTransform_double_3_3=matrix(AffineTransform_double_3_3,12,1),fixed=as.matrix(fixed)))
        .CallOctave("save", "-v4", file, "AffineTransform_double_3_3", "fixed") 
    }

}



rot2quaternion <- function(gamma) {
    qr <- 0.5*sqrt(sum(c(1,diag(gamma))))
    qi <- (gamma[3,2]-gamma[2,3])/(4*qr)
    qj <- (gamma[1,3]-gamma[3,1])/(4*qr)
    qk <- (gamma[2,1]-gamma[1,2])/(4*qr)
    return(c(qi,qj,qk))
}
