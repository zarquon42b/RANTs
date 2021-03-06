#' Create an itk transform from landmarks (rigid, similarity and affine)
#'
#' Create an itk transform from landmarks (rigid, similarity and affine) to transform an image
#' @param lmFix fix landmarks
#' @param lmMoving moving landmarks
#' @param type of transform (can be "rigid","similarity" or "affine")
#' @param IJK2RAS transform from point to image space.
#' @param file filename to save the transform to.
#' @param ... additional parameters (such as landmark weights) to be passed to \code{\link{rotonto}}, in case type != affine.
#' @details
#' lmMoving are landmarks placed on the image to be transformed to the position of lmFix
#' @return
#' writes a transform to a tempfile and returns the path to that transform
#' @importFrom Morpho rotonto applyTransform computeTransform getTrafo4x4
#' @export
landmarkTransform <- function(lmFix,lmMoving,type=c("rigid","similarity","affine"),IJK2RAS=diag(c(-1,-1,1,1)[seq_len(ncol(lmFix)+1)]),file=NULL,...) {
    m <- ncol(lmFix)
    if (is.null(file)) {
        file <- (tempfile(pattern = "transform"))
        file <- paste0(file,".mat")
    }
    lmFix <- applyTransform(lmFix,IJK2RAS)
    lmMoving <- applyTransform(lmMoving,IJK2RAS)
    scale <- FALSE
    type = match.arg(type[1],c("rigid","affine","similarity"))
    #afftrans <- computeTransform(lmFix,lmMoving,type=type)
    if (type == "affine") {
        afftrans0 <- computeTransform(lmMoving,lmFix,type=type)
    } else {
        scale <- ifelse(type == "rigid",FALSE,TRUE)
        if (!exists("reflection"))
            reflection <- FALSE
        rot <- rotonto(lmMoving,lmFix,scale=scale,reflection=reflection,...)
        afftrans0 <- getTrafo4x4(rot)
    }
    transform2mat(afftrans0,file,diag(rep(1,m+1)))
    return (file)
}

#' write an affine transform matrix to ITK mat format
#'
#' write an affine transform matrix to ITK mat format
#' @param affinemat a 4x4 (3D case) or 3x3 (2D case) homogenous transform matrix
#' @param file filename
#' @param IJK2RAS transform from point to image space.
#' @return
#' returns the character of the filename
#' @export
transform2mat <- function(affinemat,file,IJK2RAS=diag(c(-1,-1,1,1)[seq_len(nrow(affinemat))])) {
   transout <- transform2antsTrafo(affinemat=affinemat,IJK2RAS = IJK2RAS)
   writeAntsrTransform(transout,filename=file)
}

#' convert an affine transform from matrix to antsrTransform
#'
#' convert an affine transform from matrix to antsrTransform (implementation of itk::Transform
#' 
#' @param affinemat a 4x4 (3D case) or 3x3 (2D case) homogenous transform matrix
#' @param file filename
#' @param IJK2RAS transform from point to image space.
#' @return
#' returns object of class antsrTransform
#' @export
transform2antsTrafo <- function(affinemat,IJK2RAS=diag(c(-1,-1,1,1)[seq_len(nrow(affinemat))])) {
    affinemat <- IJK2RAS%*%affinemat%*%IJK2RAS
    m <- nrow(affinemat)-1
    affinemat <- affinemat[1:m,]
    AffineTransform_double_3_3 <- c(t(affinemat[1:m,1:m]))
    AffineTransform_double_3_3 <- c(AffineTransform_double_3_3,affinemat[1:m,m+1])
    fixed <- rep(0,m)
    trafo <- createAntsrTransform(dimension=m,precision="double")
    if (m == 2) {
        setAntsrTransformParameters(trafo,matrix(AffineTransform_double_3_3,6,1)) 
       
    } else {
        setAntsrTransformParameters(trafo,matrix(AffineTransform_double_3_3,12,1)) 
    }
    return(trafo)

}

