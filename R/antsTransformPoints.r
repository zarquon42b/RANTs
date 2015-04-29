#' a handy wrapper for antsApplyTransformsToPoints
#'
#'a handy wrapper for antsApplyTransformsToPoints
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
