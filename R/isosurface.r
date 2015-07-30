#' Create isosurface from 3D images
#'
#' Create isosurface from 3D images
#' @param x image of class antsImage or SimpleITK (_p_itk__simple__Image).
#' @param threshold threshold to create isosurface
#' @param as.int logical: treat image grey values as integers.
#' @param ... not used.
#' @return returns a triangulated mesh of class mesh3d.
#' @importFrom Rvcg vcgIsosurface
#' @rdname isosurface
#' @export
isosurface <- function(x,threshold,as.int=TRUE,...) UseMethod("isosurface")

#' @rdname isosurface
#' @export
isosurface.antsImage <- function(x,threshold,as.int=TRUE,...) {
    segsurf <- vcgIsosurface(ANTsR::as.array(x),threshold=threshold,spacing = antsGetSpacing(x),origin = antsGetOrigin(x),direction = antsGetDirection(x),as.int = as.int)
    return(segsurf)
}

#' @rdname isosurface
#' @export
isosurface._p_itk__simple__Image <- function(x,threshold,as.int=TRUE,...) {
    segsurf <- vcgIsosurface(SimpleITK::as.array(x),threshold=threshold,spacing = x$GetSpacing(),origin = x$GetOrigin(),direction=matrix(x$GetDirection(),3,3),as.int=as.int)
    return(segsurf)
}
