#' Create a (signed) distance image from surface mesh
#'
#' Create a (signed) distance image from surface mesh
#'
#' @param x surface mesh
#' @param spacing spacing of grid
#' @param margin percentage to add to grid
#' @param sign logical: if TRUE signed distances are computed
#' @param k k- closest faces to consider during closest point search
#' @param IJK2RAS transformation between physical and image space
#' @return
#' returns a distance image of class antsImage
#' @useDynLib RANTs
#' @importFrom Rvcg vcgClostKD
#' @importFrom Morpho applyTransform
#' @export
DistanceImage <- function(x,spacing=rep(1,3),margin=0.05,sign=TRUE,k=1,IJK2RAS=diag(c(-1,-1,1,1))) {
    x <- applyTransform(x,IJK2RAS)
    pts <- vert2points(x)
    ranges <- apply(pts,2,range)
    ranges <- apply(ranges,2,extendrange,f=margin)
    grid <- lapply(1:3,function(x) seq(from=ranges[1,x],to=ranges[2,x],by=spacing[x]))
    mygrid <- as.matrix(expand.grid(grid[[1]],grid[[2]],grid[[3]]))
    arrdims <- sapply(grid,length)
    myarr <- array(NA,dim=arrdims)
    indices <- as.matrix(expand.grid(1L:arrdims[1],1L:arrdims[2],1L:arrdims[3]))
    dist <- vcgClostKD(mygrid,x,k=k,sign=sign)$quality
    myarr <- fillArray(indices-1,dist,arrdims)
    origin <- mygrid[1,]
    outimage <- as.antsImage(myarr, spacing = spacing, origin = origin)
    return(outimage)
    



}
#' Create a (optionally) binary image from surface mesh
#'
#'Create a (optionally) binary image from surface mesh
#'
#' @param x surface mesh
#' @param spacing spacing of grid
#' @param margin percentage to add to grid
#' @param sign logical: if TRUE signed distances are computed
#' @param k k- closest faces to consider during closest point search
#' @param IJK2RAS transformation between physical and image space
#' @param inval inside value
#' @param outval outside value (only considered if binary = TRUE) otherwise outval = thresh.
#' @param thresh maximal distance to consider
#' @param binary logical: if TRUE the output will be a binary image
#' @return
#' returns a distance image of class antsImage
#' @useDynLib RANTs
#' @importFrom Rvcg vcgClostKD
#' @importFrom Morpho applyTransform
#' @export
Mesh2Image <- function(x,spacing=rep(1,3),margin=0.05,sign=FALSE,k=1,IJK2RAS=diag(c(-1,-1,1,1)),inval=255,outval=0,thresh=min(spacing),binary=FALSE) {
    x <- applyTransform(x,IJK2RAS)
    pts <- vert2points(x)
    ranges <- apply(pts,2,range)
    ranges <- apply(ranges,2,extendrange,f=margin)
    grid <- lapply(1:3,function(x) seq(from=ranges[1,x],to=ranges[2,x],by=spacing[x]))
    mygrid <- as.matrix(expand.grid(grid[[1]],grid[[2]],grid[[3]]))
    arrdims <- sapply(grid,length)
    myarr <- array(NA,dim=arrdims)
    indices <- as.matrix(expand.grid(1L:arrdims[1],1L:arrdims[2],1L:arrdims[3]))
    dist <- vcgClostKD(mygrid,x,k=k,sign=sign)$quality
    good <- which(dist <= thresh)
    if (!binary)
        outval <- thresh
    newdist <- dist*0+outval
    if (binary)
        newdist[good] <- inval
    else
        newdist[good] <- dist[good]
    myarr <- RANTs:::fillArray(indices-1,newdist,arrdims)
    origin <- mygrid[1,]
    outimage <- as.antsImage(myarr, spacing = spacing, origin = origin)
    return(outimage)
}

fillArray <- function(Ind, yr, dims) {
    out <- .Call("fillArr",Ind, yr,dims)
    return(out)
}

