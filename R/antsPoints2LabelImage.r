#' create a label image from landmark coordinates
#'
#' create a label image from landmark coordinates
#' @param pts k x 2 or k x 3 matrix containing landmark coordinates
#' @param neighbours integer: number of closest points in the image grid to label
#' @param spacing spacing of the output image
#' @param margin margin around the landmarks (determined by the extend of the bounding box) to be part of the image
#' @param labelcolors vector of length nrow(pts) with a color for each landmark label
#' @param IJK2RAS transform from point to image space.
#' @return returns an antsImage
#' @importFrom Rvcg vcgKDtree
#' @importFrom Morpho applyTransform
#' @export 
points2LabelImage <- function(pts, neighbours=16,spacing=rep(1,ncol(pts)),margin=0.1,labelcolors=seq_len(nrow(pts)),IJK2RAS=diag(c(-1,-1,1,1)[seq_len(ncol(pts)+1)])) {
    pts <- applyTransform(pts,IJK2RAS)
    if (nrow(pts) != length(labelcolors))
        stop("number of labelcolors and number of landmarks must correspond")
    m <- ncol(pts)
    ranges <- apply(pts,2,range)
    ranges <- apply(ranges,2,extendrange,f=margin)
    grid <- lapply(1:m,function(x) seq(from=ranges[1,x],to=ranges[2,x],by=spacing[x]))
    mygrid <- as.matrix(expand.grid(grid))
    arrdims <- sapply(grid,length)
    myarr <- array(0,dim=arrdims)
    gridlist <- lapply(1:m,function(x) x <- 1L:arrdims[x])
    indices <- as.matrix(expand.grid(gridlist))
    nhs <- vcgKDtree(mygrid,pts,k=neighbours)$index
    dists <- rep(0,nrow(mygrid))
    for(i in 1:nrow(nhs))
        myarr[indices[nhs[i,],]] <- i
    #myarr <- RANTs:::fillArray(indices-1,dists,arrdims)
    origin <- mygrid[1,]
    outimage <- as.antsImage(myarr, spacing = spacing, origin = origin)
    return(outimage)
}
    
