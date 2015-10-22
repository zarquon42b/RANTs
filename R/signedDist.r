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
#' @param cores integer: amount of cores used in parallel execution (be careful, as the process gets forked and then uses cores times as much RAM).
#' @param useMP logical: if TRUE the parallel computation will be during kdtree search, otherwise the mclapply is used.
#' @param maxind integer: maximum amount of grid points to evaluate at once (for memory reasons)
#' @return
#' returns a distance image of class antsImage
#' @useDynLib RANTs
#' @importFrom Rvcg vcgClostKD
#' @importFrom Morpho applyTransform
#' @importFrom parallel mclapply
#' @export
DistanceImage <- function(x,spacing=rep(1,3),margin=0.05,sign=TRUE,k=1,IJK2RAS=diag(c(-1,-1,1,1)),cores=2,useMP=TRUE,maxind=ifelse(useMP,1e6,1e5)) {
    grid <- getGrid(x=x,spacing=spacing,margin=margin,sign=sign,IJK2RAS=IJK2RAS,cores=cores,maxind=maxind,useMP = useMP)
                                        #dist <- vcgClostKD(mygrid,x,k=k,sign=sign)$quality
                                        #dist <- RvtkStatismo:::vtkPolyDistance(mygrid,x)
    myarr <- fillArray(grid$indices-1,grid$dist,grid$arrdims)
    outimage <- as.antsImage(myarr, spacing = spacing, origin = grid$origin)
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
#' @param cores integer: amount of cores used in parallel execution (be careful, as the process gets forked and then uses cores times as much RAM).
#' @param useMP logical: if TRUE the parallel computation will be during kdtree search, otherwise the mclapply is used.
#' @param maxind integer: maximum amount of grid points to evaluate at once (for memory reasons)
#' @return
#' returns a distance image of class antsImage
#' @useDynLib RANTs
#' @importFrom Rvcg vcgClostKD
#' @importFrom Morpho applyTransform
#' @export
Mesh2Image <- function(x,spacing=rep(1,3),margin=0.05,sign=FALSE,k=1,IJK2RAS=diag(c(-1,-1,1,1)),inval=255,outval=0,thresh=min(spacing),binary=FALSE,cores=2,useMP=TRUE,maxind=ifelse(useMP,1e6,1e5)) {
    grid <- getGrid(x=x,spacing=spacing,margin=margin,sign=sign,IJK2RAS=IJK2RAS,cores=cores,maxind=maxind,useMP = useMP)
    good <- which(grid$dist <= thresh)
    if (!binary)
        outval <- thresh
    newdist <- grid$dist*0+outval
    if (binary)
        newdist[good] <- inval
    else
        newdist[good] <- grid$dist[good]
    myarr <- RANTs:::fillArray(grid$indices-1L,newdist,grid$arrdims)
    
    outimage <- as.antsImage(myarr, spacing = spacing, origin = grid$origin)
    return(outimage)
}

fillArray <- function(Ind, yr, dims) {
    out <- .Call("fillArr",Ind, yr,dims)
    return(out)
}

getGrid <- function(x,spacing=rep(1,3),margin=0.05,sign=FALSE,k=1,IJK2RAS=diag(c(-1,-1,1,1)),cores=2,maxind=1e5,useMP=TRUE) {
    x <- applyTransform(x,IJK2RAS)
    if (useMP) {
        threads <- cores
        cores <- 1
    } else
        threads <- 1
    
    pts <- vert2points(x)
    ranges <- apply(pts,2,range)
    ranges <- apply(ranges,2,extendrange,f=margin)
    grid <- lapply(1:3,function(x) seq(from=ranges[1,x],to=ranges[2,x],by=spacing[x]))
    mygrid <- as.matrix(expand.grid(grid[[1]],grid[[2]],grid[[3]]))
    arrdims <- sapply(grid,length)
    myarr <- array(NA,dim=arrdims)
    indices <- as.matrix(expand.grid(1L:arrdims[1],1L:arrdims[2],1L:arrdims[3]))
    partlist <- list()
    partitions <- floor(nrow(mygrid)/maxind)
    partitions <- ifelse(partitions == 0,1,partitions)
    if (partitions > 1) {
        for (i in 1:partitions) {
            partlist[[i]] <- mygrid[(1:maxind)+(i-1)*maxind,]
        }
    } else {
        partlist[[1]] <- mygrid[1:min(maxind,nrow(mygrid)),]
    }
    if (partitions*maxind < nrow(mygrid))
        partlist[[partitions+1]] <- mygrid[(partitions*maxind+1):nrow(mygrid),]
    
                                        #dist <- lapply(partlist,vcgClostKD,mesh=x,k=k,sign=sign)$quality
    dist <- unlist(parallel::mclapply(partlist,function(y) return(vcgClostKD(y,mesh=x,k=k,sign=sign,threads=threads)$quality),mc.cores = cores))
    return(list(dist=dist,indices=indices,origin=mygrid[1,],arrdims=arrdims))
}
