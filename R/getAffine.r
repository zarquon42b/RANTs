#' convert an affine transform file, generated by ANTsR to a 4x4 transformation matrix
#'
#' convert an affine transform file, generated by ANTsR to a 4x4 transformation matrix
#' @param file affine transform file, generated by ANTsR
#' @param IJK2RAS IJK to RAS transform
#' @return 4x4 matrix
#' @export
getAffineMat <- function(file,IJK2RAS=diag(c(-1,-1,1,1))) {
    if (!require(R.matlab))
        stop("please install package R.matlab")
    
    aff <- readMat(file)
    affine <- t(matrix(aff$AffineTransform.float.3.3[1:9],3,3,byrow = T))
    trans <- aff$AffineTransform.float.3.3[10:12]
    hgamm <- rbind(cbind(affine, 0), 0)
    hgamm[m + 1, m + 1] <- 1
    htrans <- diag(m + 1)
    htrans[1:m, m + 1] <- c(-aff$fixed)
    htrans2 <- diag(m + 1)
    htrans2[1:m, m + 1] <- c(trans)
    hall <- IJK2RAS%*%htrans2%*%solve(htrans)%*%t(hgamm) %*%htrans%*%IJK2RAS
    return(hall)
}
