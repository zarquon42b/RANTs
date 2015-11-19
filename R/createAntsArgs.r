#' interface to generate arguments for antsRegistration command
#'
#' interface to generate arguments for antsRegistration command
#'
#' @param reference path to moving image
#' @param target path to fixed image
#' @param setting name of the setting to be attached to output
#' @param percent subsampling for affine transforms
#' @param affine character vector of non-elastic transform to use
#' @param affinereach  character vector assigning the bins for Mattes metric used for affine transforms.
#' @param affineconvergence character: convergence values for affine transforms MxNxO [MxNxO,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]
#' @param elastic character vector of elastic procedures (including parameters). E.g "SyN[0.2,3,0]"
#' @param elasticconvergence character: convergence valzes MxNxO [MxNxO,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]
#' @param elasticmetrics character vector containing the metrics to be used in elastic matching
#' @param metricweights  character vector of same length as \code{elasticmetrics}, assigning weights to each metric. If NULL, all selected metrics will be weighted equally.
#' @param metricreach  character vector of same length as \code{elasticmetrics}, assigning the radius/bins for each metric to consider.
#' @param dims dimensionality of data
#' @param elasticS smoothing sigmas (see command args of antsRegistration)
#' @param elasticF shrink factors (see command args of antsRegistration)
#' @param initTransform run initial coarse alignment or provide a transform file. These features include using the geometric center of the images (=0), the image intensities (=1), or the origin of the images (=2) or a list containing a character containing the path to the transform and a logical indicating if the transform is to be inverted. E.g list("mytransfrom.mat",TRUE)
#' @param itkthreads integer: specify number of threads to be used
#' @param folder character defining the path where to store the transforms
#' @param transimg logical: if TRUE the *diff and *inv will be created.
#' @param verbose logcal: if TRUE output are written to console
#' @param referenceLM k x 3 matrix for landmarks to be used in metrics PSE and JHCT
#' @param targetLM k x 3 matrix for landmarks to be used in metrics PSE and JHCT
#' @param IJK2RAS 4x4 matrix defining the transform of the landmarks into image space.
#' @return returns a list with
#' \item{antsargs}{arguments passed to antsRegistration}
#' \item{outname}{basename of outputname}
#' \item{transforms}{list containing characters naming the transforms}
#' @export
createAntsArgs <- function(reference,target,setting="custom",percent=0.1, affine=c("trans","rigid","similarity","affine"),affinereach=32,affineconverge="[10000x10000x10000,1e-8,20]",elastic="SyN[0.2,3,0]",elasticpercent=0.1,elasticconverge ="[100x10x1,0,5]",elasticmetrics=c("mattes","cc"),metricweights=NULL, metricreach=c(32,4),dims=3,elasticS="3x1x0vox",elasticF="3x2x1",initTransform=NULL,itkthreads=parallel::detectCores(),masks=NULL,folder=NULL,transimg=FALSE,verbose=FALSE,referenceLM=NULL,targetLM=NULL,IJK2RAS=diag(c(-1,-1,1,1))) {
    if (nargs() == 0) {
        antsRegistration()
        return()
    }
    if (!is.null(folder))
        folder <- paste0(folder,"/")
    Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=itkthreads)

    hasLM <- FALSE
    if (!is.null(referenceLM) && !is.null(targetLM)){
        hasLM <- TRUE
        refimage <- antsImageRead(target)
        print(dim(referenceLM))
        randstring <- paste(sample(c(0:9, letters, LETTERS), 5, replace=TRUE),collapse="")
        tlmimage <- paste0(tempdir(),"/",randstring,"_targetLM.nii.gz")
        rlmimage <- paste0(tempdir(),"/",randstring,"_referenceLM.nii.gz")
        antsImageWrite(points2LabelImage(targetLM,RAS2IJK=IJK2RAS,refimage = refimage,neighbours = 100,spacing=antsGetSpacing(refimage)),tlmimage)
        antsImageWrite(points2LabelImage(referenceLM,RAS2IJK=IJK2RAS,refimage=refimage,neighbours = 100,spacing=antsGetSpacing(refimage)),rlmimage)
    }
    
    nm <- paste0(folder,basename(target),"_fixed_",basename(reference),"_moving_",setting[1])
    if (is.null(metricweights))
        metricweights <- rep(1,length(elasticmetrics))
    mweights <- metricweights
    metrics <- elasticmetrics
    mweights <- mweights/sum(mweights)
    
    if (length(metrics) != length(mweights) || length(metrics) != length(metricreach))
        stop("elasticmetrics, metricweights and metricreach must be of same length")
    antsargs <- list(d=dims)
    affinedef <- c("translation","rigid","similarity","affine")
    if (!is.null(initTransform)) {
        if (is.numeric(initTransform))
            antsargs <- append(antsargs,list(r=paste0("[",target,",",reference,",",initTransform,"]")))
        else if (is.list(initTransform)) {
            tmparg <- paste0("[",initTransform[[1]],", ",as.integer(initTransform[[2]]),"]")
            antsargs <- append(antsargs,list(r=tmparg))
        }
    }
    affineprefix <- list(m=paste0("mattes[",target,",",reference,",1,",affinereach,",regular,",percent,"]"))
    affinesuffix <- list(c=affineconverge, s="4x2x1vox", f="6x4x2",l="1")
    for (i in affine) {
        curaff <- match.arg(i,affinedef)
        curargs <- append(affineprefix,list(t=paste0(curaff,"[0.1]")))
        curargs <- append(curargs,affinesuffix)
        antsargs <- append(antsargs,curargs)
                          
    }
    ## setup arguments for elastic matching
    if (length(elastic)) {
        elasticargs <- list()
        for( i in 1 : length(mweights)) {
           
            if (mweights[i] != 0) {
                if (!metrics[i] %in% c("PSE")) {
                    elasticargs <- append(elasticargs,list(m=paste0(metrics[i],"[",target,",",reference,",",mweights[i],",", metricreach[i],",regular,",elasticpercent,"]")))
               } else {
                    if (hasLM)
                        elasticargs <- append(elasticargs,list(m=paste0(metrics[i],"[",tlmimage,",",rlmimage,",",mweights[i],",",elasticpercent,",0,1,",metricreach[i],"]")))
                }
            }
        }
        elasticsuffix <- list(c=elasticconverge,s=elasticS,f=elasticF)
        for (i in elastic) {
            elmid <- list(t=i)
            antsargs <- append(antsargs,unlist(list(elasticargs,elmid,elasticsuffix)))
        }
    }
    if (length(masks))
        antsargs <- append(antsargs,list(x=masks))
    if (transimg)
        output <- list(o=paste0("[",nm,",",nm,"_diff.nii.gz,",nm,"_inv.nii.gz]"))
    else
        output <- list(o=paste0("[",nm,"]"))
    antsargs <- append(antsargs,list(u="1",l="1",z="1",float="1"))
    if (verbose)
         antsargs <- append(antsargs,list(verbose="1"))
    antsargs <- append(antsargs,output)

    transforms <- list()
    if (length(elastic)) {
        if (!length(affine) && is.null(initTransform))
            prenumber <- 0
        else
            prenumber <- 1
        transforms$warpfwd <- paste0(nm,prenumber,"Warp.nii.gz")
        transforms$warpinv <- paste0(nm,prenumber,"InverseWarp.nii.gz")
    }
    if (length(affine) ||  !is.null(initTransform))
        transforms$affine <-  paste0(nm,"0GenericAffine.mat")
    return(list(antsargs=antsargs,outname=basename(nm),transforms=transforms))
    
}




cmdargs <- function(antsargs) {
    x <- antsargs$antsargs
    out <- paste0("-",names(x)," ",x)
    out <- gsub("-float","--float",out)
    cat(out)
    cat("\n")
}
