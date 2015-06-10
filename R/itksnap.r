#' call itksnap to view an antsImage
#'
#' call itksnap to view an antsImage
#' @param image antsImage
#' @param overlay antsImage as overlay
#' @export
itksnap <- function(image, overlay=NA) {
      tempname1 <- paste(tempfile(), '.nii.gz', sep='')
        antsImageWrite(image, tempname1)
        command <- paste('itksnap -g', tempname1)
        if (!is.na(overlay)) {
                tempname2 <- paste(tempfile(), '.nii.gz', sep='')
                    antsImageWrite(overlay, tempname2)
                    command <- paste(command, '-o', tempname2)
            }
        system(command,wait=FALSE)
  }
