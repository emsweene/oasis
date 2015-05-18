#' @title Voxel Selection Procedure
#' @description This function does the voxel selection 
#' procedure for OASIS.
#' @param flair flair volume of class nifti
#' @param brain_mask a brain mask for the flair 
#' @param cutoff the percentile cutoff for the thresholding, 
#' passed to \code{\link{quantile}}
#' @return Object of class nifti 
#' @examples \dontrun{
#' ## when you write the dontrun, it doesn't run but will
#' ## include it in the Rd file (documentation)
#' }
#' @export
voxel_selection <- function(flair, ##the flair volume
                            brain_mask, ## a brain mask for the flair 
                            cutoff ## the percentile cutoff for the thresholding 
                            ){
  stopifnot(length(cutoff) == 1)
  ##find the value to threshold the flair volume at 
  cutpoint <- quantile(flair[brain_mask == 1], cutoff)
  outmask = flair > cutpoint
  ##return the binary mask mask of the flair values above the cutpoint 
  return(outmask) 
}