#' @title Voxel Selection Procedure
#' @description This function creates a binary mask for the voxel selection 
#' procedure for OASIS.
#' @param flair flair volume of class nifti
#' @param brain_mask brain mask of class nifti
#' @param cutoff the percentile cutoff for the thresholding, 
#' passed to \code{\link{quantile}}
#' @return Returns the voxel slection mask as an object of class nifti.  
#' @examples \dontrun{
#' flair <- readNIfTI('path/to/flair', reorient = FALSE) 
#' brain_mask <- readNIfTI('path/to/brain_mask', reorient = FALSE) 
#' voxel_selection_mask <- oasis_preproc(flair, brain_mask, cutoff = .85) }
#' @import oro.nifti
#' @export
voxel_selection <- function(flair, ##the flair volume
                            brain_mask, ## a brain mask for the flair 
                            cutoff ## the percentile cutoff for the thresholding 
                            ){
  stopifnot(length(cutoff) == 1)
  ##find the value to threshold the flair volume at 
  cutpoint <- quantile(flair[brain_mask == 1], cutoff)
  outmask <-  flair > cutpoint
  outmask <- as.nifti(outmask * brain_mask) 
  outmask <- datatype(outmask, trybyte = TRUE)
  ##return the binary mask of the flair values above the cutpoint 
  return(outmask) 
}