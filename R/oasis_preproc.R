#' @title OASIS Image Preprocessing 
#' @description This function does the required preprocessing for OASIS for the FLAIR, T2, 
#' T1, and PD volumes using fsl through fslr.  The preprcoessing steps are (1) inhomogenity
#' correct using fsl fast and  (2) rigid registration using fsl flirt to the T1 space.  The function returns 
#' a list with the inhomogenity corrected volumes. 
#' @param flair flair volume of class nifti
#' @param t1 t1 volume of class nifti
#' @param t2 t2 volume of class nifti
#' @param pd pd volume of class nifti
#' @return Returns a list of objects of class nifti, namely the preprocessed FLAIR, T1, T2, and PD registered to the space of the T1 volume.  
#' @examples \dontrun{
#' library(oro.nifti)
#' flair <- readNIfTI('path/to/flair', reorient = FALSE) 
#' t1 <- readNIfTI('path/to/t1', reorient = FALSE) 
#' t2 <- readNIfTI('path/to/t2', reorient = FALSE) 
#' pd <- readNIfTI('path/to/pd', reorient = FALSE)
#' oasis_preprocessed_data <- oasis_preproc(flair, t1, t2, pd) }
#' @export 
oasis_preproc <- function(flair, #flair volume of class nifti
                          t1, # t1 volume of class nifti
                          t2, # t2 volume of class nifti
                          pd # pd volume of class nifti
                          ){
  
  study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  
  ## inhomogenity correction for all four modalities using fsl bias correct
  study_inhomo <- lapply(study, function(x) fsl_biascorrect(x, retimg = TRUE))

  ##rigidly register to the flair, t2, and pd to the t1 using fsl flirt 
  study_inhomo_reg <- lapply(study[c(1,3:4)], function(x)  flirt(infile = x, omat = 'temp',
                                                             reffile =   study$t1, 
                                                             retimg = TRUE,  dof = 6))
  
  ##return a list with the preprocessed images and a brain mask 
  return(list(flair = study_inhomo_reg[[1]], t1 = study_inhomo[[2]],  t2 = study_inhomo_reg[[2]], pd = study_inhomo_reg[[3]]))

}