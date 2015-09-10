#' @title OASIS Image Preprocessing 
#' @description This function does the required preprocessing for OASIS for the FLAIR, T2, 
#' T1, and PD volumes using fsl through fslr.  The preprcoessing steps are (1) inhomogenity
#' correct using fsl fast and  (2) rigid registration using fsl flirt to the T1 space.  The function returns 
#' a list with the inhomogenity corrected volumes. 
#' @param flair flair volume of class nifti
#' @param t1 t1 volume of class nifti
#' @param t2 t2 volume of class nifti
#' @param pd pd volume of class nifti
#' @param verbose a logical value for printing diagnostic output 
#' @param cores numeric indicating the number of cores to be used
#' @import parallel
#' @import fslr
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
                          pd, # pd volume of class nifti
                          verbose = TRUE,
                          cores = 1
){
  
  study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  study = check_nifti(study)
  
  if (verbose){
    message("Rigidly Registering Data to T1 Space\n")
  } 
  
  ##rigidly register to the flair, t2, and pd to the t1 using fsl flirt 
  study[c("flair","t2", "pd")] <- mclapply(study[c("flair","t2", "pd")], function(x)  {
    tfile = tempfile(fileext = ".mat")
    flirt(infile = x, omat = tfile,
          reffile =   study$t1, 
          retimg = TRUE,  dof = 6)
          }, mc.cores = cores)
 
  if (verbose){
    message("Running Brain Extraction Tool\n")
  }
  
  brain_mask <- fslbet(infile = study$t1, retimg = TRUE)
  brain_mask <- brain_mask > 0
  brain_mask <- datatyper(brain_mask, trybyte= TRUE)
  
  study <- mclapply(study, function(x) x*brain_mask) 
  
  ## inhomogenity correction for all four modalities using fsl bias correct
  if (verbose){
    message("Running Inhomogeneity Correction\n")
  }
  study <- mclapply(study, function(x) fsl_biascorrect(x, 
                           retimg = TRUE,
                           verbose = verbose), 
                           mc.cores = cores)

  study[c(length(study) + 1)] <- brain_mask
  
  ##return a list with the preprocessed images and a brain mask 
  return(study)
  
  }