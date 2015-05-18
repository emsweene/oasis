#' @import fslr
#' 
oasis_preproc <- function(flair, t1, t2, pd){
  
  study = list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  
  ## inhomogenity correction for all four modalities using fsl bias correct
  study_inhomo <- lapply(study, function(x) fsl_biascorrect(x, retimg = TRUE))

  ##rigidly register to the flair, t2, and pd to the t1 using fsl flirt 
  study_inhomo_reg <- lapply(study, function(x)  flirt(infile = x, omat = 'temp',
                                                             reffile =   study$t1, 
                                                             retimg = TRUE,  dof = 6))
  
  ##create a brain mask using fsl BET 
  brain_mask <- fslbet(infile = t1, retimg = TRUE, outfile = outfile)
  brain_mask <- brain_mask > 0
  study_inhomo_reg$brain_mask <- datatyper(brain_mask, trybyte= TRUE)
  
  return(study_inhomo_reg)

}