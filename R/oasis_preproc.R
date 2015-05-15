oasis_preproc <- function(flair, t1, t2, pd){
  
  ##inhomogenity correction for all four modalities using fsl bias correct
##  for(i in 1:length(modalities)){
##    modal <- get(modalities[i])
##    assign(modalities[i], fsl_biascorrect(modal, retimg = TRUE)) 
##  }
  
  ##rigidly register to the flair, t2, and pd to the t1 using fsl flirt 
  flair <- flirt(infile = flair, reffile = t1, retimg = TRUE,  dof = 6) 
  t2 <- flirt(infile = t2, reffile = t1, retimg = TRUE,  dof = 6) 
  pd <- flirt(infile = pd, reffile = t1, retimg = TRUE,  dof = 6) 
  
  ##create a brain mask using fsl BET 
  brain_mask <- fslbet(infile = t1, retimg = TRUE)
  brain_mask <-   niftiarr(brain_mask, brain_mask[brain_mask > 0] <- 1)
  
  ##return the preprcoessed images and the brain mask
  setClass(Class="study",
           representation(
             flair ="nifti",
             t1 ="nifti",
             t2 ="nifti",
             pd ="nifti",
             brain_mask = 'nifti')
  )
  
  return(new("study",
             flair = flair,
             t1 = t1,
             t2 = t2,
             pd = pd,
             brain_mask = brain_mask))
  
  return(study)
  
}