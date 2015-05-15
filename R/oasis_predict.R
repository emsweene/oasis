#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @export
oasis_predict <- (flair, ##the flair volume
                  t1, ##the t1 volume  
                  t2, ##the t2 volume
                  pd, ##the pd volume
                  brain_mask = NULL, ## the brain mask
                  preproc = FALSE,
                  normalize = TRUE,
                  model = default_model 
  )
{
  setClass(Class="study",
           representation(
             flair ="nifti",
             t1 ="nifti",
             t2 ="nifti",
             pd ="nifti",
             brain_mask = 'nifti')
  )
  oasis_study <- (new("study",
             flair = flair,
             t1 = t1,
             t2 = t2,
             pd = pd,
             brain_mask = brain_mask))
  
  modalities <- c("flair", "t1", "t2", "pd")
  
  ## the image preproceesing 
  if(preproc == TRUE){
    images <- oasis_preproc(flair, t1, t2, pd)
  }

  
  ## the image normalization 
  if(normalize == TRUE){
    for(i in 1:length(modalities)){
      assign(modalities[i], zscore_img(get(modalities[i]), 
                                       mask = brain_mask, margin = NULL))  
    }
  }
  
  ## create a dataframe to store the oasis covariates 
  oasis.data <-  data.frame(matrix(NA, length(c(flair)), 12, 
                                   dimnames=list(c(), c(modalities,  
                                                        paste0(modalities, "_10"), 
                                                        paste0(modalities, "_20")))))
  ## create oasis covariates and add these to the dataframe
  for(i in 1:length(modalities)){
    modal <- get(modalities[i])
    oasis.data[,i] <- c(modal)
    oasis.data[,c(i + 4)] <- c(fslsmooth(modal, sigma = 10, mask = brain_mask))
    oasis.data[,c(i + 8)] <- c(fslsmooth(modal, sigma = 20, mask = brain_mask))
  }
  
  ## make the model predictions 
  predictions <- 
  predictions <- niftiarr(flair,
  
  ##create and apply the voxel selection mask 
  top_voxels <- voxel_selection(flair = flair,
                                brain_mask = mask, 
                                cutoff = .85)
  predictions <- niftiarr(predictions, predictions[top_voxels != 1] <- 0)
  
  ##smooth the probability map 
  sigma.smooth<-diag(3,3)
  k.size<- 5
  prob_map<-niftiarr(predictions, GaussSmoothArray(predictions,sigma=sigma.smooth,
                                                   ksize=k.size,mask=brain_mask))
  
  return(predictions)
  
}