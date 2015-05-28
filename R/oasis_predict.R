#' @title OASIS Prediction
#' @description This function creates the OASIS probability maps. 
#' @param flair flair volume of class nifti
#' @param t1 flair volume of class nifti
#' @param t2 flair volume of class nifti
#' @param pd flair volume of class nifti
#' @param brain_mask brain mask of class nifti, if NULL a brain mask will be created using fsl BET
#' @param preproc calls the oasis_preproc function and performs the necessary preprocessing steps for OASIS
#' @param normalize option to perform z-score normalization of the image, should be TRUE unless you train model
#' using an alternative normalization 
#' @param model an object of class glm used to make the OASIS predictions 
#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @import fslr
#' @export 
oasis_predict <- function(flair, ##flair volume of class nifti
                  t1, ##t1 volume of class nifti
                  t2, ##t2 volume of class nifti
                  pd, ##pd volume of class nifti
                  brain_mask = NULL, ##brain mask of class nifti
                  preproc = FALSE, ##option to preprocess the data
                  normalize = TRUE, ##option to normalize 
                  model = NULL ##an OASIS model of class glm
  ) {
  oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  
  ## the image preproceesing 
  if(preproc == TRUE){
    images <- oasis_preproc(oasis_study$flair, oasis_study$t1, oasis_study$t2, oasis_study$pd)
  }
 
  ## create a brain mask  
  if(is.null(brain_mask) == TRUE){
    brain_mask <- fslbet(infile = oasis_study$t1, retimg = TRUE)
    brain_mask <- brain_mask > 0
    brain_mask <- datatyper(brain_mask, trybyte= TRUE)
  } 
  
  
  ##adjust brain mask for OASIS 
  brain_mask <- fslerode(brain_mask, kopts = "-kernel box 5x5x5", retimg = TRUE)
  cutpoint <- quantile(flair[brain_mask == 1], .15)
  brain_mask[flair <= cutpoint] <- 0 
  
  ## the image normalization 
  if(normalize == TRUE){
    oasis_study <- lapply(oasis_study, function (x) zscore_img(x, brain_mask, margin = NULL))  
    }
  
  ## smooth the images using fslsmooth from the fslr package 
  oasis_study <- append(oasis_study, lapply(oasis_study, function(x) fslsmooth(x, sigma = 10, mask = brain_mask)))
  oasis_study <- append(oasis_study, lapply(oasis_study[1:4], function(x) fslsmooth(x, sigma = 20, mask = brain_mask)))
  
  ##create and apply the voxel selection mask 
  top_voxels <- voxel_selection(flair = flair,
                                brain_mask = brain_mask, 
                                cutoff = .85)
  
  ## create a dataframe to make oasis predictions on    
  oasis_study <- lapply(oasis_study, function(x) x[top_voxels == 1])
  oasis_dataframe <- do.call(cbind.data.frame, oasis_study)

  names <- c("FLAIR", "T1", "T2", "PD")
  colnames(oasis_dataframe) <-  c(names, paste0(names, "_10"),  paste0(names, "_20"))
  
  ## make the model predictions 
  if(is.null(model) == TRUE){
  predictions <- predict(oasis_model, newdata = oasis_dataframe, type = 'response')    
  } else { 
  predictions <- predict(model, newdata = oasis_dataframe, type = 'response')    
  }
  
  ##put the predictions onto the brain 
  predictions_nifti <- niftiarr(flair, 0) 
  predictions_nifti[top_voxels == TRUE] <- predictions

  
  ##smooth the probability map 
  sigma.smooth<-diag(3,3)
  k.size<- 5
  prob_map<-niftiarr(predictions_nifti, GaussSmoothArray(predictions_nifti,sigma=sigma.smooth,
                                                   ksize=k.size,mask=brain_mask))
  ##return the probability map 
  return(predictions_nifti)
  
}


