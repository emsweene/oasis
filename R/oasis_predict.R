#' @title OASIS Prediction
#' @description
#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @param flair flair volume of class nifti
#' @param t1 flair volume of class nifti
#' @param t2 flair volume of class nifti
#' @param pd flair volume of class nifti
#' @param brain_mask brain mask of class nifti, if NULL a brain mask will be created using fsl BET
#' @param preproc calls the oasis_preproc function and performs the necessary preprocessing steps for OASIS
#' @param normalize option to perform z-score normalization of the image, should be TRUE unless you train model
#' using an alternative normalization 
#' @param model an object of class glm used to make the OASIS predictions 
#' @export 
oasis_predict <- (flair, ##flair volume of class nifti
                  t1, ##t1 volume of class nifti
                  t2, ##t2 volume of class nifti
                  pd, ##pd volume of class nifti
                  brain_mask = NULL, ##brain mask of class nifti
                  preproc = FALSE, ##option to preprocess the data
                  normalize = TRUE, ##option to normalize 
                  model = default_model ##OASIS model of class glm
  ) {
  oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  
  ## the image preproceesing 
  if(preproc == TRUE){
    images <- oasis_preproc(oasis_study$flair, oasis_study$t1, oasis_study$t2, oasis_study$pd)
  }
 
  ## create a brain mask  
  if(brain_mask == NULL){
    brain_mask <- fslbet(infile = oasis_study$t1, retimg = TRUE)
    brain_mask <- brain_mask > 0
    oasis_study$brain_mask <- datatyper(brain_mask, trybyte= TRUE)
  } else {
    oasis_study$brain_mask <- brain_mask
  }

  ## the image normalization 
  if(normalize == TRUE){
    oasis_study[1:4] <- lapply(oasis_study[1:4], function (x) zscore_img(x,  oasis_study$brain_mask, margin = NULL))  
    }

  ## smooth the images using fslsmooth from the fslr package 
    oasis_smooth_10 <- lapply(oasis_study[1:4], function(x) fslsmooth(x, sigma = 10, mask = brain_mask))
    oasis_smooth_20 <- lapply(oasis_study[1:4], function(x) fslsmooth(x, sigma = 20, mask = brain_mask))
  
  do.call(rbind.data.frame,oasis_study[1:4], oasis_smooth_10,  oasis_smooth_20)
  ## create a dataframe to store the oasis covariates 
  oasis.data <-  data.frame(matrix(NA, length(c(flair)), 12, 
                                   dimnames=list(c(), c(modalities,  
                                                        paste0(modalities, "_10"), 
                                                        paste0(modalities, "_20")))))

  
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