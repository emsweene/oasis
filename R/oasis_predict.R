#' @title OASIS Prediction
#' @description This function creates the OASIS probability map from a single MRI study with FLAIR, T1, T2, and PD volumes. 
#' @param flair flair volume of class nifti
#' @param t1 t1 volume of class nifti
#' @param t2 t2 volume of class nifti
#' @param pd pd volume of class nifti
#' @param brain_mask brain mask of class nifti, if NULL a brain mask will be created using fsl BET.  Note that provided brain masks should be in the same space as the T1 volume if preproc = TRUE, as all volumes will be registered to this space 
#' @param preproc is a logical value that determines whether to call the oasis_preproc function and performs the necessary preprocessing steps for OASIS
#' @param normalize is a logical value that determines whether to perform z-score normalization of the image over the brain mask, should be TRUE unless you train model
#' using an alternative normalization 
#' @param model an object of class glm used to make the OASIS predictions 
#' @param return_preproc is a logical value that indicates whether the preprcoessed images should be returned, if NULL then the model from the OASIS paper will be used 
#' @import fslr
#' @return If return_preproc = FALSE the function reutrns a volume of class nifti containing the OASIS probability for each voxel. 
#' Otherwise, the function returns a list of volumes: the OASIS probability map, the FLAIR volume, the T1 volume, the T2 volume,
#' the PD volume, the brain mask for the subject, and the voxel selection mask. 
#' @examples \dontrun{
#' flair <- readNIfTI('path/to/flair', reorient = FALSE) 
#' t2 <- readNIfTI('path/to/t2', reorient = FALSE) 
#' t1 <- readNIfTI('path/to/t1', reorient = FALSE) 
#' pd <- readNIfTI('path/to/pd', reorient = FALSE) 
#' oasis_map <- oasis_predict(flair = flair, t1 = t1, t2 = t2, pd = pd) }
#' @export 
oasis_predict <- function(flair, ##flair volume of class nifti
                  t1, ##t1 volume of class nifti
                  t2, ##t2 volume of class nifti
                  pd, ##pd volume of class nifti
                  brain_mask = NULL, ##brain mask of class nifti
                  preproc = FALSE, ##option to preprocess the data
                  normalize = TRUE, ##option to normalize 
                  model = NULL, ##an OASIS model of class glm
                  return_preproc = FALSE ##option to return the preprocessed data
  ) {
  flair = check_nifti(flair)
  t1 = check_nifti(t1)
  t2 = check_nifti(t2)
  pd = check_nifti(pd)
  
  ##correct image dimmension
  flair <- correct_image_dim(flair)
  t1 <- correct_image_dim(t1)
  t2 <- correct_image_dim(t2)
  pd <- correct_image_dim(pd)
  
  
  ##image preproceesing 
  if(preproc == TRUE){
    ## the image preproceesing 
    oasis_study <- oasis_preproc(flair = flair, t1 = t1, t2 = t2, pd = pd)
  }else{
    ## no preprocessing  
    oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  }
  
  if (is.null(brain_mask)){
      ## create a brain mask if not supplied
      brain_mask <- fslbet(infile = t1, retimg = TRUE)
      brain_mask <- brain_mask > 0
      brain_mask <- datatyper(brain_mask, trybyte = TRUE)
    } 

  
  ##adjust brain mask for OASIS 
  brain_mask <- correct_image_dim(brain_mask)
  brain_mask <- fslerode(brain_mask, kopts = "-kernel box 5x5x5", retimg = TRUE)
  cutpoint <- quantile(oasis_study$flair[brain_mask == 1], .15)
  brain_mask[oasis_study$flair <= cutpoint] <- 0 
  
  ## the image normalization 
  if(normalize == TRUE){
    oasis_study <- lapply(oasis_study, function (x) zscore_img(x, brain_mask, margin = NULL))  
  }
  
  ## smooth the images using fslsmooth from the fslr package 
  oasis_study <- append(oasis_study, lapply(oasis_study, function(x) fslsmooth(x, sigma = 10, mask = brain_mask, smooth_mask = TRUE)))
  oasis_study <- append(oasis_study, lapply(oasis_study[1:4], function(x) fslsmooth(x, sigma = 20, mask = brain_mask, smooth_mask = TRUE)))
  
  ##create and apply the voxel selection mask 
  top_voxels <- voxel_selection(flair = oasis_study$flair,
                                brain_mask = brain_mask, 
                                cutoff = .85)
  
  ## create a dataframe to make oasis predictions on    
  oasis_study <- lapply(oasis_study, function(x) x[top_voxels == 1])
  oasis_dataframe <- do.call(cbind.data.frame, oasis_study)

  names <- c("FLAIR", "T1", "T2", "PD")
  colnames(oasis_dataframe) <-  c(names, paste0(names, "_10"),  paste0(names, "_20"))
  
  ## make the model predictions 
  if(is.null(model)){
  predictions <- predict(oasis_model, newdata = oasis_dataframe, type = 'response')    
  } else { 
  predictions <- predict(model, newdata = oasis_dataframe, type = 'response')    
  }
  
  ##put the predictions onto the brain 
  predictions_nifti <- niftiarr(flair, 0) 
  predictions_nifti[top_voxels == 1] <- predictions

  
  ##smooth the probability map 
#   sigma.smooth<-diag(3,3)
#   k.size<- 5
#   prob_map<-niftiarr(predictions_nifti, GaussSmoothArray(predictions_nifti,sigma=sigma.smooth,
#                                                    ksize=k.size,mask=brain_mask))
  # sigma = 1.25
  ##return the data
  
  if(return_preproc == TRUE){
  return(list(oasis_map = predictions_nifti, flair = flair, t1 = t1, t2 = t2,
              pd = pd, brain_mask = brain_mask, voxel_selection = top_voxels))
  } else{
  return(predictions_nifti)
  }
}


