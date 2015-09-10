#' @title OASIS Training Data Frame
#' @description This function creates the training vectors from a single MRI study that has FLAIR, T1, T2, and PD volumes 
#' as well as binary masks of lesions. The function can create a brian mask for the data (or the user can supply a brain mask), 
#' can preprocess the data, and the user may supply already normalized data if they wish to use an alternative normalization method.  
#' @param flair flair volume of class nifti
#' @param t1 t1 volume of class nifti
#' @param t2 t2 volume of class nifti
#' @param pd pd volume of class nifti
#' @param gold_standard gold standard lesion segmentation mask of class nifti
#' @param brain_mask brain mask of class nifti, if NULL a brain mask will be created using fsl BET through fslr
#' @param preproc is a logical value that determines whether to call the oasis_preproc function and performs the necessary preprocessing steps for OASIS
#' @param normalize is a logical value that determines whether to perform z-score normalization of the image over the brain mask, should be TRUE unless you train model
#' using an alternative normalization 
#' @param slices vector of desired slices to train on, if NULL then train over the entire brain mask 
#' @param orientation string value telling which oreintation the training slices are specified in, can take the values of  "axial", "sagittal", or "coronal"
#' @param return_preproc is a logical value that indicates whether the preprcoessed images should be returned 
#' @param cores numeric indicating the number of cores to be used (no more than 4 is useful for this software implementation)
#' @import fslr
#' @import parallel
#' @return If return_preproc = FALSE the function reutrns a dataframe for use with the oasis_training function. 
#' Otherwise, the function returns a list containing: a dataframe for use with the oasis_training function, the FLAIR volume, the T1 volume, the T2 volume,
#' the PD volume, the brain mask for the subject, and the voxel selection mask. 
#' @export 
oasis_train_dataframe <- function(flair, ##flair volume of class nifti
                          t1, ##t1 volume of class nifti
                          t2, ##t2 volume of class nifti
                          pd, ##pd volume of class nifti
                          gold_standard, ##gold standard mask of class nifti
                          brain_mask = NULL, ##brain mask of class nifti
                          preproc = FALSE, ##option to preprocess the data
                          normalize = TRUE, ##option to normalize 
                          slices = NULL, #slice vector
                          orientation = "axial", #slice direction
                          return_preproc = FALSE,
                          cores = 1
                          ) 
  { 
  
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
    preprocess <- oasis_preproc(flair = flair, t1 = t1, t2 = t2, pd = pd, cores = cores)
    oasis_study <- preprocess[c("flair","t1", "t2", "pd")]
    brain_mask <- preprocess[[5]]
  }else{
    ## no preprocessing  
    oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  }
  if (is.null(brain_mask) & !preproc){
    ## create a brain mask if not supplied
    brain_mask <- fslbet(infile = oasis_study$t1, retimg = TRUE)
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
    oasis_study <- mclapply(oasis_study, function (x) zscore_img(x, brain_mask, margin = NULL), mc.cores = cores)
  }
  
  ## smooth the images using fslsmooth from the fslr package 
  oasis_study <- append(oasis_study, mclapply(oasis_study, function(x) fslsmooth(x, sigma = 10, mask = brain_mask), mc.cores = cores))
  oasis_study <- append(oasis_study, mclapply(oasis_study[1:4], function(x) fslsmooth(x, sigma = 20, mask = brain_mask), mc.cores = cores))
  
  ##create and apply the voxel selection mask 
  top_voxels <- voxel_selection(flair = oasis_study$flair,
                                brain_mask = brain_mask, 
                                cutoff = .85)

  
  oasis_study[[c(length(oasis_study) + 1)]] <- gold_standard

  if(is.null(slices)){
    oasis_study <- lapply(oasis_study, function(x) x[top_voxels == 1])
  } else {
    if(orientation == "axial"){
      oasis_study <- lapply(oasis_study, function(x) x[,,slices])
      oasis_study <- lapply(oasis_study, function(x) x[top_voxels[,,slices] == 1])
      
    }    
    if(orientation == "coronal"){
      oasis_study <- lapply(oasis_study, function(x) x[,slices,])
      oasis_study <- lapply(oasis_study, function(x) x[top_voxels[,slices,] == 1])
    } 
    if(orientation == "sagittal"){
      oasis_study <- lapply(oasis_study, function(x) x[slices,,])
      oasis_study <- lapply(oasis_study, function(x) x[top_voxels[slices,,] == 1])
    }    
    
  }
  
  oasis_dataframe <- do.call(cbind.data.frame, oasis_study)
  names <- c("FLAIR", "T1", "T2", "PD")
  colnames(oasis_dataframe) <-  c(names, paste0(names, "_10"),  paste0(names, "_20"),  "GoldStandard")
  
  if(return_preproc == TRUE){
    return(list(oasis_dataframe = oasis_dataframe, flair = flair, t1 = t1, t2 = t2, 
                pd = pd, brain_mask = brain_mask, voxel_selection = top_voxels))
  } else{
    return(oasis_dataframe)
  }
  
}


