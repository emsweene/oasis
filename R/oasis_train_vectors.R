#' @title OASIS Training Data Frame
#' @description This function creates the training vectors from a single MRI study that has FLAIR, T1, T2, and PD volumes 
#' as well as binary masks of lesions. The function can create a brian mask for the data (or the user can supply a brain mask), 
#' can preprocess the data, and the user may supply already normalized data if they wish to use an alternative normalization method.  
#' @param flair flair volume of class nifti
#' @param t1 t1 volume of class nifti
#' @param t2 t2 volume of class nifti
#' @param pd pd volume of class nifti
#' @param brain_mask brain mask of class nifti, if NULL a brain mask will be created using fsl BET through fslr
#' @param preproc calls the oasis_preproc function and performs the necessary preprocessing steps for OASIS using FSL through fslr
#' @param normalize option to perform z-score normalization of the image, should be TRUE unless you are training model
#' using data alternative normalization 
#' @param model an object of class glm used to make the OASIS predictions 
#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @import fslr
#' @return oasis_dataframe dataframe for use with the oasis_training function 
#' @export 
oasis_train_vectors <- function(flair, ##flair volume of class nifti
                          t1, ##t1 volume of class nifti
                          t2, ##t2 volume of class nifti
                          pd, ##pd volume of class nifti
                          gold_standard, ##gold standard mask of class nifti
                          brain_mask = NULL, ##brain mask of class nifti
                          preproc = FALSE, ##option to preprocess the data
                          normalize = TRUE, ##option to normalize 
                          slices = NULL, #slice vector
                          orientation = "axial" #slice direction
                          ) 
  { 
  
  if(preproc == TRUE){
    ## the image preproceesing 
    images <- oasis_preproc(flair, t1, t2, pd)
    oasis_study <- list(flair = images[[1]], t1 = t1, t2 = images[[2]], pd = images[[3]])
    brain_mask <- images[[4]]
  } else{ 
    if(is.null(brain_mask) == TRUE){
      ## create a brain mask  
      brain_mask <- fslbet(infile = t1, retimg = TRUE)
      brain_mask <- brain_mask > 0
      brain_mask <- datatyper(brain_mask, trybyte= TRUE)
    } 
    else{
      ## no preprocessing  
      oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
    }
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
  
  oasis_study[[length(oasis_study) + 1]] <-  gold_standard

  if(is.null(slices) == TRUE){
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
  
  return(oasis_dataframe)
  
}


