#' @title OASIS Prediction
#' @description This function creates the OASIS probability map from a single MRI study with FLAIR, T1, T2, and PD volumes.
#' @param flair flair volume of class \code{\link{nifti}}
#' @param t1 t1 volume of class \code{\link{nifti}}
#' @param t2 t2 volume of class \code{\link{nifti}}
#' @param pd pd volume of class \code{\link{nifti}}
#' @param brain_mask brain mask of class \code{\link{nifti}},
#' if \code{NULL} a brain mask will be created using \code{\link{fslbet}}.
#' Note that provided brain masks should be in the same space as the T1 volume
#' if \code{preproc = TRUE}, as all volumes will be registered to this space
#' @param model an object of class \code{\link{glm}} used to make the OASIS predictions
#' @param return_preproc is a logical value that indicates whether the
#' preprocessed images should be returned, if \code{NULL}
#' then the model from the OASIS paper will be used
#' @param binary logical indicating whether a binary map
#' should be returned by thresholding the probability map
#' @param threshold numeric indicating the threshold value
#' for the probability map, with default of 0.16 for the OASIS paper
#' @param verbose print diagnostic messages
#' @param ... options passed to \code{\link{oasis_train_dataframe}}
#'
#' @return A list of volumes:
#' the OASIS probability map, the preprocessed volumes (if \code{return_preproc = TRUE}),
#' the brain mask for the subject, the voxel selection mask, and a thresholded,
#' binary mask (if \code{binary = TRUE}) .
#' @examples \dontrun{
#' library(neurobase)
#' flair <- readnii('path/to/flair', reorient = FALSE)
#' t2 <- readnii('path/to/t2', reorient = FALSE)
#' t1 <- readnii('path/to/t1', reorient = FALSE)
#' pd <- readnii('path/to/pd', reorient = FALSE)
#' oasis_map <- oasis_predict(flair = flair, t1 = t1, t2 = t2, pd = pd) }
#' @export
#' @importFrom stats predict
#' @importFrom fslr fslsmooth
oasis_predict <- function(flair, ##flair volume of class nifti
                          t1, ##t1 volume of class nifti
                          t2, ##t2 volume of class nifti
                          pd = NULL, ##pd volume of class nifti
                          brain_mask = NULL, ##brain mask of class nifti
                          model = NULL, ##an OASIS model of class glm
                          return_preproc = FALSE, ##option to return the preprocessed data
                          binary = FALSE,
                          threshold = 0.16,
                          verbose = TRUE,
                          ...
) {
  L = oasis_train_dataframe(flair = flair,
                            t1 = t1,
                            t2 = t2,
                            pd = pd,
                            brain_mask = brain_mask,
                            verbose = verbose,
                            return_preproc = TRUE,
                            ...)

  oasis_dataframe = L$oasis_dataframe
  brain_mask = L$oasis_dataframe
  top_voxels = L$top_voxels
  preproc = L$preproc
  rm(list = "L")


  if (verbose) {
    message("Model Prediction")
  }
  ## make the model predictions
  if (is.null(model)) {
    predictions <- predict( oasis::oasis_model,
                            newdata = oasis_dataframe,
                            type = 'response')
  } else {
    predictions <- predict(model,
                           newdata = oasis_dataframe,
                           type = 'response')
  }


  ##put the predictions onto the brain
  predictions_nifti <- niftiarr(brain_mask, 0)
  predictions_nifti[top_voxels == 1] <- predictions

  if (verbose) {
    message("Smoothing Prediction")
  }
  ##smooth the probability map
  prob_map <- fslsmooth(predictions_nifti, sigma = 1.25,
                        mask = brain_mask, retimg = TRUE,
                        smooth_mask = TRUE)

  binary_map = NULL
  if (binary == TRUE) {
    if (verbose) {
      message("Thresholding Smoothed Prediction")
    }
    binary_map <- prob_map
    binary_map[prob_map > threshold] <- 1
    binary_map[prob_map <= threshold] <- 0
  }
  L = list(oasis_map = prob_map,
           binary_map = binary_map,
           unsmoothed_map = predictions_nifti,
           flair = preproc$flair,
           t1 = preproc$t1, t2 = preproc$t2,
           pd = preproc$pd,
           brain_mask = brain_mask,
           voxel_selection = top_voxels)
  if (!return_preproc) {
    L$flair = L$t1 = L$t2 = L$pd = NULL
  }

  if (!binary) {
    L$binary_map = NULL
  }

  return(L)

}


