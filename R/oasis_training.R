#' @title OASIS Training
#' @description This function trains the OASIS model from a data frame produced by the function oasis_train_vectors 
#' @param oasis_dataframe dataframe from the oasis_train_vectors function
#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @import fslr
#' @export 
oasis_training <- function(oasis_dataframe ##dataframe of OASIS data
){
  
  ##fit the oasis model 
  oasis_model <- glm(formula = GoldStandard ~ FLAIR_10 *FLAIR  +
                       FLAIR_20*FLAIR + PD_10 *PD  + PD_20 *PD 
                     + T2_10 *T2 +  T2_20 *T2 + T1_10 *T1 
                     + T1_20 *T1 , data = oasis_dataframe, family=binomial)
  
  ##clean up the oasis model 
  oasis_model$y = c()
  oasis_model$model = c()
  oasis_model$residuals = c()
  oasis_model$fitted.values = c()
  oasis_model$effects = c()
  oasis_model$qr$qr = c()  
  oasis_model$linear.predictors = c()
  oasis_model$weights = c()
  oasis_model$prior.weights = c()
  oasis_model$data = c()
  oasis_model$family$variance = c()
  oasis_model$family$dev.resids = c()
  oasis_model$family$aic = c()
  oasis_model$family$validmu = c()
  oasis_model$family$simulate = c()
  attr(oasis_model$terms,".Environment") = c()
  attr(oasis_model$formula,".Environment") = c()
  
  return(oasis_model)
}