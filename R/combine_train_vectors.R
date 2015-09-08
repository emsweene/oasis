#' @title Combine OASIS training vectors 
#' @description This function combines multiple oasis training dataframes produced by the oasis_trian_vectors function and produces a larger dataframe 
#' that can be used in the oasis_training function.  The multiple dataframes should be passed into the function in list form. 
#' @param list_of_train_vectors list of outputs from the oasis_trian_vectors function 
#' @param preproc logical indicating whether the output of the oasis_train_vectors contains the preprocessed images as well as the dataframe. 
#' @return This function reutrns a dataframe for use with the oasis_training function. 
#' @examples \dontrun{
#' study_1_dataframe <- oasis_train_vectors(flair = FLAIR_1 ,t1 = T1_1 , t2 = T2_1, pd = PD_1, gold_standard = GOLD_STANDARD_1)
#' study_2_dataframe <- oasis_train_vectors(flair = FLAIR_2 ,t1 = T1_2 , t2 = T2_2, pd = PD_2, gold_standard = GOLD_STANDARD_2)
#' multi_study_dataframe <- combine_train_vectors(list(study_1_dataframe, study_2_dataframe)) }
#' @export 
combine_train_vectors <- function(list_of_train_vectors, preproc = FALSE){
  if(preproc == TRUE && length(ist_of_train_vectors) > 1){
    train_vectors_multi <- do.call(rbind, lapply(list_of_train_vectors, function(x) x[[1]]))
  }
  if(preproc == TRUE && length(ist_of_train_vectors) == 1){
    train_vectors_multi <- list_of_train_vectors[[1]]
  }
  else{
    train_vectors_multi <- do.call(rbind, list_of_train_vectors)  
  }
  return(train_vectors_multi)
}