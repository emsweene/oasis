#' @title Combine OASIS training vectors 
#' @description This function combines multiple oasis training dataframes produced by the oasis_trian_vectors function and produces a larger dataframe 
#' that can be used in the oasis_training function.  
#' @param list_of_train_vectors a list of outputs from the oasis_trian_vectors 
#' @return This function reutrns a dataframe for use with the oasis_training function. 
#' @export 
combine_train_vectors <- function(list_of_train_vectors){
  if(is.list(list_of_train_vectors == TRUE){
    train_vectors_multi <- do.call(rbind, lapply(list_of_train_vectors, function(x) x[[1]]))
  }else{
    train_vectors_multi <- do.call(rbind, list_of_train_vectors)  
  }
  return(train_vectors_multi)
}