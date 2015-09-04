#' @title Image Dimmension Correction
#' @description This function takes an image and drops dimmensions until the volume is a user specified dimmension. 
#' @param image volume of class nifti
#' @param dim scalar value of desired image dimmension
#' @import fslr
#' @return Reutrns a volume of class nifti of desired dimmension.  
#' @examples \dontrun{
#' flair <- readNIfTI('path/to/flair', reorient = FALSE) 
#' correct_image_dim <- function(flair, dim = 3) } 
#' @export 
correct_image_dim <- function(image, dim = 3){
  out.img<-image
  while(length(dim(out.img)) > dim){
    out.img <- drop_img_dim(out.img)
  }
  return(out.img)
}