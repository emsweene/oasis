% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oasis_predict.R
\name{oasis_predict}
\alias{oasis_predict}
\title{OASIS Prediction}
\usage{
oasis_predict(flair, t1, t2, pd = NULL, brain_mask = NULL, model = NULL,
  return_preproc = FALSE, binary = FALSE, threshold = 0.16,
  verbose = TRUE, ...)
}
\arguments{
\item{flair}{flair volume of class \code{\link{nifti}}}

\item{t1}{t1 volume of class \code{\link{nifti}}}

\item{t2}{t2 volume of class \code{\link{nifti}}}

\item{pd}{pd volume of class \code{\link{nifti}}}

\item{brain_mask}{brain mask of class \code{\link{nifti}},
if \code{NULL} a brain mask will be created using \code{\link{fslbet}}.
Note that provided brain masks should be in the same space as the T1 volume
if \code{preproc = TRUE}, as all volumes will be registered to this space}

\item{model}{an object of class \code{\link{glm}} used to make the 
OASIS predictions}

\item{return_preproc}{is a logical value that indicates whether the
preprocessed images should be returned, if \code{NULL}
then the model from the OASIS paper will be used}

\item{binary}{logical indicating whether a binary map
should be returned by thresholding the probability map}

\item{threshold}{numeric indicating the threshold value
for the probability map, with default of 0.16 for the OASIS paper}

\item{verbose}{print diagnostic messages}

\item{...}{options passed to \code{\link{oasis_train_dataframe}}}
}
\value{
A list of volumes:
the OASIS probability map, the preprocessed volumes (if 
\code{return_preproc = TRUE}),
the brain mask for the subject, the voxel selection mask, and a thresholded,
binary mask (if \code{binary = TRUE}) .
}
\description{
This function creates the OASIS probability map from a single 
MRI study with FLAIR, T1, T2, and PD volumes.
}
\examples{
library(ROCR)
p = predict( oasis::oasis_model,
    newdata = example_oasis_df,
    type = 'response')
nopd_p = predict( oasis::nopd_oasis_model,
    newdata = example_oasis_df,
    type = 'response')    
y =  example_oasis_df$GOLD_Lesions
pred = ROCR::prediction(p, y)
perf = ROCR::performance(pred, "tpr", "fpr")
plot(perf)


library(neurobase)
dl_file = function(url) {
   tfile = tempfile(fileext = ".nii.gz")
   req <- httr::GET(url,
   httr::write_disk(path = tfile))
   httr::stop_for_status(req)
   tfile
}
in_ci <- function() {
 nzchar(Sys.getenv("CI"))
}
on_cran = function() {
 identical(Sys.getenv("NOT_CRAN"), "false")
} 
if (in_ci() || on_cran()) {
  if (fslr::have.fsl() && require(httr)) {
    mods = c("FLAIR", "T1W", "T2W", "consensus_gt", "brainmask")
    base_url = file.path(
      "https://raw.githubusercontent.com/muschellij2/open_ms_data", 
      "master/cross_sectional/coregistered/patient01/")
    files = paste0(base_url, mods, ".nii.gz")
    files = sapply(files, dl_file)
    names(files) = mods

    flair <- readnii(files["FLAIR"])
    t1 <- readnii(files["T1W"])
    t2 <- readnii(files["T2W"])
    brain_mask <- readnii(files["brainmask"])
    gold_standard = readnii(files["consensus_gt"])
    oasis_preprocessed_data <- oasis_predict(flair, t1, t2, 
      brain_mask = brain_mask, preproc = TRUE)
  } 
}
}
