library(fslr)
options(fsl.path="/Applications/fsl")

datadir <- '/Users/elizabethsweeney/Dropbox/Elizabeth_Sweeney_Documents/Packages/OASIS_Test/'

flair <- readNIfTI(paste0(datadir, "2000310_FLAIR_matchToT1.nii"), reorient = FALSE)
t2 <- readNIfTI(paste0(datadir, "2000310t2_N3Corrected_clone_reg.nii"), reorient = FALSE)
t1 <- readNIfTI(paste0(datadir, "2000310_T1_matchToT1.nii"), reorient = FALSE)
pd <- readNIfTI(paste0(datadir, "2000310pd_N3Corrected_clone_reg.nii"), reorient = FALSE)

gold_standard <- readNIfTI(paste0(datadir, "2000310_mod.nii.gz"), reorient = FALSE)




