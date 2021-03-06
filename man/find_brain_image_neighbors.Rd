% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFFscreening.R
\name{find_brain_image_neighbors}
\alias{find_brain_image_neighbors}
\title{Find the neighbors for each voxel in images}
\usage{
find_brain_image_neighbors(img, mask = NULL, radius = 1)
}
\arguments{
\item{img}{an nifti object.}

\item{mask}{an nifti object.
\code{mask > 0} specifies which voxels are on the mask.
Default value is NULL indicating all voxels are considered.}

\item{radius}{the size (voxel) of neighborhood.}
}
\value{
the neighbor indices of each voxel in each row.
}
\description{
Find the neighbors for each voxel in images
}
\examples{
maskfile <- file.path(system.file("nifti", package="BFFscreening"),"AAL_MNI_2mm.nii")
mask <- oro.nifti::readNIfTI(maskfile)
imgfile <- file.path(system.file("nifti", package="BFFscreening"),"VBM_example.nii.gz")
img <- oro.nifti::readNIfTI(imgfile)
nb <- find_brain_image_neighbors(img, mask,radius=1)
}
\author{
Jian Kang <jiankang@umich.edu>
}
