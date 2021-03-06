% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blur_mask.R
\name{blur_mask}
\alias{blur_mask}
\title{Blur Mask and Get Indices for Patches}
\usage{
blur_mask(mask)
}
\arguments{
\item{mask}{binary 3D array or \code{nifti} image}
}
\value{
The blurred mask
}
\description{
Blur Mask and Get Indices for Patches
}
\examples{
user = Sys.getenv("USER")
if (user \%in\% c("johnmuschelli", "travis") &
Sys.info()["sysname"] == "Darwin"){
reticulate::use_python(paste0(
"/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
} else {
python = system("which python", intern = TRUE)
print(python)
reticulate::use_python(python)
}
dims = rep(10, 3)
n = prod(dims)
mask = rbinom(n, size = 1, prob = 0.85)
mask = array(mask, dim = dims)
res = blur_mask(mask)

}
