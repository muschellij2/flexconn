#' Blur Mask and Get Indices for Patches
#'
#' @param mask binary 3D array or \code{nifti} image
#'
#' @return The blurred mask
#' @export
#'
#' @examples
#' dims = rep(10, 3)
#' n = prod(dims)
#' mask = rbinom(n, size = 1, prob = 0.85)
#' mask = array(mask, dim = dims)
#' res = blur_mask(mask)
#'
blur_mask  = function(mask) {
  sc <- reticulate::import("scipy")
  np <- reticulate::import("numpy")

  py_mask = as.array(mask)
  py_mask = reticulate::r_to_py(py_mask)
  py_mask = py_mask$astype(np$float32)


  blurmask <-
    sc$ndimage$filters$gaussian_filter(py_mask, sigma = c(1, 1, 1))
  blurmask[blurmask < 0.0001] <- 0
  blurmask <- blurmask * 100
  return(blurmask)
}


