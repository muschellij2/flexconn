#' Blur Mask and Get Indices for Patches
#'
#' @param mask binary 3D array or \code{nifti} image
#' @param verbose print diagnostic messages
#'
#' @return A list of indices and the blurred mask
#' @export
#'
#' @examples
#' dims = rep(10, 3)
#' n = prod(dims)
#' mask = rbinom(n, size = 1, prob = 0.85)
#' mask = array(mask, dim = dims)
#' res = blur_mask(mask)
#'
blur_mask  = function(mask, verbose = TRUE) {
  sc <- reticulate::import("scipy")
  np <- reticulate::import("numpy")

  py_mask = as.array(mask)
  py_mask = reticulate::r_to_py(py_mask)
  py_mask = py_mask$astype(np$float32)

  indx <- which(mask != 0, arr.ind = TRUE)
  indx = t(indx)
  num_patches <- ncol(indx)
  if (verbose) {
    message("Number of patches used: ", num_patches, "\n")
  }
  randindx <- sample(1:num_patches, num_patches, replace = FALSE)
  newindx <- matrix(0, nrow = 3, ncol = num_patches)
  for (i in 1:num_patches) {
    for (j in 1:3) {
      newindx[j, i] <- indx[j, randindx[i]]
    }
  }

  blurmask <-
    sc$ndimage$filters$gaussian_filter(py_mask, sigma = c(1, 1, 1))
  blurmask[blurmask < 0.0001] <- 0
  blurmask <- blurmask * 100
  return(list(indices = newindx, blurred_mask = blurmask))
}


