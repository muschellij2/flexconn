#' Get Mask Indices
#'
#' @param mask binary 3D array or \code{nifti} image
#' @param seed Seed for random sampling of indices.  If \code{NULL},
#' no sampling is done
#' @param verbose print diagnostic messages
#'
#' @return A 3xV matrix with V number of patches
#' @export
#'
#' @examples
#' dims = rep(10,3)
#' mask = array(rbinom(prod(dims), size = 1, prob = 0.2),
#' dim = dims)
mask_indices = function(mask, seed = NULL, verbose = TRUE) {
  mask = check_nifti(mask, allow.array = TRUE)

  dmask = dim(mask)
  indx <- which(mask != 0, arr.ind = TRUE)
  indx = t(indx)
  attr(indx, "original_dimension") = dmask
  if (is.null(seed)) {
    return(indx)
  }
  set.seed(seed)

  num_patches <- ncol(indx)
  if (verbose) {
    message("Number of patches used: ", num_patches, "\n")
  }
  # random permutation
  randindx <- sample.int(num_patches,
                         size = num_patches, replace = FALSE)
  newindx <- indx[, randindx]

  attr(newindx, "original_dimension") = dmask
  return(newindx)
}
