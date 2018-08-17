#' Volume to Patches
#'
#' @param vol 3D array or \code{nifti} image
#' @param patchsize Vector of length 2 (or more)
#' @param indices A matrix of 3 by \code{num_patches} indices for the
#' data to extract
#' @param verbose print diagnostic messages
#'
#' @return A list of image and mask Patches
#' @export
volume_to_patches = function(
  vol, indices, patchsize,
  verbose = TRUE) {

  ndim = length(patchsize)
  dsize <- floor(patchsize / 2)
  dsize = c(dsize, rep(0, length = 5 - length(dsize)))
  num_patches = ncol(indices)
  matsize = get_matsize(num_patches, patchsize = patchsize)

  if (verbose) {
    message("Size matrix to create: ",
            paste(matsize, collapse = "x"))
  }
  patches <- array(0, dim = matsize)
  for (i in 1:num_patches) {
    x <- indices[1, i]
    y <- indices[2, i]
    z <- indices[3, i]
    zindex = (z - dsize[3]):(z + dsize[3])
    vals = vol[
      (x - dsize[1]):(x + dsize[1]),
      (y - dsize[2]):(y + dsize[2]),
      zindex
      ]
    if (ndim == 2) {
      patches[i, , , 1] <- vals
    }
    if (ndim == 3) {
      patches[i, , , , 1] <- vals
    }
  }
  return(patches)
}
