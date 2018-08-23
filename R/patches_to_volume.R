#' Convert Patches to a Volume
#'
#' @param patches A matrix of patches, generally from
#' \code{\link{volume_to_patches}}
#' @param indices A matrix of 3 by \code{num_patches} indices for the
#' data to extract
#' @param verbose print diagnostic messages
#' @param patchsize Vector of length 2 (or more)
#' @param orig_dim Original dimensions of volume
#'
#' @return A 3D volume
#' @export
patches_to_volume = function(
  patches,
  indices,
  patchsize,
  verbose = TRUE,
  orig_dim = NULL) {

  num_patches = ncol(indices)
  if (is.null(orig_dim)) {
    orig_dim = attributes(indices)$original_dimension
  }

  tindx = t(indices)
  if (any(duplicated(tindx))) {
    warning("Duplicated indices! May be wrong")
    indices = t(unique(tindx))
  }
  matsize = get_matsize(num_patches, patchsize = patchsize)


  if (is.null(orig_dim)) {
    stop("original dimensions cannot be estimated from this data!")
  }

  ndim = length(patchsize)
  index = (patchsize + 1) / 2
  stopifnot( abs(index - as.integer(index)) < 1e-5)

  dsize <- floor(patchsize / 2)
  dsize = c(dsize, rep(0, length = 5 - length(dsize)))

  if (verbose) {
    message("Size matrix to create: ",
            paste(matsize, collapse = "x"))
  }
  vol <- array(dim = orig_dim)

  for (i in 1:num_patches) {
    x <- indices[1, i]
    y <- indices[2, i]
    z <- indices[3, i]
    if (ndim == 2) {
      vol[x,y,z] = patches[i, index, index, 1]
    }
    if (ndim == 3) {
      vol[x,y,z] = patches[i, index, index, index, 1]
    }
  }
  return(vol)
}


#' @rdname patches_to_volume
#' @param invert_pad Run \code{\link{invert_pad_image}} on the image
#' before returning.
#' @export
get_volume_from_patch <- function(
  patches,
  indices,
  patchsize,
  verbose = TRUE,
  invert_pad = TRUE,
  orig_dim = NULL) {

  ndim = length(patchsize)
  if (!ndim %in% c(2, 3)) {
    stop("Patch Size is not length 2 or 3!")
  }

  padsize = patchsize_to_padsize(patchsize)

  vol = patches_to_volume(
    patches = patches,
    indices = indices,
    patchsize = patchsize,
    verbose = verbose,
    orig_dim = orig_dim)
  if (invert_pad) {
    vol = invert_pad_image(vol = vol, padsize = padsize)
  }
  return(vol)
}
