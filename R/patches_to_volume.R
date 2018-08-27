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
#' @examples
#' user = Sys.getenv("USER")
#' if (user == "johnmuschelli") {
#' reticulate::use_python(paste0(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#' }
#' library(neurobase)
#' fname = system.file("extdata", "MPRAGE.nii.gz", package = "flexconn")
#' t1 = readnii(fname)
#' fname = system.file("extdata", "FLAIR.nii.gz", package = "flexconn")
#' flair = readnii(fname)
#' patchsize = c(3, 3)
#' verbose = TRUE
#' patch = get_patches(t1, flair, mask = NULL, patchsize = patchsize,
#' only_patches = FALSE, normalize = FALSE)
#' out_mask = patch$mask
#' masked_t1 = mask_img(t1, out_mask)
#' t1_remade = get_volume_from_patches(patch$t1_patches,
#' indices = patch$indices,
#' patchsize = patchsize,
#' invert_pad = patch$pad)
#' stopifnot(all(dim(t1) == dim(t1_remade)))
#' t1vals = mask_vals(t1, out_mask)
#' t1_remade_vals = mask_vals(t1_remade, out_mask)
#' testthat::expect_equal(t1vals, t1_remade_vals)
#'
patches_to_volume = function(
  patches,
  indices,
  patchsize,
  verbose = TRUE,
  orig_dim = NULL) {

  num_patches = ncol(indices)
  matsize = get_matsize(num_patches, patchsize = patchsize)
  dpatches = dim(patches)
  stopifnot(all(dpatches == matsize))

  if (is.null(orig_dim)) {
    orig_dim = attributes(indices)$original_dimension
  }

  tindx = t(indices)
  if (any(duplicated(tindx))) {
    warning("Duplicated indices! May be wrong")
    indices = t(unique(tindx))
  }

  if (is.null(orig_dim)) {
    stop("original dimensions cannot be estimated from this data!")
  }

  ndim = length(patchsize)
  index = (patchsize + 1) / 2
  stopifnot( abs(index - as.integer(index)) < 1e-5)

  dsize <- floor(patchsize / 2)
  dsize = c(dsize, rep(0, length = 5 - length(dsize)))

  if (verbose) {
    message("Volume array to create: ",
            paste(orig_dim, collapse = "x"))
  }
  vol <- array(dim = orig_dim)

  for (i in 1:num_patches) {
    x <- indices[1, i]
    y <- indices[2, i]
    z <- indices[3, i]
    if (ndim == 2) {
      vol[x,y,z] = patches[i, index[1], index[2], 1]
    }
    if (ndim == 3) {
      vol[x,y,z] = patches[i, index[1], index[2], index[3], 1]
    }
  }
  return(vol)
}


#' @rdname patches_to_volume
#' @param invert_pad Run \code{\link{invert_pad_image}} on the image
#' before returning.
#' @export
get_volume_from_patches <- function(
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
