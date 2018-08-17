#' Get Patches from T1 and FLAIR image
#'
#' @param t1 3D array or \code{nifti} image
#' @param flair 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param pad Run \code{\link{pad_image}} on the image
#' before getting patches (pads then normalizes if \code{normalize = TRUE})
#' @param normalize Run \code{\link{normalize_image}} on the image
#' before getting patches
#' @param only_patches Only return the patches, not additional information.
#'
#' @note If \code{mask = NULL}, a mask will be created based on
#' voxels greater than the 75th percentile of the FLAIR image.
#'
#' @return A list of T1, FLAIR, and Mask Patches
#' @export
#'
#' @importFrom reticulate import
#' @examples
#' library(neurobase)
#' fname = system.file("extdata", "MPRAGE.nii.gz", package = "flexconn")
#' t1 = readnii(fname)
#' fname = system.file("extdata", "FLAIR.nii.gz", package = "flexconn")
#' flair = readnii(fname)
#' mask = t1 > 0
#' patchsize = c(3, 3)
#' verbose = TRUE
#' patch = get_patches(t1, flair, mask, patchsize = patchsize)
#'
#' rm(patch)
get_patches <- function(
  t1, flair, mask = NULL, patchsize,
  pad = TRUE,
  normalize = TRUE,
  verbose = TRUE,
  only_patches = TRUE) {

  if (is.null(mask)) {
    flair = check_nifti(flair, allow.array = TRUE)
    qflair = quantile(flair[ flair != 0], probs = 0.75)
    mask = flair >= qflair
  }

  t1_patches = get_patch_from_volume(
    t1, mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    normalize = normalize,
    contrast = "T1",
    pad = pad)
  t1_patches = t1_patches$image_patches
  fl_patches = get_patch_from_volume(
    flair, mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    normalize = normalize,
    contrast = "FLAIR",
    pad = pad)

  L = list(t1_patches = t1_patches,
       fl_patches = fl_patches$image_patches,
       mask_patches = fl_patches$mask_patches
       )
  if (!only_patches) {
       L$blurred_mask = fl_patches$blurred_mask
       L$indices = fl_patches$indices
       L$mask = mask
       L$patchsize = patchsize
  }
  return(L)
}

#' Get Patches from 3D Volume
#'
#' @param vol 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image, for the brain usually
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param normalize Run \code{\link{normalize_image}} on the image
#' before getting patches
#' @param contrast What imaging sequence of MRI is this volume, passed
#' to \code{\link{normalize_image}}
#' @param pad Run \code{\link{pad_image}} on the image
#' before getting patches (pads then normalizes if \code{normalize = TRUE})
#' @param ... not used
#'
#' @return A list of image and mask Patches
#' @export
#'
#' @importFrom reticulate import
#' @importFrom neurobase check_nifti
#' @examples
#' library(neurobase)
#' fname = system.file("extdata", "MPRAGE.nii.gz", package = "flexconn")
#' vol = readnii(fname)
#' mask = vol > 0
#' patchsize = c(5, 5)
#' verbose = TRUE
#' patch = get_patch_from_volume(vol, mask, patchsize = patchsize,
#' contrast = "T1")
#'
#' rm(patch)
#'
#' patch = get_patch_from_volume(vol, mask, patchsize = c(1,1,1),
#' contrast = "T1")
#' rm(patch)
#'
get_patch_from_volume <- function(
  ..., patchsize) {

  ndim = length(patchsize)

  func = switch(as.character(ndim),
                "2" = get_2d_patch_from_volume,
                "3" = get_3d_patch_from_volume)
  if (is.null(func)) {
    stop("Patch Size is not length 2 or 3!")
  }
  args = list(...)
  args$patchsize = patchsize
  res = do.call(func, args = args)
  return(res)
}

#' @rdname get_patch_from_volume
#' @export
get_num_patches = function(mask) {
  mask = check_nifti(mask, allow.array = TRUE)
  num_patches = sum(mask != 0)
  return(num_patches)
}

#' @rdname get_patch_from_volume
#' @export
get_2d_patch_from_volume <- function(
  vol, mask = NULL, patchsize, verbose = TRUE,
  pad = TRUE,
  normalize = TRUE, contrast) {

  res = norm_pad(
    vol = vol, mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    pad = pad,
    normalize = normalize, contrast = contrast)
  mask = res$mask
  vol = res$vol

  num_patches = sum(mask != 0)

  bmask = blur_mask(mask, verbose = verbose)
  blurmask = bmask$blurred_mask
  newindx = bmask$indices
  dsize <- floor(patchsize / 2)

  matsize = get_matsize(num_patches, patchsize = patchsize)
  if (verbose) {
    message("Size matrix to create: ", paste(matsize, collapse = "x"))
  }
  t1_patches <- array(0, dim = matsize)
  mask_patches <- array(0, dim = matsize)
  for (i in 1:num_patches) {
    x <- newindx[1, i]
    y <- newindx[2, i]
    z <- newindx[3, i]
    t1_patches[i, , , 1] <-
      vol[(x - dsize[1]):(x + dsize[1]), (y - dsize[2]):(y + dsize[2]), z]
    mask_patches[i, , , 1] <-
      blurmask[(x - dsize[1]):(x + dsize[1]), (y - dsize[2]):(y + dsize[2]), z]
  }
  list(image_patches = t1_patches, mask_patches = mask_patches,
       blurred_mask = blurmask,
       indices = newindx,
       padded_mask = mask)
}

#' @rdname get_patch_from_volume
#' @export
get_3d_patch_from_volume <- function(
  vol, mask = NULL, patchsize, verbose = TRUE,
  pad = TRUE,
  normalize = TRUE, contrast) {

  res = norm_pad(
    vol = vol, mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    pad = pad,
    normalize = normalize, contrast = contrast)
  mask = res$mask
  vol = res$vol

  num_patches = sum(mask != 0)

  bmask = blur_mask(mask, verbose = verbose)
  blurmask = bmask$blurred_mask
  newindx = bmask$indices
  dsize <- floor(patchsize / 2)

  matsize = get_matsize(num_patches, patchsize = patchsize)
  if (verbose) {
    message("Size matrix to create: ", paste(matsize, collapse = "x"))
  }
  t1_patches <- array(0, dim = matsize)
  mask_patches <- array(0, dim = matsize)
  for (i in 1:num_patches) {
    x <- newindx[1, i]
    y <- newindx[2, i]
    z <- newindx[3, i]
    t1_patches[i, , , , 1] <-
      vol[(x - dsize[1]):(x + dsize[1]),
          (y - dsize[2]):(y + dsize[2]),
          (z - dsize[3]):(z + dsize[3])
          ]
    mask_patches[i, , , , 1] <-
      blurmask[
        (x - dsize[1]):(x + dsize[1]),
        (y - dsize[2]):(y + dsize[2]),
        (z - dsize[3]):(z + dsize[3])
        ]
  }
  list(image_patches = t1_patches, mask_patches = mask_patches,
       blurred_mask = blurmask,
       indices = newindx,
       padded_mask = mask)
}



#' @rdname get_patch_from_volume
#' @export
norm_pad = function(
  vol, mask = NULL, patchsize, verbose = TRUE,
  pad = TRUE,
  normalize = TRUE, contrast) {
  vol = check_nifti(vol, allow.array = TRUE)
  if (normalize) {
    vol = normalize_image(vol = vol, contrast = contrast, verbose = verbose)
  }
  if (is.null(mask)) {
    mask = vol != 0
  }
  mask = check_nifti(mask, allow.array = TRUE)
  if (pad) {
    padsize = patchsize_to_padsize(patchsize)
    vol <- pad_image(vol, padsize = padsize)
    mask <- pad_image(mask, padsize = padsize)
  }
  return(list(vol = vol, mask = mask))
}
