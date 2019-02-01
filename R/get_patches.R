#' Get Patches from T1 and FLAIR image
#'
#' @param t1 3D array or \code{nifti} image
#' @param flair 3D array or \code{nifti} image
#' @param t2 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param pad Run \code{\link{pad_image}} on the image
#' before getting patches (pads then normalizes if \code{normalize = TRUE})
#' @param normalize Run \code{\link{normalize_image}} on the image
#' before getting patches
#' @param only_patches Only return the patches, not additional information.
#' @param seed Seed for random sampling of indices.  If \code{NULL},
#' no sampling is done, passed to \code{\link{mask_indices}}
#' @param ... additional arguments to pass to
#' \code{\link{get_patch_from_volume}}
#'
#' @note If \code{mask = NULL}, a mask will be created based on
#' voxels greater than the 75th percentile of the FLAIR image.
#'
#' @return A list of T1, FLAIR, and Mask Patches
#' @export
#'
#' @importFrom reticulate import
#' @examples
#' user = Sys.getenv("USER")
#' if (user %in% c("johnmuschelli", "travis") &
#' Sys.info()["sysname"] == "Darwin"){
#' reticulate::use_python(paste0(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#' } else {
#' python = system("which python", intern = TRUE)
#' reticulate::use_python(python)
#' }
#' library(neurobase)
#' fname = system.file("extdata", "MPRAGE.nii.gz", package = "flexconn")
#' t1 = readnii(fname)
#' fname = system.file("extdata", "FLAIR.nii.gz", package = "flexconn")
#' flair = readnii(fname)
#' patchsize = c(3, 3)
#' verbose = TRUE
#' patch = get_patches(t1, flair, mask = NULL, patchsize = patchsize)
#'
#' rm(patch)
get_patches <- function(
  t1, flair = NULL, t2 = NULL,
  mask = NULL, patchsize,
  pad = TRUE,
  normalize = TRUE,
  verbose = TRUE,
  only_patches = FALSE,
  seed = NULL,
  ...) {

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
    pad = pad,
    seed = seed,
    run_mask_patches = TRUE,
    ...)
  mask_patches = t1_patches$mask_patches
  t1_patches = t1_patches$image_patches
  fl_patches = NULL
  if (!is.null(flair)) {
    fl_patches = get_patch_from_volume(
      flair, mask = mask,
      patchsize = patchsize,
      verbose = verbose,
      normalize = normalize,
      contrast = "FLAIR",
      pad = pad,
      seed = seed,
      run_mask_patches = FALSE,
      ...)
  }

  t2_patches = NULL
  if (!is.null(t2)) {
    t2_patches = get_patch_from_volume(
      t2, mask = mask,
      patchsize = patchsize,
      verbose = verbose,
      normalize = normalize,
      contrast = "T2",
      pad = pad,
      seed = seed,
      run_mask_patches = FALSE,
      ...)
  }
  L = list(
    t1_patches = t1_patches
  )
  L$fl_patches = fl_patches$image_patches
  L$mask_patches = mask_patches
  L$t2_patches = t2_patches

  if (!only_patches) {
    L$blurred_mask = fl_patches$blurred_mask
    L$indices = fl_patches$indices
    L$mask = mask
    L$patchsize = patchsize
    L$pad = pad
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
#' @param seed Seed for random sampling of indices.  If \code{NULL},
#' no sampling is done, passed to \code{\link{mask_indices}}
#' @param run_mask_patches Should the patches for the mask be
#' made as well?
#'
#' @param ... not used
#'
#' @return A list of image and mask Patches
#' @export
#'
#' @importFrom reticulate import
#' @importFrom neurobase check_nifti
#' @examples
#' user = Sys.getenv("USER")
#' if (user == "johnmuschelli") {
#' reticulate::use_python(paste0(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#' mpl = reticulate::import("matplotlib")
#' mpl$use('TkAgg')
#' }
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
get_patch_from_volume <- function(
  vol, mask = NULL, patchsize, verbose = TRUE,
  pad = TRUE,
  normalize = TRUE, contrast,
  seed = NULL,
  run_mask_patches = TRUE,
  peak_estimator = c("Python", "R")
) {

  ndim = length(patchsize)
  if (!ndim %in% c(2, 3)) {
    stop("Patch Size is not length 2 or 3!")
  }

  res = norm_pad(
    vol = vol, mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    pad = pad,
    normalize = normalize, contrast = contrast,
    peak_estimator = peak_estimator)
  mask = res$mask
  vol = res$vol

  # num_patches = get_num_patches(mask)

  indices = mask_indices(mask, seed = seed, verbose = verbose)

  # blurmask = bmask$blurred_mask
  # indices = bmask$indices
  # dsize <- floor(patchsize / 2)

  t1_patches = volume_to_patches(
    vol = vol,
    indices = indices,
    patchsize = patchsize,
    verbose = verbose)

  mask_patches = NULL
  blurmask = NULL
  if (run_mask_patches) {
    blurmask = blur_mask(mask)
    mask_patches = volume_to_patches(
      vol = blurmask,
      indices = indices,
      patchsize = patchsize,
      verbose = verbose)
  }

  L = list(image_patches = t1_patches)
  L$mask_patches = mask_patches
  L$blurred_mask = blurmask
  L$indices = indices
  L$padded_mask = mask
  L$padded_vol = vol
  return(L)
}




#' @rdname get_patch_from_volume
#' @export
get_num_patches = function(mask) {
  mask = check_nifti(mask, allow.array = TRUE)
  num_patches = sum(mask != 0)
  return(num_patches)
}



#' @rdname get_patch_from_volume
#' @param peak_estimator Which functions to use to estimate peak,
#' either wrapper Python code or native R code.
#' @export
norm_pad = function(
  vol,
  mask = NULL,
  patchsize,
  verbose = TRUE,
  pad = TRUE,
  normalize = TRUE,
  peak_estimator = c("Python", "R"),
  contrast) {
  vol = check_nifti(vol, allow.array = TRUE)
  if (normalize) {
    vol = normalize_image(vol = vol, contrast = contrast,
                          verbose = verbose,
                          peak_estimator = peak_estimator)
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


#' @rdname get_patch_from_volume
#' @export
get_mask_patches = function(mask, patchsize,
                            pad = TRUE, seed = NULL,
                            verbose = TRUE) {

  ndim = length(patchsize)
  if (!ndim %in% c(2, 3)) {
    stop("Patch Size is not length 2 or 3!")
  }

  mask = check_nifti(mask, allow.array = TRUE)
  if (pad) {
    padsize = patchsize_to_padsize(patchsize)
    mask <- pad_image(mask, padsize = padsize)
  }
  # num_patches = get_num_patches(mask)

  blurmask = blur_mask(mask)
  indices = mask_indices(mask, seed = seed, verbose = verbose)

  mask_patches = volume_to_patches(
    vol = blurmask,
    indices = indices,
    patchsize = patchsize,
    verbose = verbose)


  list(mask_patches = mask_patches,
       blurred_mask = blurmask,
       indices = indices,
       padded_mask = mask)
}
