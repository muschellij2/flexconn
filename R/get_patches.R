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
get_patches <- function(t1, flair, mask, patchsize,
                        pad = TRUE,
                        normalize = TRUE,
                        verbose = TRUE) {

  if (pad) {
    padsize = patchsize_to_padsize(patchsize)
    t1 <- pad_image(t1, padsize = padsize)
    flair <- pad_image(flair, padsize = padsize)
    mask <- pad_image(mask, padsize = padsize)
  }

  t1_patches = get_patch_from_volume(t1, mask,
                                     patchsize,
                                     verbose = verbose,
                                     normalize = normalize,
                                     contrast = "T1")
  fl_patches = get_patch_from_volume(flair, mask,
                                     patchsize,
                                     verbose = verbose,
                                     normalize = normalize,
                                     contrast = "FLAIR")

  list(t1_patches$image_patches,
       fl_patches$image_patches,
       fl_patches$mask_patches)
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
get_patch_from_volume <- function(vol, mask, patchsize, verbose = TRUE,
                                  normalize = TRUE, contrast) {

  vol = check_nifti(vol, allow.array = TRUE)
  if (normalize) {
    vol = normalize_image(vol = vol, contrast = contrast, verbose = verbose)
  }
  mask = check_nifti(mask, allow.array = TRUE)
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
  list(image_patches = t1_patches, mask_patches = mask_patches)
}
