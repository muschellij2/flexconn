#' Get Patches from T1 and FLAIR image
#'
#' @param t1 3D array or \code{nifti} image
#' @param flair 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#'
#' @return A list of T1, FLAIR, and Mask Patches
#' @export
#'
#' @importFrom reticulate import
get_patches <- function(t1, flair, mask, patchsize, verbose = TRUE) {


  t1_patches = get_patch_from_volume(t1, mask,
                                     patchsize, verbose)
  fl_patches = get_patch_from_volume(flair, mask,
                                     patchsize, verbose)

  list(t1_patches$image_patches,
       fl_patches$image_patches,
       fl_patches$mask_patches)
}

#' Get Patches from 3D Volume
#'
#' @param vol 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#'
#' @return A list of image and mask Patches
#' @export
#'
#' @importFrom reticulate import
#' @importFrom neurobase check_nifti
get_patch_from_volume <- function(vol, mask, patchsize, verbose = TRUE) {

  vol = check_nifti(vol, allow.array = TRUE)
  mask = check_nifti(mask, allow.array = TRUE)
  num_patches = sum(mask != 0)

  bmask = blur_mask(mask, verbose = verbose)
  blurmask = bmask$blurred_mask
  newindx = bmask$indices
  dsize <- floor(patchsize / 2)

  matsize = get_matsize(num_patches, patchsize = patchsize)
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
