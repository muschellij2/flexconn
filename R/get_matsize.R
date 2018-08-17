#' Get Matrix Size to Input to Model
#'
#' @param num_patches Number of Patches
#' @param patchsize Patch Size
#'
#' @return A vector
#' @export
#'
#' @examples
#' get_matsize(5, c(35, 35))
get_matsize = function(num_patches, patchsize) {
  if (all(patchsize %% 2 != 1)) {
    stop("Patch Sizes must be odd!")
  }
  matsize <- c(num_patches, patchsize, 1)
  return(matsize)
}
