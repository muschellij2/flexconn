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
  matsize <- c(num_patches, patchsize, 1)
  return(matsize)
}