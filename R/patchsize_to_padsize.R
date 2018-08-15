#' Patch Size to Pad Size
#'
#' @param patchsize A vector of patch sizes
#'
#' @return A scalar number

#' @export
#' @examples
#' patchsize_to_padsize(c(10, 5))
#' patchsize_to_padsize(c(5, 5))
patchsize_to_padsize = function(patchsize) {
  padsize <- max(patchsize + 1) / 2
  return(padsize)
}
