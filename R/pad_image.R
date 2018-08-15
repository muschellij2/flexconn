#' Pad Image
#'
#' @param vol 3D Array or matrix
#' @param padsize Scalar number
#'
#' @return Array with same number of dimensions as \code{vol} but increased
#' dimensions.
#' @export
#'
#' @examples
#' vol = array(rnorm(10^3), dim = rep(10, 3))
#' dim(vol)
#' padded = pad_image(vol, padsize = 2)
#' dim(padded)
pad_image <- function(vol, padsize) {
  vol = as.array(vol)
  dim_orig <- dim(vol)
  dim_pad <- dim_orig + 2 * padsize
  padded <- array(0, dim = dim_pad)
  padded[(padsize + 1):(dim_orig[1] + padsize), (padsize + 1):(dim_orig[2] + padsize),
         (padsize + 1):(dim_orig[3] + padsize)] <- vol
  padded
}
