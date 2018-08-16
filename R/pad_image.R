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
#' orig = invert_pad_image(padded, padsize = 2)
#' testthat::expect_equal(orig, vol)
pad_image <- function(vol, padsize) {
  vol = as.array(vol)
  dim_orig <- dim(vol)
  dim_pad <- dim_orig + 2 * padsize
  padded <- array(0, dim = dim_pad)
  padded[(padsize + 1):(dim_orig[1] + padsize), (padsize + 1):(dim_orig[2] + padsize),
         (padsize + 1):(dim_orig[3] + padsize)] <- vol
  padded
}

#' @rdname pad_image
#' @export
invert_pad_image <- function(vol, padsize) {
  dim_pad <- dim(vol)
  dim_orig <- dim_pad - 2 * padsize
  vol = as.array(vol)
  unpadded = vol[(padsize + 1):(dim_orig[1] + padsize),
                 (padsize + 1):(dim_orig[2] + padsize),
                 (padsize + 1):(dim_orig[3] + padsize)]
  unpadded
}
