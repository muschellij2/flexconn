#' FLEXCONN Convolution Chain
#'
#' @param object Model or layer object
#' @param ds Some fixed number???
#' @param num_filters Nubmer of filters
#' @param kernel_size_1 Size of first kernel
#' @param kernel_size_2 Size of second kernel
#' @param prefix output prefix for each layer
#' @param ndim Number of dimensions for convolution chain
#'
#' @return Output model
#' @export
#'
#' @examples
#' library(tensorflow)
#'  t1_input <- layer_input(shape = shape(NULL, NULL, 1))
#'  t1 <- t1_input %>%
#'  flexconn_conv_chain(
#'    prefix = "t1"
#'  )
#'
#'  t1_3d_input <- layer_input(shape = shape(NULL, NULL, NULL, 1))
#'  t1_3d <- t1_input %>%
#'  flexconn_conv_chain(
#'    prefix = "t1"
#'  )
#' @importFrom keras layer_conv_2d layer_conv_3d %>%
flexconn_conv_chain <- function(
  object,
  ds = 2,
  num_filters = 128,
  kernel_size_1 = 3,
  kernel_size_2 = 5,
  prefix = NULL,
  ndim ) {

  if (missing(ndim)) {
    ndim = object$shape$ndims - 2
  }

  ndim = ndim[1]
  func = switch(as.character(ndim),
                "2" = layer_conv_2d,
                "3" = layer_conv_3d)
  if (is.null(func)) {
    stop(paste0("ndim, number of dimensions: ", ndim, "not supported"))
  }


  object %>%
    func(
      filters = num_filters,
      kernel_size = kernel_size_1,
      activation = "relu",
      padding = "same",
      name = paste0(prefix, "_conv1")
    ) %>%
    func(
      filters = num_filters / ds,
      kernel_size = kernel_size_2,
      activation = "relu",
      padding = "same",
      name = paste0(prefix, "_conv2")
    ) %>%
    func(
      filters = num_filters / (ds * 2),
      kernel_size = kernel_size_1,
      activation = "relu",
      padding = "same",
      name = paste0(prefix, "_conv3")
    ) %>%
    func(
      filters = num_filters / (ds ^ 3),
      kernel_size = kernel_size_2,
      activation = "relu",
      padding = "same",
      name = paste0(prefix, "_conv4")
    ) %>%
    func(
      filters = num_filters / (ds ^ 4),
      kernel_size = kernel_size_1,
      activation = "relu",
      padding = "same",
      name = paste0(prefix, "_conv5")
    )

}
