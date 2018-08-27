#' FLEXCONN Model
#'
#' @param ndim Number of dimensions for convolution chain (2 or 3)
#' @param ds Some fixed number???
#' @param num_filters Nubmer of filters
#' @param kernel_size_1 Size of first kernel
#' @param kernel_size_2 Size of second kernel
#' @param prefixes The prefixes for each of the layers of the model.
#'
#' @return A keras model object
#' @export
#'
#' @importFrom keras layer_concatenate layer_input keras_model
#' @examples
#' user = Sys.getenv("USER")
#' if (user == "johnmuschelli") {
#' reticulate::use_python(paste0(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#' }
#' res = flexconn_model()
#' res = flexconn_model(prefixes = c("t1", "fl", "t2"))
flexconn_model = function(
  prefixes = c("t1", "fl"),
  ndim = 2,
  ds = 2,
  num_filters = 128,
  kernel_size_1 = 3,
  kernel_size_2 = 5
  ) {

  args = lapply(seq(ndim), function(x) NULL)
  args = c(args, 1)
  shaper = do.call(tensorflow::shape, args = args)

  ndim = ndim[1]
  func = switch(as.character(ndim),
                "2" = layer_conv_2d,
                "3" = layer_conv_3d)
  if (is.null(func)) {
    stop(paste0("ndim, number of dimensions: ", ndim, "not supported"))
  }

  outs = lapply(prefixes, function(prefix) {
    t1_input <- layer_input(shape = shaper)
    t1 <- t1_input %>%
      flexconn_conv_chain(
        ds = ds,
        num_filters = num_filters,
        kernel_size_1 = kernel_size_1,
        kernel_size_2 = kernel_size_2,
        prefix = prefix,
        ndim = ndim
      )
    return(list(input = t1_input, mod = t1))
  })

  # t1_input <- layer_input(shape = shaper)
  # t1 <- t1_input %>%
  #   flexconn_conv_chain(
  #     ds = ds,
  #     num_filters = num_filters,
  #     kernel_size_1 = kernel_size_1,
  #     kernel_size_2 = kernel_size_2,
  #     prefix = "t1",
  #     ndim = ndim
  #   )
  #
  # fl_input <- layer_input(shape = shaper)
  # fl <- fl_input %>%
  #   flexconn_conv_chain(
  #     ds = ds,
  #     num_filters = num_filters,
  #     kernel_size_1 = kernel_size_1,
  #     kernel_size_2 = kernel_size_2,
  #     prefix = "fl",
  #     ndim = ndim
  #   )
  # concat <- layer_concatenate(list(t1, fl), axis = -1)

  cat_list = lapply(outs, function(x) {
    x$mod
  })

  concat <- layer_concatenate(cat_list, axis = -1)

  combined <- concat %>%
    flexconn_conv_chain(
      ds = ds,
      num_filters = num_filters,
      kernel_size_1 = kernel_size_1,
      kernel_size_2 = kernel_size_2,
      prefix = "combined",
      ndim = ndim
    ) %>%
    func(
      filters = 1,
      kernel_size = 3,
      activation = "relu",
      padding = "same",
      name = "conv_final"
    )
  # inputs = list(t1_input, fl_input)
  inputs = lapply(outs, function(x) {
    x$input
  })
  model <- keras_model(inputs = inputs, outputs = combined)
  return(model)
}
