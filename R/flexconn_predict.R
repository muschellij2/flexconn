#' Predict from FLEXCONN model
#'
#' @param model A keras model object trained
#' @param t1 3D array or \code{nifti} image
#' @param flair 3D array or \code{nifti} image
#' @param mask binary 3D array or \code{nifti} image, but of
#' the whole area of prediction, not the ground truth
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param batch_size Size of batches for prediction.
#' Integer. Passed to
#' \code{\link[keras]{predict.keras.engine.training.Model}}.
#' @param ... additional arguments to \code{\link{get_patches}}
#'
#' @note If \code{mask = NULL}, a mask will be generated for
#' \code{t1 > 0}.
#' @return A vector of predictions, based on the indices of the
#' mask
#' @export
#'
#' @examples
flexconn_predict = function(
  model, t1, flair, mask = NULL,
  patchsize, verbose = TRUE, ..., batch_size = 1) {

  patches = get_patches(
    t1 = t1, flair = flair, patchsize = patchsize,
    mask = mask,
    ...)
  t1_test <- patches$t1_patches
  fl_test <- patches$fl_patches

  preds <- model %>% predict(
    list(t1_test,
         fl_test),
    verbose = verbose,
    batch_size = batch_size)
  return(preds)
}
