#' Predict from FLEXCONN model
#'
#' @param model A keras model object trained
#' @param t1 3D array or \code{nifti} image
#' @param flair 3D array or \code{nifti} image
#' @param t2 3D array or \code{nifti} image, optional if model has this
#'
#' @param mask binary 3D array or \code{nifti} image, but of
#' the whole area of prediction, not the ground truth
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param batch_size Size of batches for prediction.
#' Integer. Passed to
#' \code{\link[keras]{predict.keras.engine.training.Model}}.
#' @param ... additional arguments to \code{\link{get_patches}}
#' @param type type of prediction to use, patch-based or slice/volume
#' based
#' @note If \code{mask = NULL}, a mask will be generated for
#' \code{t1 > 0}.
#' @return A vector of predictions, based on the indices of the
#' mask
#' @export
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict
#' @rdname flexconn_predict
flexconn_predict_patch = function(
  model, t1, flair, t2 = NULL,
  mask = NULL,
  patchsize, verbose = TRUE, ..., batch_size = 1) {

  preds = flexconn_predict_with_patches(
    model = model,
    t1 = t1,
    flair = flair,
    t2 = t2,
    mask = mask,
    patchsize = patchsize,
    verbose = verbose,
    ... = ...,
    batch_size = batch_size
  )
  preds = preds$preds
  return(preds)
}


#' @rdname flexconn_predict
#' @export
flexconn_predict = function(
  model, t1, flair, t2 = NULL,
  mask = NULL,
  type = c("volume", "patch"),
  patchsize, verbose = TRUE,
  normalize = TRUE, ..., batch_size = 1) {

  n_images = tryCatch({
    length(model$input)
  })
  if (inherits(n_images, "try-error")) {
    n_images = length(model$input_layers)
  }
  if ( n_images == 2 & !is.null(t2)) {
    stop("T2 can't be specified - only 2 images as inputs")
  }
  type = match.arg(type)
  if (type == "patch") {
    res = flexconn_predict_with_patches(
      model = model,
      t1 = t1,
      flair = flair,
      t2 = t2,
      mask = mask,
      patchsize = patchsize,
      verbose = verbose,
      normalize = normalize,
      ... = ...,
      batch_size = batch_size
    )
    preds = res$preds

    if (verbose) {
      message("Creating Volume from Patch Prediction")
    }
    vol = get_volume_from_patches(
      patches = preds,
      indices = res$indices,
      patchsize = patchsize,
      verbose = verbose,
      invert_pad = res$pad
    )
  }

  if (type == "volume") {
    vol = flexconn_predict_with_volume(
      model = model, t1 = t1, flair = flair, t2 = t2,
      verbose = verbose, normalize = normalize,
      ... = ..., batch_size = batch_size)
  }
  return(vol)

}

#' @rdname flexconn_predict
#' @export
flexconn_predict_with_patches = function(
  model, t1, flair, t2 = NULL,
  mask = NULL,
  patchsize, verbose = TRUE, ..., batch_size = 1) {

  if (verbose) {
    message("Creating Patches")
  }
  patches = get_patches(
    t1 = t1, flair = flair,
    t2 = t2, patchsize = patchsize,
    mask = mask,
    only_patches = FALSE,
    ...)
  t1_test <- patches$t1_patches
  fl_test <- patches$fl_patches
  t2_test <- patches$t2_patches

  L = list(t1_test,
           fl_test)
  L$t2_test = t2_test
  if (verbose) {
    message("Predicting from Patches")
  }
  preds <- model %>%
    predict(L,
            verbose = verbose,
            batch_size = batch_size)
  patches$preds = preds
  return(patches)
}


#' @rdname flexconn_predict
#' @importFrom neurobase niftiarr
#' @param normalize Run \code{\link{normalize_image}} on the image
#' before prediction
#' @export
flexconn_predict_with_volume = function(
  model, t1, flair, t2 = NULL,
  verbose = TRUE, normalize = TRUE,
  ..., batch_size = 1) {

  n_images = tryCatch({
    length(model$input)
  })
  if (inherits(n_images, "try-error")) {
    n_images = length(model$input_layers)
  }
  if ( n_images == 2 & !is.null(t2)) {
    stop("T2 can't be specified - only 2 images as inputs")
  }

  # patches and last channel
  ndim = length(model$input_shape[[1]]) - 2
  t1 = check_nifti(t1)
  flair = check_nifti(flair)
  if (!is.null(t2)) {
    t2 = check_nifti(t2)
  }
  if (normalize) {
    if (verbose) {
      message("Normalizing Images")
    }
    t1 = normalize_image(
      vol = t1, contrast = "T1",
      verbose = verbose,
      ...)
    flair = normalize_image(
      vol = flair, contrast = "FLAIR",
      verbose = verbose,
      ...)
    if (!is.null(t2)) {
      t2 = normalize_image(
        vol = t2, contrast = "T2",
        verbose = verbose,
        ...)
    }
  }

  dims = dim(t1)
  arr_dim = c(1, dims[seq(ndim)], 1)

  slice_2d_t1 = array(0, dim = arr_dim)
  slice_2d_t2 = slice_2d_fl = slice_2d_t1

  output_image = array(dim = dims)

  if (verbose) {
    message("Predicting from model")
  }
  if (ndim == 2) {
    n_slices = dims[3]
    k = n_slices
    if (verbose) {
      pb = txtProgressBar(min = 0, max = n_slices)
    }
    for (k in seq(n_slices)) {
      slice_2d_t1[1, , , 1] = t1[, , k]
      slice_2d_fl[1, , , 1] = flair[,, k]
      slice_2d_t2 = NULL
      if (!is.null(t2)) {
        slice_2d_t2[1, , , 1] = t2[,, k]
      }
      L = list(t1 = slice_2d_t1,
               flair = slice_2d_fl)
      L$t2 = slice_2d_t2
      names(L) = NULL
      pred = model %>%
        predict(L,
                verbose = FALSE, # pb does it
                batch_size = batch_size)
      if (verbose) {
        setTxtProgressBar(pb, value = k)
      }
      output_image[, , k] = pred[1, , , 1]
    }
    if (verbose) {
      close(pb)
    }
  }

  if (ndim == 3) {
    slice_2d_t1[1, , , , 1] = t1
    slice_2d_fl[1, , , , 1] = flair
    slice_2d_t2 = NULL
    if (!is.null(t2)) {
      slice_2d_t2[1, , , , 1] = t2
    }
    L = list(t1 = slice_2d_t1,
             flair = slice_2d_fl)
    L$t2 = slice_2d_t2
    names(L) = NULL
    pred = model %>%
      predict(L,
              verbose = verbose,
              batch_size = batch_size)
    output_image[, , k] = pred[1, , , , 1]
  }

  output_image = niftiarr(img = t1, arr = output_image)
  return(output_image)
}
