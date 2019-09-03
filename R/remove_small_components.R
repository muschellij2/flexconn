
#' Remove connected objects with volume <27 voxels
#'
#' @param img A probability image from \code{\link{flexconn_predict_with_volume}}
#' or \code{\link{flexconn_predict}}, should be
#' object of class \code{nifti} or an array
#' @param threshold threshold to apply to image to create a mask (
#' less than or equal)
#' @param voxel_threshold Threshold of number of voxels to remove if smaller
#' than this number (strict inequality)
#'
#' @return An array or object of class \code{nifti}
#' @export
#'
#' @examples
#' user = Sys.getenv("USER")
#' if (user %in% c("johnmuschelli", "travis") &
#' Sys.info()["sysname"] == "Darwin"){
#' reticulate::use_python(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3")
#' } else {
#' python = system("which python", intern = TRUE)
#' print(python)
#' reticulate::use_python(python)
#' }
#' set.seed(100)
#' d = rep(50, 3)
#' img = array(runif(prod(d), max = 1), dim = d)
#' cc = remove_small_components(img, threshold = 0.9)
remove_small_components = function(img, threshold = 0.34,
                                   voxel_threshold = 27) {
  sc = reticulate::import("scipy")
  np <- reticulate::import("numpy")
  seg = as.array(img >= threshold)

  seg = reticulate::r_to_py(seg)
  seg = seg$astype(np$float32)

  se = sc$ndimage$morphology$generate_binary_structure(3L, 1L)
  lab = sc$ndimage$label(seg, structure = se)
  label = lab[[1]]
  ncomp = lab[[2]]
  tab = table(label)
  remove = as.numeric(names(tab)[ tab < voxel_threshold])
  label[ label %in% remove] = 0
  label[label > 0] = 1
  if (inherits(img, "nifti")) {
    label = niftiarr(img = img, label)
  }
  return(label)
}
