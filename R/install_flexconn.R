#' Install FLEXCONN Requirements
#'
#' @param ... arguments to pass to \code{\link[reticulate]{py_install}}
#' @param extra_packages Additional packages to install
#'
#' @return Nothing
#' @export
#'
#' @examples \dontrun{
#' install_flexconn()
#' }
install_flexconn <- function(..., extra_packages = NULL) {
  packages = c(
    "numpy",
    "tensorflow",
    "keras",
    "scipy",
    "statsmodels",
    "tqdm",
    "nibabel",
    "h5py",
    "scikit-learn",
    "termcolor",
    "travis",
    "h5py",
    "pyyaml",
    "Pillow"
  )
  packages = c(packages, extra_packages)

  reticulate::py_install(packages = packages, ...)
}
