
#' Python Numpy Percentile
#'
#' @param x a numeric vector
#' @param probs Quantile to estimate
#'
#' @return A numeric
#' @export
#'
#' @examples
#' set.seed(11)
#' n = 45
#' x =rnorm(n)
#' x = sort(x)
#' q = numpy_percentile(x)
#' user = Sys.getenv("USER")
#' if (user == "johnmuschelli") {
#' reticulate::use_python(paste0(
#' "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#' }
#' np = reticulate::import("numpy")
#' q2 = np$percentile(x, 99L)
#' q2
#' q2 - q
#' quant = quantile(x, probs = 0.99)
#' q2 - quant
#'
#' np$percentile(sort(as.double(x)), 99L)
#' np$percentile(sort(x), 99L)
#'
numpy_percentile = function(x, probs = 0.99) {
  x = sort(x)
  # see
  # https://stackoverflow.com/a/48799350/2549133
  index = length(x) * probs
  frac = index %% 1
  i = floor(index)
  j = ceiling(index)
  q = x[i] + (x[j] - x[i]) * frac
  return(q)
}
