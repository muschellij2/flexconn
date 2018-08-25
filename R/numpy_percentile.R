
#' Python Numpy Percentile
#'
#' @param x a numeric vector
#' @param probs Quantile to estimate
#'
#' @return A numeric
#' @export
#'
#' @examples
#' numpy_percentile(rnorm(100))
numpy_percentile = function(x, probs = 0.99) {
  x = sort(x)
  index = length(x) * probs
  frac = index %% 1
  i = floor(index)
  j = ceiling(index)
  q = x[i] + (x[j] - x[i]) * frac
  return(q)
}
