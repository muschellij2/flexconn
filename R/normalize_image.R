
#' Normalize Image Based on the Peak
#'
#' @param t1 3D array or \code{nifti} image
#' @param contrast What imaging sequence of MRI is this volume
#' @param verbose print diagnostic messages
#'
#' @return The peak of the image is returned for \code{image_peak}
#' and the normalized image in \code{normalize_image}
#' @export
#' @importFrom pracma findpeaks
#' @importFrom KernSmooth bkde
#' @importFrom stats quantile
#'
#' @examples
#' set.seed(1)
#' vol = array(rpois(10^3, lambda = 10), dim = rep(10, 3))
#' vol[ vol < 6] = 0
#' normed = normalize_image(vol, "T1")
#'
normalize_image <- function(vol,
                            contrast = c("T1", "PD", "T2", "FLAIR", "FL"),
                            verbose = TRUE) {
  vol = vol / normalize_image(vol, contrast)
}

#' @rdname normalize_image
#' @export
image_peak <- function(vol,
                            contrast = c("T1", "PD", "T2", "FLAIR", "FL"),
                            verbose = TRUE) {

  contrast = toupper(contrast)
  contrast = match.arg(contrast)
  contrast = substr(contrast, 1, 2)

  temp <- vol[vol != 0]
  q <- quantile(temp, probs = .99)
  temp <- temp[temp <= q]
  bw <- q / 80
  gridsize <- 128 ## scipy rounds up from 80 to 128
  if (verbose) {
    message(paste0("99th quantile is",
                   q,
                   ", bandwidth =",
                   bw,
                   ", gridsize =",
                   gridsize,
                   "\n"))
  }
  kde <- KernSmooth::bkde(temp, bandwidth = bw, gridsize = gridsize)
  x <- kde$y
  y <- kde$x
  indx <- pracma::findpeaks(x)[, 2]
  heights <- x[indx]
  peaks <- y[indx]
  cat(length(peaks), "peaks found.", "\n")
  if (contrast == "T1") {
    peak <- peaks[length(peaks)]
  } else if (contrast %in% c("T2", "PD", "FL")) {
    # peak_height <- which.max(heights)
    peak <- peaks[which.max(heights)]
  }
  if (verbose) {
    message(paste0("Peak found at", peak, "for", contrast, "\n"))
  }
  peak
}

