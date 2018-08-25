
#' Normalize Image Based on the Peak
#'
#' @param vol 3D array or \code{nifti} image
#' @param contrast What imaging sequence of MRI is this volume
#' @param verbose print diagnostic messages
#' @param gridsize Grid size to use kernel density estiamtion.
#' Passed to \code{\link[KernSmooth]{bkde}}
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
                            verbose = TRUE,
                            gridsize = 128) {
  vol = check_nifti(vol, allow.array = TRUE)
  vol = vol / image_peak(vol, contrast, gridsize = gridsize)$peak
}

#' @rdname normalize_image
#' @export
image_peak <- function(vol,
                       contrast = c("T1", "PD", "T2", "FLAIR", "FL"),
                       verbose = TRUE,
                       gridsize = 128) ## scipy rounds up from 80 to 128)
{

  contrast = toupper(contrast)
  contrast = match.arg(contrast)
  contrast = substr(contrast, 1, 2)

  vol = check_nifti(vol, allow.array = TRUE)

  temp <- vol[vol != 0]

  # temp = sort(temp)
  # index = length(temp) * 0.99
  # frac = index %% 1
  # i = floor(index)
  # j = ceiling(index)
  # q = temp[i] + (temp[j] - temp[i]) * frac

  q <- quantile(temp, probs = .99)
  temp <- temp[temp <= q]
  bw <- q / 80
  if (verbose) {
    message(paste0("99th quantile is ",
                   round(q, 3),
                   ", bandwidth =",
                   round(bw, 3),
                   ", gridsize =",
                   round(gridsize, 3)
    ))
  }
  kde <- KernSmooth::bkde(temp, bandwidth = bw, gridsize = gridsize)
  x <- kde$y
  y <- kde$x
  indx <- pracma::findpeaks(x)[, 2]
  heights <- x[indx]
  peaks <- y[indx]
  if (verbose) {
    message(paste0(length(peaks), " peaks found.", "\n"))
  }
  if (contrast == "T1") {
    peak <- peaks[length(peaks)]
  } else if (contrast %in% c("T2", "PD", "FL")) {
    # peak_height <- which.max(heights)
    peak <- peaks[which.max(heights)]
  }
  if (verbose) {
    message(paste("Peak found at", peak, "for", contrast))
  }
  L = list(peak = peak,
           q99 = q,
           bandwidth = bw)
  return(L)
}


#' @rdname normalize_image
#' @export
flexconn_image_peak <- function(vol,
                                contrast = c("T1", "PD", "T2", "FLAIR", "FL"),
                                verbose = TRUE,
                                gridsize = 80L) ## scipy rounds up from 80 to 128)
{
  contrast = toupper(contrast)
  contrast = match.arg(contrast)
  contrast = substr(contrast, 1, 2)

  sm = reticulate::import("statsmodels.api")
  np = reticulate::import("numpy")
  ss = reticulate::import("scipy.signal")

  vol = check_nifti(vol, allow.array = TRUE)

  temp <- vol[vol != 0]

  temp = sort(temp)
  index = length(temp) * 0.99
  frac = index %% 1
  i = floor(index)
  j = ceiling(index)
  q = temp[i] + (temp[j] - temp[i]) * frac

  temp <- temp[temp <= q]
  bw <- q / 80
  if (verbose) {
    message(paste0("99th quantile is ",
                   round(q, 3),
                   ", bandwidth =",
                   round(bw, 3),
                   ", gridsize =",
                   round(gridsize, 3)
    ))
  }

  kde = sm$nonparametric$KDEUnivariate(temp)

  kde$fit(kernel='gau', bw=bw, gridsize=gridsize, fft=TRUE)
  x_mat = 100.0 * kde$density
  y_mat = kde$support

  indx = ss$argrelextrema(x_mat, np$greater)
  indx = np$asarray(indx, dtype=np$int)
  heights = x_mat[indx]
  peaks = y_mat[indx]

  if (verbose) {
    message(paste0(length(peaks), " peaks found.", "\n"))
  }
  if (contrast == "T1") {
    peak <- peaks[length(peaks)]
  } else if (contrast %in% c("T2", "PD", "FL")) {
    # peak_height <- which.max(heights)
    peak <- peaks[which.max(heights)]
  }
  if (verbose) {
    message(paste("Peak found at", peak, "for", contrast))
  }
  L = list(peak = peak,
           q99 = q,
           bandwidth = bw)
  return(L)

}
