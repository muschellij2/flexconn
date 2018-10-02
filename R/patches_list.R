#' Make a patches list for a population of images
#'
#' @param t1 a character vector of filenames for T1 images
#' @param flair  a character vector of filenames for FLAIR images
#' @param mask  a character vector of filenames for mask images
#' @param t2 a character vector of filenames for T2 images
#' @param patchsize Vector of length 2 (or more)
#' @param verbose print diagnostic messages
#' @param outfile Output file to write arrays, an \code{.rds} file.
#' If this exists, it will be read in versus making patches again
#' @param write_file Should an output file be written
#' @param ... Additional arguments to pass to \code{\link{get_patches}}
#'
#' @importFrom RNifti readNifti
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @return A list of named arrays
#' @export
#'
#' @examples
#' user = Sys.getenv("USER")
#' if (user == "johnmuschelli") {
#'   setwd("/Volumes/DATA_LOCAL/Projects/ms_lesion_challenge/atlases/none")
#'   reticulate::use_python(paste0(
#'    "/Library/Frameworks/Python.framework/Versions/3.5/bin/python3"))
#'    i = 1:5
#' t1 <- paste0("atlas", i, "_T1.nii.gz")
#' flair = paste0("atlas", i, "_FL.nii.gz")
#' t2 = NULL
#' mask = paste0("atlas", i, "_mask.nii.gz")
#' verbose = TRUE
#' patchsize = c(35, 35)
#' outfile = NULL
#' write_file = TRUE
#' }
patches_list = function(
  t1, flair, mask,
  t2 = NULL,
  patchsize = c(35, 35),
  verbose = TRUE,
  outfile = NULL,
  write_file = TRUE,
  ...) {

  if (!is.null(outfile)) {
    if (file.exists(outfile)) {
      all_arrays = readr::read_rds(outfile)
      attr(all_arrays, "outfile") = outfile
      return(all_arrays)
    }
  } else {
    if (write_file) {
      outfile = tempfile(fileext = ".rds")
    }
  }


  ndim = length(patchsize)
  L = list()
  L$mask = mask
  L$t1 = t1
  L$fl = flair
  L$t2 = t2

  chars = sapply(L, is.character)
  if (!all(chars)) {
    stop("files need to be passed!")
  }
  # need all the same length
  lengths = sapply(L, length)
  stopifnot(length(unique(lengths)) == 1)


  masks = lapply(L$mask, RNifti::readNifti)
  mask_num_patches = sapply(masks, sum)
  num_patches = sum(mask_num_patches)

  rm(masks);
  for (i in 1:10) gc()
  if (verbose) {
    message(paste0("Total number of lesion patches = " , num_patches))
  }

  matsize = get_matsize(num_patches, patchsize)

  all_arrays = lapply(L, function(x) {
    array(0, dim = matsize)
  })
  n_array = names(all_arrays)
  count2 <- 1
  count1 <- 1

  i = 1
  if (verbose) {
    pb = txtProgressBar(min = 1, max = n_atlas)
  }
  for (i in seq(n_atlas)) {

    fnames = sapply(L, function(x) x[i])
    fnames = unlist(fnames)

    na_null = function(x) {
      if (is.na(x)) {
        return(NULL)
      }
      x
    }


    res = get_patches(
      t1 = na_null(fnames["t1"]),
      flair = na_null(fnames["fl"]),
      mask = na_null(fnames["mask"]),
      t2 = na_null(fnames["t2"]),
      patchsize = patchsize,
      only_patches = TRUE,
      ...
    )

    if (verbose) {
      message(paste0("Dim of T1 patches: ",
                     paste(dim(res$t1_patches), collapse = "x")
      ))
    }

    pdim <- nrow(res$t1_patches)
    if (pdim > 0) {

      count2 <- count1 + pdim - 1
      # if (verbose) {
      #   message(paste("Atlas", i, "indices:", count1, count2))
      # }


      ind = count1:count2
      if (ndim == 2) {
        for (iname in n_array) {
          all_arrays[[iname]][ ind, , , ] = res[[ paste0(iname, "_patches")]]
        }
      }

      if (ndim == 3) {
        for (iname in n_array) {
          all_arrays[[iname]][ ind, , , , ] = res[[ paste0(iname, "_patches")]]
        }
      }

      count1 <- count1 + pdim
    }
    if (verbose) {
      setTxtProgressBar(pb, value = i)
    }
  }
  if (verbose) {
    close(pb)
  }

  if (verbose) {
    message(paste0("Total number of patches collected = ", count2))
  }


  if (!is.null(outfile)) {
    readr::write_rds(
      all_arrays,
      path = outfile,
      compress = "xz")
    attr(all_arrays, "outfile") = outfile
  }
  return(all_arrays)
}
