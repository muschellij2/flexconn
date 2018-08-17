library(keras)
library(neurobase)
library(flexconn)

ndim = 2
run_model = flexconn_model(ndim = ndim)


run_model %>% compile(
  optimizer = optimizer_adam(lr =  0.0001),
  loss = "mean_squared_error",
  metrics = c("mean_squared_error")
)


# Configuration -----------------------------------------------------------

atlas_dir <- "lesion_challenge"
n_atlas <- 1
psize <- 35
model_dir <- "saved_models"

patchsize <- rep(psize, ndim)
padsize = patchsize_to_padsize(patchsize = patchsize)


# Read data ---------------------------------------------------------------

num_patches <- 0
for (i in 1:n_atlas) {
  p <-
    readnii(file.path(atlas_dir, paste0("atlas", i, "_mask.nii.gz")))
  num_patches <- num_patches + img_data(p) %>% sum()
}

cat("Total number of lesion patches =" , num_patches, "\n")

matsize = get_matsize(num_patches, patchsize)

t1_patches <- array(0, dim = matsize)
fl_patches <- array(0, dim = matsize)
mask_patches <- array(0, dim = matsize)

count2 <- 1
count1 <- 1

for (i in 1:n_atlas) {

  t1name <- file.path(atlas_dir, paste0("atlas", i, "_T1.nii.gz"))
  cat("Reading", t1name, "\n")
  t1 <- readnii(t1name) %>% img_data()
  flname <- file.path(atlas_dir, paste0("atlas", i, "_FL.nii.gz"))
  cat("Reading", flname, "\n")
  fl <- readnii(flname) %>% img_data()
  maskname <-
    file.path(atlas_dir, paste0("atlas", i, "_mask.nii.gz"))
  cat("Reading", maskname, "\n")
  mask <- readnii(maskname) %>% img_data()

  t1 <- normalize_image(t1, 'T1')
  fl <- normalize_image(fl, 'FL')

  padded_t1 <- pad_image(t1, padsize)
  padded_fl <- pad_image(fl, padsize)
  padded_mask <- pad_image(mask, padsize)

  cat("T1 orig dim: ", dim(t1), "\n")
  cat("T1 padded to dim: ", dim(padded_t1), "\n")

  c(t1_patches_a, fl_patches_a, mask_patches_a) %<-%
    get_patches(padded_t1, padded_fl, padded_mask, patchsize)

  cat("Dim of T1 patches:", dim(t1_patches_a), "\n")

  pdim <- nrow(t1_patches_a)
  if (pdim > 0) {
    count2 <- count1 + pdim - 1
    cat("Atlas", i, "indices:", count1, count2, "\n")

    if (ndim == 2) {
      t1_patches[  count1:count2, , , ] <- t1_patches_a
      fl_patches[  count1:count2, , , ] <- fl_patches_a
      mask_patches[count1:count2, , , ] <- mask_patches_a
    }

    if (ndim == 3) {
      t1_patches[  count1:count2, , , , ] <- t1_patches_a
      fl_patches[  count1:count2, , , , ] <- fl_patches_a
      mask_patches[count1:count2, , , , ] <- mask_patches_a
    }

    count1 <- count1 + pdim
  }

}

cat("Total number of patches collected = ", count2, "\n")
cat("Size of the input matrix is ", dim(mask_patches), "\n")




history <- model %>% fit(
  x = list(t1_patches, fl_patches),
  y = mask_patches,
  batch_size = batch_size,
  epochs = 10,
  validation_split = 0.3
)

plot(history, metrics = "loss")

preds <-
  model %>% predict(list(t1_patches[1:10, , , , drop = FALSE],
                         fl_patches[1:10, , , , drop = FALSE]))
dim(preds)
true <- mask_patches[1:10, , , , drop = FALSE]
cor(true, preds)
preds[1, , , ]
preds[1, , , ] %>% image(col = grey.colors(n=10))
true[1, , , ]
true[1, , , ] %>% image(col = grey.colors(n=10))

