library(keras)
library(neurobase)
library(flexconn)
library(RNifti)

ndim = 2
run_model = flexconn_model(ndim = ndim)


run_model %>% compile(
  optimizer = optimizer_adam(lr =  0.0001),
  loss = "mean_squared_error",
  metrics = c("mean_squared_error")
)


# Configuration -----------------------------------------------------------

atlas_dir <- "lesion_challenge/atlas_with_mask1"
n_atlas <- 5

psize <- 35
patchsize <- rep(psize, ndim)
padsize = patchsize_to_padsize(patchsize = patchsize)

model_exists <- TRUE

model_dir <- "saved_models"
# change this to use a model of your choice
restore_path <- file.path(model_dir, "weights.05-51.56.hdf5")
tmp_data_dir <- "tmp_data"

# model configuration
batch_size <- 128


# Read data ---------------------------------------------------------------

num_patches <- 0
for (i in 1:n_atlas) {
  p <- readNifti(file.path(atlas_dir, paste0("atlas", i, "_mask.nii.gz")))
  num_patches <- num_patches + sum(p)
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

  c(t1_patches_a, fl_patches_a, mask_patches_a) %<-%
    get_patches(
      t1 = t1,
      fl = fl,
      mask = mask,
      patchsize = patchsize,
      normalize = TRUE,
      pad = TRUE,
      only_patches = TRUE
    )

  cat("Dim of T1 patches:", dim(t1_patches_a), "\n")

  pdim <- nrow(t1_patches_a)
  if (pdim > 0) {
    count2 <- count1 + pdim - 1
    cat("Atlas", i, "indices:", count1, count2, "\n")

    if (ndim == 2) {
      t1_patches[count1:count2, , ,] <- t1_patches_a
      fl_patches[count1:count2, , ,] <- fl_patches_a
      mask_patches[count1:count2, , ,] <- mask_patches_a
    }

    if (ndim == 3) {
      t1_patches[count1:count2, , , ,] <- t1_patches_a
      fl_patches[count1:count2, , , ,] <- fl_patches_a
      mask_patches[count1:count2, , , ,] <- mask_patches_a
    }

    count1 <- count1 + pdim
  }

}

cat("Total number of patches collected = ", count2, "\n")
cat("Size of the input matrix is ", dim(mask_patches), "\n")


# Optional train-test split --------------------------------------------------------

train_indx <- sample(1:num_patches, num_patches * 0.7)
# optionally save indices
save_split <- FALSE
if (save_split) {
  train_indx %>% saveRDS(file.path(tmp_data_dir, "train_indx.rds"))
} else {
  train_indx <- readRDS(file.path(tmp_data_dir, "train_indx.rds"))
}

c(t1_train, t1_test) %<-% list(t1_patches[train_indx, , , , drop = FALSE], t1_patches[-train_indx, , ,  , drop = FALSE])
c(fl_train, fl_test) %<-% list(fl_patches[train_indx, , ,  , drop = FALSE], fl_patches[-train_indx, , , , drop = FALSE])
c(mask_train, mask_test) %<-% list(mask_patches[train_indx, , , , drop = FALSE], mask_patches[-train_indx, , ,  , drop = FALSE])



# Train or load model -----------------------------------------------------


if (!model_exists) {
  model <- flexconn_model()
  model %>% compile(
    optimizer = optimizer_adam(lr =  0.0001),
    loss = "mean_squared_error",
    metrics = c("mean_squared_error")
  )
  history <- model %>% fit(
    x = list(t1_train, fl_train),
    y = mask_train,
    batch_size = batch_size,
    epochs = 10,
    validation_split = 0.2,
    callbacks = list(
      callback_model_checkpoint(filepath = "weights.{epoch:02d}-{val_loss:.2f}.hdf5"),
      callback_early_stopping(patience = 1)
    )
  )
  saveRDS(history, file.path(tmp_data_dir, "history.rds"))
} else {
  model <- load_model_hdf5(restore_path)
  history <- readRDS(file.path(tmp_data_dir, "history.rds"))
}


# Evaluate model ----------------------------------------------------------


plot(history, metrics = "loss")

# model %>% evaluate(x = list(t1_test,
#                             fl_test),
#                    y = mask_test,
#                    batch_size = 1)



# Get predictions ---------------------------------------------------------

test_fl <- system.file("extdata/FLAIR.nii.gz", package = "flexconn")
test_t1 <-
  system.file("extdata/MPRAGE.nii.gz", package = "flexconn")

predicted <- model %>% flexconn_predict(t1 = test_t1,
                                        flair = test_fl,
                                        patchsize = c(35, 35))


image(predicted, z = 100, plot.type = "single")
