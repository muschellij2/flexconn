library(keras)
library(neurobase)
library(flexconn)


# Configuration -----------------------------------------------------------

ndim <- 2
psize <- 35
patchsize <- rep(psize, ndim)
padsize = patchsize_to_padsize(patchsize = patchsize)

model_dir <- "saved_models"
# change this to use a model of your choice
restore_path <- file.path(model_dir, "weights.05-51.56.hdf5")

# model configuration
batch_size <- 1


# Get predictions ---------------------------------------------------------

model <- load_model_hdf5(restore_path)

test_fl <- system.file("extdata/FLAIR.nii.gz", package = "flexconn")
test_t1 <-
  system.file("extdata/MPRAGE.nii.gz", package = "flexconn")

predicted <- model %>% flexconn_predict(t1 = test_t1,
                                        flair = test_fl,
                                        patchsize = patchsize)

mask <- remove_small_components(predicted)

image(predicted, z = 100, plot.type = "single")
image(mask, z = 100, plot.type = "single")
