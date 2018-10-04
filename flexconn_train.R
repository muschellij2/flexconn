reticulate::use_python(
  "/software/apps/anaconda/5.2/python/3.6/bin/python")

library(keras)
library(neurobase)
library(flexconn)
library(RNifti)


# Configuration -----------------------------------------------------------

atlas_dir <- "lesion_challenge/atlas_with_mask1"
n_atlas <- 5

ndim <- 2
psize <- 35
patchsize <- rep(psize, ndim)
padsize = patchsize_to_padsize(patchsize = patchsize)

model_dir <- "saved_models"
tmp_data_dir <- "tmp_data"

# model configuration
batch_size <- 128


t1 = list.files("_T1", path = atlas_dir,
  full.names = TRUE)
flair = list.files("_FL", path = atlas_dir,
  full.names = TRUE)
mask = list.files("_mask", path = atlas_dir,
  full.names = TRUE)

outfile = file.path(tmp_data_dir, "patches.rds")
train = patches_list(
  t1 = t1,
  flair = flair,
  mask = mask,
  outfile = outfile
  )
num_patches = nrow(train$mask)

train_indx <- sample(1:num_patches,
  floor(num_patches * 0.7))

# Read data ---------------------------------------------------------------

# num_patches <- 0
# for (i in 1:n_atlas) {
#   p <- readNifti(file.path(atlas_dir, paste0("atlas", i, "_mask.nii.gz")))
#   num_patches <- num_patches + sum(p)
# }

# cat("Total number of lesion patches =" , num_patches, "\n")

# matsize = get_matsize(num_patches, patchsize)

# t1_patches <- array(0, dim = matsize)
# fl_patches <- array(0, dim = matsize)
# mask_patches <- array(0, dim = matsize)

# count2 <- 1
# count1 <- 1

# for (i in 1:n_atlas) {
#   t1name <- file.path(atlas_dir, paste0("atlas", i, "_T1.nii.gz"))
#   cat("Reading", t1name, "\n")
#   t1 <- readnii(t1name) %>% img_data()
#   flname <- file.path(atlas_dir, paste0("atlas", i, "_FL.nii.gz"))
#   cat("Reading", flname, "\n")
#   fl <- readnii(flname) %>% img_data()
#   maskname <-
#     file.path(atlas_dir, paste0("atlas", i, "_mask.nii.gz"))
#   cat("Reading", maskname, "\n")
#   mask <- readnii(maskname) %>% img_data()

#   c(t1_patches_a, fl_patches_a, mask_patches_a) %<-%
#     get_patches(
#       t1 = t1,
#       fl = fl,
#       mask = mask,
#       patchsize = patchsize,
#       normalize = TRUE,
#       pad = TRUE,
#       only_patches = TRUE
#     )

#   cat("Dim of T1 patches:", dim(t1_patches_a), "\n")

#   pdim <- nrow(t1_patches_a)
#   if (pdim > 0) {
#     count2 <- count1 + pdim - 1
#     cat("Atlas", i, "indices:", count1, count2, "\n")

#     if (ndim == 2) {
#       t1_patches[count1:count2, , ,] <- t1_patches_a
#       fl_patches[count1:count2, , ,] <- fl_patches_a
#       mask_patches[count1:count2, , ,] <- mask_patches_a
#     }

#     if (ndim == 3) {
#       t1_patches[count1:count2, , , ,] <- t1_patches_a
#       fl_patches[count1:count2, , , ,] <- fl_patches_a
#       mask_patches[count1:count2, , , ,] <- mask_patches_a
#     }

#     count1 <- count1 + pdim
#   }

# }

# cat("Total number of patches collected = ", count2, "\n")
# cat("Size of the input matrix is ", dim(mask_patches), "\n")


# Optional train-test split --------------------------------------------------------

# train_indx <- sample(1:num_patches,
#   floor(num_patches * 0.7))
# # optionally save indices
# save_split <- FALSE
# if (save_split) {
#   train_indx %>% saveRDS(file.path(tmp_data_dir, "train_indx.rds"))
# } else {
#   train_indx <- readRDS(file.path(tmp_data_dir, "train_indx.rds"))
# }

# c(t1_train, t1_test) %<-% list(t1_patches[train_indx, , , , drop = FALSE], t1_patches[-train_indx, , ,  , drop = FALSE])
# c(fl_train, fl_test) %<-% list(fl_patches[train_indx, , ,  , drop = FALSE], fl_patches[-train_indx, , , , drop = FALSE])
# c(mask_train, mask_test) %<-% list(mask_patches[train_indx, , , , drop = FALSE], mask_patches[-train_indx, , ,  , drop = FALSE])

test = lapply(train, function(x) {
  x[-train_indx, , , , drop = FALSE]
})

train = lapply(train, function(x) {
  x[train_indx, , , , drop = FALSE]
})



model <- flexconn_model()
model %>% compile(
  optimizer = optimizer_adam(lr =  0.0001),
  loss = "mean_squared_error",
  metrics = c("mean_squared_error")
)

history <- model %>% fit(
  x = unname(train[ c("t1", "fl")]),
  y = train$mask,
  batch_size = batch_size,
  epochs = 10,
  validation_split = 0.2,
  callbacks = list(
    callback_model_checkpoint(
      filepath = file.path(model_dir, "weights.{epoch:02d}-{val_loss:.2f}.hdf5")
    ),
    callback_early_stopping(patience = 1)
  )
)

# Train model -----------------------------------------------------

# ds <- 2
# num_filters <- 128
# kernel_size_1 <- 3
# kernel_size_2 <- 5
#  batch_size <- 128
#  conv_chain <- function(prev_layer,
#                        ds,
#                        num_filters,
#                        kernel_size_1,
#                        kernel_size_2,
#                        prefix = NULL) {
#   prev_layer %>%
#     layer_conv_2d(
#       filters = num_filters,
#       kernel_size = kernel_size_1,
#       activation = "relu",
#       padding = "same",
#       name = paste0(prefix, "_conv1")
#     ) %>%
#     layer_conv_2d(
#       filters = num_filters / ds,
#       kernel_size = kernel_size_2,
#       activation = "relu",
#       padding = "same",
#       name = paste0(prefix, "_conv2")
#     ) %>%
#     layer_conv_2d(
#       filters = num_filters / (ds * 2),
#       kernel_size = kernel_size_1,
#       activation = "relu",
#       padding = "same",
#       name = paste0(prefix, "_conv3")
#     ) %>%
#     layer_conv_2d(
#       filters = num_filters / (ds ^ 3),
#       kernel_size = kernel_size_2,
#       activation = "relu",
#       padding = "same",
#       name = paste0(prefix, "_conv4")
#     ) %>%
#     layer_conv_2d(
#       filters = num_filters / (ds ^ 4),
#       kernel_size = kernel_size_1,
#       activation = "relu",
#       padding = "same",
#       name = paste0(prefix, "_conv5")
#     )
#  }
# t1_input <- layer_input(
#   shape = shape(NULL, NULL, 1))
# t1 <- t1_input %>%
#   conv_chain(
#     ds = ds,
#     num_filters = num_filters,
#     kernel_size_1 = kernel_size_1,
#     kernel_size_2 = kernel_size_2,
#     prefix = "t1"
#   )
#  fl_input <- layer_input(shape = shape(NULL, NULL, 1))
# fl <- fl_input %>%
#   conv_chain(
#     ds = ds,
#     num_filters = num_filters,
#     kernel_size_1 = kernel_size_1,
#     kernel_size_2 = kernel_size_2,
#     prefix = "fl"
#   )
#  concat <- layer_concatenate(list(t1, fl), axis = -1)
#  combined <- concat %>%
#   conv_chain(
#     ds = ds,
#     num_filters = num_filters,
#     kernel_size_1 = kernel_size_1,
#     kernel_size_2 = kernel_size_2,
#     prefix = "combined"
#   ) %>%
#   layer_conv_2d(
#     filters = 1,
#     kernel_size = 3,
#     activation = "relu",
#     padding = "same",
#     name = "conv_final"
#   )
# model <-
#   keras_model(inputs = list(t1_input, fl_input),
#     outputs = combined)
# model %>% compile(
#   optimizer = optimizer_adam(lr =  0.0001),
#   loss = "mean_squared_error",
#   metrics = c("mean_squared_error")
# )
#  history <- model %>% fit(
#   x = unname(train[ c("t1", "flair")]),
#   y = train$mask,
#   batch_size = batch_size,
#   epochs = 10,
#   validation_split = 0.3,
#   verbose = 2
# )

saveRDS(history, file.path(tmp_data_dir, "history.rds"))
model %>% save_model_hdf5(file.path(model_dir, "flexconn_final.hdf5"))


# Evaluate model ----------------------------------------------------------

plot(history, metrics = "loss")

model %>% evaluate(x = list(t1_test,
                            fl_test),
                   y = mask_test,
                   batch_size = 1)
