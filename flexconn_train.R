# reticulate::use_python(
#   "/software/apps/anaconda/5.2/python/3.6/bin/python")

library(keras)
library(neurobase)
library(flexconn)
library(RNifti)


# Configuration -----------------------

atlas_dir <- "lesion_challenge/atlas_with_mask1"
n_atlas <- 21

ndim <- 2
psize <- 35
patchsize <- rep(psize, ndim)
padsize = patchsize_to_padsize(patchsize = patchsize)

model_dir <- "saved_models"
tmp_data_dir <- "tmp_data"

# model configuration
batch_size <- 512


# Read data -------------------------------
t1 = list.files("_T1", path = atlas_dir,
  full.names = TRUE)
flair = list.files("_FL", path = atlas_dir,
  full.names = TRUE)
mask = list.files("_mask", path = atlas_dir,
  full.names = TRUE)

t1 = t1[seq(n_atlas)]
flair = flair[seq(n_atlas)]
mask = mask[seq(n_atlas)]

outfile = file.path(tmp_data_dir,
  paste0("patches_", n_atlas, ".rds"))
train = patches_list(
  t1 = t1,
  flair = flair,
  mask = mask,
  outfile = outfile
  )
num_patches = nrow(train$mask)


# Optional train-test split --------------------------
train_indx <- sample(1:num_patches,
  floor(num_patches * 0.7))


test = lapply(train, function(x) {
  x[-train_indx, , , , drop = FALSE]
})

train = lapply(train, function(x) {
  x[train_indx, , , , drop = FALSE]
})

# Train model --------------------------------------

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
      filepath = file.path(model_dir, 
        "weights.{epoch:02d}-{val_loss:.2f}.hdf5")
    ),
    callback_early_stopping(patience = 1)
  )
)


saveRDS(history, 
  file.path(tmp_data_dir, "history.rds"))
model %>%
  save_model_hdf5(
    file.path(model_dir,
      "flexconn_final.hdf5"))


# Evaluate model -----------------------

png("model_history.png")
plot(history, metrics = "loss")
dev.off()

model %>% 
evaluate(x = unname(test[ c("t1", "fl")]),
         y = test$mask,
         batch_size = 1)
