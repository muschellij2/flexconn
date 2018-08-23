context("Trying to blur a mask")

set.seed(1)
dims = rep(80, 3)
n = prod(dims)
mask = rbinom(n, size = 1, prob = 0.85)
mask = array(mask, dim = dims)

test_that("Blurring gives list", {
  res = blur_mask(mask)
  # expect_equal(names(res), c("indices", "blurred_mask"))
  expect_equal(sum(res), 43549300.0063181)
  expect_equal(class(res), "array")
})
