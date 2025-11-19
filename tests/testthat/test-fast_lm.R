# tests/testthat/test-fast_lm.R

test_that("fast_lm (formula interface) matches lm coefficients", {
  set.seed(1)
  n  <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 2 * x1 - 3 * x2 + rnorm(n, sd = 0.5)
  dat <- data.frame(y, x1, x2)

  fit_fast <- fast_lm(y ~ x1 + x2, data = dat)
  fit_lm   <- lm(y ~ x1 + x2, data = dat)

  expect_true(
    isTRUE(all.equal(coef(fit_fast), coef(fit_lm)))
  )
})

test_that("fast_lm (x/y interface) matches formula interface", {
  set.seed(2)
  n  <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- -0.5 + 0.8 * x1 + 1.2 * x2 + rnorm(n)
  X  <- cbind(x1, x2)
  dat <- data.frame(y, x1, x2)

  fit_formula <- fast_lm(y ~ x1 + x2, data = dat)
  fit_xy      <- fast_lm(x = X, y = y, intercept = TRUE)

  expect_true(
    isTRUE(all.equal(coef(fit_formula), coef(fit_xy)))
  )
})

test_that("fast_lm returns object with expected structure", {
  set.seed(3)
  x1 <- rnorm(50)
  x2 <- rnorm(50)
  y  <- x1 - x2 + rnorm(50)
  dat <- data.frame(y, x1, x2)

  fit <- fast_lm(y ~ x1 + x2, data = dat)

  expect_s3_class(fit, "fastlm")

  expect_true(all(c(
    "coefficients", "fitted.values", "residuals",
    "sigma2", "call", "model.matrix", "response"
  ) %in% names(fit)))

  n <- nrow(dat)
  p <- length(coef(fit))

  expect_length(fit$fitted.values, n)
  expect_length(fit$residuals, n)
  expect_equal(dim(fit$model.matrix), c(n, p))
})

test_that("summary.fastlm works and returns sensible output", {
  set.seed(4)
  x1 <- rnorm(80)
  x2 <- rnorm(80)
  y  <- 2 + 1.5 * x1 + rnorm(80)
  dat <- data.frame(y, x1, x2)

  fit <- fast_lm(y ~ x1 + x2, data = dat)
  s   <- summary(fit)

  expect_s3_class(s, "summary.fastlm")
  expect_true(all(c("call", "coefficients", "sigma", "df") %in% names(s)))

  expect_equal(rownames(s$coefficients), names(coef(fit)))
})
