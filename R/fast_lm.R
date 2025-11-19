#' Fast Linear Regression (pure R implementation)
#'
#' @description
#' `fast_lm()` fits an ordinary least squares (OLS) linear model
#' using explicit matrix formulas. The interface is similar to `lm()`,
#' and the estimated coefficients are identical to those from `lm()`
#' (up to numerical precision).
#'
#' You can use either a formula interface (`y ~ x1 + x2`) with a
#' data frame, or an `x`/`y` interface with numeric matrices/vectors.
#'
#' @param formula A model formula (e.g. `y ~ x1 + x2`). Optional if
#'   `x` and `y` are supplied.
#' @param data Optional data frame containing the variables in the formula.
#' @param x Optional numeric matrix (or object coercible to matrix) of predictors.
#' @param y Optional numeric response vector.
#' @param intercept Logical; for the `x`/`y` interface only, whether to
#'   include an intercept column (default `TRUE`).
#'
#' @return An object of class `"fastlm"` with components:
#' \itemize{
#'   \item \code{coefficients} Named numeric vector of regression coefficients.
#'   \item \code{fitted.values} Numeric vector of fitted values.
#'   \item \code{residuals} Numeric vector of residuals.
#'   \item \code{sigma2} Residual variance estimate.
#'   \item \code{call} The matched call.
#'   \item \code{formula} The model formula (if used).
#'   \item \code{terms} The terms object (if formula interface used).
#'   \item \code{model.matrix} The model matrix used for fitting.
#'   \item \code{response} The response vector.
#' }
#'
#' @examples
#' set.seed(1)
#' n  <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y  <- 1 + 2 * x1 - 3 * x2 + rnorm(n)
#' dat <- data.frame(y, x1, x2)
#'
#' fit_fast <- fast_lm(y ~ x1 + x2, data = dat)
#' fit_lm   <- lm(y ~ x1 + x2, data = dat)
#'
#' all.equal(coef(fit_fast), coef(fit_lm))
#'
#' @export
fast_lm <- function(formula, data, x, y, intercept = TRUE) {
  call <- match.call()

  # ------------------------------------------------------------------
  # x / y interface
  # ------------------------------------------------------------------
  if (!missing(x) && !missing(y)) {
    X <- as.matrix(x)
    y_vec <- as.numeric(y)

    if (intercept) {
      X <- cbind("(Intercept)" = 1, X)
    }

    # OLS: beta = (X'X)^(-1) X'y
    XtX <- crossprod(X)          # t(X) %*% X
    Xty <- crossprod(X, y_vec)   # t(X) %*% y
    beta <- solve(XtX, Xty)      # p x 1

    beta <- drop(beta)
    names(beta) <- colnames(X)

    fitted <- as.vector(X %*% beta)
    resid  <- y_vec - fitted

    n <- nrow(X)
    p <- ncol(X)
    sigma2 <- sum(resid^2) / (n - p)

    out <- list(
      coefficients  = beta,
      fitted.values = fitted,
      residuals     = resid,
      sigma2        = sigma2,
      call          = call,
      formula       = NULL,
      terms         = NULL,
      model.matrix  = X,
      response      = y_vec
    )
    class(out) <- "fastlm"
    return(out)
  }

  # ------------------------------------------------------------------
  # formula interface
  # ------------------------------------------------------------------
  if (missing(formula)) {
    stop("Provide either a formula or both x and y.", call. = FALSE)
  }

  mf <- model.frame(formula, data = data)
  terms <- attr(mf, "terms")
  y_vec <- model.response(mf)
  X <- model.matrix(terms, data = mf)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y_vec)
  beta <- solve(XtX, Xty)

  beta <- drop(beta)
  names(beta) <- colnames(X)

  fitted <- as.vector(X %*% beta)
  resid  <- y_vec - fitted

  n <- nrow(X)
  p <- ncol(X)
  sigma2 <- sum(resid^2) / (n - p)

  out <- list(
    coefficients  = beta,
    fitted.values = fitted,
    residuals     = resid,
    sigma2        = sigma2,
    call          = call,
    formula       = formula,
    terms         = terms,
    model.matrix  = X,
    response      = y_vec
  )
  class(out) <- "fastlm"
  out
}

#' @export
print.fastlm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
summary.fastlm <- function(object, ...) {
  X <- object$model.matrix
  y <- object$response
  n <- nrow(X)
  p <- ncol(X)

  coef <- object$coefficients
  resid <- object$residuals
  rss <- sum(resid^2)

  sigma2 <- rss / (n - p)
  XtX_inv <- solve(crossprod(X))

  se <- sqrt(diag(XtX_inv) * sigma2)
  tval <- coef / se
  pval <- 2 * pt(-abs(tval), df = n - p)

  coef_tab <- cbind(
    Estimate = coef,
    "Std. Error" = se,
    "t value" = tval,
    "Pr(>|t|)" = pval
  )

  ans <- list(
    call = object$call,
    coefficients = coef_tab,
    sigma = sqrt(sigma2),
    df = c(p, n - p)
  )
  class(ans) <- "summary.fastlm"
  ans
}

#' @export
print.summary.fastlm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat(
    "\nResidual standard error:",
    format(signif(x$sigma, 4)),
    "on", x$df[2], "degrees of freedom\n"
  )
  invisible(x)
}
