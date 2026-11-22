# S3 methods for classes "nbreg" and "zinbreg"
# ---------------------------------------------------------------------------

# --- shared helpers ---------------------------------------------------------

coef_table <- function(est, se) {
  z <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(tab) <- names(est)
  tab
}

print_call <- function(x) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
}

print_convergence <- function(x) {
  if (!isTRUE(x$converged))
    cat("\nWarning: the algorithm did not converge.\n")
}

# Offset for prediction on new data: offset() terms in the formula plus an
# `offset =` argument used in the original call (evaluated in `newdata`).
pred_offset <- function(object, tt, mf, newdata) {
  off <- rep(0, nrow(mf))
  o <- model.offset(mf)
  if (!is.null(o)) off <- off + o
  if (!is.null(object$call$offset))
    off <- off + eval(object$call$offset, newdata, environment(tt))
  off
}

# ===========================================================================
# nbreg
# ===========================================================================

#' Methods for `nbreg` objects
#'
#' Standard S3 methods for fits from [nbreg()].
#'
#' @param x,object An object of class `"nbreg"` (or `"summary.nbreg"` for the
#'   print method of a summary).
#' @param digits Number of significant digits to print.
#' @param type For `residuals`: `"pearson"` (default) or `"response"`. For
#'   `predict`: `"response"` (the mean) or `"link"` (the linear predictor).
#' @param newdata Optional data frame for `predict`.
#' @param ... Further arguments (ignored).
#'
#' Wald confidence intervals are available through the default `confint()`
#' method, which uses `coef()` and `vcov()`.
#'
#' @return `coef`, `vcov`, `fitted`, `residuals` and `predict` return numeric
#'   vectors or matrices; `logLik` an object of class `"logLik"` with a `df`
#'   attribute so that `AIC()` and `BIC()` work; `nobs` an integer; `summary` an
#'   object of class `"summary.nbreg"`.
#' @name nbreg-methods
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rnbinom(200, mu = exp(0.5 + 0.4 * d$x), size = 2)
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' summary(fit)
#' AIC(fit)
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print_call(x)
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$se.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 2L),
      " on ", x$df, " Df\n", sep = "")
  print_convergence(x)
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) object$nobs

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r <- object$y - mu
  if (type == "pearson") r <- r / sqrt(mu + mu^2 / object$theta)
  r
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata))
    return(switch(type, response = object$fitted.values,
                  link = object$linear.predictors))
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
  eta <- drop(X %*% object$coefficients) + pred_offset(object, tt, mf, newdata)
  names(eta) <- rownames(X)
  switch(type, response = exp(eta), link = eta)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  structure(
    list(call = object$call,
         coefficients = coef_table(object$coefficients, se),
         theta = object$theta, se.theta = object$se.theta,
         loglik = object$loglik, df = object$df,
         AIC = -2 * object$loglik + 2 * object$df,
         iter = object$iter, converged = object$converged,
         nobs = object$nobs),
    class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  print_call(x)
  cat("Coefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.legend = TRUE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$se.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 2L),
      " on ", x$df, " Df;  AIC: ", format(x$AIC, digits = digits + 2L), "\n",
      sep = "")
  cat("Number of observations: ", x$nobs, ";  outer iterations: ", x$iter,
      "\n", sep = "")
  print_convergence(x)
  invisible(x)
}

# ===========================================================================
# zinbreg
# ===========================================================================

#' Methods for `zinbreg` objects
#'
#' Standard S3 methods for fits from [zinbreg()].
#'
#' `coef` and `vcov` return the coefficients (and their covariance matrix) of
#' the count part, the zero part, or both (`"full"`, with names `count_*` and
#' `zero_*`, as in `pscl::zeroinfl()`). The covariance of the full fit
#' including `Log(theta)` is stored in the fitted object as `vcov`.
#'
#' For `predict`, `type = "response"` gives the marginal mean
#' \eqn{(1-\pi)\mu}, `"count"` the NB mean \eqn{\mu}, and `"zero"` the
#' structural-zero probability \eqn{\pi}. Pearson residuals use the ZINB
#' variance \eqn{(1-\pi)\mu\{1 + \mu(\pi + 1/\theta)\}}.
#'
#' @param x,object An object of class `"zinbreg"` (or `"summary.zinbreg"` for
#'   the print method of a summary).
#' @param model Which part of the model: `"full"` (default), `"count"` or
#'   `"zero"`.
#' @param digits Number of significant digits to print.
#' @param type For `residuals`: `"pearson"` (default) or `"response"`. For
#'   `predict`: `"response"`, `"count"` or `"zero"`.
#' @param newdata Optional data frame for `predict`.
#' @param ... Further arguments (ignored).
#'
#' @return See [nbreg-methods] for the common methods. `summary` returns an
#'   object of class `"summary.zinbreg"` with coefficient tables for the count
#'   and zero parts.
#' @name zinbreg-methods
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' pz <- plogis(-1 + 0.8 * d$z)
#' d$y <- ifelse(runif(n) < pz, 0, rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2))
#' fit <- zinbreg(y ~ x | z, data = d)
#' summary(fit)
#' coef(fit, model = "zero")
#' head(predict(fit, type = "zero"))
#' predict(fit, newdata = data.frame(x = 0, z = 0), type = "response")
NULL

zinb_block <- function(object, model) {
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  switch(model, full = seq_len(p + q), count = seq_len(p), zero = p + seq_len(q))
}

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print_call(x)
  cat("Count model coefficients (negative binomial with log link):\n")
  print.default(format(x$coefficients$count, digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$se.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 2L),
      " on ", x$df, " Df\n", sep = "")
  print_convergence(x)
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  switch(model,
         full = c(stats::setNames(object$coefficients$count,
                                  paste0("count_", names(object$coefficients$count))),
                  stats::setNames(object$coefficients$zero,
                                  paste0("zero_", names(object$coefficients$zero)))),
         count = object$coefficients$count,
         zero = object$coefficients$zero)
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  idx <- zinb_block(object, model)
  object$vcov[idx, idx, drop = FALSE]
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) object$nobs

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  m <- object$fitted.values
  r <- object$y - m
  if (type == "pearson") {
    mu <- object$mu
    r <- r / sqrt(m * (1 + mu * (object$pi + 1 / object$theta)))
  }
  r
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata))
    return(switch(type, response = object$fitted.values, count = object$mu,
                  zero = object$pi))
  ttX <- stats::delete.response(object$terms$count)
  mfX <- stats::model.frame(ttX, newdata, na.action = stats::na.pass,
                            xlev = object$xlevels$count)
  X <- model.matrix(ttX, mfX, contrasts.arg = object$contrasts$count)
  mu <- exp(drop(X %*% object$coefficients$count) +
              pred_offset(object, ttX, mfX, newdata))
  ttZ <- object$terms$zero
  mfZ <- stats::model.frame(ttZ, newdata, na.action = stats::na.pass,
                            xlev = object$xlevels$zero)
  Z <- model.matrix(ttZ, mfZ, contrasts.arg = object$contrasts$zero)
  pi_ <- stats::plogis(drop(Z %*% object$coefficients$zero))
  out <- switch(type, response = (1 - pi_) * mu, count = mu, zero = pi_)
  names(out) <- rownames(X)
  out
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  structure(
    list(call = object$call,
         count = coef_table(object$coefficients$count, se[seq_len(p)]),
         zero = coef_table(object$coefficients$zero, se[p + seq_len(q)]),
         theta = object$theta, se.theta = object$se.theta,
         logtheta = object$logtheta, se.logtheta = object$se.logtheta,
         loglik = object$loglik, df = object$df,
         AIC = -2 * object$loglik + 2 * object$df,
         iter = object$iter, converged = object$converged,
         nobs = object$nobs),
    class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  print_call(x)
  cat("Count model coefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$count, digits = digits, signif.legend = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$zero, digits = digits, signif.legend = TRUE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$se.theta, digits = digits), ")\n", sep = "")
  cat("Log(theta): ", format(x$logtheta, digits = digits),
      "  (SE ", format(x$se.logtheta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 2L),
      " on ", x$df, " Df;  AIC: ", format(x$AIC, digits = digits + 2L), "\n",
      sep = "")
  cat("Number of observations: ", x$nobs, ";  EM iterations: ", x$iter, "\n",
      sep = "")
  print_convergence(x)
  invisible(x)
}
