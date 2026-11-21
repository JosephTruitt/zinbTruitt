#' Control parameters for negative binomial regression
#'
#' Tolerances and iteration limits for [nbreg()].
#'
#' @param epsilon Convergence tolerance. The outer (beta, theta) loop stops when
#'   the relative change in the log-likelihood is below `epsilon`.
#' @param maxit Maximum number of outer iterations (IRLS for beta alternated
#'   with Newton-Raphson for theta).
#' @param irls_maxit Maximum number of IRLS iterations in each beta update.
#' @param theta_maxit Maximum number of Newton-Raphson iterations in each theta
#'   update.
#' @param theta_max Upper bound for theta. If theta reaches this value the model
#'   is essentially Poisson and a warning is issued.
#' @param trace Logical; print the log-likelihood and theta at each outer
#'   iteration?
#'
#' @return A list of control settings.
#' @export
#' @examples
#' nbreg_control(epsilon = 1e-12, trace = TRUE)
nbreg_control <- function(epsilon = 1e-13, maxit = 100L, irls_maxit = 50L,
                          theta_maxit = 50L, theta_max = 1e6, trace = FALSE) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0)
    stop("'epsilon' must be a single positive number")
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1)
    stop("'maxit' must be a positive integer")
  if (!is.numeric(irls_maxit) || irls_maxit < 1)
    stop("'irls_maxit' must be a positive integer")
  if (!is.numeric(theta_maxit) || theta_maxit < 1)
    stop("'theta_maxit' must be a positive integer")
  if (!is.numeric(theta_max) || theta_max <= 1)
    stop("'theta_max' must be a number greater than 1")
  list(epsilon = epsilon, maxit = as.integer(maxit),
       irls_maxit = as.integer(irls_maxit),
       theta_maxit = as.integer(theta_maxit),
       theta_max = theta_max, trace = isTRUE(trace))
}

#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' Fits the NB2 regression model with log link,
#' \eqn{Var(Y) = \mu + \mu^2/\theta}, without any general-purpose optimizer.
#' The coefficients \eqn{\beta} are updated by iteratively reweighted least
#' squares (Fisher scoring) with \eqn{\theta} held fixed, and \eqn{\theta} is
#' updated by Newton-Raphson on \eqn{\log\theta} with \eqn{\beta} held fixed.
#' The two updates are alternated until the log-likelihood converges.
#'
#' The parameterization matches `MASS::glm.nb()`. Standard errors of
#' \eqn{\hat\beta} are \eqn{(X^\top \hat W X)^{-1}}, treating \eqn{\theta} as
#' known, and the standard error of \eqn{\hat\theta} is
#' \eqn{\{-H(\hat\theta)\}^{-1/2}}.
#'
#' @param formula A model formula, e.g. `y ~ x1 + x2`. An `offset()` term is
#'   allowed.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional numeric offset vector (e.g. `log(exposure)`),
#'   evaluated in `data`.
#' @param control A list of control settings from [nbreg_control()].
#'
#' @return An object of class `"nbreg"`, a list with components including
#'   `coefficients`, `vcov`, `theta`, `se.theta`, `loglik`, `df`, `iter`,
#'   `converged`, `fitted.values`, `linear.predictors`, `y`, `offset`, `nobs`
#'   and `call`.
#' @export
#' @importFrom stats model.frame model.matrix model.response model.offset
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rpois(200, exp(0.5 + 0.4 * d$x) * rgamma(200, shape = 2, rate = 2))
#' fit <- nbreg(y ~ x, data = d)
#' fit$coefficients
#' fit$theta
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf, "numeric")
  if (is.null(y)) stop("a response variable is required in 'formula'")
  if (!is.numeric(y)) stop("the response must be numeric")
  if (length(y) == 0L) stop("no observations left after removing missing values")
  if (any(y < 0)) stop("the response must be non-negative counts")
  if (any(abs(y - round(y)) > 1e-8)) stop("the response must be integer-valued counts")

  X <- model.matrix(mt, mf)
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, length(y))
  if (qr(X)$rank < ncol(X))
    stop("the model matrix is rank deficient (collinear covariates)")
  if (nrow(X) <= ncol(X))
    stop("need more observations than coefficients")

  fit <- nb_fit(X, y, off, a = rep(1, length(y)), control = control)

  fit$call <- cl
  fit$terms <- mt
  fit$xlevels <- stats::.getXlevels(mt, mf)
  fit$contrasts <- attr(X, "contrasts")
  fit$na.action <- attr(mf, "na.action")
  class(fit) <- "nbreg"
  fit
}

# ---------------------------------------------------------------------------
# Internal workhorse. Works on matrices, accepts prior weights `a` so that it
# can be reused for the count part of the ZINB EM algorithm (a = 1 - tau).
# ---------------------------------------------------------------------------
nb_fit <- function(X, y, off, a, control, beta_start = NULL,
                   theta_start = NULL) {
  n <- length(y)
  p <- ncol(X)

  # Start: Poisson IRLS (theta = Inf)
  if (is.null(beta_start)) {
    mu0 <- y + 0.1
    beta <- nb_irls(X, y, off, a, theta = Inf,
                    beta = nb_wls(X, log(mu0) - off, a * mu0),
                    control = control)$beta
  } else {
    beta <- beta_start
  }
  mu <- nb_mu(X, beta, off)

  if (is.null(theta_start)) {
    # moment estimate from Var = mu + mu^2/theta
    den <- sum(a * ((y - mu)^2 - mu))
    theta <- if (den > 0) sum(a * mu^2) / den else 1
    theta <- min(max(theta, 0.05), 100)
  } else {
    theta <- theta_start
  }

  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  theta_capped <- FALSE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    beta <- nb_irls(X, y, off, a, theta, beta, control)$beta
    mu <- nb_mu(X, beta, off)
    th <- nb_theta_nr(y, mu, a, theta, control)
    theta <- th$theta
    theta_capped <- th$capped
    ll_new <- nb_loglik(y, mu, theta, a)
    if (control$trace)
      cat(sprintf("iter %3d  loglik = %.10f  theta = %.6f\n", iter, ll_new, theta))
    done <- abs(ll_new - ll) < control$epsilon * (abs(ll_new) + 0.1)
    ll <- ll_new
    if (done) { converged <- TRUE; break }
  }
  if (!converged)
    warning("nbreg did not converge in ", control$maxit,
            " iterations; consider increasing 'maxit'", call. = FALSE)
  if (theta_capped)
    warning("theta reached the upper bound (", format(control$theta_max),
            "); the data show no overdispersion and the Poisson model may be ",
            "more appropriate", call. = FALSE)

  eta <- drop(X %*% beta) + off
  mu <- exp(eta)
  w <- nb_weights(a, mu, theta)
  xtwx <- crossprod(X * sqrt(w))
  vc <- solve(xtwx)
  dimnames(vc) <- list(colnames(X), colnames(X))

  hth <- nb_hess_theta(y, mu, a, theta)
  se_theta <- if (is.finite(hth) && hth < 0) sqrt(-1 / hth) else NA_real_
  se_logtheta <- if (is.finite(se_theta)) se_theta / theta else NA_real_

  names(beta) <- colnames(X)
  list(coefficients = beta, vcov = vc, theta = theta, se.theta = se_theta,
       logtheta = log(theta), se.logtheta = se_logtheta,
       loglik = ll, df = p + 1L, iter = iter, converged = converged,
       fitted.values = mu, linear.predictors = eta, y = y, offset = off,
       weights = a, nobs = sum(a > 0), n = n, control = control)
}

nb_mu <- function(X, beta, off) {
  eta <- drop(X %*% beta) + off
  exp(pmin(eta, 700))
}

nb_loglik <- function(y, mu, theta, a) {
  if (is.infinite(theta)) return(sum(a * stats::dpois(y, mu, log = TRUE)))
  sum(a * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
}

# IRLS working weights: a * (dmu/deta)^2 / V(mu) = a * theta * mu / (theta + mu)
nb_weights <- function(a, mu, theta) {
  if (is.infinite(theta)) return(a * mu)
  a * theta * mu / (theta + mu)
}

# Weighted least squares via QR: solves min sum w (z - X b)^2
nb_wls <- function(X, z, w) {
  sw <- sqrt(w)
  qr.coef(qr(X * sw), z * sw)
}

# IRLS for beta with theta fixed (Fisher scoring) with step halving safeguard
nb_irls <- function(X, y, off, a, theta, beta, control) {
  mu <- nb_mu(X, beta, off)
  ll <- nb_loglik(y, mu, theta, a)
  for (it in seq_len(control$irls_maxit)) {
    w <- nb_weights(a, mu, theta)
    z <- log(mu) - off + (y - mu) / mu
    beta_new <- nb_wls(X, z, w)
    mu_new <- nb_mu(X, beta_new, off)
    ll_new <- nb_loglik(y, mu_new, theta, a)
    half <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-8 * abs(ll)) && half < 20L) {
      beta_new <- (beta + beta_new) / 2
      mu_new <- nb_mu(X, beta_new, off)
      ll_new <- nb_loglik(y, mu_new, theta, a)
      half <- half + 1L
    }
    done <- abs(ll_new - ll) < control$epsilon * (abs(ll_new) + 0.1)
    beta <- beta_new
    mu <- mu_new
    ll <- ll_new
    if (done) break
  }
  list(beta = beta, mu = mu, loglik = ll)
}

# Score and second derivative of the log-likelihood in theta
nb_score_theta <- function(y, mu, a, theta) {
  sum(a * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
             log(theta + mu) - (y + theta) / (theta + mu)))
}

nb_hess_theta <- function(y, mu, a, theta) {
  sum(a * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
             2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

# Newton-Raphson for phi = log(theta) with beta (mu) fixed
nb_theta_nr <- function(y, mu, a, theta, control) {
  phi <- log(theta)
  phi_max <- log(control$theta_max)
  capped <- FALSE
  ll <- nb_loglik(y, mu, exp(phi), a)
  for (it in seq_len(control$theta_maxit)) {
    th <- exp(phi)
    U <- nb_score_theta(y, mu, a, th)
    H <- nb_hess_theta(y, mu, a, th)
    Uphi <- th * U
    Hphi <- th^2 * H + th * U
    # Newton step; if curvature is not negative, take a bounded ascent step
    step <- if (is.finite(Hphi) && Hphi < 0) -Uphi / Hphi else sign(Uphi) * min(abs(Uphi), 1)
    step <- max(min(step, 5), -5)
    phi_new <- min(phi + step, phi_max)
    ll_new <- nb_loglik(y, mu, exp(phi_new), a)
    half <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-10 * abs(ll)) && half < 30L) {
      step <- step / 2
      phi_new <- min(phi + step, phi_max)
      ll_new <- nb_loglik(y, mu, exp(phi_new), a)
      half <- half + 1L
    }
    delta <- abs(phi_new - phi)
    phi <- phi_new
    ll <- ll_new
    if (delta < 1e-10) break
  }
  if (phi >= phi_max - 1e-8) capped <- TRUE
  list(theta = exp(phi), capped = capped)
}
