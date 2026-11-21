#' Control parameters for zero-inflated negative binomial regression
#'
#' Tolerances and iteration limits for [zinbreg()].
#'
#' @param epsilon Convergence tolerance for the EM algorithm: stop when the
#'   relative change in the observed-data log-likelihood is below `epsilon`.
#' @param maxit Maximum number of EM iterations.
#' @param irls_maxit Maximum number of IRLS iterations in the zero-part M-step.
#' @param theta_max Upper bound for theta (Poisson limit); a warning is issued
#'   if it is reached.
#' @param nb_epsilon Convergence tolerance used by the count-part (NB) M-step.
#' @param trace Logical; print the log-likelihood at each EM iteration?
#'
#' @return A list of control settings.
#' @export
#' @examples
#' zinbreg_control(epsilon = 1e-10, maxit = 500)
zinbreg_control <- function(epsilon = 1e-12, maxit = 2000L, irls_maxit = 50L,
                            theta_max = 1e6, nb_epsilon = 1e-10,
                            trace = FALSE) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0)
    stop("'epsilon' must be a single positive number")
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1)
    stop("'maxit' must be a positive integer")
  if (!is.numeric(irls_maxit) || irls_maxit < 1)
    stop("'irls_maxit' must be a positive integer")
  if (!is.numeric(theta_max) || theta_max <= 1)
    stop("'theta_max' must be a number greater than 1")
  if (!is.numeric(nb_epsilon) || nb_epsilon <= 0)
    stop("'nb_epsilon' must be a positive number")
  list(epsilon = epsilon, maxit = as.integer(maxit),
       irls_maxit = as.integer(irls_maxit), theta_max = theta_max,
       nb_epsilon = nb_epsilon, trace = isTRUE(trace))
}

#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' Fits the ZINB model: with probability \eqn{\pi_i} an observation is a
#' structural zero, otherwise it follows an NB(\eqn{\mu_i}, \eqn{\theta})
#' distribution, with \eqn{\log \mu_i = x_i^\top\beta + o_i} and
#' \eqn{\mathrm{logit}\, \pi_i = z_i^\top\gamma}. The likelihood is maximized
#' by the EM algorithm. The E-step computes the posterior probability
#' \eqn{\tau_i} of being a structural zero. The zero-part M-step is a logistic
#' regression with fractional response \eqn{\tau_i} fitted by IRLS, and the
#' count-part M-step is the NB fit of [nbreg()] with prior weights
#' \eqn{1 - \tau_i}. No general-purpose optimizer is used.
#'
#' Standard errors come from the inverse of the observed information matrix of
#' the observed-data log-likelihood with respect to
#' \eqn{(\beta, \gamma, \log\theta)}, computed from the analytic Hessian.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before the bar, zero-model covariates after it. If there is no
#'   `|`, the same covariates are used in both parts. An `offset()` term in the
#'   count part is allowed.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional numeric offset for the count model, evaluated in
#'   `data`.
#' @param control A list of control settings from [zinbreg_control()].
#'
#' @return An object of class `"zinbreg"`, a list with components including
#'   `coefficients` (a list with elements `count` and `zero`), `vcov` (for
#'   `count_*`, `zero_*` and `Log(theta)`), `theta`, `se.theta`, `logtheta`,
#'   `se.logtheta`, `loglik`, `df`, `iter`, `converged`, `loglik.path`,
#'   `fitted.values`, `mu`, `pi`, `tau`, `y` and `call`.
#' @export
#' @importFrom stats model.frame model.matrix model.response model.offset
#' @importFrom stats terms delete.response
#' @importFrom Formula as.Formula
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' pz <- plogis(-1 + 0.8 * d$z)
#' d$y <- ifelse(runif(n) < pz, 0, rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2))
#' fit <- zinbreg(y ~ x | z, data = d)
#' fit$coefficients
#' fit$theta
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  formula <- Formula::as.Formula(formula)
  one_part <- length(formula)[2L] < 2L
  if (one_part) {
    formula <- Formula::as.Formula(stats::formula(formula),
                                   stats::formula(formula, lhs = 0L, rhs = 1L))
  } else if (length(formula)[2L] > 2L) {
    stop("'formula' may have at most two parts separated by '|'")
  }
  mtZ_full <- stats::terms(formula, data = data, rhs = 2L)
  if (!one_part && !is.null(attr(mtZ_full, "offset")))
    stop("offsets are only supported in the count part of the model")

  mf$formula <- formula
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  y <- model.response(mf, "numeric")
  if (is.null(y)) stop("a response variable is required in 'formula'")
  if (!is.numeric(y)) stop("the response must be numeric")
  if (length(y) == 0L) stop("no observations left after removing missing values")
  if (any(y < 0)) stop("the response must be non-negative counts")
  if (any(abs(y - round(y)) > 1e-8)) stop("the response must be integer-valued counts")
  if (!any(y == 0)) stop("the response contains no zeros; zero inflation cannot be estimated")

  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, length(y))
  if (qr(X)$rank < ncol(X))
    stop("the count-model matrix is rank deficient (collinear covariates)")
  if (qr(Z)$rank < ncol(Z))
    stop("the zero-model matrix is rank deficient (collinear covariates)")
  if (nrow(X) <= ncol(X) + ncol(Z) + 1L)
    stop("need more observations than parameters")

  fit <- zinb_fit(X, Z, y, off, control)

  fit$call <- cl
  fit$formula <- formula
  fit$terms <- list(count = mtX, zero = mtZ, full = attr(mf, "terms"))
  fit$xlevels <- list(count = stats::.getXlevels(mtX, mf),
                      zero = stats::.getXlevels(mtZ, mf))
  fit$contrasts <- list(count = attr(X, "contrasts"), zero = attr(Z, "contrasts"))
  fit$na.action <- attr(mf, "na.action")
  class(fit) <- "zinbreg"
  fit
}

# ---------------------------------------------------------------------------
# EM algorithm (matrix interface)
# ---------------------------------------------------------------------------
zinb_fit <- function(X, Z, y, off, control) {
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  nbctl <- nbreg_control(epsilon = control$nb_epsilon,
                         theta_max = control$theta_max)

  # Starting values: NB fit for (beta, theta); logistic IRLS of I(y = 0) on Z
  start <- suppressWarnings(nb_fit(X, y, off, a = rep(1, n), control = nbctl))
  beta <- start$coefficients
  theta <- start$theta
  gamma <- zi_irls(Z, as.numeric(y == 0), rep(0, q), control)
  mu <- nb_mu(X, beta, off)
  zeta <- drop(Z %*% gamma)

  ll <- zinb_loglik(y, mu, zeta, theta)
  path <- ll
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    # E-step: posterior probability of a structural zero
    logf0 <- stats::dnbinom(0, size = theta, mu = mu, log = TRUE)
    tau <- ifelse(y == 0, stats::plogis(zeta - logf0), 0)

    # M-step, zero part: logistic IRLS with fractional response tau
    gamma <- zi_irls(Z, tau, gamma, control)
    zeta <- drop(Z %*% gamma)

    # M-step, count part: NB fit with prior weights 1 - tau
    nbf <- suppressWarnings(
      nb_fit(X, y, off, a = 1 - tau, control = nbctl,
             beta_start = beta, theta_start = theta))
    beta <- nbf$coefficients
    theta <- nbf$theta
    mu <- nb_mu(X, beta, off)

    ll_new <- zinb_loglik(y, mu, zeta, theta)
    path <- c(path, ll_new)
    if (control$trace)
      cat(sprintf("iter %4d  loglik = %.10f  theta = %.6f\n", iter, ll_new, theta))
    done <- abs(ll_new - ll) < control$epsilon * (abs(ll_new) + 0.1)
    ll <- ll_new
    if (done) { converged <- TRUE; break }
  }
  if (!converged)
    warning("zinbreg did not converge in ", control$maxit,
            " EM iterations; consider increasing 'maxit'", call. = FALSE)
  if (theta >= control$theta_max * 0.999)
    warning("theta reached the upper bound (", format(control$theta_max),
            "); the count part shows no overdispersion and a zero-inflated ",
            "Poisson model may be more appropriate", call. = FALSE)
  if (any(diff(path) < -1e-6 * (abs(ll) + 0.1)))
    warning("the log-likelihood decreased during EM; results may be unreliable",
            call. = FALSE)

  # Posterior probabilities and fitted quantities at the final estimates
  logf0 <- stats::dnbinom(0, size = theta, mu = mu, log = TRUE)
  tau <- ifelse(y == 0, stats::plogis(zeta - logf0), 0)
  pi_hat <- stats::plogis(zeta)

  # Standard errors: inverse observed information, analytic Hessian
  H <- zinb_hessian(X, Z, y, off, beta, gamma, theta)
  nm <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
          "Log(theta)")
  dimnames(H) <- list(nm, nm)
  vc <- tryCatch(solve(-H), error = function(e) {
    warning("the observed information matrix is singular; standard errors ",
            "are not available", call. = FALSE)
    matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
  })
  se_logtheta <- sqrt(vc[p + q + 1L, p + q + 1L])

  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)
  list(coefficients = list(count = beta, zero = gamma), vcov = vc,
       theta = theta, se.theta = theta * se_logtheta,
       logtheta = log(theta), se.logtheta = se_logtheta,
       loglik = ll, df = p + q + 1L, iter = iter, converged = converged,
       loglik.path = path, fitted.values = (1 - pi_hat) * mu,
       mu = mu, pi = pi_hat, tau = tau, y = y, offset = off,
       nobs = n, control = control)
}

# Observed-data log-likelihood, stable via log-sum-exp
zinb_loglik <- function(y, mu, zeta, theta) {
  lpi <- stats::plogis(zeta, log.p = TRUE)                     # log pi
  l1p <- stats::plogis(zeta, lower.tail = FALSE, log.p = TRUE) # log(1 - pi)
  lnb <- stats::dnbinom(y, size = theta, mu = mu, log = TRUE)  # log f_NB
  a <- lpi
  b <- l1p + lnb
  m <- pmax(a, b)
  l0 <- m + log(exp(a - m) + exp(b - m))                       # y = 0
  sum(ifelse(y == 0, l0, b))
}

# Logistic regression with fractional response tau by IRLS:
# w = pi (1 - pi), working response z = eta + (tau - pi) / (pi (1 - pi))
zi_irls <- function(Z, tau, gamma, control) {
  ll_fun <- function(eta)
    sum(tau * stats::plogis(eta, log.p = TRUE) +
          (1 - tau) * stats::plogis(eta, lower.tail = FALSE, log.p = TRUE))
  eta <- drop(Z %*% gamma)
  ll <- ll_fun(eta)
  for (it in seq_len(control$irls_maxit)) {
    pi_ <- stats::plogis(eta)
    v <- pmax(pi_ * (1 - pi_), 1e-12)
    zz <- eta + (tau - pi_) / v
    gamma_new <- nb_wls(Z, zz, v)
    eta_new <- drop(Z %*% gamma_new)
    ll_new <- ll_fun(eta_new)
    half <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-8 * abs(ll)) && half < 20L) {
      gamma_new <- (gamma + gamma_new) / 2
      eta_new <- drop(Z %*% gamma_new)
      ll_new <- ll_fun(eta_new)
      half <- half + 1L
    }
    done <- abs(ll_new - ll) < control$nb_epsilon * (abs(ll_new) + 0.1)
    gamma <- gamma_new
    eta <- eta_new
    ll <- ll_new
    if (done) break
  }
  gamma
}

# Analytic Hessian of the observed-data log-likelihood with respect to
# (beta, gamma, phi = log theta)
zinb_hessian <- function(X, Z, y, off, beta, gamma, theta) {
  p <- ncol(X)
  q <- ncol(Z)
  mu <- nb_mu(X, beta, off)
  zeta <- drop(Z %*% gamma)
  pi_ <- stats::plogis(zeta)
  th <- theta
  s <- th + mu

  # Derivatives of log f_NB(y; mu, theta) w.r.t. eta = log mu and phi = log theta.
  # For y = 0 these are also the derivatives of log f0 = theta * log(theta/s).
  g <- digamma(y + th) - digamma(th) + log(th) + 1 - log(s) - (y + th) / s
  gp <- trigamma(y + th) - trigamma(th) + 1 / th - 2 / s + (y + th) / s^2
  e1 <- th * (y - mu) / s
  ee <- -th * (y + th) * mu / s^2
  ep <- th * (y - mu) * mu / s^2
  p1 <- th * g
  pp <- th * (g + th * gp)

  hzz <- hze <- hzp <- hee <- hep <- hpp <- numeric(length(y))

  pos <- y > 0
  hzz[pos] <- -(pi_ * (1 - pi_))[pos]
  hee[pos] <- ee[pos]
  hep[pos] <- ep[pos]
  hpp[pos] <- pp[pos]

  z0 <- !pos
  f0 <- exp(th * (log(th) - log(s)))
  P <- pi_ + (1 - pi_) * f0
  pq <- pi_ * (1 - pi_)
  omp <- 1 - pi_
  Pz <- pq * (1 - f0)
  Pzz <- pq * (1 - 2 * pi_) * (1 - f0)
  Pe <- omp * f0 * e1
  Pee <- omp * f0 * (e1^2 + ee)
  Pp <- omp * f0 * p1
  Ppp <- omp * f0 * (p1^2 + pp)
  Pze <- -pq * f0 * e1
  Pzp <- -pq * f0 * p1
  Pep <- omp * f0 * (e1 * p1 + ep)

  hzz[z0] <- (Pzz / P - Pz^2 / P^2)[z0]
  hze[z0] <- (Pze / P - Pz * Pe / P^2)[z0]
  hzp[z0] <- (Pzp / P - Pz * Pp / P^2)[z0]
  hee[z0] <- (Pee / P - Pe^2 / P^2)[z0]
  hep[z0] <- (Pep / P - Pe * Pp / P^2)[z0]
  hpp[z0] <- (Ppp / P - Pp^2 / P^2)[z0]

  H <- matrix(0, p + q + 1L, p + q + 1L)
  ib <- seq_len(p)
  ig <- p + seq_len(q)
  ip <- p + q + 1L
  H[ib, ib] <- crossprod(X, X * hee)
  H[ig, ig] <- crossprod(Z, Z * hzz)
  H[ib, ig] <- crossprod(X, Z * hze)
  H[ig, ib] <- t(H[ib, ig])
  H[ib, ip] <- crossprod(X, hep)
  H[ip, ib] <- H[ib, ip]
  H[ig, ip] <- crossprod(Z, hzp)
  H[ip, ig] <- H[ig, ip]
  H[ip, ip] <- sum(hpp)
  H
}
