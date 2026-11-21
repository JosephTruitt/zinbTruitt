zi_ref <- function(formula, data) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12, maxit = 10000))
}

# Table 4 tolerances for ZINB vs zeroinfl(): all 1e-3
expect_matches_zeroinfl <- function(fit, ref, tol = 1e-3) {
  s <- sqrt(diag(vcov(ref)))
  est <- c(fit$coefficients$count, fit$coefficients$zero)
  se <- sqrt(diag(fit$vcov))[seq_along(est)]
  expect_true(fit$converged)
  expect_lt(max(abs(est - c(coef(ref, "count"), coef(ref, "zero")))), tol)
  expect_lt(max(abs(se - s)), tol)
  expect_lt(abs(fit$logtheta - log(ref$theta)), tol)
  expect_lt(abs(fit$se.logtheta - ref$SE.logtheta), tol)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), tol)
}

test_that("zinbreg matches pscl::zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  expect_matches_zeroinfl(zinbreg(f, data = bioChemists), zi_ref(f, bioChemists))
})

test_that("zinbreg matches pscl::zeroinfl on NMES1988", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  expect_matches_zeroinfl(zinbreg(f, data = NMES1988), zi_ref(f, NMES1988))
})

test_that("zinbreg handles different count and zero covariates", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | ment + kid5
  expect_matches_zeroinfl(zinbreg(f, data = bioChemists), zi_ref(f, bioChemists))
})

test_that("zinbreg handles an intercept-only zero part and a single covariate", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f1 <- visits ~ hospital + chronic | 1
  expect_matches_zeroinfl(zinbreg(f1, data = NMES1988), zi_ref(f1, NMES1988))
  f2 <- visits ~ chronic
  expect_matches_zeroinfl(zinbreg(f2, data = NMES1988), zi_ref(f2, NMES1988))
})

test_that("zinbreg handles offsets (formula and argument)", {
  skip_if_not_installed("pscl")
  set.seed(11)
  n <- 600
  d <- data.frame(x = rnorm(n), z = rnorm(n), e = runif(n, 1, 5))
  pz <- plogis(-0.5 + 0.8 * d$z)
  d$y <- ifelse(runif(n) < pz, 0,
                rnbinom(n, mu = d$e * exp(0.5 + 0.4 * d$x), size = 2))
  f <- y ~ x + offset(log(e)) | z
  ref <- zi_ref(f, d)
  expect_matches_zeroinfl(zinbreg(f, data = d), ref)
  expect_matches_zeroinfl(zinbreg(y ~ x | z, data = d, offset = log(e)), ref)
})

test_that("EM never decreases the observed-data log-likelihood", {
  data("bioChemists", package = "pscl", envir = environment())
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  expect_true(all(diff(fit$loglik.path) > -1e-8))
  expect_equal(fit$loglik, tail(fit$loglik.path, 1))
})

test_that("analytic Hessian agrees with numDeriv::hessian", {
  skip_if_not_installed("numDeriv")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | ment + kid5
  fit <- zinbreg(f, data = bioChemists)
  X <- model.matrix(~ fem + mar + kid5 + phd + ment, bioChemists)
  Z <- model.matrix(~ ment + kid5, bioChemists)
  y <- bioChemists$art
  p <- ncol(X); q <- ncol(Z)
  th <- c(fit$coefficients$count, fit$coefficients$zero, fit$logtheta)
  llf <- function(par)
    zinbTruitt:::zinb_loglik(y, exp(drop(X %*% par[1:p])),
                             drop(Z %*% par[p + 1:q]), exp(par[p + q + 1]))
  Hn <- numDeriv::hessian(llf, th)
  Ha <- zinbTruitt:::zinb_hessian(X, Z, y, rep(0, length(y)), th[1:p],
                                  th[p + 1:q], exp(th[p + q + 1]))
  expect_lt(max(abs(Hn - Ha)) / max(abs(Ha)), 1e-5)
})

test_that("zinbreg recovers known parameters on a large simulated sample", {
  set.seed(42)
  n <- 20000
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  pz <- plogis(-1 + 0.8 * d$z)
  d$y <- ifelse(runif(n) < pz, 0, rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2))
  fit <- zinbreg(y ~ x | z, data = d)
  expect_equal(unname(fit$coefficients$count), c(1, 0.5), tolerance = 0.05)
  expect_equal(unname(fit$coefficients$zero), c(-1, 0.8), tolerance = 0.1)
  expect_equal(fit$theta, 2, tolerance = 0.1)
})

test_that("zinbreg warns when maxit is reached", {
  data("bioChemists", package = "pscl", envir = environment())
  expect_warning(
    zinbreg(art ~ fem + ment, data = bioChemists,
            control = zinbreg_control(maxit = 3)),
    "did not converge")
})

test_that("zinbreg gives informative errors for bad input", {
  d <- data.frame(x = rnorm(30), z = rnorm(30), y = rpois(30, 3) + 1)
  expect_error(zinbreg(y ~ x | z, data = d), "no zeros")
  d$y <- c(0, 0, 0, rpois(27, 3))
  d$yneg <- d$y - 1
  expect_error(zinbreg(yneg ~ x | z, data = d), "non-negative")
  d$yfrac <- d$y + 0.5
  expect_error(zinbreg(yfrac ~ x | z, data = d), "integer")
  d$xdup <- d$x
  expect_error(zinbreg(y ~ x + xdup | z, data = d), "rank deficient")
  expect_error(zinbreg(y ~ x | z + offset(x), data = d), "offsets")
  expect_error(zinbreg(y ~ x | z | x, data = d), "at most two")
  expect_error(zinbreg_control(epsilon = 0), "epsilon")
})
