nb_ref <- function(formula, data) {
  MASS::glm.nb(formula, data = data,
               control = glm.control(epsilon = 1e-10, maxit = 100))
}

expect_matches_glmnb <- function(fit, ref) {
  s <- summary(ref)
  expect_true(fit$converged)
  expect_lt(max(abs(fit$coefficients - coef(ref))), 1e-5)
  expect_lt(max(abs(sqrt(diag(fit$vcov)) - s$coefficients[, 2])), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
}

test_that("nbreg matches MASS::glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  expect_matches_glmnb(nbreg(f, data = bioChemists), nb_ref(f, bioChemists))
})

test_that("nbreg matches MASS::glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  expect_matches_glmnb(nbreg(f, data = NMES1988), nb_ref(f, NMES1988))
})

test_that("nbreg handles offsets (in formula and as argument)", {
  skip_if_not_installed("MASS")
  set.seed(2)
  d <- data.frame(x = rnorm(300), e = runif(300, 1, 5))
  d$y <- rnbinom(300, mu = d$e * exp(0.3 + 0.5 * d$x), size = 1.5)
  ref <- nb_ref(y ~ x + offset(log(e)), d)
  expect_matches_glmnb(nbreg(y ~ x + offset(log(e)), data = d), ref)
  expect_matches_glmnb(nbreg(y ~ x, data = d, offset = log(e)), ref)
})

test_that("nbreg works with a single covariate and intercept only", {
  skip_if_not_installed("MASS")
  set.seed(5)
  d <- data.frame(x = rnorm(200))
  d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
  expect_matches_glmnb(nbreg(y ~ x, data = d), nb_ref(y ~ x, d))
  expect_matches_glmnb(nbreg(y ~ 1, data = d), nb_ref(y ~ 1, d))
})

test_that("nbreg warns when theta hits the Poisson limit", {
  set.seed(3)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(1 + 0.2 * d$x))
  expect_warning(nbreg(y ~ x, data = d), "theta reached the upper bound")
})

test_that("nbreg warns when maxit is reached", {
  set.seed(5)
  d <- data.frame(x = rnorm(200))
  d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
  expect_warning(nbreg(y ~ x, data = d, control = nbreg_control(maxit = 1)),
                 "did not converge")
})

test_that("nbreg gives informative errors for bad input", {
  d <- data.frame(x = rnorm(20), y = rpois(20, 3))
  expect_error(nbreg(~ x, data = d), "response")
  d$yneg <- d$y - 10
  expect_error(nbreg(yneg ~ x, data = d), "non-negative")
  d$yfrac <- d$y + 0.5
  expect_error(nbreg(yfrac ~ x, data = d), "integer")
  d$xdup <- d$x
  expect_error(nbreg(y ~ x + xdup, data = d), "rank deficient")
  expect_error(nbreg_control(epsilon = -1), "epsilon")
})

test_that("nbreg drops missing values", {
  set.seed(8)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, mu = exp(1 + 0.3 * d$x), size = 2)
  d$x[1:5] <- NA
  fit <- nbreg(y ~ x, data = d)
  expect_equal(fit$nobs, 95)
})
