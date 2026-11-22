# ---- nbreg methods ---------------------------------------------------------
nb_setup <- function() {
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  list(d = bioChemists, f = f, fit = nbreg(f, data = bioChemists),
       ref = MASS::glm.nb(f, data = bioChemists,
                          control = glm.control(epsilon = 1e-10, maxit = 100)))
}

test_that("nbreg methods agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  s <- nb_setup()
  expect_equal(coef(s$fit), coef(s$ref), tolerance = 1e-5)
  expect_equal(unname(vcov(s$fit)), unname(vcov(s$ref)), tolerance = 1e-4)
  expect_equal(nobs(s$fit), nobs(s$ref))
  expect_equal(unname(fitted(s$fit)), unname(fitted(s$ref)), tolerance = 1e-5)
  expect_equal(unname(residuals(s$fit, "pearson")),
               unname(residuals(s$ref, "pearson")), tolerance = 1e-4)
  expect_equal(unname(residuals(s$fit, "response")),
               unname(residuals(s$ref, "response")), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(s$fit)), as.numeric(logLik(s$ref)),
               tolerance = 1e-6)
  expect_equal(attr(logLik(s$fit), "df"), attr(logLik(s$ref), "df"))
  expect_equal(AIC(s$fit), AIC(s$ref), tolerance = 1e-6)
  expect_equal(BIC(s$fit), BIC(s$ref), tolerance = 1e-6)
})

test_that("nbreg predict works with and without newdata", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  s <- nb_setup()
  expect_equal(predict(s$fit), fitted(s$fit))
  expect_equal(predict(s$fit, type = "link"), s$fit$linear.predictors)
  nd <- s$d[c(1, 50, 300), ]
  expect_equal(unname(predict(s$fit, newdata = nd)),
               unname(predict(s$ref, newdata = nd, type = "response")),
               tolerance = 1e-5)
  expect_equal(unname(predict(s$fit, newdata = nd, type = "link")),
               unname(predict(s$ref, newdata = nd, type = "link")),
               tolerance = 1e-5)
  # fitted values reproduced by predict on the full data
  expect_equal(unname(predict(s$fit, newdata = s$d)), unname(fitted(s$fit)),
               tolerance = 1e-8)
})

test_that("nbreg predict handles offsets given in the formula or as argument", {
  set.seed(2)
  d <- data.frame(x = rnorm(300), e = runif(300, 1, 5))
  d$y <- rnbinom(300, mu = d$e * exp(0.3 + 0.5 * d$x), size = 1.5)
  f1 <- nbreg(y ~ x + offset(log(e)), data = d)
  f2 <- nbreg(y ~ x, data = d, offset = log(e))
  nd <- data.frame(x = c(-1, 0, 1), e = c(1, 2, 3))
  expected <- nd$e * exp(coef(f1)[1] + coef(f1)[2] * nd$x)
  expect_equal(unname(predict(f1, nd)), unname(expected), tolerance = 1e-8)
  expect_equal(unname(predict(f2, nd)), unname(predict(f1, nd)), tolerance = 1e-6)
})

test_that("nbreg print and summary work", {
  skip_if_not_installed("pscl")
  s <- nb_setup()
  expect_output(print(s$fit), "Theta")
  sm <- summary(s$fit)
  expect_s3_class(sm, "summary.nbreg")
  expect_equal(colnames(sm$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(unname(sm$coefficients[, 2]), unname(summary(s$ref)$coefficients[, 2]),
               tolerance = 1e-5)
  expect_equal(unname(sm$coefficients[, 4]), unname(summary(s$ref)$coefficients[, 4]),
               tolerance = 1e-4)
  expect_output(print(sm), "Coefficients")
  ci <- confint(s$fit)
  expect_equal(dim(ci), c(6L, 2L))
  expect_true(all(ci[, 1] < coef(s$fit) & coef(s$fit) < ci[, 2]))
})

# ---- zinbreg methods ---------------------------------------------------------
zi_setup <- function() {
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | fem + kid5 + ment
  list(d = bioChemists, fit = zinbreg(f, data = bioChemists),
       ref = pscl::zeroinfl(f, data = bioChemists, dist = "negbin",
                            control = pscl::zeroinfl.control(reltol = 1e-12,
                                                             maxit = 10000)))
}

test_that("zinbreg methods agree with zeroinfl", {
  skip_if_not_installed("pscl")
  s <- zi_setup()
  tol <- 1e-3
  expect_equal(coef(s$fit), coef(s$ref), tolerance = tol)
  expect_equal(coef(s$fit, "count"), coef(s$ref, "count"), tolerance = tol)
  expect_equal(coef(s$fit, "zero"), coef(s$ref, "zero"), tolerance = tol)
  expect_equal(unname(vcov(s$fit)), unname(vcov(s$ref)), tolerance = 1e-2)
  expect_equal(dim(vcov(s$fit, "zero")), c(4L, 4L))
  expect_equal(nobs(s$fit), s$ref$n)
  expect_equal(unname(fitted(s$fit)), unname(fitted(s$ref)), tolerance = tol)
  expect_equal(unname(residuals(s$fit, "response")),
               unname(residuals(s$ref, "response")), tolerance = tol)
  expect_equal(unname(residuals(s$fit, "pearson")),
               unname(residuals(s$ref, "pearson")), tolerance = tol)
  expect_equal(as.numeric(logLik(s$fit)), as.numeric(logLik(s$ref)),
               tolerance = 1e-6)
  expect_equal(attr(logLik(s$fit), "df"), attr(logLik(s$ref), "df"))
  expect_equal(AIC(s$fit), AIC(s$ref), tolerance = 1e-6)
  expect_equal(BIC(s$fit), BIC(s$ref), tolerance = 1e-6)
})

test_that("zinbreg predict supports response, count and zero types", {
  skip_if_not_installed("pscl")
  s <- zi_setup()
  tol <- 1e-3
  for (ty in c("response", "count", "zero")) {
    expect_equal(unname(predict(s$fit, type = ty)),
                 unname(predict(s$ref, type = ty)), tolerance = tol)
  }
  nd <- s$d[c(1, 50, 300), ]
  for (ty in c("response", "count", "zero")) {
    expect_equal(unname(predict(s$fit, newdata = nd, type = ty)),
                 unname(predict(s$ref, newdata = nd, type = ty)), tolerance = tol)
  }
  # relationship between the three types
  expect_equal(predict(s$fit, nd, "response"),
               (1 - predict(s$fit, nd, "zero")) * predict(s$fit, nd, "count"))
  # fitted values reproduced by predict on the full data
  expect_equal(unname(predict(s$fit, newdata = s$d)), unname(fitted(s$fit)),
               tolerance = 1e-8)
})

test_that("zinbreg predict handles offsets", {
  set.seed(11)
  n <- 600
  d <- data.frame(x = rnorm(n), z = rnorm(n), e = runif(n, 1, 5))
  pz <- plogis(-0.5 + 0.8 * d$z)
  d$y <- ifelse(runif(n) < pz, 0,
                rnbinom(n, mu = d$e * exp(0.5 + 0.4 * d$x), size = 2))
  f1 <- zinbreg(y ~ x + offset(log(e)) | z, data = d)
  f2 <- zinbreg(y ~ x | z, data = d, offset = log(e))
  nd <- data.frame(x = c(-1, 0, 1), z = c(0, 1, -1), e = c(1, 2, 3))
  expect_equal(unname(predict(f1, nd, "count")),
               unname(nd$e * exp(coef(f1, "count")[1] + coef(f1, "count")[2] * nd$x)),
               tolerance = 1e-8)
  expect_equal(unname(predict(f2, nd)), unname(predict(f1, nd)), tolerance = 1e-5)
  expect_equal(unname(predict(f1, d)), unname(fitted(f1)), tolerance = 1e-8)
})

test_that("zinbreg print and summary work", {
  skip_if_not_installed("pscl")
  s <- zi_setup()
  expect_output(print(s$fit), "Zero-inflation")
  sm <- summary(s$fit)
  expect_s3_class(sm, "summary.zinbreg")
  expect_equal(dim(sm$count), c(6L, 4L))
  expect_equal(dim(sm$zero), c(4L, 4L))
  expect_equal(unname(sm$count[, 2]),
               unname(summary(s$ref)$coefficients$count[1:6, 2]), tolerance = 1e-3)
  expect_equal(unname(sm$zero[, 2]),
               unname(summary(s$ref)$coefficients$zero[, 2]), tolerance = 1e-3)
  expect_output(print(sm), "Log\\(theta\\)")
  ci <- confint(s$fit)
  expect_equal(dim(ci), c(10L, 2L))
  expect_equal(rownames(ci), names(coef(s$fit)))
})

test_that("predict works with factor covariates and new factor levels are rejected", {
  set.seed(4)
  n <- 500
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)), x = rnorm(n))
  d$y <- ifelse(runif(n) < 0.2, 0, rnbinom(n, mu = exp(1 + 0.3 * d$x + (d$g == "b")), size = 2))
  fit <- zinbreg(y ~ g + x | x, data = d)
  nd <- data.frame(g = factor(c("a", "c"), levels = c("a", "b", "c")), x = c(0, 1))
  expect_length(predict(fit, nd), 2)
  expect_error(predict(fit, data.frame(g = "z", x = 0)))
})
