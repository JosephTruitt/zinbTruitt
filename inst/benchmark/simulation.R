# Simulation study: does zinbreg() recover known ZINB parameters?
#
# Run from the package root:  source("inst/benchmark/simulation.R")
# (requires the zinbTruitt package to be installed)
#
# Design: R = 200 data sets of size n = 500 from a ZINB model with one
# covariate in the count part (x) and one in the zero part (z).
#   count:  log(mu) = 1.0 + 0.5 x          (beta0, beta1)
#   zero:   logit(pi) = -1.0 + 0.8 z       (gamma0, gamma1)
#   theta = 2   (log theta = log 2)
# For each parameter we report bias, mean estimated SE vs Monte Carlo SD, and
# the empirical coverage of 95% Wald intervals.

suppressPackageStartupMessages(library(zinbTruitt))

out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

R <- 200
n <- 500
truth <- c(count_Intercept = 1.0, count_x = 0.5,
           zero_Intercept = -1.0, zero_z = 0.8,
           Log_theta = log(2))

simulate_data <- function(n) {
  x <- rnorm(n)
  z <- rnorm(n)
  pi_ <- plogis(truth["zero_Intercept"] + truth["zero_z"] * z)
  mu <- exp(truth["count_Intercept"] + truth["count_x"] * x)
  y <- ifelse(runif(n) < pi_, 0, rnbinom(n, mu = mu, size = exp(truth["Log_theta"])))
  data.frame(y = y, x = x, z = z)
}

set.seed(651)  # reproducible
est <- se <- matrix(NA_real_, R, length(truth), dimnames = list(NULL, names(truth)))
conv <- logical(R)
n_warn <- 0L

for (r in seq_len(R)) {
  d <- simulate_data(n)
  fit <- withCallingHandlers(
    tryCatch(zinbreg(y ~ x | z, data = d), error = function(e) NULL),
    warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
  if (is.null(fit)) next
  conv[r] <- fit$converged
  est[r, ] <- c(coef(fit, "count"), coef(fit, "zero"), fit$logtheta)
  se[r, ] <- sqrt(c(diag(vcov(fit)), fit$se.logtheta^2))
}

ok <- conv & stats::complete.cases(est) & stats::complete.cases(se)
cat("Replicates:", R, " usable (converged, finite):", sum(ok),
    " warnings:", n_warn, "\n\n")

z975 <- qnorm(0.975)
true_mat <- matrix(truth, sum(ok), length(truth), byrow = TRUE)
cover <- (est[ok, ] - z975 * se[ok, ] <= true_mat) &
  (true_mat <= est[ok, ] + z975 * se[ok, ])

results <- data.frame(
  parameter = names(truth),
  truth = unname(truth),
  mean_estimate = colMeans(est[ok, ]),
  bias = colMeans(est[ok, ]) - truth,
  mc_sd = apply(est[ok, ], 2, sd),
  mean_se = colMeans(se[ok, ]),
  coverage_95 = colMeans(cover),
  row.names = NULL)
results$mc_se_of_bias <- results$mc_sd / sqrt(sum(ok))

cat("=== ZINB simulation: R =", R, ", n =", n, "===\n")
print(results, digits = 3, row.names = FALSE)

saveRDS(list(results = results, est = est, se = se, converged = conv,
             truth = truth, R = R, n = n),
        file.path(out_dir, "simulation_results.rds"))
write.csv(results, file.path(out_dir, "simulation_summary.csv"), row.names = FALSE)
cat("\nResults saved to", out_dir, "\n")
