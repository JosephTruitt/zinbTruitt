# Benchmark: zinbTruitt vs MASS::glm.nb() and pscl::zeroinfl()
#
# Run from the package root:  source("inst/benchmark/benchmark.R")
# (requires the zinbTruitt package to be installed, plus MASS, pscl, AER)
#
# Output: agreement tables (printed and saved to inst/benchmark/results/).

suppressPackageStartupMessages({
  library(zinbTruitt)
  library(MASS)
  library(pscl)
})
data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Table 4 tolerances (maximum absolute difference)
tol_nb   <- c(coef = 1e-5, se = 1e-5, theta = 1e-4, loglik = 1e-4)
tol_zinb <- c(coef = 1e-3, se = 1e-3, theta = 1e-3, loglik = 1e-3)

# Reference fits with tight tolerances
ref_nb <- function(f, d)
  glm.nb(f, data = d, control = glm.control(epsilon = 1e-10, maxit = 100))
ref_zinb <- function(f, d)
  zeroinfl(f, data = d, dist = "negbin",
           control = zeroinfl.control(reltol = 1e-12, maxit = 10000))

# Median elapsed time over `reps` fits
time_it <- function(expr_fun, reps = 5) {
  median(vapply(seq_len(reps), function(i)
    unname(system.time(expr_fun())["elapsed"]), numeric(1)))
}

datasets <- list(
  bioChemists = list(
    data = bioChemists,
    nb_formula = art ~ fem + mar + kid5 + phd + ment,
    zinb_formula = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment),
  NMES1988 = list(
    data = NMES1988,
    nb_formula = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb_formula = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance)
)

est_tables <- list()   # side-by-side estimates
summ_rows <- list()    # max-abs-difference summary
time_rows <- list()

for (dn in names(datasets)) {
  d <- datasets[[dn]]$data

  # ---------------- NB ----------------
  f <- datasets[[dn]]$nb_formula
  mine <- nbreg(f, data = d)
  ref <- ref_nb(f, d)
  s_ref <- summary(ref)

  tab <- data.frame(
    dataset = dn, model = "NB", term = names(coef(mine)),
    ours = unname(coef(mine)), reference = unname(coef(ref)),
    se_ours = sqrt(diag(vcov(mine))),
    se_reference = unname(s_ref$coefficients[, 2]))
  tab$abs_diff_coef <- abs(tab$ours - tab$reference)
  tab$abs_diff_se <- abs(tab$se_ours - tab$se_reference)
  tab <- rbind(tab, data.frame(
    dataset = dn, model = "NB", term = "theta",
    ours = mine$theta, reference = ref$theta,
    se_ours = mine$se.theta, se_reference = ref$SE.theta,
    abs_diff_coef = abs(mine$theta - ref$theta),
    abs_diff_se = abs(mine$se.theta - ref$SE.theta)))
  est_tables[[paste(dn, "NB")]] <- tab

  summ_rows[[paste(dn, "NB")]] <- data.frame(
    dataset = dn, model = "NB",
    coef = max(tab$abs_diff_coef[tab$term != "theta"]),
    se = max(tab$abs_diff_se[tab$term != "theta"]),
    theta = abs(mine$theta - ref$theta),
    loglik = abs(mine$loglik - as.numeric(logLik(ref))),
    loglik_ours = mine$loglik, loglik_reference = as.numeric(logLik(ref)))
  time_rows[[paste(dn, "NB")]] <- data.frame(
    dataset = dn, model = "NB",
    ours_sec = time_it(function() nbreg(f, data = d)),
    reference_sec = time_it(function() ref_nb(f, d)),
    iterations = mine$iter)

  # ---------------- ZINB ----------------
  f <- datasets[[dn]]$zinb_formula
  mine <- zinbreg(f, data = d)
  ref <- ref_zinb(f, d)
  se_ref <- sqrt(diag(vcov(ref)))
  est_ref <- coef(ref)
  est_mine <- coef(mine)
  se_mine <- sqrt(diag(vcov(mine)))

  tab <- data.frame(
    dataset = dn, model = "ZINB", term = names(est_mine),
    ours = unname(est_mine), reference = unname(est_ref[names(est_mine)]),
    se_ours = unname(se_mine), se_reference = unname(se_ref[names(est_mine)]))
  tab$abs_diff_coef <- abs(tab$ours - tab$reference)
  tab$abs_diff_se <- abs(tab$se_ours - tab$se_reference)
  tab <- rbind(tab, data.frame(
    dataset = dn, model = "ZINB", term = "Log(theta)",
    ours = mine$logtheta, reference = log(ref$theta),
    se_ours = mine$se.logtheta, se_reference = ref$SE.logtheta,
    abs_diff_coef = abs(mine$logtheta - log(ref$theta)),
    abs_diff_se = abs(mine$se.logtheta - ref$SE.logtheta)))
  est_tables[[paste(dn, "ZINB")]] <- tab

  summ_rows[[paste(dn, "ZINB")]] <- data.frame(
    dataset = dn, model = "ZINB",
    coef = max(tab$abs_diff_coef[tab$term != "Log(theta)"]),
    se = max(tab$abs_diff_se[tab$term != "Log(theta)"]),
    theta = abs(mine$logtheta - log(ref$theta)),
    loglik = abs(mine$loglik - as.numeric(logLik(ref))),
    loglik_ours = mine$loglik, loglik_reference = as.numeric(logLik(ref)))
  time_rows[[paste(dn, "ZINB")]] <- data.frame(
    dataset = dn, model = "ZINB",
    ours_sec = time_it(function() zinbreg(f, data = d)),
    reference_sec = time_it(function() ref_zinb(f, d)),
    iterations = mine$iter)
}

estimates <- do.call(rbind, est_tables); rownames(estimates) <- NULL
summary_tab <- do.call(rbind, summ_rows); rownames(summary_tab) <- NULL
timing <- do.call(rbind, time_rows); rownames(timing) <- NULL

# Pass/fail against the Table 4 tolerances
summary_tab$pass <- with(summary_tab, ifelse(
  model == "NB",
  coef <= tol_nb["coef"] & se <= tol_nb["se"] &
    theta <= tol_nb["theta"] & loglik <= tol_nb["loglik"],
  coef <= tol_zinb["coef"] & se <= tol_zinb["se"] &
    theta <= tol_zinb["theta"] & loglik <= tol_zinb["loglik"]))

cat("\n=== Maximum absolute differences (ours vs reference) ===\n")
cat("NB tolerances: coef/se 1e-5, theta 1e-4, loglik 1e-4\n")
cat("ZINB tolerances: all 1e-3 (theta column is log(theta) for ZINB)\n\n")
print(summary_tab, digits = 3, row.names = FALSE)

cat("\n=== Computation time (median of 5 fits, seconds) ===\n\n")
print(timing, digits = 3, row.names = FALSE)

cat("\n=== Side-by-side estimates ===\n")
for (nm in names(est_tables)) {
  cat("\n---", nm, "---\n")
  print(est_tables[[nm]][, -(1:2)], digits = 6, row.names = FALSE)
}

saveRDS(list(estimates = estimates, summary = summary_tab, timing = timing),
        file.path(out_dir, "benchmark_results.rds"))
write.csv(estimates, file.path(out_dir, "benchmark_estimates.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)
write.csv(timing, file.path(out_dir, "benchmark_timing.csv"), row.names = FALSE)
cat("\nResults saved to", out_dir, "\n")
