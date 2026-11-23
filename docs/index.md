# zinbTruitt

**zinbTruitt** fits count-data regression models without any
general-purpose optimizer
([`optim()`](https://rdrr.io/r/stats/optim.html) and friends are never
called):

- [`nbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg.md):
  negative binomial (NB2) regression, using IRLS (Fisher scoring) for
  the coefficients alternated with Newton-Raphson for the dispersion
  `theta`. Comparable to
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).
- [`zinbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/zinbreg.md):
  zero-inflated negative binomial regression, using an EM algorithm
  whose M-steps are IRLS fits. Standard errors come from the observed
  information matrix. Comparable to
  `pscl::zeroinfl(..., dist = "negbin")`.

Developed for the STA 651 Final Project (Fall 2026), Miami University.

**Website and report:** <https://JosephTruitt.github.io/zinbTruitt/>

## Installation

``` r

# install.packages("devtools")
devtools::install_github("JosephTruitt/zinbTruitt")
```

## Usage

``` r

library(zinbTruitt)
data("bioChemists", package = "pscl")

# Negative binomial regression
fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_nb)

# Zero-inflated negative binomial: count covariates | zero covariates
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment,
                    data = bioChemists)
summary(fit_zinb)

AIC(fit_nb, fit_zinb)
predict(fit_zinb, type = "response")   # also "count" and "zero"
```

## Methods provided

Both `nbreg` and `zinbreg` objects support `print`, `summary`, `coef`,
`vcov`, `logLik`, `fitted`, `residuals` (response and Pearson), `nobs`,
and `predict` (with `newdata`).

## Agreement with the reference implementations

[`nbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg.md)
reproduces [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)
and
[`zinbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/zinbreg.md)
reproduces `pscl::zeroinfl(..., dist = "negbin")` on
[`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html) and
[`AER::NMES1988`](https://rdrr.io/pkg/AER/man/NMES1988.html), well
inside the required tolerances (coefficients and standard errors to
about 1e-8 for NB and 1e-4 for ZINB). The benchmark and simulation
scripts are in `inst/benchmark/` with saved results in
`inst/benchmark/results/`.

## Report

The full report (methods, benchmarking, simulation study and the
NMES1988 analysis) is a Quarto document,
`vignettes/articles/report.qmd`. It is rendered to `report/report.pdf`
and appears as an article on the [package
website](https://JosephTruitt.github.io/zinbTruitt/).

## R CMD check notes

*(To be updated after the final
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
on the submission machine.)*

## No general-purpose optimizers

No function in `R/` calls
[`optim()`](https://rdrr.io/r/stats/optim.html),
[`optimize()`](https://rdrr.io/r/stats/optimize.html),
[`nlm()`](https://rdrr.io/r/stats/nlm.html),
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html),
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html),
[`glm()`](https://rdrr.io/r/stats/glm.html),
[`lm()`](https://rdrr.io/r/stats/lm.html) or any other prohibited
function. `numDeriv` is used only in the tests to verify the analytic
Hessian.
