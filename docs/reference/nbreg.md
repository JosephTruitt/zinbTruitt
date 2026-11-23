# Negative binomial regression by IRLS and Newton-Raphson

Fits the NB2 regression model with log link, \\Var(Y) = \mu +
\mu^2/\theta\\, without any general-purpose optimizer. The coefficients
\\\beta\\ are updated by iteratively reweighted least squares (Fisher
scoring) with \\\theta\\ held fixed, and \\\theta\\ is updated by
Newton-Raphson on \\\log\theta\\ with \\\beta\\ held fixed. The two
updates are alternated until the log-likelihood converges.

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term is allowed.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional numeric offset vector (e.g. `log(exposure)`), evaluated in
  `data`.

- control:

  A list of control settings from
  [`nbreg_control()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`, a list with components including
`coefficients`, `vcov`, `theta`, `se.theta`, `loglik`, `df`, `iter`,
`converged`, `fitted.values`, `linear.predictors`, `y`, `offset`, `nobs`
and `call`.

## Details

The parameterization matches
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). Standard
errors of \\\hat\beta\\ are \\(X^\top \hat W X)^{-1}\\, treating
\\\theta\\ as known, and the standard error of \\\hat\theta\\ is
\\\\-H(\hat\theta)\\^{-1/2}\\.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- rpois(200, exp(0.5 + 0.4 * d$x) * rgamma(200, shape = 2, rate = 2))
fit <- nbreg(y ~ x, data = d)
fit$coefficients
#> (Intercept)           x 
#>   0.4152361   0.3895846 
fit$theta
#> [1] 1.942204
```
