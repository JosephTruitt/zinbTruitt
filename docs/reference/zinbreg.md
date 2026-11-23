# Zero-inflated negative binomial regression by EM with IRLS M-steps

Fits the ZINB model: with probability \\\pi_i\\ an observation is a
structural zero, otherwise it follows an NB(\\\mu_i\\, \\\theta\\)
distribution, with \\\log \mu_i = x_i^\top\beta + o_i\\ and
\\\mathrm{logit}\\ \pi_i = z_i^\top\gamma\\. The likelihood is maximized
by the EM algorithm. The E-step computes the posterior probability
\\\tau_i\\ of being a structural zero. The zero-part M-step is a
logistic regression with fractional response \\\tau_i\\ fitted by IRLS,
and the count-part M-step is the NB fit of
[`nbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg.md)
with prior weights \\1 - \tau_i\\. No general-purpose optimizer is used.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before the bar, zero-model covariates after it. If there is no `|`,
  the same covariates are used in both parts. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the count
  part is allowed.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional numeric offset for the count model, evaluated in `data`.

- control:

  A list of control settings from
  [`zinbreg_control()`](https://JosephTruitt.github.io/zinbTruitt/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`, a list with components including
`coefficients` (a list with elements `count` and `zero`), `vcov` (for
`count_*`, `zero_*` and `Log(theta)`), `theta`, `se.theta`, `logtheta`,
`se.logtheta`, `loglik`, `df`, `iter`, `converged`, `loglik.path`,
`fitted.values`, `mu`, `pi`, `tau`, `y` and `call`.

## Details

Standard errors come from the inverse of the observed information matrix
of the observed-data log-likelihood with respect to \\(\beta, \gamma,
\log\theta)\\, computed from the analytic Hessian.

## Examples

``` r
set.seed(1)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
pz <- plogis(-1 + 0.8 * d$z)
d$y <- ifelse(runif(n) < pz, 0, rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2))
fit <- zinbreg(y ~ x | z, data = d)
fit$coefficients
#> $count
#> (Intercept)           x 
#>   1.0657977   0.4689358 
#> 
#> $zero
#> (Intercept)           z 
#>  -1.1858462   0.9021801 
#> 
fit$theta
#> [1] 1.523283
```
