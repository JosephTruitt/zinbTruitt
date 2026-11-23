# Methods for `zinbreg` objects

Standard S3 methods for fits from
[`zinbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'zinbreg'
predict(object, newdata, type = c("response", "count", "zero"), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x, object:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the print
  method of a summary).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (ignored).

- model:

  Which part of the model: `"full"` (default), `"count"` or `"zero"`.

- type:

  For `residuals`: `"pearson"` (default) or `"response"`. For `predict`:
  `"response"`, `"count"` or `"zero"`.

- newdata:

  Optional data frame for `predict`.

## Value

See
[nbreg-methods](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg-methods.md)
for the common methods. `summary` returns an object of class
`"summary.zinbreg"` with coefficient tables for the count and zero
parts.

## Details

`coef` and `vcov` return the coefficients (and their covariance matrix)
of the count part, the zero part, or both (`"full"`, with names
`count_*` and `zero_*`, as in
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)). The
covariance of the full fit including `Log(theta)` is stored in the
fitted object as `vcov`.

For `predict`, `type = "response"` gives the marginal mean
\\(1-\pi)\mu\\, `"count"` the NB mean \\\mu\\, and `"zero"` the
structural-zero probability \\\pi\\. Pearson residuals use the ZINB
variance \\(1-\pi)\mu\\1 + \mu(\pi + 1/\theta)\\\\.

## Examples

``` r
set.seed(1)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
pz <- plogis(-1 + 0.8 * d$z)
d$y <- ifelse(runif(n) < pz, 0, rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2))
fit <- zinbreg(y ~ x | z, data = d)
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Count model coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.06580    0.07757  13.740  < 2e-16 ***
#> x            0.46894    0.06501   7.213 5.47e-13 ***
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.1858     0.2936  -4.039 5.38e-05 ***
#> z             0.9022     0.2132   4.231 2.33e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.523  (SE 0.327)
#> Log(theta): 0.4209  (SE 0.2146)
#> Log-likelihood: -768.513 on 5 Df;  AIC: 1547.03
#> Number of observations: 400;  EM iterations: 61
coef(fit, model = "zero")
#> (Intercept)           z 
#>  -1.1858462   0.9021801 
head(predict(fit, type = "zero"))
#>         1         2         3         4         5         6 
#> 0.4460837 0.6281704 0.1506057 0.1767564 0.1734525 0.1787621 
predict(fit, newdata = data.frame(x = 0, z = 0), type = "response")
#>        1 
#> 2.223808 
```
