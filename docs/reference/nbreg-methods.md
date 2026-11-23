# Methods for `nbreg` objects

Standard S3 methods for fits from
[`nbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("response", "link"), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x, object:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the print
  method of a summary).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (ignored).

  Wald confidence intervals are available through the default
  [`confint()`](https://rdrr.io/r/stats/confint.html) method, which uses
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html).

- type:

  For `residuals`: `"pearson"` (default) or `"response"`. For `predict`:
  `"response"` (the mean) or `"link"` (the linear predictor).

- newdata:

  Optional data frame for `predict`.

## Value

`coef`, `vcov`, `fitted`, `residuals` and `predict` return numeric
vectors or matrices; `logLik` an object of class `"logLik"` with a `df`
attribute so that [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work; `nobs` an integer;
`summary` an object of class `"summary.nbreg"`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- rnbinom(200, mu = exp(0.5 + 0.4 * d$x), size = 2)
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.4990       0.3995  
#> 
#> Theta: 2.572  (SE 0.6704)
#> Log-likelihood: -346.51 on 3 Df
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.49896    0.07224   6.907 4.94e-12 ***
#> x            0.39946    0.07522   5.311 1.09e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.572  (SE 0.6704)
#> Log-likelihood: -346.51 on 3 Df;  AIC: 699.02
#> Number of observations: 200;  outer iterations: 3
AIC(fit)
#> [1] 699.0201
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>        1        2        3 
#> 1.104624 1.647012 2.455719 
```
