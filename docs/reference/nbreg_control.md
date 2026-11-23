# Control parameters for negative binomial regression

Tolerances and iteration limits for
[`nbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/nbreg.md).

## Usage

``` r
nbreg_control(
  epsilon = 1e-13,
  maxit = 100L,
  irls_maxit = 50L,
  theta_maxit = 50L,
  theta_max = 1e+06,
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance (default 1e-13). The outer (beta, theta) loop
  stops when the relative change in the log-likelihood is below
  `epsilon`.

- maxit:

  Maximum number of outer iterations (IRLS for beta alternated with
  Newton-Raphson for theta).

- irls_maxit:

  Maximum number of IRLS iterations in each beta update.

- theta_maxit:

  Maximum number of Newton-Raphson iterations in each theta update.

- theta_max:

  Upper bound for theta. If theta reaches this value the model is
  essentially Poisson and a warning is issued.

- trace:

  Logical; print the log-likelihood and theta at each outer iteration?

## Value

A list of control settings.

## Examples

``` r
nbreg_control(epsilon = 1e-12, trace = TRUE)
#> $epsilon
#> [1] 1e-12
#> 
#> $maxit
#> [1] 100
#> 
#> $irls_maxit
#> [1] 50
#> 
#> $theta_maxit
#> [1] 50
#> 
#> $theta_max
#> [1] 1e+06
#> 
#> $trace
#> [1] TRUE
#> 
```
