# Control parameters for zero-inflated negative binomial regression

Tolerances and iteration limits for
[`zinbreg()`](https://JosephTruitt.github.io/zinbTruitt/reference/zinbreg.md).

## Usage

``` r
zinbreg_control(
  epsilon = 1e-12,
  maxit = 2000L,
  irls_maxit = 50L,
  theta_max = 1e+06,
  nb_epsilon = 1e-10,
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance for the EM algorithm: stop when the relative
  change in the observed-data log-likelihood is below `epsilon`.

- maxit:

  Maximum number of EM iterations.

- irls_maxit:

  Maximum number of IRLS iterations in the zero-part M-step.

- theta_max:

  Upper bound for theta (Poisson limit); a warning is issued if it is
  reached.

- nb_epsilon:

  Convergence tolerance used by the count-part (NB) M-step.

- trace:

  Logical; print the log-likelihood at each EM iteration?

## Value

A list of control settings.

## Examples

``` r
zinbreg_control(epsilon = 1e-10, maxit = 500)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 500
#> 
#> $irls_maxit
#> [1] 50
#> 
#> $theta_max
#> [1] 1e+06
#> 
#> $nb_epsilon
#> [1] 1e-10
#> 
#> $trace
#> [1] FALSE
#> 
```
