# 1 Introduction

# 1 Introduction

Count outcomes such as the number of articles an author publishes or the
number of physician visits a patient makes are routinely analyzed with
Poisson regression. The Poisson model forces the variance to equal the
mean, an assumption that real count data usually violate in two ways.
First, the data are often *overdispersed*: the variance exceeds the mean
because of unobserved heterogeneity. The negative binomial (NB) model
handles this by adding a dispersion parameter $`\theta`$. Second, the
data often contain *more zeros* than even an NB model predicts, because
some units are never at risk of a positive count. The zero-inflated
negative binomial (ZINB) model handles this by mixing the NB
distribution with a point mass at zero (Lambert, 1992; Cameron and
Trivedi, 2013).

The standard R implementations, `MASS::glm.nb()` and `pscl::zeroinfl()`,
are reliable but treat the estimation as a black box. The goal of this
project was to build an R package, `zinbTruitt`, that fits both models
from first principles and reproduces the reference implementations. The
estimation uses only closed-form building blocks: iteratively reweighted
least squares (IRLS) for regression coefficients, Newtonâ€“Raphson for
the dispersion parameter, and the EM algorithm (Dempster, Laird and
Rubin, 1977) for the zero-inflation mixture. No call to `optim()` or any
other general-purpose optimizer is made.

The report describes the algorithms (Section 2), the package (Section
3), agreement with the reference functions on two data sets (Section 4),
a simulation study (Section 5), an analysis of physician office visits
in the 1987â€“88 National Medical Expenditure Survey (Section 6), and
limitations (Section 7). All numbers, tables, and figures are produced
by code in this document that calls the installed package.

# 2 Methods

## 2.1 Negative binomial regression

Let $`Y_i`$ be a count with covariate vector $`x_i`$ and optional offset
$`o_i`$. The NB2 model has probability function
``` math
P(Y_i = y) = \frac{\Gamma(y + \theta)}{\Gamma(\theta)\,y!}\left(\frac{\theta}{\theta + \mu_i}\right)^{\theta}\left(\frac{\mu_i}{\theta + \mu_i}\right)^{y}, \qquad \log \mu_i = x_i^\top\beta + o_i,
```
so that $`E(Y_i) = \mu_i`$ and
$`\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta`$ with $`\theta > 0`$. As
$`\theta \to \infty`$ the model reduces to Poisson. With prior weights
$`a_i`$ (all equal to one for an ordinary fit) the log-likelihood is
``` math
\ell(\beta, \theta) = \sum_i a_i\Big[\log\Gamma(y_i + \theta) - \log\Gamma(\theta) - \log y_i! + \theta\log\theta + y_i\log\mu_i - (\theta + y_i)\log(\theta + \mu_i)\Big].
```
This is the parameterization of `MASS::glm.nb()`.

**IRLS for $`\beta`$ (with $`\theta`$ fixed).** The log link gives
$`\partial\mu_i/\partial\eta_i = \mu_i`$ and the variance function is
$`V(\mu_i) = \mu_i + \mu_i^2/\theta`$. Fisher scoring is then equivalent
to repeated weighted least squares with working weights and working
response
``` math
w_i = a_i\,\frac{(\partial\mu_i/\partial\eta_i)^2}{V(\mu_i)} = a_i\,\frac{\theta\mu_i}{\theta + \mu_i}, \qquad u_i = (\eta_i - o_i) + \frac{y_i - \mu_i}{\mu_i},
```
``` math
\beta^{(t+1)} = (X^\top W X)^{-1}X^\top W u, \qquad W = \mathrm{diag}(w_1, \dots, w_n).
```
The working response $`u_i`$ is a first-order Taylor expansion of the
link function at the observed $`y_i`$, and its approximate variance is
$`1/w_i = 1/\mu_i + 1/\theta`$, the sum of a Poisson component and an
overdispersion component. The weighted least squares step is solved by a
QR decomposition. Because $`u_i`$ and $`w_i`$ depend on $`\mu_i`$, and
hence on $`\beta`$, the step is repeated until the log-likelihood stops
changing.

**Newtonâ€“Raphson for $`\theta`$ (with $`\beta`$ fixed).** The score
and second derivative of $`\ell`$ in $`\theta`$ are
``` math
U(\theta) = \sum_i a_i\Big[\psi(y_i + \theta) - \psi(\theta) + \log\theta + 1 - \log(\theta + \mu_i) - \tfrac{y_i + \theta}{\theta + \mu_i}\Big],
```
``` math
H(\theta) = \sum_i a_i\Big[\psi'(y_i + \theta) - \psi'(\theta) + \tfrac{1}{\theta} - \tfrac{2}{\theta + \mu_i} + \tfrac{y_i + \theta}{(\theta + \mu_i)^2}\Big],
```
where $`\psi`$ and $`\psi'`$ are the digamma and trigamma functions. To
keep $`\theta`$ positive the update is made on $`\phi = \log\theta`$
using $`U_\phi = \theta U`$ and $`H_\phi = \theta^2 H + \theta U`$:
``` math
\phi \leftarrow \phi - U_\phi / H_\phi .
```
Steps are bounded and halved if the log-likelihood does not increase,
and $`\theta`$ is capped at a large value (default $`10^6`$), at which a
warning reports the Poisson limit.

**Outer loop.** The fit starts from a Poisson IRLS fit and a moment
estimate of $`\theta`$ based on
$`\mathrm{Var}(Y) = \mu + \mu^2/\theta`$. The IRLS and Newtonâ€“Raphson
updates are then alternated until the relative change in the
log-likelihood falls below $`10^{-13}`$.

## 2.2 Zero-inflated negative binomial regression

In the ZINB model an observation is a structural zero with probability
$`\pi_i`$ and otherwise follows NB($`\mu_i, \theta`$), where
$`\log\mu_i = x_i^\top\beta + o_i`$ and
$`\mathrm{logit}\,\pi_i = z_i^\top\gamma`$. The covariates $`z_i`$ may
differ from $`x_i`$. Writing
$`f_{NB}(0;\mu,\theta) = \{\theta/(\theta + \mu)\}^\theta`$,
``` math
P(Y_i = 0) = \pi_i + (1 - \pi_i)f_{NB}(0;\mu_i,\theta), \qquad P(Y_i = y) = (1 - \pi_i)f_{NB}(y;\mu_i,\theta),\ y \ge 1,
```
and $`E(Y_i) = (1 - \pi_i)\mu_i`$. The observed-data log-likelihood
contains the logarithm of a sum for the zeros, so it is not the
log-likelihood of an exponential-family GLM and IRLS cannot be applied
to it directly. Had the latent indicators $`D_i`$ (1 for a structural
zero) been observed, the problem would separate into a logistic
regression for $`\gamma`$ and an NB fit for $`(\beta, \theta)`$ on the
non-structural observations. The EM algorithm exploits this.

**E-step.** Given current $`(\beta, \gamma, \theta)`$, compute the
posterior probability of a structural zero,
``` math
\tau_i = P(D_i = 1 \mid y_i) = \begin{cases} \dfrac{\pi_i}{\pi_i + (1 - \pi_i)f_{NB}(0;\mu_i,\theta)}, & y_i = 0,\\[2ex] 0, & y_i > 0. \end{cases}
```
In the code this is evaluated as
$`\tau_i = \mathrm{plogis}\{z_i^\top\gamma - \log f_{NB}(0;\mu_i,\theta)\}`$
for zeros, which is numerically stable.

**M-step, zero part.** Maximize
$`\sum_i \{\tau_i\log\pi_i + (1 - \tau_i)\log(1 - \pi_i)\}`$ over
$`\gamma`$. This is a logistic regression with fractional response
$`\tau_i`$, fitted by IRLS with
``` math
w_i = \pi_i(1 - \pi_i), \qquad u_i = z_i^\top\gamma + \frac{\tau_i - \pi_i}{\pi_i(1 - \pi_i)} .
```

**M-step, count part.** Maximize
$`\sum_i (1 - \tau_i)\log f_{NB}(y_i;\mu_i,\theta)`$ over
$`(\beta, \theta)`$. This is exactly the NB fit above with prior weights
$`a_i = 1 - \tau_i`$, so the same code is reused, warm-started at the
current values.

**Starting values and convergence.** $`(\beta, \theta)`$ start at the
plain NB fit and $`\gamma`$ at a logistic IRLS fit of $`I(y_i = 0)`$ on
$`z_i`$. The algorithm stops when the observed-data log-likelihood
changes by a relative amount below $`10^{-12}`$. EM never decreases the
observed-data log-likelihood, and the package stores the log-likelihood
at every iteration so this can be checked directly (Section 4).

## 2.3 Standard errors

For NB regression,
$`\widehat{\mathrm{Var}}(\hat\beta) = (X^\top\hat W X)^{-1}`$, treating
$`\theta`$ as known as `glm.nb()` does, and
$`\mathrm{SE}(\hat\theta) = \{-H(\hat\theta)\}^{-1/2}`$.

EM does not produce standard errors, so for ZINB the observed
information matrix of the observed-data log-likelihood is computed at
the maximum likelihood estimate with respect to
$`(\beta, \gamma, \log\theta)`$, and its inverse gives the covariance
matrix. The Hessian is derived analytically. For each observation let
$`\eta_i = x_i^\top\beta + o_i`$, $`\zeta_i = z_i^\top\gamma`$, and
$`\phi = \log\theta`$. For $`y_i > 0`$ the contribution is
$`\log(1 - \pi_i) + \log f_{NB}`$, whose second derivatives are those of
the NB log-density in $`(\eta, \phi)`$ together with
$`-\pi_i(1 - \pi_i)`$ for $`\zeta`$. For $`y_i = 0`$ the contribution is
$`\log P_i`$ with $`P_i = \pi_i + (1 - \pi_i)f_0`$, so each second
derivative is $`P_{ab}/P - P_aP_b/P^2`$ with the partial derivatives of
$`P_i`$ available in closed form. The per-observation derivatives are
assembled into the block Hessian
``` math
H = \begin{pmatrix} X^\top D_{\eta\eta}X & X^\top D_{\eta\zeta}Z & X^\top d_{\eta\phi}\\ \cdot & Z^\top D_{\zeta\zeta}Z & Z^\top d_{\zeta\phi}\\ \cdot & \cdot & \sum_i h_{\phi\phi,i}\end{pmatrix}.
```
The standard error of $`\hat\theta`$ follows from the delta method,
$`\mathrm{SE}(\hat\theta) = \hat\theta\,\mathrm{SE}(\log\hat\theta)`$.
As a check, the analytic Hessian is compared with `numDeriv::hessian()`
applied to the packageâ€™s own log-likelihood function in Section 4.

# 3 The R package

The package is available at <https://github.com/JosephTruitt/zinbTruitt>
and installs with

``` r
# install.packages("devtools")
devtools::install_github("JosephTruitt/zinbTruitt")
```

Its main functions are summarized in [Table 1](#tbl-funs).

| Function | Purpose |
|:---|:---|
| nbreg(formula, data, offset, control) | NB regression by IRLS and Newtonâ€“Raphson |
| zinbreg(formula, data, offset, control) | ZINB regression by EM with IRLS M-steps; formula y ~ x1 + x2 \| z1 + z2 |
| nbreg_control(), zinbreg_control() | Tolerances, maximum iterations, trace option |

Table 1: Exported functions of the zinbTruitt package.

Both fitting functions return S3 objects (classes `"nbreg"` and
`"zinbreg"`) that store the coefficients, covariance matrix,
$`\hat\theta`$ and its standard error, log-likelihood, iteration count,
convergence flag, fitted values, and the matched call. The methods
`print`, `summary`, `coef`, `vcov`, `logLik` (with a `df` attribute, so
`AIC()` and `BIC()` work), `fitted`, `residuals` (Pearson and response),
`nobs`, and `predict` (with `newdata`) are provided for both classes.
For `zinbreg`, `predict` accepts `type = "response"` for
$`(1 - \pi)\mu`$, `"count"` for $`\mu`$, and `"zero"` for $`\pi`$. Wald
intervals are available through `confint()`. Without a `|` in the
formula, `zinbreg` uses the same covariates in both parts, as
`zeroinfl()` does. The package is documented with roxygen2, tested with
testthat (including comparisons against `glm.nb()` and `zeroinfl()`),
and checked on Windows, macOS, and Linux by GitHub Actions. A short
usage example:

``` r
fit_nb   <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | fem + kid5 + ment,
                    data = bioChemists)
AIC(fit_nb, fit_zinb)
```

             df      AIC
    fit_nb    7 3135.917
    fit_zinb 11 3125.049

``` r
predict(fit_zinb, newdata = bioChemists[1:3, ], type = "zero")
```

               1            2            3 
    0.0008138696 0.0066822590 0.0066822590 

# 4 Benchmarking

The fits from `nbreg()` and `zinbreg()` were compared with
`MASS::glm.nb()` and `pscl::zeroinfl(..., dist = "negbin")` on
`pscl::bioChemists` (915 biochemistry PhD students; response `art`) and
`AER::NMES1988` (4,406 individuals; response `visits`). The reference
functions were run with tight tolerances (`epsilon = 1e-10` for
`glm.nb()` and `reltol = 1e-12` for `zeroinfl()`). The NB models use all
covariates in the count part; the ZINB models use the same covariates in
the count part and, for bioChemists, in the zero part, and for NMES1988
the covariates of Section 6.

[Table 2](#tbl-agree) reports the maximum absolute differences between
the two implementations, and the full side-by-side estimates are in
Appendix A. All differences are well inside the required tolerances (the
â€œThetaâ€? column is $`\theta`$ for NB and $`\log\theta`$ for ZINB).
The remaining ZINB differences of order $`10^{-5}`$ reflect the stopping
rules of the two optimizers: on bioChemists, when both are run to very
tight tolerances, the coefficients agree to about $`2 \times 10^{-6}`$.

| Data        | Model | Coefficients | Std. errors | Theta   | Log-lik. |
|:------------|:------|:-------------|:------------|:--------|:---------|
| bioChemists | NB    | 7.8e-09      | 8.2e-11     | 3.0e-09 | 0.0e+00  |
| bioChemists | ZINB  | 5.8e-05      | 3.5e-05     | 4.9e-06 | 6.2e-09  |
| NMES1988    | NB    | 6.2e-08      | 1.7e-10     | 1.3e-09 | 1.8e-12  |
| NMES1988    | ZINB  | 3.3e-05      | 2.2e-05     | 5.2e-06 | 2.4e-08  |

Table 2: Maximum absolute difference between zinbTruitt and the
reference functions. Tolerances (coefficients / standard errors / theta
/ log-likelihood): NB 1e-5, 1e-5, 1e-4, 1e-4; ZINB 1e-3 for all.

[Table 3](#tbl-timing) compares median computation times over five fits.
The NB fit is about twice as fast as `glm.nb()`; the ZINB fit is slower
than `zeroinfl()` because EM is a first-order method that needs tens to
over a hundred iterations, while BFGS exploits curvature information.
The package is not designed to be faster.

| Data        | Model | Ours  | Reference | Iterations |
|:------------|:------|:------|:----------|-----------:|
| bioChemists | NB    | 0.030 | 0.040     |          4 |
| bioChemists | ZINB  | 1.140 | 0.170     |        121 |
| NMES1988    | NB    | 0.090 | 0.170     |          3 |
| NMES1988    | ZINB  | 2.150 | 0.890     |         52 |

Table 3: Median elapsed time (seconds) over five fits, and number of
iterations of the packageâ€™s algorithm (outer iterations for NB, EM
iterations for ZINB).

Two further checks support the implementation. First, the observed-data
log-likelihood never decreased across EM iterations in any ZINB fit
(monotone for all fits: TRUE), as theory requires. Second, the analytic
Hessian agrees with the numerical Hessian from `numDeriv` to a maximum
relative difference of 1.1e-09 on a ZINB fit to bioChemists with
different count and zero covariates.

# 5 Simulation study

To check that `zinbreg()` recovers known parameters, $`R = 200`$ data
sets of size $`n = 500`$ were simulated from a ZINB model with one
covariate in each part: $`\log\mu_i = 1 + 0.5x_i`$,
$`\mathrm{logit}\,\pi_i = -1 + 0.8z_i`$, and $`\theta = 2`$, with
$`x_i, z_i \sim N(0,1)`$ independent. The seed was fixed (651). For each
parameter, [Table 4](#tbl-sim) reports the mean estimate, the bias, the
Monte Carlo standard deviation of the estimates, the mean of the
estimated standard errors, and the empirical coverage of 95% Wald
intervals. The parameter $`\log\theta`$ is reported because the Wald
interval is constructed on that scale.

| Parameter       | Truth  | Mean   | Bias   | MC SD | Mean SE | Coverage |
|:----------------|:-------|:-------|:-------|:------|:--------|:---------|
| Count intercept | 1.000  | 0.998  | -0.002 | 0.068 | 0.065   | 0.915    |
| Count x         | 0.500  | 0.502  | 0.002  | 0.057 | 0.054   | 0.945    |
| Zero intercept  | -1.000 | -1.025 | -0.025 | 0.235 | 0.221   | 0.960    |
| Zero z          | 0.800  | 0.831  | 0.031  | 0.198 | 0.182   | 0.940    |
| Log theta       | 0.693  | 0.728  | 0.034  | 0.214 | 0.214   | 0.945    |

Table 4: Simulation results (R = 200, n = 500): bias, Monte Carlo SD of
the estimates, mean estimated standard error, and coverage of 95% Wald
intervals.

All 200 of the 200 fits converged. Biases are small relative to their
Monte Carlo standard errors (the largest absolute bias is 0.034), and
the mean estimated standard errors are close to the Monte Carlo standard
deviations, though slightly smaller for some parameters, as is typical
of Wald standard errors at moderate sample sizes. Coverage ranges from
0.915 to 0.960. With 200 replicates the Monte Carlo error of a coverage
proportion is about 0.015, so the lowest value is within sampling error
of the nominal level, although it hints that the intervals for the count
intercept are a little liberal.

# 6 Data analysis: physician office visits in NMES1988

## 6.1 Data and models

The `NMES1988` data contain the number of physician office visits
(`visits`) for 4,406 individuals aged 66 and over in the 1987â€“88 US
National Medical Expenditure Survey (Deb and Trivedi, 1997; Zeileis,
Kleiber and Jackman, 2008). The covariates are the number of hospital
stays (`hospital`), self-perceived health (`health`: poor, average,
excellent; average is the reference), the number of chronic conditions
(`chronic`), `gender`, years of education (`school`), and private
insurance (`insurance`). A Poisson model (fitted with `glm()` for
reference), an NB model, and a ZINB model were fitted. All three use
`hospital`, `health`, `chronic`, `gender`, `school`, and `insurance` in
the count part; the ZINB zero part uses `hospital`, `chronic`, `gender`,
`school`, and `insurance`, the minimum set suggested in the project
instructions (`health` was left out of the zero part to keep it to the
suggested set).

## 6.2 Overdispersion and excess zeros

The mean number of visits is 5.77 while the variance is 45.69, a
variance-to-mean ratio of 7.9. After adjusting for the covariates the
overdispersion remains: the Pearson chi-square statistic divided by the
residual degrees of freedom from the Poisson model is 6.71, far above
the value of one expected under a correct Poisson model. Zeros are also
more common than the Poisson model predicts: 683 of the 4406 individuals
(15.5%) had no visits, whereas the Poisson model expects only 47.
[Figure 1](#fig-freq) compares the observed frequencies of small counts
with those expected from the three models.

![](report_files/figure-html/fig-freq-1.png)

Figure 1: Observed and expected frequencies of physician visits (0 to
15) under the Poisson, NB, and ZINB models.

## 6.3 Model comparison

| Model   | Parameters | Log-lik. | AIC     | BIC     | Expected zeros |
|:--------|-----------:|:---------|:--------|:--------|:---------------|
| Poisson |          8 | -17971.6 | 35959.2 | 36010.4 | 46.7           |
| NB      |          9 | -12170.6 | 24359.1 | 24416.6 | 608.0          |
| ZINB    |         15 | -12090.7 | 24211.4 | 24307.3 | 708.9          |

Table 5: Model comparison for NMES1988: log-likelihood, AIC, BIC, and
expected number of zeros (observed: 683).

[Table 5](#tbl-ic) shows that the NB model improves on the Poisson model
by thousands of log-likelihood units, confirming the overdispersion. The
ZINB model improves on the NB model by a further 79.8 log-likelihood
units for six additional parameters; AIC falls by 147.7 and BIC by
109.3, so both criteria favor ZINB. The NB model underpredicts the
number of zeros (608 expected against 683 observed), and the ZINB model
slightly overpredicts it (709), but is much closer to the observed value
in relative terms (an error of 26 zeros against 75 for NB).

## 6.4 Estimates and interpretation

[Table 6](#tbl-irr) reports incidence rate ratios (IRRs,
$`e^{\hat\beta}`$) from the NB model and from the count part of the ZINB
model, and [Table 7](#tbl-orr) reports odds ratios (ORs,
$`e^{\hat\gamma}`$) for the zero part of the ZINB model, all with 95%
Wald confidence intervals.

| Covariate | NB IRR (95% CI) | ZINB count IRR (95% CI) |
|:---|:---|:---|
| Intercept | 2.53 (2.28, 2.82) | 3.30 (2.95, 3.69) |
| Hospital stays (per stay) | 1.24 (1.20, 1.29) | 1.22 (1.18, 1.27) |
| Health: poor vs. average | 1.36 (1.23, 1.49) | 1.33 (1.22, 1.45) |
| Health: excellent vs. average | 0.71 (0.63, 0.80) | 0.73 (0.65, 0.82) |
| Chronic conditions (per condition) | 1.19 (1.16, 1.22) | 1.14 (1.11, 1.16) |
| Male vs. female | 0.88 (0.83, 0.94) | 0.92 (0.87, 0.98) |
| Education (per year) | 1.03 (1.02, 1.04) | 1.02 (1.01, 1.03) |
| Private insurance: yes vs. no | 1.25 (1.16, 1.35) | 1.13 (1.05, 1.23) |

Table 6: Incidence rate ratios (95% Wald CI) for the NB model and the
count part of the ZINB model.

| Covariate                          | OR (95% CI)       | p-value  |
|:-----------------------------------|:------------------|:---------|
| Intercept                          | 0.95 (0.56, 1.62) | 0.8614   |
| Hospital stays (per stay)          | 0.45 (0.20, 1.02) | 0.0571   |
| Chronic conditions (per condition) | 0.29 (0.20, 0.41) | \<0.0001 |
| Male vs. female                    | 1.91 (1.29, 2.83) | 0.0012   |
| Education (per year)               | 0.92 (0.87, 0.97) | 0.0014   |
| Private insurance: yes vs. no      | 0.31 (0.20, 0.48) | \<0.0001 |

Table 7: Odds ratios (95% Wald CI) of being a structural zero, ZINB zero
part.

**Count part (NB and ZINB).** The two models agree qualitatively.
Holding the other covariates fixed, each additional chronic condition is
associated with 19% more expected visits in the NB model (IRR 1.19) and
14% more in the count part of the ZINB model (IRR 1.14). Each additional
hospital stay is associated with about 22% more visits (ZINB IRR 1.22).
People in poor health are expected to make 33% more visits than people
in average health (IRR 1.33), and people in excellent health 27% fewer
(IRR 0.73). Private insurance is associated with 13% more visits (IRR
1.13), and each year of education with about 2% more. Men are expected
to make 8% fewer visits than women in the ZINB count part (IRR 0.92), a
smaller difference than in the NB model (IRR 0.88). The IRRs of the ZINB
count part are generally closer to one than those of the NB model for
chronic conditions, gender, and insurance, because part of what the NB
model attributes to those covariates is absorbed by the zero part.

**Zero part (ZINB).** The zero part describes the odds of belonging to
the group that has no visits regardless of circumstances. Each
additional chronic condition multiplies the odds of being a structural
zero by 0.29, a reduction of 71%. Private insurance multiplies the odds
by 0.31 (69% lower), and each year of education by 0.92. Men have 1.91
times the odds of women of being structural zeros (91% higher). The
effect of hospital stays is borderline (OR 0.45, $`p =`$ 0.0571). In
short, healthier, uninsured, less educated men are the most likely to be
non-users of physician office care, and among users, sicker people visit
more often. Because these odds ratios come from a latent class, they
should be read as statements about a modeled non-user group, not as
directly observed behavior.

## 6.5 Which model do we prefer?

The ZINB model is the better description of these data. It has much
lower AIC and BIC than the NB model, it reproduces the number of zeros
much better, and its zero-part effects (chronic conditions, insurance,
gender, education) are strong and plausible rather than weakly
identified. The estimated dispersion parameter also rises from
$`\hat\theta =`$ 1.21 to 1.48 (larger $`\theta`$ means less
overdispersion), which means that some of the apparent overdispersion in
the NB fit was really excess zeros. The cost is complexity: the model
has six more parameters and its structural-zero interpretation is a
modeling device that the data alone cannot confirm. The slight
overprediction of zeros by the ZINB model and the large remaining
dispersion ($`\hat\theta`$ is still small) suggest that a hurdle model
or a model with additional covariates in the zero part could be
explored, but within the NB family the ZINB model is clearly preferred
here.

# 7 Discussion

The package reproduces `glm.nb()` and `zeroinfl()` to well within the
required tolerances using only IRLS, Newtonâ€“Raphson, and EM, and the
simulation shows that it recovers known parameters. Several limitations
are worth noting. EM converges linearly, so ZINB fits take more
iterations and more time than the quasi-Newton optimizer used by
`zeroinfl()`, particularly when the zero part is weakly identified. When
the data support no zero inflation, the maximum likelihood estimate lies
on the boundary $`\pi \to 0`$; the zero-part intercept then diverges to
$`-\infty`$ and its standard error is meaningless, as happens for
`art ~ ment | 1` on bioChemists, where the package and `zeroinfl()`
agree on everything except the arbitrary value reached by that
intercept. The likelihood of a mixture can have several local maxima,
and the package uses a single starting strategy. Standard errors are
Wald-type and treat the model as correctly specified, and the NB
standard errors treat $`\theta`$ as known. The package supports offsets
only in the count part, has no user-facing prior weights, and uses only
the log link for the count model and the logit link for the zero model.
Natural extensions are zero-inflated Poisson and hurdle models,
likelihood-ratio or Vuong tests for the zero part, profile-likelihood or
bootstrap intervals, observation weights, and a compiled implementation
of the inner loops.

# 8 AI-use statement

I used Claude (Anthropic) as a programming and writing assistant for
this project. I asked it to implement `nbreg()` and `zinbreg()`
according to the equations in Sections 2.1â€“3.2 of the project
description (with the rule that `optim()`, `glm()` and the other
functions in Table 2 are prohibited, recorded in `CLAUDE.md`), to write
the S3 methods, roxygen documentation, testthat tests, benchmark and
simulation scripts, and to draft this report and the CI workflow. I also
asked Claude to explain the algorithms in detail (the IRLS weights and
working response, the Newtonâ€“Raphson update for $`\log\theta`$, the EM
E- and M-steps, and the fractional-response logistic regression), and I
used these explanations to prepare for the oral walkthrough. Claude
additionally derived the analytic Hessian used for the ZINB standard
errors.

# References

- Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of
  Count Data* (2nd ed.). Cambridge University Press.
- Deb, P. and Trivedi, P. K. (1997). Demand for medical care by the
  elderly: A finite mixture approach. *Journal of Applied Econometrics*,
  12(3), 313â€“336.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum
  likelihood from incomplete data via the EM algorithm. *Journal of the
  Royal Statistical Society B*, 39(1), 1â€“38.
- Lambert, D. (1992). Zero-inflated Poisson regression, with an
  application to defects in manufacturing. *Technometrics*, 34(1),
  1â€“14.
- Long, J. S. (1997). *Regression Models for Categorical and Limited
  Dependent Variables*. Sage.
- McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*
  (2nd ed.). Chapman & Hall.
- Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
  with S* (4th ed.). Springer.
- Wickham, H. and Bryan, J. (2023). *R Packages* (2nd ed.). Oâ€™Reilly.
  <https://r-pkgs.org>
- Zeileis, A., Kleiber, C. and Jackman, S. (2008). Regression models for
  count data in R. *Journal of Statistical Software*, 27(8), 1â€“25.

# Appendix A: Side-by-side estimates

For each data set and model, the estimates and standard errors from
`zinbTruitt` and from the reference function, with absolute differences.
For NB the last row is $`\theta`$; for ZINB it is $`\log\theta`$.

| Term        | Ours     | Reference | Abs. diff. | SE ours | SE ref. | SE abs. diff. |
|:------------|:---------|:----------|:-----------|:--------|:--------|:--------------|
| (Intercept) | 0.25614  | 0.25614   | 7.8e-09    | 0.13735 | 0.13735 | 8.2e-11       |
| femWomen    | -0.21642 | -0.21642  | 1.6e-09    | 0.07264 | 0.07264 | 2.7e-11       |
| marMarried  | 0.15049  | 0.15049   | 1.0e-09    | 0.08210 | 0.08210 | 4.1e-11       |
| kid5        | -0.17642 | -0.17642  | 1.8e-09    | 0.05281 | 0.05281 | 6.5e-12       |
| phd         | 0.01527  | 0.01527   | 1.9e-10    | 0.03587 | 0.03587 | 1.4e-11       |
| ment        | 0.02908  | 0.02908   | 6.8e-10    | 0.00321 | 0.00321 | 6.1e-12       |
| theta       | 2.26439  | 2.26439   | 3.0e-09    | 0.27118 | 0.27118 | 1.1e-06       |

bioChemists NB {.caption-top}

| Term               | Ours     | Reference | Abs. diff. | SE ours | SE ref. | SE abs. diff. |
|:-------------------|:---------|:----------|:-----------|:--------|:--------|:--------------|
| count\_(Intercept) | 0.41675  | 0.41675   | 9.3e-07    | 0.14360 | 0.14360 | 3.3e-07       |
| count_femWomen     | -0.19551 | -0.19551  | 8.7e-07    | 0.07559 | 0.07559 | 4.8e-08       |
| count_marMarried   | 0.09758  | 0.09758   | 5.1e-07    | 0.08445 | 0.08445 | 9.7e-09       |
| count_kid5         | -0.15173 | -0.15173  | 6.4e-07    | 0.05421 | 0.05421 | 8.9e-10       |
| count_phd          | -0.00070 | -0.00070  | 2.4e-07    | 0.03627 | 0.03627 | 7.8e-08       |
| count_ment         | 0.02479  | 0.02479   | 1.2e-07    | 0.00349 | 0.00349 | 2.5e-07       |
| zero\_(Intercept)  | -0.19174 | -0.19169  | 5.6e-05    | 1.32278 | 1.32282 | 3.5e-05       |
| zero_femWomen      | 0.63594  | 0.63593   | 4.4e-06    | 0.84889 | 0.84892 | 2.5e-05       |
| zero_marMarried    | -1.49941 | -1.49947  | 5.8e-05    | 0.93864 | 0.93867 | 3.2e-05       |
| zero_kid5          | 0.62842  | 0.62843   | 8.0e-06    | 0.44277 | 0.44278 | 9.4e-06       |
| zero_phd           | -0.03770 | -0.03771  | 1.2e-05    | 0.30800 | 0.30801 | 7.9e-06       |
| zero_ment          | -0.88226 | -0.88229  | 2.9e-05    | 0.31621 | 0.31623 | 1.4e-05       |
| Log(theta)         | 0.97636  | 0.97636   | 4.9e-06    | 0.13547 | 0.13547 | 1.1e-06       |

bioChemists ZINB {.caption-top}

| Term            | Ours     | Reference | Abs. diff. | SE ours | SE ref. | SE abs. diff. |
|:----------------|:---------|:----------|:-----------|:--------|:--------|:--------------|
| (Intercept)     | 0.92926  | 0.92926   | 6.2e-08    | 0.05459 | 0.05459 | 1.5e-10       |
| hospital        | 0.21777  | 0.21777   | 3.5e-08    | 0.02018 | 0.02018 | 6.5e-11       |
| healthpoor      | 0.30501  | 0.30501   | 3.7e-08    | 0.04851 | 0.04851 | 8.6e-11       |
| healthexcellent | -0.34181 | -0.34181  | 9.4e-09    | 0.06092 | 0.06092 | 1.7e-10       |
| chronic         | 0.17492  | 0.17492   | 1.7e-08    | 0.01209 | 0.01209 | 1.3e-12       |
| gendermale      | -0.12649 | -0.12649  | 5.3e-09    | 0.03122 | 0.03122 | 4.5e-11       |
| school          | 0.02682  | 0.02682   | 9.6e-10    | 0.00439 | 0.00439 | 7.1e-12       |
| insuranceyes    | 0.22440  | 0.22440   | 1.5e-08    | 0.03946 | 0.03946 | 8.2e-11       |
| theta           | 1.20660  | 1.20660   | 1.3e-09    | 0.03356 | 0.03355 | 4.2e-07       |

NMES1988 NB {.caption-top}

| Term | Ours | Reference | Abs. diff. | SE ours | SE ref. | SE abs. diff. |
|:---|:---|:---|:---|:---|:---|:---|
| count\_(Intercept) | 1.19372 | 1.19371 | 3.8e-06 | 0.05666 | 0.05666 | 4.1e-08 |
| count_hospital | 0.20148 | 0.20148 | 2.8e-08 | 0.02036 | 0.02036 | 1.7e-08 |
| count_healthpoor | 0.28513 | 0.28513 | 2.3e-08 | 0.04509 | 0.04509 | 5.2e-08 |
| count_healthexcellent | -0.31934 | -0.31934 | 3.8e-08 | 0.06040 | 0.06040 | 5.1e-08 |
| count_chronic | 0.12900 | 0.12900 | 5.3e-07 | 0.01193 | 0.01193 | 1.7e-08 |
| count_gendermale | -0.08028 | -0.08028 | 7.0e-07 | 0.03102 | 0.03102 | 4.9e-08 |
| count_school | 0.02142 | 0.02142 | 9.5e-09 | 0.00436 | 0.00436 | 1.6e-08 |
| count_insuranceyes | 0.12586 | 0.12586 | 1.9e-06 | 0.04159 | 0.04159 | 6.5e-08 |
| zero\_(Intercept) | -0.04689 | -0.04688 | 9.9e-06 | 0.26856 | 0.26856 | 1.2e-07 |
| zero_hospital | -0.80045 | -0.80048 | 3.3e-05 | 0.42079 | 0.42082 | 2.2e-05 |
| zero_chronic | -1.24791 | -1.24794 | 3.1e-05 | 0.17831 | 0.17832 | 7.3e-06 |
| zero_gendermale | 0.64768 | 0.64768 | 3.0e-06 | 0.20011 | 0.20011 | 3.6e-06 |
| zero_school | -0.08377 | -0.08377 | 1.5e-06 | 0.02625 | 0.02625 | 3.6e-07 |
| zero_insuranceyes | -1.17561 | -1.17562 | 4.8e-06 | 0.22012 | 0.22012 | 5.1e-06 |
| Log(theta) | 0.39415 | 0.39414 | 5.2e-06 | 0.03504 | 0.03504 | 2.8e-08 |

NMES1988 ZINB {.caption-top}

# Appendix B: Session information

    R version 4.5.2 (2025-10-31 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8 
    [2] LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] knitr_1.50            pscl_1.5.9            MASS_7.3-65          
    [4] zinbTruitt_0.0.0.9000

    loaded via a namespace (and not attached):
     [1] digest_0.6.39       Formula_1.2-5       numDeriv_2016.8-1.1
     [4] fastmap_1.2.0       xfun_0.52           htmltools_0.5.8.1  
     [7] rmarkdown_2.29      cli_3.6.5           compiler_4.5.2     
    [10] rstudioapi_0.18.0   tools_4.5.2         evaluate_1.0.4     
    [13] yaml_2.3.10         rlang_1.3.0         jsonlite_2.0.0     
