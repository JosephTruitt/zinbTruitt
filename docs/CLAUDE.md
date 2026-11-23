# Project rules (STA 651 Final Project: zinbTruitt)

R package implementing NB (`nbreg`) and ZINB (`zinbreg`) regression.

## Hard rules

- NEVER use these inside `R/`: optim, optimize, optimise, nlm, nlminb,
  uniroot, glm, glm.fit, lm, lm.fit, lm.wfit, MASS::glm.nb,
  MASS::theta.ml, pscl::zeroinfl, countreg, glmmTMB, gamlss, VGAM,
  optimx, maxLik, stats4::mle, bbmle.
- Allowed: solve, qr, chol, crossprod, %\*%, lgamma, digamma, trigamma,
  dnbinom, plogis, qlogis, model.frame, model.matrix, model.response,
  model.offset, Formula package, numDeriv::hessian (standard errors
  only).
- Reference functions (glm, glm.nb, zeroinfl) are allowed ONLY in
  tests/, inst/benchmark/, and the report.
- Estimation: IRLS for beta, Newton-Raphson on log(theta), EM with IRLS
  M-steps for ZINB.

## Workflow

- Small, frequent commits with clear messages.
- Check before submitting: grep R/ for prohibited functions;
  devtools::check() must have 0 errors and 0 warnings.
