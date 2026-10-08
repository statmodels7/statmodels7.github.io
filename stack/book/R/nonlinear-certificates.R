# Certificates for Part III, chapter 12 (nonlinear models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   nl        the Michaelis-Menten fit is nls()'s, estimates and standard
#             errors; the two sigmas are close and not equal; each
#             log-likelihood is a sum of dnorm, ours at our sigma and nls()'s
#             at the maximum likelihood sigma;
#   links     the log-link fit has the same log-likelihood and its
#             coefficients are the logarithms of the identity fit's; the
#             interval for K carried back is positive and not symmetric;
#             the gamma fit is glm()'s with the inverse link, phi is the
#             maximum likelihood estimate (maximized here by hand), and the
#             two log-likelihoods agree to four decimals;
#   covariates under ml() the log-likelihoods are nls()'s, for both models;
#             the coefficients map onto nls()'s; the likelihood ratio is
#             recomputed;
#   random    the restricted criterion of both packages is the Laplace
#             approximation with the Gauss-Newton matrix, written out here:
#             at our estimates it is our value and at nlme's it is nlme's,
#             and ours is the larger; under ML the two maxima agree; the
#             typical tree is the population curve written out; for
#             theophylline ours is again the larger criterion, the estimates
#             are close, and the average curve is a Gauss-Hermite rule
#             written out here;
#   lasso     no level is missing, the levels at zero are the ones the text
#             names, C and F are kept, the coefficients at zero satisfy the
#             lasso's optimality condition (certificate ratio at most 1), and
#             the random effect gives every enzyme a departure;
#   breakpoint the change sits between day 17 and day 18, its size is within
#             0.1 of the simulated 0.5, the fit certifies converged, the AIC
#             prefers it, the log-likelihood is a sum of dpois, and the sharp
#             jump is rejected with a message naming `smoothed`;
#   function  the function, the formula and nls() give the same estimates,
#             with and without the gradient.

assert_nonlinear_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part III chapter 12: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  # an empty or mismatched comparison is a failure, never a vacuous pass
  rel <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (!length(b) || length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
    max(abs(a - b)) / max(1, max(abs(b)))
  }
  ll <- function(f, type = "conditional") as.numeric(logLik(f, type = type))
  state <- function(f) statmod_certificate(f)$state

  # --- nl --------------------------------------------------------------------
  pu <- data$puro
  ns <- ref$nls_mm
  cf <- coef(fits$mm)$mu
  ci <- confint(fits$mm)
  if (rel(cf[c("nl.Vm", "nl.K")] / stats::coef(ns), c(1, 1)) > 1e-5 ||
      rel(ci[c("mu:nl.Vm", "mu:nl.K"), "se"] /
            summary(ns)$coefficients[, "Std. Error"], c(1, 1)) > 1e-4) {
    fail("the Michaelis-Menten fit against nls()")
  }
  sg <- exp(coef(fits$mm)$sigma[["(Intercept)"]])
  mu_h <- predict(fits$mm, what = "mu", newdata = pu)
  rss <- sum(stats::residuals(ns)^2)
  # both are sqrt(RSS / (n - p)): nls() in closed form at its coefficients,
  # statmod() to the tolerance of its search, which the text puts in the
  # fourth decimal place
  rss_s <- sum((pu$rate - mu_h)^2)
  if (abs(sqrt(rss_s / (nrow(pu) - 2)) / stats::sigma(ns) - 1) > 1e-7 ||
      abs(sg - stats::sigma(ns)) >= 1e-3 || abs(sg - stats::sigma(ns)) < 1e-5) {
    fail("the two estimates of sigma")
  }
  # the start: K searched, Vm the least squares solution at that K
  bt <- modelterms7::term_build(modelterms7::nl(~ Vm * conc / (K + conc)), pu)
  s0 <- modelterms7::term_coef_start(bt, target = pu$rate)
  z0 <- pu$conc / (s0[[2]] + pu$conc)
  if (abs(s0[[1]] - sum(z0 * pu$rate) / sum(z0^2)) > 1e-8 * abs(s0[[1]]) ||
      abs(s0[[2]] / cf[["nl.K"]] - 1) > 0.5) {
    fail("the starting values of nl()")
  }
  if (abs(ll(fits$mm) - sum(stats::dnorm(pu$rate, mu_h, sg, log = TRUE))) > 1e-8 ||
      abs(as.numeric(stats::logLik(ns)) -
            sum(stats::dnorm(pu$rate, stats::fitted(ns), sqrt(rss / nrow(pu)),
                             log = TRUE))) > 1e-8) {
    fail("the two log-likelihoods of the Michaelis-Menten fit")
  }

  # --- links -----------------------------------------------------------------
  cl <- coef(fits$mm_log)$mu
  if (rel(exp(cl) / cf, c(1, 1)) > 1e-5 ||
      abs(ll(fits$mm_log) - ll(fits$mm)) > 1e-7) {
    fail("the log-link fit against the identity fit")
  }
  ck <- exp(unlist(confint(fits$mm_log)["mu:nl.K", c("estimate", "lower", "upper")]))
  if (any(ck <= 0) || abs((ck[3] - ck[1]) - (ck[1] - ck[2])) < 1e-6) {
    fail("the interval for K on the log scale")
  }
  bg <- stats::coef(ref$glm_mm)
  cg <- coef(fits$mm_gam)$mu
  if (rel(c(cg[["nl.Vm"]], exp(cg[["nl.K"]])) / c(1 / bg[[1]], bg[[2]] / bg[[1]]),
          c(1, 1)) > 1e-5) {
    fail("the gamma fit against glm()")
  }
  mg <- predict(fits$mm_gam, what = "mu", newdata = pu)
  nll <- function(lp) -sum(stats::dgamma(pu$rate, shape = 1 / exp(lp),
                                         rate = 1 / (exp(lp) * mg), log = TRUE))
  op <- stats::optimize(nll, c(-10, 2), tol = 1e-12)
  if (rel(exp(coef(fits$mm_gam)$phi[["(Intercept)"]]), exp(op$minimum)) > 1e-5 ||
      abs(ll(fits$mm_gam) + op$objective) > 1e-8 ||
      abs(ll(fits$mm_gam) - as.numeric(stats::logLik(ref$glm_mm))) >= 5e-5) {
    fail("the maximum likelihood dispersion of the gamma fit")
  }

  # --- covariates ------------------------------------------------------------
  P <- data$Puromycin
  ns0 <- stats::nls(rate ~ Vm[state] * conc / (K + conc), data = P,
                    start = list(Vm = c(200, 160), K = 0.05))
  cs <- coef(fits$st)$mu
  if (abs(ll(fits$st) - as.numeric(stats::logLik(ref$nls_st))) > 1e-6 ||
      abs(ll(fits$st0) - as.numeric(stats::logLik(ns0))) > 1e-6 ||
      rel(c(cs[["nl.Vm.(Intercept)"]],
            cs[["nl.Vm.(Intercept)"]] + cs[["nl.Vm.stateuntreated"]],
            cs[["nl.K.(Intercept)"]],
            cs[["nl.K.(Intercept)"]] + cs[["nl.K.stateuntreated"]]) /
            stats::coef(ref$nls_st), rep(1, 4)) > 1e-5 ||
      rel(claims$lr_st, 2 * (ll(fits$st) - ll(fits$st0))) > 1e-12) {
    fail("the two states against nls()")
  }

  # --- random: Orange ----------------------------------------------------------
  og <- data$orange
  y <- og$circumference
  tr <- as.integer(og$Tree)
  mu_of <- function(beta, b) (beta[1] + b[tr]) / (1 + exp((beta[2] - og$age) / beta[3]))
  # the Laplace approximation with the Gauss-Newton matrix, the fixed
  # coefficients integrated with a flat prior (the restricted criterion)
  laplace <- function(tau, sig) {
    pl <- function(p) -(sum(stats::dnorm(y, mu_of(p[1:3], p[4:8]), sig, log = TRUE)) +
                          sum(stats::dnorm(p[4:8], 0, tau, log = TRUE)))
    o <- stats::optim(c(190, 720, 345, rep(0, 5)), pl, method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
    o <- stats::optim(o$par, pl, method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
    J <- numDeriv::jacobian(function(q) mu_of(q[1:3], q[4:8]), o$par)
    H <- crossprod(J) / sig^2 + diag(c(0, 0, 0, rep(1 / tau^2, 5)))
    -o$value - 0.5 * as.numeric(determinant(H)$modulus) + 4 * log(2 * pi)
  }
  vc <- as.numeric(nlme::VarCorr(ref$nlme_or)[, "StdDev"])
  at_ours <- laplace(hyper(fits$or)$estimate,
                     exp(coef(fits$or)$sigma[["(Intercept)"]]))
  at_nlme <- laplace(vc[1], vc[2])
  if (abs(at_ours - ll(fits$or, "marginal")) > 1e-3 ||
      abs(at_nlme - as.numeric(stats::logLik(ref$nlme_or))) > 2e-3 ||
      !(ll(fits$or, "marginal") > as.numeric(stats::logLik(ref$nlme_or)))) {
    fail("the restricted criterion of the orange trees")
  }
  l_ml <- ll(fits$or_ml, "marginal")
  l_nm <- as.numeric(stats::logLik(ref$nlme_or_ml))
  com <- coef(fits$or_ml)$mu[c("nl.Asym.(Intercept)", "nl.xmid", "nl.scal")]
  if (abs(l_ml - l_nm) >= 1e-3 || l_ml < l_nm ||
      max(abs(com / nlme::fixef(ref$nlme_or_ml) - 1)) >= 5e-3) {
    fail("the maximum likelihood fit of the orange trees against nlme()")
  }
  co <- coef(fits$or)$mu
  nd <- data.frame(age = c(500, 1000, 1500), Tree = factor(c("1", "1", "9")))
  pop <- co[["nl.Asym.(Intercept)"]] / (1 + exp((co[["nl.xmid"]] - nd$age) / co[["nl.scal"]]))
  if (rel(predict(fits$or, what = "mu", newdata = nd, random = "zero"), pop) > 1e-10) {
    fail("the curve of the typical tree")
  }

  # --- random: theophylline ----------------------------------------------------
  ct <- coef(fits$th)$mu
  if (!(ll(fits$th, "marginal") > as.numeric(stats::logLik(ref$nlme_th))) ||
      max(abs(ct[c("nl.lKe", "nl.lKa.(Intercept)", "nl.lCl.(Intercept)")] -
                nlme::fixef(ref$nlme_th))) > 0.05) {
    fail("the theophylline fit against nlme()")
  }
  tau <- hyper(fits$th)$estimate
  fol <- function(t, d, ua, uc) {
    d * exp(ct[["nl.lKe"]] + ct[["nl.lKa.(Intercept)"]] + ua - ct[["nl.lCl.(Intercept)"]] - uc) *
      (exp(-exp(ct[["nl.lKe"]]) * t) - exp(-exp(ct[["nl.lKa.(Intercept)"]] + ua) * t)) /
      (exp(ct[["nl.lKa.(Intercept)"]] + ua) - exp(ct[["nl.lKe"]]))
  }
  k <- 40L
  J <- matrix(0, k, k)
  i <- seq_len(k - 1L)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  x <- e$values * sqrt(2)
  w <- e$vectors[1L, ]^2
  nt <- data.frame(Time = c(1, 5, 15), Dose = 4.5, Subject = factor("1"))
  gh <- vapply(nt$Time, function(t)
    sum(outer(w, w) * outer(x * tau[1], x * tau[2],
                            function(a, b) fol(t, 4.5, a, b))), numeric(1))
  if (rel(predict(fits$th, what = "mu", newdata = nt, random = "marginal"), gh) > 1e-8 ||
      rel(predict(fits$th, what = "mu", newdata = nt, random = "zero"),
          fol(nt$Time, 4.5, 0, 0)) > 1e-10) {
    fail("the typical and the average theophylline curves")
  }

  # --- lasso -----------------------------------------------------------------
  dl <- coef(fits$las)$mu[paste0("nl.Vm.lasso.enzyme", LETTERS[1:8])]
  dr <- coef(fits$ran)$mu[paste0("nl.Vm.random.", LETTERS[1:8])]
  zs <- statmod_certificate(fits$las)$zeros
  if (anyNA(dl) || length(fits$las@aliased) ||
      !identical(LETTERS[1:8][dl == 0], claims$zero_las) ||
      any(dl[c(3L, 6L)] == 0) ||
      # C and F close to the simulated departures, the others kept small
      any(abs(dl[c(3L, 6L)] - (claims$vm_true[c(3L, 6L)] - 200)) > 15) ||
      any(abs(dl[-c(3L, 6L)]) > 10) ||
      !(zs[["ratio"]] <= 1) || zs[["n"]] != length(claims$zero_las) ||
      any(dr == 0)) {
    fail("the lasso on the levels of Vm")
  }

  # --- breakpoint ------------------------------------------------------------
  dc <- coef(fits$dec)$mu
  dd <- data$decay
  if (!(dc[["nl.r.jump.psi1"]] > 17 && dc[["nl.r.jump.psi1"]] < 18) ||
      abs(dc[["nl.r.jump.delta1"]] - 0.5) > 0.1 ||
      !identical(state(fits$dec), "converged") ||
      !(AIC(fits$dec) < AIC(fits$dec0)) ||
      abs(ll(fits$dec) - sum(stats::dpois(dd$y, predict(fits$dec, what = "mu",
                                                        newdata = dd), log = TRUE))) > 1e-8) {
    fail("the change of the rate of decay")
  }
  msg <- tryCatch({
    statmod(y ~ 0 + nl(~ a * exp(-r * x), r ~ jump(day),
                       links = list(a = log_link(), r = log_link()),
                       start = list(a = 50, r = 1)),
            distrib = poisson_distrib(link_mu = identity_link()), data = dd)
    ""
  }, error = function(e) conditionMessage(e))
  if (!grepl("smoothed", msg, fixed = TRUE)) fail("the rejection of a sharp jump inside nl()")

  # --- function --------------------------------------------------------------
  if (rel(coef(fits$fn)$mu, coef(fits$fm)$mu) > 1e-6 ||
      rel(coef(fits$fn)$mu, stats::coef(ref$nls_fn)) > 1e-5 ||
      rel(coef(fits$fg)$mu, coef(fits$fn)$mu) > 1e-6) {
    fail("the function, the formula and nls()")
  }
  invisible(TRUE)
}
