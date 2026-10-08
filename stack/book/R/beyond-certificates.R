# Certificates for Part I, chapter 2 (beyond the basics).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   variances    V_b = (H + S)^{-1} and V_f = V_b H V_b with H the expected
#                information and S assembled here from penalties7; V_b and V_f
#                have a zero mu-sigma block and V_u does not; V_u - V_b has no
#                negative eigenvalue; the three predictions share the fit; the
#                frequentist standard error is below the Bayesian one at every
#                time and the unconditional one is not below it anywhere; the
#                coverage of the Bayesian intervals lies in [0.90, 0.97], above
#                the frequentist one and not above the unconditional one;
#   likelihood   the fixed effects, their standard errors, the two standard
#                deviations, the correlation and sigma of Orthodont are lme()'s;
#                the marginal log-likelihood is lme()'s with df 6 and the
#                conditional one is the sum of stats::dnorm at the fitted means;
#   tests        Wald from the estimate and se; the likelihood ratio is
#                2 (logLik(nls) - logLik(lm with K held)); the score statistic
#                is U' I^{-1} U with U by numDeriv on the gaussian log-likelihood
#                at the restricted ML fit and I the expected information there;
#                the gradient statistic is U_K (K_hat - 0.05); none rejects at
#                5%; the column r of summary(test = "lr") is the signed root;
#   invert       at the ends of the likelihood ratio interval the ML profile,
#                computed with lm() at K held, is the chi-squared quantile; the
#                Wald interval is symmetric and the likelihood ratio one is not;
#                the interval of nls() contains the likelihood ratio one, its
#                ends are where the profile t equals the t quantile, and
#                L = n log(1 + tau^2/(n - 2));
#   lr           the marginal likelihood ratio is anova()'s; the REML difference
#                is not it; statmod_test() on Sex gives a larger statistic;
#                the ML sd of the random intercepts is larger without Sex, and
#                statmod's ML sds are lme()'s;
#   random       the conditional prediction of patient 1 is exp(eta) with that
#                patient's effect; a new patient raises an error by default; the
#                marginal mean is the typical one times exp(tau^2/2); the group
#                interval is exp(eta +- z sqrt(se_eta^2 + tau^2)); the ends and
#                the median of the prediction interval are the quantiles of a
#                Poisson mixed over a gaussian predictor, by quadrature;
#   start        the fit from fit0 reaches the REML criterion of the fit from
#                the intercepts within 1e-3, with the same fitted values to
#                1e-2 g, in fewer evaluations;
#   outer        sigma is sqrt(RSS/(n - tau_mu)) under REML, sqrt(RSS/(n -
#                tau_mu + 2)) under ML and sqrt(RSS/n) at the joint mode, in
#                that order; the three lambdas agree within 10%; the held lambda
#                is 0.001, marked held and fixed, with fewer edf, and sigma is
#                still the REML one;
#   weights      frequency weights equal the expanded data; the mean equals
#                lm(weights); sigma of lm differs; the offset on sigma gives
#                lm's sigma and standard errors; the offset model is glm()'s;
#   fixed        the fixed distribution has mu and sigma only, its
#                log-likelihood is the sum of stats::dt, and its cAIC is above
#                the gaussian one.

assert_beyond_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part I chapter 2: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  z <- stats::qnorm(0.975)
  mc <- data$mcycle
  n <- nrow(mc)

  # --- variances -------------------------------------------------------------
  f1 <- fits$fit1
  sp <- f1@spec
  de <- statmod_design(sp)
  H <- -hessian(f1, expected = TRUE)
  offs <- cumsum(c(0L, vapply(de, function(e) e$npar, integer(1))))
  S <- matrix(0, nrow(H), nrow(H))
  nm <- character(0)
  for (k in seq_along(de)) {
    par <- names(de)[k]
    nm <- c(nm, paste0(par, ":", de[[k]]$coef_names))
    for (bk in names(de[[k]]$blocks)) {
      pen <- modelterms7::term_penalty(sp@terms[[par]][[bk]])
      if (is.null(pen)) next
      cols <- offs[k] + de[[k]]$blocks[[bk]]
      S[cols, cols] <- as.matrix(penalties7::penalty_hessian(
        pen, f1@coefficients[[par]][de[[k]]$blocks[[bk]]],
        list(lambda = f1@hyper[[par]][[bk]][["lambda"]])))
    }
  }
  Vb <- solve(H + S)
  Vf <- Vb %*% H %*% Vb
  dimnames(Vb) <- dimnames(Vf) <- list(nm, nm)
  vb <- stats::vcov(f1, type = "bayesian")
  vf <- stats::vcov(f1, type = "frequentist")
  vu <- stats::vcov(f1, type = "unconditional")
  if (rel(vb[nm, nm], Vb) > 1e-6) fail("the Bayesian variance")
  if (rel(vf[nm, nm], Vf) > 1e-6) fail("the frequentist variance")
  imu <- grep("^mu:", nm); isg <- grep("^sigma:", nm)
  if (any(vb[nm[imu], nm[isg]] != 0) || any(vf[nm[imu], nm[isg]] != 0) ||
      all(vu[nm[imu], nm[isg]] == 0)) {
    fail("the mu-sigma block of the three variances")
  }
  ev <- eigen(vu - vb, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8 * max(abs(ev), 1)) fail("the unconditional variance")
  pb <- claims$pb; pf <- claims$pf; pu <- claims$pu
  if (rel(pf$fit, pb$fit) > 1e-12 || rel(pu$fit, pb$fit) > 1e-12) {
    fail("the three predictions")
  }
  if (!(max(pf$se / pb$se) < 1)) fail("the frequentist standard errors")
  if (!(min(pu$se / pb$se) >= 1 - 1e-10)) fail("the unconditional standard errors")
  cv <- claims$cv
  if (cv[["bayesian"]] < 0.90 || cv[["bayesian"]] > 0.97 ||
      !(cv[["frequentist"]] < cv[["bayesian"]]) ||
      cv[["unconditional"]] < cv[["bayesian"]]) {
    fail("the coverage of the three intervals")
  }

  # --- likelihood --------------------------------------------------------------
  r1 <- fits$reml1; l1 <- ref$lme1
  if (rel(fits$reml1@coefficients$mu[1:2], nlme::fixef(l1)) > 1e-6) {
    fail("the fixed effects of Orthodont")
  }
  ci <- confint(r1, parm = c("mu:(Intercept)", "mu:I(age - 11)"))
  if (rel(ci$se, sqrt(diag(stats::vcov(l1)))) > 1e-5) {
    fail("the standard errors of the fixed effects")
  }
  sb <- summary(r1)@tables$mu[[2]]$table
  vc <- nlme::VarCorr(l1)
  lme_sd <- as.numeric(vc[1:2, "StdDev"])
  lme_cor <- as.numeric(vc[2, 3])
  if (rel(sb$estimate[1:2], lme_sd) > 1e-4 || abs(sb$estimate[3] - lme_cor) > 1e-3) {
    fail("the standard deviations and the correlation of Orthodont")
  }
  if (rel(exp(r1@coefficients$sigma), l1$sigma) > 1e-5) fail("sigma of Orthodont")
  lm_ <- logLik(r1, type = "marginal")
  if (rel(as.numeric(lm_), as.numeric(stats::logLik(l1))) > 1e-6 ||
      attr(lm_, "df") != 6) {
    fail("the marginal log-likelihood of Orthodont")
  }
  mu_o <- stats::fitted(r1)
  sg_o <- exp(r1@coefficients$sigma)
  lc <- sum(stats::dnorm(nlme::getResponse(l1), mu_o, sg_o, log = TRUE))
  if (rel(lc, as.numeric(logLik(r1))) > 1e-8) fail("the conditional log-likelihood")

  # --- tests -----------------------------------------------------------------
  pu_d <- data$puro
  np <- nrow(pu_d)
  nls_p <- ref$nls_p
  th_hat <- stats::coef(nls_p)
  tests <- claims$tests
  ciw <- claims$ciw
  if (rel(tests$wald@statistic, (ciw$estimate - 0.05)^2 / ciw$se^2) > 1e-8) {
    fail("the Wald statistic")
  }
  prof <- function(K) {
    X <- pu_d$conc / (K + pu_d$conc)
    m <- stats::lm(pu_d$rate ~ 0 + X)
    s2 <- sum(stats::resid(m)^2) / np
    list(ll = sum(stats::dnorm(pu_d$rate, stats::fitted(m), sqrt(s2), log = TRUE)),
         vm = unname(stats::coef(m)), s = sqrt(s2))
  }
  llmax <- as.numeric(stats::logLik(nls_p))
  if (rel(tests$lr@statistic, 2 * (llmax - prof(0.05)$ll)) > 1e-6) {
    fail("the likelihood ratio for K")
  }
  ll <- function(th) sum(stats::dnorm(pu_d$rate, th[1] * pu_d$conc / (th[2] + pu_d$conc),
                                      exp(th[3]), log = TRUE))
  pr05 <- prof(0.05)
  thr <- c(pr05$vm, 0.05, log(pr05$s))
  U <- numDeriv::grad(ll, thr)
  J <- cbind(pu_d$conc / (0.05 + pu_d$conc), -thr[1] * pu_d$conc / (0.05 + pu_d$conc)^2)
  I <- matrix(0, 3, 3)
  I[1:2, 1:2] <- crossprod(J) / exp(2 * thr[3]); I[3, 3] <- 2 * np
  if (rel(tests$score@statistic, drop(t(U) %*% solve(I, U))) > 1e-5) {
    fail("the score statistic")
  }
  if (rel(tests$gradient@statistic, U[2] * (th_hat[["K"]] - 0.05)) > 1e-3) {
    fail("the gradient statistic")
  }
  if (!all(vapply(tests, function(t) t@p.value > 0.05, logical(1)))) {
    fail("the four p-values")
  }
  sl <- summary(fits$puro, test = "lr")@tables$mu[[1]]$table
  t0 <- statmod_test(fits$puro, param = "mu", coefname = "nl.K", value = 0,
                     test = "lr")
  rK <- sl[grepl("nl.K", sl$name), ]
  sc <- setdiff(names(rK), c("name", "estimate", "se", "p_value", "lower", "upper"))
  if (!length(sc) || abs(abs(as.numeric(rK[[sc[1]]])) - sqrt(t0@statistic)) > 1e-6) {
    fail("the signed root of summary(test = \"lr\")")
  }

  # --- invert -------------------------------------------------------------------
  cil <- claims$cil
  kh <- ciw$estimate
  q <- stats::qchisq(0.95, 1)
  for (e in c(cil$lower, cil$upper)) {
    if (abs(2 * (llmax - prof(e)$ll) - q) > 1e-4) fail("the likelihood-ratio interval for K")
  }
  ws <- c(kh - ciw$lower, ciw$upper - kh)
  ls <- c(kh - cil$lower, cil$upper - kh)
  if (abs(ws[1] - ws[2]) > 1e-10 || abs(ls[1] - ls[2]) < 1e-3 * ls[1]) {
    fail("the symmetry of the two intervals")
  }
  cn <- suppressMessages(stats::confint(nls_p))["K", ]
  if (!(cn[1] < cil$lower && cn[2] > cil$upper)) fail("the interval of nls()")
  # the profile t of nls() at its own ends is the t quantile, and L is
  # n log(1 + tau^2 / (n - 2)) at every K
  s1 <- sum(stats::resid(nls_p)^2)
  tau_at <- function(K) {
    X <- pu_d$conc / (K + pu_d$conc)
    sqrt((sum(stats::resid(stats::lm(pu_d$rate ~ 0 + X))^2) - s1) / (s1 / (np - 2)))
  }
  qt_ <- stats::qt(0.975, np - 2)
  if (abs(tau_at(cn[1]) - qt_) > 1e-2 || abs(tau_at(cn[2]) - qt_) > 1e-2) {
    fail("the profile t of nls()")
  }
  for (K in c(0.05, cil$lower, cil$upper)) {
    L <- 2 * (llmax - prof(K)$ll)
    if (abs(L - np * log(1 + tau_at(K)^2 / (np - 2))) > 1e-6) {
      fail("the relation between L and the profile t")
    }
  }

  # --- lr ----------------------------------------------------------------------
  an <- stats::anova(ref$lme1_ml, ref$lme2_ml)
  if (abs(claims$lr - an$L.Ratio[2]) > 1e-4) fail("the likelihood ratio for Sex")
  dreml <- 2 * (as.numeric(logLik(fits$reml2, type = "marginal")) -
                  as.numeric(logLik(fits$reml1, type = "marginal")))
  if (abs(dreml - claims$lr) < 0.1) fail("the REML difference")
  if (!(claims$t_sex@statistic > claims$lr)) fail("the test on Sex")
  sd_int <- function(m) as.numeric(nlme::VarCorr(m)[1, "StdDev"])
  s1ml <- summary(fits$ml1)@tables$mu[[2]]$table$estimate[1]
  s2ml <- summary(fits$ml2)@tables$mu[[2]]$table$estimate[1]
  if (!(sd_int(ref$lme1_ml) > sd_int(ref$lme2_ml)) ||
      rel(c(s1ml, s2ml), c(sd_int(ref$lme1_ml), sd_int(ref$lme2_ml))) > 1e-3) {
    fail("the random-intercept sd under ML with and without Sex")
  }

  # --- random --------------------------------------------------------------------
  fe <- fits$epil
  tau <- hyper(fe)$estimate
  new <- data$new
  err <- tryCatch({ predict(fe, what = "mu", newdata = new); FALSE },
                  error = function(e) TRUE)
  if (!err) fail("the error for a new patient")
  cf <- fe@coefficients$mu
  cn_e <- statmod_design(fe@spec)$mu$coef_names
  lev1 <- match(MASS::epil$subject[1], sort(unique(MASS::epil$subject)))
  b1 <- cf[[match(paste0("random.", lev1), cn_e)]]
  rows <- MASS::epil[1:4, ]
  pk <- predict(fe, what = "mu", newdata = rows)
  rows$subject <- 0L
  p0 <- predict(fe, what = "mu", newdata = rows, random = "zero")
  if (rel(pk, p0 * exp(b1)) > 1e-8) fail("the prediction for a known patient")
  pz <- claims$pz; pm <- claims$pm; pg <- claims$pg; pr <- claims$pr
  if (rel(pm$fit, pz$fit * exp(tau^2 / 2)) > 1e-8) fail("the marginal mean")
  eta <- log(pz$fit); se_eta <- pz$se / pz$fit
  s_g <- sqrt(se_eta^2 + tau^2)
  if (rel(c(pg$lower, pg$upper), exp(eta + c(-1, 1) * z * s_g)) > 1e-6) {
    fail("the group interval")
  }
  cdf <- function(k) stats::integrate(function(u)
    stats::ppois(k, exp(u)) * stats::dnorm(u, eta, s_g), -Inf, Inf,
    rel.tol = 1e-10)$value
  qk <- function(p) { k <- 0; while (cdf(k) < p) k <- k + 1; k }
  if (!identical(as.numeric(c(pr$lower, pr$fit, pr$upper)),
                 as.numeric(c(qk(0.025), qk(0.5), qk(0.975))))) {
    fail("the prediction interval of the seizure count")
  }

  # --- start -----------------------------------------------------------------
  fs <- fits$fit1_s
  if (abs(fs@criterion - f1@criterion) > 1e-3 ||
      rel(as.numeric(logLik(fs)), as.numeric(logLik(f1))) > 1e-4 ||
      max(abs(stats::fitted(fs) - stats::fitted(f1))) > 1e-2 ||
      !(nrow(fs@history$outer) < nrow(f1@history$outer))) {
    fail("the fit started from fit0")
  }

  # --- outer -----------------------------------------------------------------
  s_of <- function(f) exp(f@coefficients$sigma[[1]])
  rss_of <- function(f) sum((mc$accel - stats::fitted(f))^2)
  tmu_of <- function(f) sum(f@edf$edf[f@edf$parameter == "mu"])
  f0 <- fits$fit0; fn <- fits$fit0_none; fm <- fits$fit0_ml; fh <- fits$fit0_held
  if (rel(s_of(f0), sqrt(rss_of(f0) / (n - tmu_of(f0)))) > 1e-5 ||
      rel(s_of(fm), sqrt(rss_of(fm) / (n - tmu_of(fm) + 2))) > 1e-5 ||
      rel(s_of(fn), sqrt(rss_of(fn) / n)) > 1e-6) {
    fail("the three estimates of sigma")
  }
  if (!(s_of(fn) < s_of(fm) && s_of(fm) < s_of(f0))) fail("the order of the three sigmas")
  lam <- c(hyper(f0)$estimate, hyper(fn)$estimate, hyper(fm)$estimate)
  if (max(lam) / min(lam) > 1.1) fail("the three smoothing parameters")
  hh <- hyper(fh)
  if (!isTRUE(hh$held) || hh$source != "fixed" || hh$estimate != 1e-3 ||
      !(fh@edf$edf[2] < f0@edf$edf[2]) ||
      rel(s_of(fh), sqrt(rss_of(fh) / (n - tmu_of(fh)))) > 1e-5) {
    fail("the held smoothing parameter")
  }

  # --- weights ---------------------------------------------------------------
  fw <- fits$fit_w; fb <- fits$fit_big; lw <- ref$lm_w; fp <- fits$fit_prec
  if (rel(unlist(fw@coefficients), unlist(fb@coefficients)) > 1e-6) {
    fail("frequency weights against the expanded data")
  }
  if (rel(fw@coefficients$mu, stats::coef(lw)) > 1e-8) fail("the weighted mean")
  w <- data$w
  sf <- sqrt(sum(w * stats::resid(lw)^2) / (sum(w) - 2))
  if (rel(s_of(fw), sf) > 1e-6 || rel(s_of(fw), summary(lw)$sigma) < 0.1) {
    fail("sigma under frequency weights")
  }
  cp <- confint(fp, parm = c("mu:(Intercept)", "mu:times"))
  if (rel(exp(fp@coefficients$sigma[[1]]), summary(lw)$sigma) > 1e-6 ||
      rel(cp$se, sqrt(diag(stats::vcov(lw)))) > 1e-6) {
    fail("precision weights as an offset")
  }

  # --- offset ----------------------------------------------------------------
  if (rel(fits$fit_i@coefficients$mu, stats::coef(ref$glm_i)) > 1e-6) {
    fail("the offset model")
  }

  # --- fixed -----------------------------------------------------------------
  ft <- fits$fit_t4
  if (!identical(ft@spec@distrib@params, c("mu", "sigma"))) fail("the fixed distribution")
  mu_t <- predict(ft, what = "mu", newdata = mc)
  sg_t <- predict(ft, what = "sigma", newdata = mc)
  lt <- sum(stats::dt((mc$accel - mu_t) / sg_t, df = 4, log = TRUE) - log(sg_t))
  if (rel(lt, as.numeric(logLik(ft))) > 1e-8) fail("the log-likelihood of the t model")
  if (!(AIC(ft) > AIC(f1))) fail("the cAIC of the t model")

  invisible(TRUE)
}
