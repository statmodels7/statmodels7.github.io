# Certificates for Part II, chapter 4 (generalized linear models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   counts      the Poisson fit of `ships` is glm()'s: coefficients, variance
#               matrix, log-likelihood and AIC; the rate ratios are the
#               exponentials of the intervals; the Pearson statistic is the
#               sum of glm()'s squared Pearson residuals and its ratio is
#               above one;
#   binary      the binomial fit of `budworm`, written with cbind(), is
#               glm()'s, and so is each additive fit under the three links;
#               the probit link has the largest log-likelihood, within one
#               unit of the others; the interaction is not significant; the
#               Bernoulli fit of `birthwt` is glm()'s under four links, the
#               odds ratios are the exponentials of the intervals and the
#               four log-likelihoods lie within one unit;
#   positive    the gamma fit of GAGurine has glm()'s mean; the ml() estimate
#               of phi is the reciprocal of MASS::gamma.shape(); the default
#               estimate is that times n/(n - 2) to three digits; the ratio
#               of the two standard errors is the square root of the ratio of
#               the two dispersions; logLik(glm) is the gamma density at
#               deviance / n; the inverse gaussian AIC exceeds the gamma's
#               and the smooth's is below it;
#   counts 2    the Pearson ratio of the Poisson fit is glm()'s; the ML fit of
#               the negative binomial is glm.nb()'s and the REML estimate of
#               theta is the smaller; the negative binomial has the smallest
#               AIC of the five; the expected classes sum to the sample
#               size and the negative binomial expects more zeros than the
#               Poisson; the beta-binomial fit recovers log sigma, its
#               standard errors exceed the binomial's by about the square
#               root of the binomial Pearson ratio and its AIC is smaller;
#   dispersion  the likelihood ratios are the twice differences of the
#               marginal log-likelihoods, the gamma dispersion depends on the
#               age and the negative binomial one does not on Eth and Lrn;
#   links       the square-root fit is glm()'s, with a larger AIC than the
#               log link; the user-written link inverts, has no analytic
#               derivative, passes check_link() and fits.

assert_glm_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 4: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  mu_coef <- function(f) unname(coef(f)$mu)

  # --- counts: ships --------------------------------------------------------
  sh <- data$ships
  g <- ref$ships
  f <- fits$ships
  if (rel(mu_coef(f), unname(stats::coef(g))) > 1e-6) {
    fail("the Poisson coefficients of ships")
  }
  if (rel(unname(stats::vcov(f)), unname(stats::vcov(g))) > 1e-3) {
    fail("the variance matrix of ships")
  }
  if (abs(as.numeric(logLik(f)) - as.numeric(stats::logLik(g))) > 1e-6 ||
      abs(AIC(f) - stats::AIC(g)) > 1e-5) fail("the log-likelihood of ships")
  r <- claims$ratios
  ci <- confint(f)
  ci <- ci[ci$coefficient != "(Intercept)", ]
  if (rel(r$ratio, exp(ci$estimate)) > 1e-10 ||
      rel(r$lower, exp(ci$lower)) > 1e-10 ||
      rel(r$upper, exp(ci$upper)) > 1e-10) fail("the rate ratios")
  if (!all(r$lower < r$ratio & r$ratio < r$upper)) fail("the ratio intervals")
  x2 <- sum(stats::residuals(g, type = "pearson")^2)
  if (abs(claims$x2_ships / x2 - 1) > 1e-4 ||
      claims$df_ships != g$df.residual ||
      !(claims$x2_ships / claims$df_ships > 1)) fail("the Pearson ratio of ships")
  if (!(r["typeC", "ratio"] == min(r$ratio[1:4])) ||
      !(r["typeE", "ratio"] == max(r$ratio[1:4]))) fail("types C and E")

  # --- binary: budworm ------------------------------------------------------
  bw <- data$budworm
  g <- ref$bw
  f <- fits$bw
  if (rel(mu_coef(f), unname(stats::coef(g))) > 1e-5 ||
      abs(as.numeric(logLik(f)) - as.numeric(stats::logLik(g))) > 1e-6) {
    fail("the binomial fit of budworm against glm()")
  }
  pint <- summary(f)@tables$mu[[1]]$table$p_value[4]
  if (!(pint > 0.05)) fail("the insignificant interaction")
  add <- c(logit = "logit", probit = "probit", cloglog = "cloglog")
  add_fit <- list(logit = fits$bw_logit, probit = fits$bw_probit,
                  cloglog = fits$bw_cloglog)
  ll_bw <- numeric(3)
  for (k in seq_along(add)) {
    gk <- stats::glm(cbind(numdead, n - numdead) ~ sex + ldose,
                     family = stats::binomial(link = add[[k]]), data = bw)
    fk <- add_fit[[k]]
    if (rel(mu_coef(fk), unname(stats::coef(gk))) > 1e-5 ||
        abs(as.numeric(logLik(fk)) - as.numeric(stats::logLik(gk))) > 1e-6) {
      fail(paste("the", names(add)[k], "fit of budworm"))
    }
    ll_bw[k] <- as.numeric(logLik(fk))
  }
  lt <- claims$link_table
  if (rel(lt$logLik, ll_bw) > 1e-8 ||
      rel(lt$AIC, -2 * ll_bw + 2 * 3) > 1e-8) fail("the link table of budworm")
  if (which.max(ll_bw) != 2L || diff(range(ll_bw)) > 1) {
    fail("the probit link and the spread of the links")
  }

  # --- binary: birthwt ------------------------------------------------------
  bt <- data$birthwt
  g <- ref$bwt
  f <- fits$bwt
  if (rel(mu_coef(f), unname(stats::coef(g))) > 1e-5 ||
      abs(as.numeric(logLik(f)) - as.numeric(stats::logLik(g))) > 1e-6) {
    fail("the Bernoulli fit of birthwt")
  }
  ci <- confint(f)
  ci <- ci[ci$coefficient != "(Intercept)", ]
  o <- claims$odds
  if (rel(o$odds_ratio, exp(ci$estimate)) > 1e-10 ||
      rel(o$lower, exp(ci$lower)) > 1e-10 ||
      rel(o$upper, exp(ci$upper)) > 1e-10) fail("the odds ratios")
  bl <- c("probit", "cloglog", "cauchit")
  bf <- list(fits$bwt_probit, fits$bwt_cloglog, fits$bwt_cauchit)
  ll_b <- as.numeric(logLik(f))
  for (k in 1:3) {
    gk <- stats::glm(low ~ age + lwt + smoke + ht,
                     family = stats::binomial(link = bl[k]), data = bt)
    if (abs(as.numeric(logLik(bf[[k]])) - as.numeric(stats::logLik(gk))) > 1e-4) {
      fail(paste("the", bl[k], "fit of birthwt"))
    }
    ll_b <- c(ll_b, as.numeric(logLik(bf[[k]])))
  }
  if (diff(range(ll_b)) > 1 || rel(claims$bwt_links$logLik, ll_b) > 1e-8) {
    fail("the links of birthwt")
  }

  # --- positive: GAGurine ---------------------------------------------------
  ga <- data$gag
  n <- nrow(ga)
  g <- ref$gam
  f <- fits$gam
  if (rel(mu_coef(f), unname(stats::coef(g))) > 1e-5) {
    fail("the gamma coefficients")
  }
  phi <- claims$phi_table
  phi_ml <- 1 / MASS::gamma.shape(g)$alpha
  if (abs(phi["statmod, ml()", "phi"] / phi_ml - 1) > 1e-4 ||
      abs(phi["MASS::gamma.shape", "phi"] / phi_ml - 1) > 1e-10) {
    fail("the maximum likelihood dispersion")
  }
  if (abs(phi["statmod, default", "phi"] /
          (phi["statmod, ml()", "phi"] * n / (n - 2)) - 1) > 1e-3) {
    fail("the REML dispersion")
  }
  if (abs(phi["glm, Pearson", "phi"] / summary(g)$dispersion - 1) > 1e-10 ||
      abs(phi["glm, deviance / n", "phi"] / (g$deviance / n) - 1) > 1e-10) {
    fail("the dispersions of glm")
  }
  se_a <- sqrt(diag(stats::vcov(f)))[["mu:Age"]]
  se_g <- summary(g)$coefficients["Age", "Std. Error"]
  if (abs(se_g / se_a / sqrt(phi["glm, Pearson", "phi"] /
                             phi["statmod, default", "phi"]) - 1) > 1e-3) {
    fail("the ratio of the standard errors")
  }
  mu_g <- stats::fitted(g)
  ll_glm <- sum(stats::dgamma(ga$GAG, shape = n / g$deviance,
                              scale = mu_g * g$deviance / n, log = TRUE))
  if (abs(ll_glm - as.numeric(stats::logLik(g))) > 1e-6) {
    fail("the log-likelihood of glm for a gamma response")
  }
  if (!(as.numeric(logLik(fits$gam_ml)) >= as.numeric(logLik(f)) - 1e-8) ||
      !(as.numeric(logLik(fits$gam_ml)) > as.numeric(stats::logLik(g)))) {
    fail("the maximized gamma log-likelihood")
  }
  if (!(AIC(fits$ig) > AIC(f)) || !(AIC(fits$gam_s) < AIC(f) - 50)) {
    fail("the AIC of the inverse gaussian and of the smooth")
  }

  # --- counts that vary too much: quine -------------------------------------
  q <- data$quine
  gp <- stats::glm(Days ~ Eth + Sex + Age + Lrn, family = stats::poisson,
                   data = q)
  qf <- claims$quine_fits
  if (abs(qf$pearson[1] / (sum(stats::residuals(gp, "pearson")^2) /
                           gp$df.residual) - 1) > 1e-4) {
    fail("the Pearson ratio of the Poisson fit of quine")
  }
  if (!(which.min(qf$AIC) == 2L) || !all(qf$pearson[4:5] < 1) ||
      abs(qf$pearson[2] - 1) > 0.1) fail("the AIC and Pearson ratios of quine")
  if (rel(qf$AIC, c(AIC(fits$q_pois), AIC(fits$q_nb2), AIC(fits$q_nb1),
                    AIC(fits$q_pig1), AIC(fits$q_pig2))) > 1e-10) {
    fail("the AIC table of quine")
  }
  gn <- ref$q_nb
  mlf <- fits$q_nb2_ml
  th_ml <- exp(coef(mlf)$theta[[1]])
  th_re <- exp(coef(fits$q_nb2)$theta[[1]])
  if (abs(th_ml / gn$theta - 1) > 1e-4 || !(th_re < th_ml) ||
      rel(mu_coef(mlf), unname(stats::coef(gn))) > 1e-5 ||
      abs(as.numeric(logLik(mlf)) - as.numeric(stats::logLik(gn))) > 1e-4) {
    fail("the negative binomial fit against glm.nb()")
  }
  if (abs(sum(claims$observed) - nrow(q)) > 1e-9 ||
      abs(sum(claims$exp_pois) - nrow(q)) > 1e-6 ||
      abs(sum(claims$exp_nb2) - nrow(q)) > 1e-6 ||
      !(claims$exp_nb2[1] > 100 * claims$exp_pois[1])) {
    fail("the expected frequencies of quine")
  }
  sb <- data$sim_bb
  fb <- fits$b_bb
  fn <- fits$b_bin
  se_bb <- sqrt(diag(stats::vcov(fb)))[c("mu:(Intercept)", "mu:x")]
  se_bi <- sqrt(diag(stats::vcov(fn)))
  gbn <- stats::glm(cbind(y, 20 - y) ~ x, family = stats::binomial, data = sb)
  pr <- sum(stats::residuals(gbn, "pearson")^2) / gbn$df.residual
  ratio2 <- (se_bb / se_bi)^2
  if (abs(coef(fb)$sigma[[1]] - log(0.5)) > 2.5 *
      sqrt(diag(stats::vcov(fb)))[["sigma:(Intercept)"]] ||
      any(ratio2 < 0.7 * pr) || any(ratio2 > 1.3 * pr) ||
      !(AIC(fb) < AIC(fn) - 100)) fail("the beta-binomial fit")

  # --- a dispersion with covariates -----------------------------------------
  lr_g <- 2 * (as.numeric(logLik(fits$gam_phi, type = "marginal")) -
                 as.numeric(logLik(fits$gam, type = "marginal")))
  lr_q <- 2 * (as.numeric(logLik(fits$q_theta, type = "marginal")) -
                 as.numeric(logLik(fits$q_nb2, type = "marginal")))
  if (abs(lr_g - claims$lr_gam) > 1e-8 || abs(lr_q - claims$lr_q) > 1e-8) {
    fail("the likelihood ratios")
  }
  if (!(stats::pchisq(lr_g, 1, lower.tail = FALSE) < 0.01) ||
      !(stats::pchisq(lr_q, 2, lower.tail = FALSE) > 0.05) ||
      !(AIC(fits$gam_phi) < AIC(fits$gam)) ||
      !(AIC(fits$q_theta) > AIC(fits$q_nb2))) fail("the dispersion tests")
  if (!(coef(fits$gam_phi)$phi[["Age"]] > 0)) fail("the sign of the age effect")

  # --- links ----------------------------------------------------------------
  gs <- ref$q_sqrt
  if (rel(mu_coef(fits$q_sqrt), unname(stats::coef(gs))) > 1e-3 ||
      abs(AIC(fits$q_sqrt) - stats::AIC(gs)) > 1e-2 ||
      !(AIC(fits$q_sqrt) > AIC(fits$q_pois))) fail("the square-root link fit")
  lk <- linkfunctions7::sqrt_link()
  ev <- c(1, 2, 3)
  if (max(abs(linkfunctions7::linkinv(lk, eta = ev) - ev^2)) > 1e-12 ||
      max(abs(linkfunctions7::dlinkinv(lk, eta = ev) - 2 * ev)) > 1e-12) {
    fail("the square-root link")
  }
  sq <- claims$sq_logit
  th <- c(0.1, 0.5, 0.9)
  if (max(abs(linkfunctions7::linkinv(sq,
        eta = linkfunctions7::linkfun(sq, theta = th)) - th)) > 1e-12) {
    fail("the inverse of the user link")
  }
  fo <- linkfunctions7::link_fallback_orders(sq)
  if (fo$forward != 0 || fo$inverse != 0) fail("the orders of the user link")
  chk <- utils::capture.output(ck <- linkfunctions7::check_link(sq))
  if (!(ck$invertibility_theta && ck$invertibility_eta && ck$monotonicity &&
        ck$inverse_theorem)) fail("the checks of the user link")
  d_num <- linkfunctions7::dlinkinv(sq, eta = c(-1, 0, 1))
  e <- c(-1, 0, 1)
  d_exact <- 2 * stats::plogis(e) * stats::dlogis(e)
  if (max(abs(d_num - d_exact)) > 1e-5) fail("the numerical derivative of the user link")
  ll_sq <- as.numeric(logLik(fits$bw_sq))
  if (!is.finite(ll_sq) || abs(ll_sq - ll_bw[1]) > 1) {
    fail("the fit with the user link")
  }

  invisible(TRUE)
}
