# Certificates for Part V, chapter 15 (writing a distribution and a link).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   genpois   the probability function written here in (theta, lambda), summed
#             over the data at the fitted parameters, is the fit's
#             log-likelihood and glmmTMB's; the coefficients are glmmTMB's;
#             the fallback mean and variance are mu and phi * mu; every check
#             passes and the derivative rows are skipped; wool B has the
#             smaller dispersion; the AIC of the two count models differ by
#             less than 4; the closed score agrees with numDeriv's gradient of
#             the log-probability and with the finite differences; the fit
#             with the score reaches the same log-likelihood in less time;
#             ml(marginal = "all") reproduces glmmTMB with a random intercept
#             and ml() differs from it slightly;
#   burr      the density is gamlss.dist's GB2 with nu = 1, so the linear fit's
#             log-likelihood is recomputed from dGB2; statmod's is above
#             gamlss's by less than one; the Burr cAIC is below the gamma's by
#             less than ten; the closed cdf is the integral of the density and
#             the closed quantile inverts it; the closed score agrees with
#             numDeriv; every check passes, and the closed-form table has the
#             score's row; the two smooth fits reach the same log-likelihood;
#             the quantile residuals have a standard deviation within 0.05 of
#             one and the closed cdf computes them faster;
#   link      alpha = 1 is the logit link and a small alpha the complementary
#             log-log link; alpha_hat is between 1 and 2 and maximizes the
#             profile; the fit is glm()'s with the hand-written link; the
#             logit and cloglog fits are glm()'s; the logit log-likelihood is
#             below the Aranda-Ordaz one; no order is analytic before
#             dlinkinv() and the first inverse order is after, and it passes.

assert_own_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Chapter 'Writing a distribution and a link': ", what,
         " no longer agrees with the package.", call. = FALSE)
  }
  ll <- function(f) as.numeric(stats::logLik(f))
  checks_pass <- function(ck) all(ck$status == "OK")

  # --- genpois -------------------------------------------------------------------
  wb <- data$warpbreaks
  gp_lpmf <- function(y, mu, phi) {
    th <- mu / sqrt(phi)
    la <- 1 - 1 / sqrt(phi)
    log(th) + (y - 1) * log(th + la * y) - th - la * y - lgamma(y + 1)
  }
  X <- stats::model.matrix(~ wool * tension, wb)
  Z <- stats::model.matrix(~ wool, wb)
  ll_gp <- tryCatch({
    mu <- exp(drop(X %*% fits$gp@coefficients$mu))
    phi <- exp(drop(Z %*% fits$gp@coefficients$phi))
    sum(gp_lpmf(wb$breaks, mu, phi))
  }, error = function(e) NA_real_)
  if (!is.finite(ll_gp)) fail("the generalized Poisson log-likelihood")
  if (abs(ll_gp - ll(fits$gp)) > 1e-6) fail("the generalized Poisson log-likelihood")
  if (abs(ll(fits$gp) - ll(ref$tmb_gp)) > 1e-4) fail("the generalized Poisson against glmmTMB")
  cf_tmb <- c(glmmTMB::fixef(ref$tmb_gp)$cond, glmmTMB::fixef(ref$tmb_gp)$disp)
  cf_sm <- c(fits$gp@coefficients$mu, fits$gp@coefficients$phi)
  if (max(abs(unname(cf_tmb) - unname(cf_sm))) > 1e-3) {
    fail("the generalized Poisson coefficients against glmmTMB")
  }
  yy <- 0:400
  pp <- exp(gp_lpmf(yy, 4, 2.5))
  m1 <- sum(yy * pp)
  v1 <- sum(yy^2 * pp) - m1^2
  if (abs(m1 - 4) > 1e-8 || abs(v1 - 10) > 1e-6) fail("the generalized Poisson moments by direct summation")
  if (abs(mean(claims$gp, theta = list(mu = 4, phi = 2.5)) - 4) > 1e-6 ||
      abs(variance(claims$gp, theta = list(mu = 4, phi = 2.5)) - 10) > 1e-5) {
    fail("the fallback mean and variance")
  }
  if (!checks_pass(claims$chk_gp)) fail("the checks of the generalized Poisson")
  sk <- attr(claims$chk_gp, "skipped")
  if (is.null(sk) || !NROW(sk)) fail("the skipped derivative rows")
  if (fits$gp@coefficients$phi[[2]] >= 0) fail("the dispersion of wool B")
  if (abs(stats::AIC(fits$gp) - stats::AIC(fits$gp_nb)) >= 4) fail("the AIC of the two count models")
  if (claims$score_gap > 1e-7) fail("the closed score against the finite differences")
  yv <- c(0, 1, 3, 7, 12)
  g_cl <- distrib_gradient(claims$gp, y = yv, theta = list(mu = 4, phi = 2.5))
  for (i in seq_along(yv)) {
    gn <- numDeriv::grad(function(p) gp_lpmf(yv[i], p[1], p[2]), c(4, 2.5))
    if (max(abs(c(g_cl$mu[i], g_cl$phi[i]) - gn)) > 1e-7) fail("the closed score against numDeriv")
  }
  if (abs(ll(fits$gp_s) - ll(fits$gp)) > 1e-6) fail("the fit with the score")
  if (claims$t[["gp_score"]] >= claims$t[["gp_density"]]) fail("the time of the fit with the score")
  tr <- ref$tmb_gp_r
  sd_tmb <- attr(glmmTMB::VarCorr(tr)$cond$tension, "stddev")
  if (abs(fits$gp_ra@coefficients$mu[[1]] - glmmTMB::fixef(tr)$cond[[1]]) > 1e-3 ||
      abs(hyper(fits$gp_ra)$estimate[[1]] - sd_tmb) > 1e-3 ||
      abs(as.numeric(stats::logLik(fits$gp_ra, type = "marginal")) - ll(tr)) > 1e-3) {
    fail("ml(marginal = \"all\") against glmmTMB")
  }
  d_int <- abs(fits$gp_r@coefficients$mu[[1]] - glmmTMB::fixef(tr)$cond[[1]])
  if (d_int < 1e-4 || d_int > 0.05) fail("the joint-mode intercept against glmmTMB")

  # --- burr ----------------------------------------------------------------------
  rt <- data$rent
  fl <- fits$burr_lin
  Xb <- stats::model.matrix(~ Fl + loc, rt)
  Zb <- stats::model.matrix(~ loc, rt)
  ll_gb2 <- tryCatch({
    mu_b <- exp(drop(Xb %*% fl@coefficients$mu))
    c_b <- exp(drop(Zb %*% fl@coefficients$c))
    k_b <- exp(fl@coefficients$k[[1]])
    sum(gamlss.dist::dGB2(rt$R, mu = mu_b, sigma = c_b, nu = 1, tau = k_b,
                          log = TRUE))
  }, error = function(e) NA_real_)
  if (!is.finite(ll_gb2) || abs(ll_gb2 - ll(fl)) > 1e-6) {
    fail("the Burr density against GB2")
  }
  dg <- ll(fl) - ll(ref$gl_burr)
  if (dg <= 0 || dg >= 1) fail("the Burr fit against gamlss")
  da <- summary(fits$rent_gamma)@aic - summary(fits$burr)@aic
  if (da <= 0 || da >= 10) fail("the cAIC of the Burr and the gamma")
  bd <- claims$bd
  th <- list(mu = 800, c = 3, k = 1.5)
  dens <- function(y) exp(log(3) + log(1.5) - log(800) + 2 * log(y / 800) -
                            2.5 * log1p((y / 800)^3))
  for (q in c(200, 800, 3000)) {
    Fi <- stats::integrate(dens, 0, q, rel.tol = 1e-12)$value
    if (abs(distrib_cdf(bd, q = q, theta = th) - Fi) > 1e-9) fail("the closed cdf")
  }
  pr <- c(0.01, 0.5, 0.99)
  if (max(abs(distrib_cdf(bd, q = distrib_quantile(bd, p = pr, theta = th),
                          theta = th) - pr)) > 1e-12) fail("the closed quantile")
  b_lpdf <- function(y, p) log(p[2]) + log(p[3]) - log(p[1]) +
    (p[2] - 1) * log(y / p[1]) - (p[3] + 1) * log1p((y / p[1])^p[2])
  yv <- c(300, 800, 2000)
  g_b <- distrib_gradient(bd, y = yv, theta = th)
  for (i in seq_along(yv)) {
    gn <- numDeriv::grad(function(p) b_lpdf(yv[i], p), c(800, 3, 1.5))
    gc <- c(g_b$mu[i], g_b$c[i], g_b$k[i])
    if (max(abs(gc - gn) / pmax(1e-3, abs(gn))) > 1e-6) fail("the closed Burr score")
  }
  if (!checks_pass(claims$chk_burr) || !checks_pass(claims$chk_burr_closed)) {
    fail("the checks of the Burr distribution")
  }
  if (!any(grepl("^gradient", claims$chk_burr_closed$check)) ||
      any(grepl("^gradient", claims$chk_burr$check))) {
    fail("the score's row in the check table")
  }
  if (abs(ll(fits$burr) - ll(fits$burr_c)) > 1e-3) fail("the two smooth Burr fits")
  if (abs(stats::sd(claims$r_burr) - 1) > 0.05 ||
      abs(stats::sd(claims$r_burr_c) - 1) > 0.05) {
    fail("the Burr quantile residuals")
  }
  if (claims$t[["resid_closed"]] >= claims$t[["resid_density"]]) fail("the time of the residuals")

  # --- link ----------------------------------------------------------------------
  eta <- seq(-4, 3, by = 0.5)
  if (max(abs(linkinv(ao_link(1), eta = eta) - stats::plogis(eta))) > 1e-14) {
    fail("alpha = 1 as the logit link")
  }
  if (max(abs(linkinv(ao_link(1e-7), eta = eta) - (1 - exp(-exp(eta))))) > 1e-5) {
    fail("a small alpha as the complementary log-log link")
  }
  ah <- claims$alpha_hat
  if (ah <= 1 || ah >= 2) fail("the estimate of alpha")
  pl <- claims$ao_loglik(c(ah))
  if (claims$ao_loglik(ah * 1.05) > pl || claims$ao_loglik(ah / 1.05) > pl) {
    fail("the maximum of the profile")
  }
  mn <- data$menarche
  if (max(abs(unname(fits$ao@coefficients$mu) - unname(stats::coef(ref$glm_ao)))) > 1e-5 ||
      abs(ll(fits$ao) - ll(ref$glm_ao)) > 1e-6) {
    fail("the Aranda-Ordaz fit against glm()")
  }
  g_lo <- stats::glm(cbind(Menarche, Total - Menarche) ~ Age,
                     family = stats::binomial(), data = mn)
  g_cl <- stats::glm(cbind(Menarche, Total - Menarche) ~ Age,
                     family = stats::binomial(link = "cloglog"), data = mn)
  if (abs(ll(fits$logit) - ll(g_lo)) > 1e-6 || abs(ll(fits$cloglog) - ll(g_cl)) > 1e-6) {
    fail("the logit and cloglog fits against glm()")
  }
  if (ll(fits$logit) >= ll(fits$ao)) fail("the logit against the Aranda-Ordaz link")
  ao0 <- attr(claims$chk_ao, "analytic_orders")
  ao1 <- attr(claims$chk_ao_d, "analytic_orders")
  if (is.null(ao0) || is.null(ao1) || any(unlist(ao0) != 0) ||
      ao1$inverse != 1 || ao1$forward != 0) {
    fail("the analytic orders of the link")
  }
  inv1 <- claims$chk_ao_d$inverse_link_derivatives
  if (is.null(inv1) || !isTRUE(inv1[["order_1"]]) || !all(is.na(inv1[-1])) ||
      !all(is.na(claims$chk_ao$inverse_link_derivatives))) {
    fail("the check of the analytic dlinkinv()")
  }
  invisible(TRUE)
}
