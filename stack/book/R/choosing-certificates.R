# Certificates for Part IV, chapter 14 (choosing).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   counts      the ml() negative binomial is glm.nb()'s (log-likelihood and
#               theta) and gamlss's NBI; the Poisson-inverse gaussian is
#               gamlss's PIG; the Pearson ratio of the Poisson model is
#               recomputed from its fitted means and is above 10, those of
#               the three others are between 0.8 and 1.4; pig1 has the
#               smallest AIC of the four; the reml() theta is below the ml()
#               one; theta ~ log_lake lowers the AIC by more than the gap
#               between pig1 and negbin2; the negbin2 quantile residuals have
#               a standard deviation within 0.1 of one and the Poisson ones
#               above 2;
#   positive    the gamma model's log-likelihood is mgcv's within 0.1; gamma
#               and generalized gamma are within 5 of cAIC and below the
#               lognormal and the inverse gaussian by more than 50; the
#               lognormal and inverse gaussian quantile residuals have a
#               third moment below -0.3; the model of phi lowers the cAIC by
#               more than 20;
#   hyper       the REML log-likelihood is gaulss's within 0.05; sigma's
#               smoothing parameter is above 1e6 under every criterion; BIC
#               spends fewer edf than REML; the corrected cAIC is the naive
#               one plus twice the change of edf, and the change is
#               positive; the simulation: the corrected cAIC chooses the
#               random intercept less often in both settings; the lasso fits
#               estimate both smoothing parameters (neither left at its start
#               of 1), cross-validation chooses a lambda more than five times
#               BIC's, neither sets an attribute to zero, B1 is reported not
#               identified and the coefficients differ by less than 0.02; the
#               t beats the gaussian by more than 50 of cAIC, nu
#               is between 8 and 13 and the two informations give the same
#               log-likelihood within 0.01;
#   convergence the five states; the LAD fit within 0.1 of quantreg, its sum
#               of absolute residuals above the minimum by less than 0.01, its
#               flag FALSE and its state unknown; the budget fit is not
#               converged; the expected-information fit is converged with the
#               age slope of log nu negative; both fits from the joint mode
#               are converged, the expected one higher than the first by more
#               than 0.1, with the age slope positive;
#   limits      every refusal carries the phrase that names its cause; the
#               Wald coverage is below the nls() coverage.

assert_choosing_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Chapter 'Choosing': ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  ll <- function(f) as.numeric(stats::logLik(f))
  caic <- function(f) summary(f)@aic

  # --- counts ------------------------------------------------------------------
  sp <- data$species
  gnb <- MASS::glm.nb(fish ~ log(lake), data = sp)
  if (abs(ll(fits$sp_nb2) - as.numeric(stats::logLik(gnb))) > 1e-4 ||
      abs(exp(fits$sp_nb2@coefficients$theta[[1]]) - gnb$theta) > 1e-3) {
    fail("the negative binomial against glm.nb()")
  }
  if (abs(ll(fits$sp_nb2) - as.numeric(stats::logLik(ref$gl_nbi))) > 1e-3) {
    fail("the negative binomial against gamlss")
  }
  if (abs(ll(fits$sp_pig) - as.numeric(stats::logLik(ref$gl_pig))) > 1e-2) {
    fail("the Poisson-inverse gaussian against gamlss")
  }
  mu_p <- exp(fits$sp_pois@coefficients$mu[[1]] +
                fits$sp_pois@coefficients$mu[[2]] * log(sp$lake))
  x2 <- sum((sp$fish - mu_p)^2 / mu_p) / (nrow(sp) - 2)
  st <- claims$sp_table
  if (abs(st$pearson[1] - x2) > 1e-6 * x2) fail("the Pearson ratio of the Poisson model")
  if (!(st$pearson[1] > 10 && all(st$pearson[2:4] > 0.8 & st$pearson[2:4] < 1.4))) {
    fail("the four Pearson ratios")
  }
  if (st$family[which.min(st$AIC)] != "pig1") fail("the smallest AIC of the four families")
  if (!(exp(fits$sp_nb2_reml@coefficients$theta[[1]]) <
        exp(fits$sp_nb2@coefficients$theta[[1]]))) {
    fail("the reml() theta against the ml() theta")
  }
  if (!(stats::AIC(fits$sp_nb2) - stats::AIC(fits$sp_nb2_th) >
        stats::AIC(fits$sp_nb2) - stats::AIC(fits$sp_pig))) {
    fail("the AIC of the modelled dispersion")
  }
  if (abs(stats::sd(claims$qr_nb2) - 1) > 0.1 || !(stats::sd(claims$qr_pois) > 2)) {
    fail("the spread of the quantile residuals")
  }

  # --- positive ----------------------------------------------------------------
  if (abs(ll(fits$rent_gamma) - as.numeric(stats::logLik(ref$g_rent))) > 0.1) {
    fail("the gamma model against mgcv")
  }
  rt <- claims$rent_table
  if (abs(rt$cAIC[1] - rt$cAIC[4]) > 5 || !(min(rt$cAIC[2:3]) > max(rt$cAIC[c(1, 4)]) + 50)) {
    fail("the cAIC of the four families")
  }
  m3 <- function(f) mean(stats::residuals(f)^3)
  if (!(m3(fits$rent_lnorm) < -0.3 && m3(fits$rent_ig) < -0.3)) {
    fail("the skewness of the lognormal and inverse gaussian residuals")
  }
  if (!(caic(fits$rent_gamma) - caic(fits$rent_phi) > 20)) fail("the cAIC of the model of phi")

  # --- hyper -------------------------------------------------------------------
  if (abs(ll(fits$ab_reml) - as.numeric(stats::logLik(ref$g_ab))) > 0.05) {
    fail("the REML fit against gaulss")
  }
  at <- claims$ab_table
  if (!all(at[, "lambda_sigma"] > 1e6)) fail("the smoothing parameter of sigma")
  if (!(at["bic", "edf"] < at["reml", "edf"])) fail("the edf of BIC against REML")
  dn <- claims$s_corr@df - claims$s_naive@df
  if (!(dn > 0) || abs(claims$s_corr@aic - claims$s_naive@aic - 2 * dn) > 1e-6) {
    fail("the corrected cAIC")
  }
  if (!(claims$none[["corrected"]] < claims$none[["naive"]] &&
        claims$some[["corrected"]] <= claims$some[["naive"]])) {
    fail("the simulation of the corrected cAIC")
  }
  for (nm in c("rent_bic", "rent_cv")) {
    h <- hyper(fits[[nm]])
    sm <- h$estimate[h$source == "reml"]
    if (length(sm) != 2L || any(abs(sm - 1) < 1e-8)) {
      fail(paste("the smoothing parameters of", nm))
    }
  }
  if (!(is.finite(claims$lam_bic) && is.finite(claims$lam_cv) &&
        claims$lam_cv / claims$lam_bic > 5)) {
    fail("the two lasso choices")
  }
  lc <- function(f) { cf <- f@coefficients$mu
    names(cf) <- statmod_design(f@spec)$mu$coef_names
    cf[startsWith(names(cf), "lasso.")] }
  cb <- lc(fits$rent_bic); cc <- lc(fits$rent_cv)
  nb <- stats::coef(fits$rent_bic)$mu
  if (!is.na(nb[["lasso.B1"]]) || !is.na(stats::coef(fits$rent_cv)$mu[["lasso.B1"]])) {
    fail("the column B1 reported as not identified")
  }
  keep <- names(cb) != "lasso.B1"
  if (any(cb[keep] == 0) || any(cc[keep] == 0) || max(abs(cb[keep] - cc[keep])) > 0.02) {
    fail("the lasso coefficients under BIC and cross-validation")
  }
  if (!(caic(fits$film_g) - caic(fits$film_t) > 50)) fail("the t against the gaussian")
  nu_t <- exp(fits$film_t@coefficients$nu[[1]])
  if (!(nu_t > 8 && nu_t < 13) || abs(ll(fits$film_t) - ll(fits$film_te)) > 0.01) {
    fail("the degrees of freedom of the t under the two informations")
  }

  # --- convergence -------------------------------------------------------------
  states <- vapply(fits[c("fit1", "fit_n", "fit_s", "fit_l", "fit_l2")],
                   function(f) statmod_certificate(f)$state, character(1))
  if (!identical(unname(states), c("converged", "boundary", "unknown", "unknown",
                                   "converged"))) {
    fail("the four states of the certificate")
  }
  sl <- data$stackloss
  rq <- quantreg::rq(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data = sl)
  if (max(abs(fits$lad@coefficients$mu - stats::coef(rq))) > 0.1) {
    fail("the Laplace fit against quantreg")
  }
  obj_rq <- sum(abs(stats::resid(rq)))
  if (abs(claims$obj_rq - obj_rq) > 1e-8 ||
      !(claims$obj_lad > obj_rq + 1e-6 && claims$obj_lad < obj_rq + 0.01)) {
    fail("the sum of absolute residuals")
  }
  if (isTRUE(fits$lad@converged) || statmod_certificate(fits$lad)$state != "unknown") {
    fail("the flag and the state of the Laplace fit")
  }
  if (claims$cert_b$state == "converged") fail("the fit with a reduced budget")
  if (statmod_certificate(fits$gag_e)$state != "converged" ||
      !(fits$gag_e@coefficients$nu[[2]] < 0)) {
    fail("the fit with the expected information")
  }
  if (statmod_certificate(fits$gag_je)$state != "converged" ||
      statmod_certificate(fits$gag_jo)$state != "converged" ||
      !(claims$gap_gag > 0.1) || !(fits$gag_je@coefficients$nu[[2]] > 0)) {
    fail("the fits from the joint mode")
  }

  # --- limits ------------------------------------------------------------------
  msg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
  d_mv <- data.frame(x = seq(0, 1, length.out = 20))
  d_mv$Y <- cbind(d_mv$x, rev(d_mv$x))
  mc <- MASS::mcycle
  checks <- list(
    "does not fit a multivariate response" =
      msg(statmod(Y ~ x, distrib = mvgaussian1_distrib(2), data = d_mv)),
    "censored likelihood" =
      msg(statmod(cens(accel, accel < -100) ~ times, distrib = gaussian1_distrib(),
                  data = mc)),
    "at most one structural term" =
      msg(statmod(r ~ gas(p = 1, q = 1) | sigma ~ gas(p = 1, q = 1),
                  distrib = gaussian1_distrib(),
                  data = data.frame(r = stats::rnorm(50)))),
    "needs 2 derivatives" =
      msg(statmod(accel ~ s(times), distrib = laplace_distrib(), data = mc)),
    "single finite number" =
      msg(fixed(gaussian1_distrib(), sigma = c(1, 2))),
    "only a formula says" =
      msg(statmod(rate ~ 0 + nl(function(conc, Vm, K) Vm * conc / (K + conc),
                                K ~ state, params = c("Vm", "K"),
                                start = list(Vm = 200, K = 0.1)),
                  distrib = gaussian1_distrib(), data = Puromycin)),
    "nothing to select" =
      msg(statmod(accel ~ s(times), distrib = gaussian1_distrib(), data = mc,
                  outer_criterion = cv())))
  for (k in names(checks)) {
    if (!grepl(k, checks[[k]], fixed = TRUE)) fail(paste0("the refusal '", k, "'"))
  }
  if (!(claims$cv_p[["wald"]] < claims$cv_p[["nls"]])) fail("the two coverages")
  invisible(TRUE)
}
