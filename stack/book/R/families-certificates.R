# Certificates for Part II, chapter 5 (distributions).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   heavy       the Student t fit agrees with gamlss's TF under ml(); the
#               interval of nu excludes the gaussian limit; the t has a much
#               smaller AIC than the gaussian; nu ~ Age has a negative slope,
#               is significant by the marginal likelihood ratio and converges
#               under the expected-information criterion; on gaussian data nu
#               reaches the boundary and the certificate names it;
#   skew        the residual skewness exceeds the skew normal's limit, the
#               limit is the closed form, the direct and the centred skew
#               normal reach the same maximum and slope, the AIC falls from the gaussian to
#               the skew normal to the skew t, the skew normal maximum is the
#               hand-written likelihood's, and the fitted nu of the skew t is
#               between 2 and 3;
#   beta        the ml() fit is gamlss's BE and glmmTMB's beta family, the
#               default precision is the smaller, the precision equation is not
#               significant;
#   zeros       the zero-inflated and hurdle fits are glmmTMB's under ml(), the
#               hurdle has a much smaller AIC, the signs of its za
#               coefficients, and a zi equation improves the zero-inflated fit;
#   wrappers    the truncated and folded fits are the hand-written
#               likelihoods', the truncated one recovers the truth and the
#               plain gaussian does not, the transformed gaussian is the
#               lognormal and the logit-normal beats the beta on the swiss data;
#   fixed       the fit with a known theta is glm() with negative.binomial();
#   own         the user family passes check_distrib() and is survreg's
#               log-logistic.

assert_families_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 5: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  ll <- function(f, ...) as.numeric(logLik(f, ...))

  # --- heavy tails ----------------------------------------------------------
  ga <- data$gag
  g <- ref$gl_ht
  f <- fits$ht_ml
  co <- c(coef(f)$mu, coef(f)$sigma, coef(f)$nu)
  cg <- c(stats::coef(g, what = "mu"), stats::coef(g, what = "sigma"),
          stats::coef(g, what = "nu"))
  if (max(abs(unname(co) - unname(cg))) > 1e-4 ||
      abs(ll(f) - as.numeric(stats::logLik(g))) > 1e-4) {
    fail("the Student t fit against gamlss")
  }
  ci <- confint(fits$ht)
  ci <- ci[ci$parameter == "nu", ]
  nu_hat <- exp(ci$estimate)
  if (!(nu_hat > 2 && nu_hat < 4) || !(exp(ci$upper) < 10)) {
    fail("the degrees of freedom of the Student t fit")
  }
  if (!(AIC(fits$hg) - AIC(fits$ht) > 50)) fail("the AIC of the t against the gaussian")
  cov_of <- function(fit, q) {
    m <- predict(fit, "mu", ga); s <- predict(fit, "sigma", ga)
    mean(ga$GAG < m + q(0.05) * s | ga$GAG > m + q(0.95) * s)
  }
  nu_h <- exp(coef(fits$ht)$nu[[1]])
  if (abs(cov_of(fits$ht, function(p) stats::qt(p, nu_h)) - 0.10) > 0.04 ||
      abs(cov_of(fits$hg, stats::qnorm) - 0.10) > 0.04) {
    fail("the coverage of the bands")
  }
  if (!(coef(fits$ht1)$nu[["Age"]] < 0)) fail("the sign of the age effect on nu")
  lr <- 2 * (ll(fits$ht1, type = "marginal") - ll(fits$ht0, type = "marginal"))
  if (abs(lr - claims$lr_ht) > 1e-8 || !(stats::pchisq(lr, 1, lower.tail = FALSE) < 0.01)) {
    fail("the likelihood ratio of nu ~ Age")
  }
  if (statmod_certificate(fits$ht1)$state != "converged") {
    fail("the certificate of the fit with nu ~ Age")
  }
  cb <- claims$cert_bt
  if (cb$state != "boundary" || !identical(cb$boundary, "nu/(Intercept)") ||
      !(exp(coef(fits$bt)$nu[[1]]) > 1e6) || !(AIC(fits$bg) < AIC(fits$bt))) {
    fail("the degrees of freedom at the boundary")
  }

  # --- asymmetry ------------------------------------------------------------
  gm <- (4 - pi) / 2 * (2 / pi)^1.5 / (1 - 2 / pi)^1.5
  if (abs(claims$gamma_max - gm) > 1e-12 || abs(gm - 0.9952717) > 1e-6 ||
      !(claims$skew_lm > gm)) fail("the skewness limit of the skew normal")
  sk <- claims$skew_table
  if (!(sk$AIC[3] < sk$AIC[2] && sk$AIC[2] < sk$AIC[1]) ||
      rel(sk$AIC, c(AIC(fits$sg), AIC(fits$ss), AIC(fits$st))) > 1e-10) {
    fail("the AIC of the asymmetric models")
  }
  # the direct and the centred skew normal are one family: the same maximum,
  # the same slope, and the direct fit converges from its moment start
  if (statmod_certificate(fits$sn1)$state != "converged" ||
      abs(ll(fits$sn1) - ll(fits$ss)) > 1e-3 ||
      abs(as.numeric(logLik(fits$sn1, type = "marginal")) -
            as.numeric(logLik(fits$ss, type = "marginal"))) > 1e-3 ||
      abs(coef(fits$sn1)$mu[["Horsepower"]] -
            coef(fits$ss)$mu[["Horsepower"]]) > 1e-4) {
    fail("the direct and the centred skew normal")
  }
  if (abs(-claims$opt_sn$value - ll(fits$ss_ml)) > 1e-3) {
    fail("the skew normal log-likelihood against the hand-written one")
  }
  nu_st <- exp(coef(fits$st)$nu[[1]])
  if (!(nu_st > 2 && nu_st < 3)) fail("the degrees of freedom of the skew t")
  hp150 <- data.frame(Horsepower = 150)
  th <- list(mu = predict(fits$st, what = "mu", newdata = hp150),
             sigma = predict(fits$st, what = "sigma", newdata = hp150),
             alpha = predict(fits$st, what = "alpha", newdata = hp150),
             nu = predict(fits$st, what = "nu", newdata = hp150))
  if (!is.finite(variance(skewt_distrib(), theta = th)) ||
      !is.na(skewness(skewt_distrib(), theta = th))) {
    fail("the moments of the skew t at the fitted degrees of freedom")
  }
  th$nu <- 6
  if (!is.finite(variance(skewt_distrib(), theta = th)) ||
      !is.finite(skewness(skewt_distrib(), theta = th))) {
    fail("the moments of the skew t at six degrees of freedom")
  }

  # --- proportions ----------------------------------------------------------
  b <- fits$b0_ml
  co_b <- coef(b)$mu
  if (rel(unname(co_b), unname(stats::coef(ref$gl_be, what = "mu"))) > 1e-3 ||
      rel(unname(co_b), unname(glmmTMB::fixef(ref$tmb_be)$cond)) > 1e-3 ||
      abs(ll(b) - as.numeric(stats::logLik(ref$gl_be))) > 1e-3 ||
      abs(ll(b) - as.numeric(stats::logLik(ref$tmb_be))) > 1e-3) {
    fail("the beta fit against gamlss and glmmTMB")
  }
  ph <- claims$phi_table
  if (rel(ph$phi[1:3], rep(ph$phi[1], 3)) > 1e-3 || !(ph$phi[4] < ph$phi[1])) {
    fail("the precision of the beta fits")
  }
  lrb <- 2 * (ll(fits$b1, type = "marginal") - ll(fits$b0, type = "marginal"))
  if (abs(lrb - claims$lr_b) > 1e-8 || !(stats::pchisq(lrb, 1, lower.tail = FALSE) > 0.05)) {
    fail("the likelihood ratio of the precision equation")
  }

  # --- zeros ----------------------------------------------------------------
  if (abs(ll(fits$zi_ml) - as.numeric(stats::logLik(ref$tmb_zi))) > 1e-2 ||
      abs(ll(fits$za_ml) - as.numeric(stats::logLik(ref$tmb_za))) > 1e-2) {
    fail("the zero models against glmmTMB")
  }
  if (!(AIC(fits$za) < AIC(fits$zn) - 50) ||
      !(abs(AIC(fits$zi) - AIC(fits$zn)) < 3) ||
      !(AIC(fits$zi2) < AIC(fits$zi))) fail("the AIC of the zero models")
  cz <- coef(fits$za)$za
  if (!all(cz[c("OpeningM", "OpeningS", "SolderThin")] < 0)) {
    fail("the signs of the hurdle coefficients")
  }

  # --- wrappers -------------------------------------------------------------
  tr <- fits$tr
  se_t <- sqrt(diag(stats::vcov(tr)))
  cf_t <- c(coef(tr)$mu, coef(tr)$sigma)
  if (any(abs(cf_t - c(0.5, 1.5, 0)) > 3 * se_t[seq_along(cf_t)])) {
    fail("the truncated fit against the truth")
  }
  if (!(abs(coef(fits$tg)$mu[["x"]] - 1.5) > abs(coef(tr)$mu[["x"]] - 1.5))) {
    fail("the bias of the gaussian fit to truncated data")
  }
  if (abs(-claims$opt_tr$value - ll(tr)) > 1e-3) fail("the truncated log-likelihood")
  if (abs(-claims$opt_fo$value - ll(fits$fo)) > 1e-3) fail("the folded log-likelihood")
  if (abs(ll(fits$tx) - ll(fits$ln)) > 1e-6 ||
      rel(unname(coef(fits$tx)$mu), unname(coef(fits$ln)$mu)) > 1e-5) {
    fail("the transformed gaussian against the lognormal")
  }
  if (!(AIC(fits$ln_b) < AIC(fits$b0))) fail("the logit-normal against the beta")

  # --- a known parameter ----------------------------------------------------
  if (!identical(unname(claims$d_fixed@params), "mu") ||
      rel(unname(coef(fits$fx)$mu), unname(stats::coef(ref$glm_fx))) > 1e-3 ||
      abs(ll(fits$fx) - as.numeric(stats::logLik(ref$glm_fx))) > 1e-3) {
    fail("the fit with a known dispersion against glm()")
  }

  # --- a user distribution --------------------------------------------------
  ck <- claims$chk_ll
  sk <- attr(ck, "skipped")
  if (nrow(ck) != 9L || !all(ck$status == "OK") || is.null(sk) ||
      nrow(sk) != 4L || !all(grepl("vs finite differences", sk$check)) ||
      !all(grepl("numerical fallback", sk$reason))) {
    fail("the checks of the user distribution")
  }
  sv <- ref$sv_ll
  co_l <- c(coef(fits$ll)$mu, coef(fits$ll)$shape)
  if (rel(unname(co_l), c(unname(stats::coef(sv)), -log(sv$scale))) > 1e-3 ||
      abs(ll(fits$ll) - as.numeric(stats::logLik(sv))) > 1e-3) {
    fail("the user distribution against survreg")
  }

  invisible(TRUE)
}
