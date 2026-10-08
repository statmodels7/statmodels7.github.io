# Certificates for Part II, chapter 8 (generalized linear mixed models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   poisson    the epil random intercept is glmmTMB's with REML = TRUE (fixed
#              effects, standard deviation, restricted log-likelihood);
#              ml(marginal = "all") is glmer()'s, and ml() at the joint mode
#              is below it and differs from it;
#   bernoulli  ml(marginal = "all") is glmer()'s; adaptive quadrature with 25
#              nodes gives a larger standard deviation than the Laplace one;
#   negbin     the random intercept and the random effects on theta are
#              glmmTMB's (mean, dispersion, both standard deviations,
#              restricted log-likelihood); the food treatment on theta and
#              the nests' effect on theta are both large; the nests' standard
#              deviation on the mean falls;
#   zero       the zero-inflated Poisson with random effects on mu and zi and
#              the zero-inflated negative binomial are glmmTMB's; the
#              negative binomial has a smaller cAIC than the Poisson with the
#              same equation for zi, and its zero inflation in the unaffected
#              streams runs to zero with a large standard error; with a random
#              effect on zi the negative binomial is not certified, glmmTMB
#              returns no log-likelihood, and the two agree on the standard
#              deviation of that effect but not on the intercept;
#   predict    the marginal probabilities are the integrals over the prior,
#              computed here from the coefficients; the Zeger approximation is
#              close; the marginal log odds ratio is smaller than the
#              conditional one by about the approximation's factor; the group
#              interval of theta is exp(eta0 +- 1.96 s) with s from the
#              confidence standard error and kappa, and its standard error is
#              the lognormal one and not the delta method; the prediction
#              intervals of a new nest and of a new stream are the smallest
#              integers with G(y) >= p, with G a double integral over the
#              random effects, and their standard deviations are the
#              closed-form ones.

assert_glmm_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 8: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  ml <- function(f) as.numeric(logLik(f, type = "marginal"))
  near <- function(a, b, tol) {
    a <- as.numeric(a); b <- as.numeric(b)
    length(a) == length(b) && all(is.finite(a)) &&
      all(abs(a - b) <= tol * pmax(1, abs(b)))
  }
  state <- function(f) statmod_certificate(f)$state
  hy <- function(f) hyper(f)$estimate
  tmb_sd <- function(m, comp) {
    vc <- glmmTMB::VarCorr(m)[[comp]]
    as.numeric(attr(vc[[1]], "stddev"))
  }
  z <- stats::qnorm(0.975)

  # --- counts with a random intercept ----------------------------------------
  fe <- coef(fits$ep)$mu[1:6]
  if (!near(fe, glmmTMB::fixef(ref$ep)$cond, 1e-5) ||
      !near(hy(fits$ep), tmb_sd(ref$ep, "cond"), 1e-4) ||
      !near(ml(fits$ep), stats::logLik(ref$ep), 1e-7)) {
    fail("the epil random intercept against glmmTMB")
  }
  gl <- ref$ep_glmer
  if (!near(coef(fits$ep_all)$mu[1:6], lme4::fixef(gl), 1e-3) ||
      !near(hy(fits$ep_all), attr(lme4::VarCorr(gl)$subject, "stddev"), 1e-3) ||
      !near(ml(fits$ep_all), stats::logLik(gl), 1e-6)) {
    fail("ml(marginal = \"all\") against glmer() on epil")
  }
  if (!(ml(fits$ep_ml) < ml(fits$ep_all) &&
        abs(coef(fits$ep_ml)$mu[["trtprogabide"]] -
              coef(fits$ep_all)$mu[["trtprogabide"]]) > 1e-3)) {
    fail("the joint-mode ml() against the Laplace maximum")
  }

  # --- binary responses --------------------------------------------------------
  gb <- ref$bc
  if (!near(coef(fits$bc_all)$mu[1:4], lme4::fixef(gb), 1e-3) ||
      !near(hy(fits$bc_all), attr(lme4::VarCorr(gb)$ID, "stddev"), 1e-3) ||
      !near(ml(fits$bc_all), stats::logLik(gb), 1e-6)) {
    fail("the bacteria model against glmer()")
  }
  sd25 <- attr(lme4::VarCorr(ref$bc_25)$ID, "stddev")
  if (!(sd25 > hy(fits$bc_all) * 1.02)) fail("the quadrature standard deviation")
  if (!(all(coef(fits$bc)$mu[c("trtdrug", "trtdrug+")] < 0) &&
        coef(fits$bc)$mu[["late"]] < 0)) {
    fail("the signs of the bacteria coefficients")
  }

  # --- negative binomial counts -------------------------------------------------
  if (!near(coef(fits$ow1)$mu[1:4], glmmTMB::fixef(ref$ow1)$cond, 1e-4) ||
      !near(exp(coef(fits$ow1)$theta[[1]]), glmmTMB::sigma(ref$ow1), 1e-4) ||
      !near(ml(fits$ow1), stats::logLik(ref$ow1), 1e-7)) {
    fail("the Owls random intercept against glmmTMB")
  }
  if (!near(coef(fits$ow2)$mu[1:4], glmmTMB::fixef(ref$ow2)$cond, 1e-4) ||
      !near(coef(fits$ow2)$theta[1:2], glmmTMB::fixef(ref$ow2)$disp, 1e-3) ||
      !near(hy(fits$ow2), c(tmb_sd(ref$ow2, "cond"), tmb_sd(ref$ow2, "disp")),
            1e-3) ||
      !near(ml(fits$ow2), stats::logLik(ref$ow2), 1e-7)) {
    fail("the random effects on theta against glmmTMB")
  }
  if (!(stats::pchisq(claims$lr_food, df = 1, lower.tail = FALSE) < 1e-6 &&
        claims$lr_nest > 10 && hy(fits$ow2)[1] < hy(fits$ow1) &&
        coef(fits$ow2)$theta[["FoodTreatmentSatiated"]] < 0)) {
    fail("the tests on the dispersion")
  }

  # --- zero inflation ---------------------------------------------------------
  fx <- glmmTMB::fixef(ref$zi)
  if (!near(coef(fits$zi)$mu[1:8], fx$cond, 1e-4) ||
      !near(coef(fits$zi)$zi[1:2], fx$zi, 1e-3) ||
      !near(hy(fits$zi), c(tmb_sd(ref$zi, "cond"), tmb_sd(ref$zi, "zi")), 1e-3) ||
      !near(ml(fits$zi), stats::logLik(ref$zi), 1e-7)) {
    fail("the zero-inflated Poisson against glmmTMB")
  }
  fz <- glmmTMB::fixef(ref$zinb)
  if (!near(coef(fits$zinb)$mu[1:8], fz$cond, 1e-3) ||
      !near(ml(fits$zinb), stats::logLik(ref$zinb), 1e-7)) {
    fail("the zero-inflated negative binomial against glmmTMB")
  }
  se_a1 <- sqrt(diag(vcov(fits$zinb)))[["zi:minedno"]]
  if (!(AIC(fits$zinb) < AIC(fits$zip2) && sum(coef(fits$zinb)$zi) < -3 &&
        se_a1 > 5)) {
    fail("the zero inflation of the unaffected streams")
  }
  # the negative binomial with a random effect on zi: not certified here and
  # not converged in glmmTMB; similar estimates of the sd on zi, different
  # zi intercepts and slopes (0.186.0: 1.924, 1.00 and -5.64 against 1.954,
  # 0.61 and -5.29), and the REML criterion -817.819 (0.179.0 read -817.957,
  # 0.178.0 stopped at -820.322; 0.182.0/0.183.0 moved it to -816.368,
  # distributions7 0.75.0 to -816.375; statmodels7 0.197.0, which repairs an
  # indefinite observed information before the expected fallback, stops at
  # -817.819 in 17.7 s against 368 s; with the levels of site in their
  # original order it stops at -817.824. The outer curvature is not readable
  # at either point, a rank-one matrix of order 1e14 with an equilibrated
  # eigenvalue of -5, and the criterion jumps by 4 between neighbouring
  # evaluations under 0.196.0, so the value is an arithmetic accident)
  zr <- claims$zinb_re
  tr <- claims$tmb_zinb_re
  sd_zr <- attr(glmmTMB::VarCorr(tr)$zi[[1]], "stddev")
  if (state(zr) %in% c("converged", "boundary") ||
      !is.na(as.numeric(stats::logLik(tr))) ||
      !near(hy(zr)[2], sd_zr, 0.05) ||
      abs(coef(zr)$zi[[1]] - glmmTMB::fixef(tr)$zi[[1]]) < 0.2 ||
      abs(as.numeric(stats::logLik(zr, type = "marginal")) + 817.819) > 1e-3) {
    fail("the random effect on the zero inflation of the negative binomial")
  }
  # and the summary gives its two standard deviations no standard error
  zt <- unlist(lapply(summary(zr)@tables, function(eq) lapply(eq, function(bl) {
    tb <- bl$table
    if (is.data.frame(tb) && "role" %in% names(tb)) tb$se[tb$role == "estimated"]
  })))
  if (length(zt) != 2L || any(is.finite(zt))) {
    fail("the random effect on the zero inflation of the negative binomial")
  }

  # --- the typical group and the average group ---------------------------------
  b <- coef(fits$bc)$mu
  eta0 <- b[["(Intercept)"]] + b[["late"]] + c(0, b[["trtdrug"]], b[["trtdrug+"]])
  tau <- hy(fits$bc)
  pbar <- vapply(eta0, function(e) stats::integrate(function(u)
    stats::plogis(e + u) * stats::dnorm(u, sd = tau), -40 * tau, 40 * tau,
    rel.tol = 1e-12)$value, 0)
  if (!near(claims$p_zero, stats::plogis(eta0), 1e-10) ||
      !near(claims$p_marg, pbar, 1e-8)) {
    fail("the marginal probabilities")
  }
  c16 <- 16 * sqrt(3) / (15 * pi)
  if (!(max(abs(claims$p_zeger - pbar)) < 0.02)) fail("the Zeger approximation")
  ratio <- claims$lor_marg / claims$lor_cond
  if (!(abs(claims$lor_marg) < abs(claims$lor_cond) &&
        abs(ratio - 1 / sqrt(1 + c16^2 * tau^2)) < 0.05)) {
    fail("the attenuation of the log odds ratio")
  }

  # --- the parameter of a new group ---------------------------------------------
  nn <- claims$new_nest
  conf <- predict(fits$ow2, what = "link:theta", newdata = nn,
                  random = "zero", se = TRUE)
  s <- sqrt(conf$se^2 + hy(fits$ow2)[2]^2)
  gt <- claims$g_theta
  sd_ln <- exp(conf$fit + s^2 / 2) * sqrt(expm1(s^2))
  if (!near(claims$g_link$se, s, 1e-8) ||
      !near(c(gt$lower, gt$upper), exp(c(conf$fit - z * s, conf$fit + z * s)),
            1e-8) ||
      !near(gt$se, sd_ln, 1e-8) || !near(gt$fit, exp(conf$fit), 1e-10)) {
    fail("the group interval of theta")
  }
  if (!all(gt$se > exp(conf$fit) * s * 1.3)) fail("the delta method being smaller")

  # --- a new observation --------------------------------------------------------
  # G(y) by a double integral over the two random effects, which shares no
  # code with predict()
  G2 <- function(cdf, y, sb, sc) {
    stats::integrate(function(cc) vapply(cc, function(ci)
      stats::integrate(function(bb) cdf(y, bb, ci) * stats::dnorm(bb, sd = sb),
                       -9 * sb, 9 * sb, rel.tol = 1e-11)$value, 0) *
        stats::dnorm(cc, sd = sc), -9 * sc, 9 * sc, rel.tol = 1e-10)$value
  }
  ends_ok <- function(cdf, pr, sb, sc) {
    ok <- TRUE
    for (k in c("lower", "fit", "upper")) {
      p <- c(lower = 0.025, fit = 0.5, upper = 0.975)[[k]]
      y <- pr[[k]]
      ok <- ok && G2(cdf, y, sb, sc) >= p &&
        (y == 0 || G2(cdf, y - 1, sb, sc) < p)
    }
    ok
  }
  em <- predict(fits$ow2, what = "link:mu", newdata = nn, random = "zero")
  et <- predict(fits$ow2, what = "link:theta", newdata = nn, random = "zero")
  hw <- hy(fits$ow2)
  for (i in seq_len(nrow(nn))) {
    cdf <- function(y, bb, cc) stats::pnbinom(y, mu = exp(em[i] + bb),
                                              size = exp(et[i] + cc))
    if (!ends_ok(cdf, claims$p_calls[i, ], hw[1], hw[2])) {
      fail("the prediction interval of a new nest")
    }
    Emu <- exp(em[i] + hw[1]^2 / 2)
    Emu2 <- exp(2 * em[i] + 2 * hw[1]^2)
    v <- Emu + Emu2 * exp(-et[i] + hw[2]^2 / 2) + Emu2 - Emu^2
    if (!near(claims$p_calls$se[i], sqrt(v), 1e-6)) {
      fail("the standard deviation of a new nest's count")
    }
  }
  ns <- claims$new_site
  sm <- predict(fits$zi, what = "link:mu", newdata = ns, random = "zero")
  sz <- predict(fits$zi, what = "link:zi", newdata = ns, random = "zero")
  hs <- hy(fits$zi)
  for (i in seq_len(nrow(ns))) {
    cdf <- function(y, bb, cc) {
      p0 <- stats::plogis(sz[i] + cc)
      p0 + (1 - p0) * stats::ppois(y, exp(sm[i] + bb))
    }
    if (!ends_ok(cdf, claims$p_sal[i, ], hs[1], hs[2])) {
      fail("the prediction interval of a new stream")
    }
    q <- stats::integrate(function(cc) (1 - stats::plogis(sz[i] + cc)) *
                            stats::dnorm(cc, sd = hs[2]), -Inf, Inf,
                          rel.tol = 1e-12)$value
    Emu <- exp(sm[i] + hs[1]^2 / 2)
    Emu2 <- exp(2 * sm[i] + 2 * hs[1]^2)
    v <- q * (Emu + Emu2) - (q * Emu)^2
    if (!near(claims$p_sal$se[i], sqrt(v), 1e-6)) {
      fail("the standard deviation of a new stream's count")
    }
  }

  if (!all(vapply(fits, state, "") %in% c("converged", "boundary"))) {
    fail("the certificates of chapter 8")
  }
  invisible(TRUE)
}
