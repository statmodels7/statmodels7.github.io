# Certificates for Part II, chapter 7 (linear mixed models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   intercept  the random intercept is lme()'s (fixed effects, the two
#              standard deviations, the restricted log-likelihood, the
#              predicted random effects, the interval of tau); sigma ~ Sex is
#              lme() with varIdent(), the girls' standard deviation is below
#              half the boys', and the likelihood ratio rejects a common one;
#   slope      the correlated and the uncorrelated random slopes are lmer()'s
#              (covariance matrix, standard deviations, restricted
#              log-likelihood); the likelihood ratio does not reject a zero
#              correlation; the mixed model shrinks the lines; the simulated
#              slope runs to the boundary, the certificate names it, lmer()
#              calls the fit singular, and the likelihood ratio against the
#              random intercept is small;
#   structure  the unstructured, compound symmetry and diagonal matrices are
#              lme()'s; the diagonal is rejected at 5% and compound symmetry
#              is not; the precision parametrization gives the unstructured
#              criterion; the AR(1) fit is glmmTMB's and recovers the truth;
#   groupings  the nested fit is lme()'s and equals compound symmetry; the
#              crossed fit is lmer()'s;
#   sigma      the random effect on sigma is glmmTMB's; the shared block has a
#              correlation whose interval covers zero, and the likelihood
#              ratio does not reject zero;
#   heavy      the gaussian prior inflates the standard deviation, the t
#              prior has few degrees of freedom and a small scale, more than
#              halves the error of the random intercepts, spends fewer degrees
#              of freedom and has the smaller AIC; with nu held at 1e6 it is
#              the gaussian fit; the exact marginal likelihood by quadrature
#              puts nu below one with a narrower scale, is above its value at
#              the Laplace point by more than one, and pulls sixteen of the
#              seventeen near groups to zero; the logistic, Cauchy,
#              pseudo-Huber and skew t priors compare as the text says; the
#              multivariate t with an AR(1) scale beats the gaussian on cAIC
#              and error, agrees on the standard deviation and correlation,
#              and estimates nu above 3 on every seed;
#   structure  (also) the diagonal keeps the unstructured standard deviations
#              and widens the machine contrasts; the precision readings;
#   predict    the conditional prediction is lmer()'s, the zero one is
#              lmer()'s with re.form = NA and equals the marginal one; the
#              group and prediction standard errors are the formulas of the
#              text, the second with sigma^2 exp(2 v); the marginal sigma is the zero one times exp(kappa^2/2);
#              a new group's interval is the closed form under the gaussian
#              distribution, has G = 0.025 and 0.975 at its ends by
#              integrate() under the t, with no standard error since nu < 2,
#              and is within 0.15 of 1e6 simulated measurements under the
#              logistic, by integrate() over the quantile nodes.

assert_lmm_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 7: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  ml <- function(f) as.numeric(logLik(f, type = "marginal"))
  near <- function(a, b, tol) {
    a <- as.numeric(a); b <- as.numeric(b)
    length(a) == length(b) && all(is.finite(a)) &&
      all(abs(a - b) <= tol * pmax(1, abs(b)))
  }
  state <- function(f) statmod_certificate(f)$state
  sd_lme <- function(m) as.numeric(nlme::VarCorr(m)[, "StdDev"])
  hy <- function(f) hyper(f)$estimate
  hyper_row <- function(f, param, k, name) {
    tb <- summary(f)@tables[[param]][[k]]$table
    tb[tb$name == name, ]
  }

  # --- a random intercept ----------------------------------------------------
  fe <- coef(fits$ri)$mu[c("(Intercept)", "I(age - 11)", "SexFemale")]
  if (!near(fe, nlme::fixef(ref$ri), 1e-6)) fail("the fixed effects of the random intercept")
  if (!near(c(hy(fits$ri), exp(coef(fits$ri)$sigma[[1]])), sd_lme(ref$ri), 1e-5)) {
    fail("the standard deviations of the random intercept")
  }
  if (!near(ml(fits$ri), logLik(ref$ri), 1e-7)) fail("the REML criterion of the random intercept")
  # the chapter compares the two vectors by position, so the order is checked
  if (!identical(names(claims$b_ri),
                 paste0("random.", rownames(nlme::ranef(ref$ri)))) ||
      max(abs(claims$b_ri - claims$b_lme)) > 1e-5) {
    fail("the predicted random intercepts")
  }
  iv <- nlme::intervals(ref$ri, which = "var-cov")$reStruct$Subject
  r1 <- hyper_row(fits$ri, "mu", 2, "sigma")
  if (!near(c(r1$lower, r1$upper), c(iv$lower, iv$upper), 1e-3)) {
    fail("the interval of the random intercept's standard deviation")
  }
  sg <- exp(cumsum(coef(fits$rs)$sigma))
  rt <- coef(ref$rs$modelStruct$varStruct, unconstrained = FALSE)[["Female"]]
  if (!near(sg, c(ref$rs$sigma, ref$rs$sigma * rt), 1e-4) ||
      !near(ml(fits$rs), logLik(ref$rs), 1e-7)) {
    fail("sigma ~ Sex against varIdent()")
  }
  if (!(sg[2] < sg[1] / 2 && claims$lr_sex > 10 &&
        pchisq(claims$lr_sex, df = 1, lower.tail = FALSE) < 1e-3)) {
    fail("the test of a common residual standard deviation")
  }
  if (!all(c(state(fits$ri), state(fits$rs)) == "converged")) fail("the certificates of 7.1")

  # --- a random slope --------------------------------------------------------
  S_sc <- param_value(log_cholesky(2), eta = hy(fits$sc))
  if (!near(S_sc, lme4::VarCorr(ref$sc)$Subject[, ], 1e-3) ||
      !near(ml(fits$sc), logLik(ref$sc), 1e-7) ||
      !near(coef(fits$sc)$mu[1:2], lme4::fixef(ref$sc), 1e-6)) {
    fail("the correlated random slope against lmer()")
  }
  vs <- lme4::VarCorr(ref$su)
  sd_su <- c(attr(vs$Subject, "stddev"), attr(vs$Subject.1, "stddev"))
  if (!near(hy(fits$su), sd_su, 1e-4) || !near(ml(fits$su), logLik(ref$su), 1e-7)) {
    fail("the uncorrelated random slope against lmer()")
  }
  if (!(pchisq(claims$lr_cor, df = 1, lower.tail = FALSE) > 0.5)) fail("the test of the correlation")
  if (!(stats::sd(claims$mix_slo) < stats::sd(claims$ols[, 2]) &&
        stats::sd(claims$mix_int) < stats::sd(claims$ols[, 1]))) {
    fail("the shrinkage of the lines")
  }
  cb <- statmod_certificate(fits$bd)
  S_bd <- param_value(log_cholesky(2), eta = hy(fits$bd))
  if (!(cb$state == "boundary" && any(grepl("sigma_log_L2", cb$boundary)) &&
        sqrt(S_bd[2, 2]) < 0.05 && lme4::isSingular(ref$bd) &&
        near(ml(fits$bd), logLik(ref$bd), 1e-6) && claims$lr_bd < 0.1)) {
    fail("the random slope at the boundary")
  }
  rb <- hyper_row(fits$bd, "mu", 2, "sd[t]")
  if (!is.na(rb$se)) fail("the missing standard error at the boundary")
  if (!all(c(state(fits$sc), state(fits$su)) == "converged")) fail("the certificates of 7.2")

  # --- structures for the covariance matrix ----------------------------------
  if (!near(c(ml(fits$un), ml(fits$cs), ml(fits$dg)),
            c(logLik(ref$un), logLik(ref$cs), logLik(ref$dg)), 1e-6)) {
    fail("the three structures against lme()")
  }
  S_cs <- param_value(compound_symmetry(3), eta = hy(fits$cs))
  vc <- nlme::VarCorr(ref$cs)
  v_lme <- as.numeric(vc[1, "Variance"])
  r_lme <- as.numeric(vc[2, 3])
  if (!near(c(S_cs[1, 1], S_cs[1, 2] / S_cs[1, 1]), c(v_lme, r_lme), 1e-3)) {
    fail("compound symmetry against lme()")
  }
  if (!(pchisq(claims$lr_dg, df = 3, lower.tail = FALSE) < 0.05 &&
        pchisq(claims$lr_cs, df = 4, lower.tail = FALSE) > 0.05)) {
    fail("the tests of the structures")
  }
  if (!near(ml(fits$pr), ml(fits$un), 1e-6)) fail("the precision parametrization")
  # the diagonal matrix keeps the unstructured standard deviations and widens
  # the standard errors of the machine contrasts
  S_un <- param_value(log_cholesky(3), eta = hy(fits$un))
  se_of <- function(f) sqrt(diag(vcov(f)))[c("mu:MachineB", "mu:MachineC")]
  if (!near(hy(fits$dg), sqrt(diag(S_un)), 1e-3) ||
      !all(se_of(fits$dg) > se_of(fits$un))) {
    fail("the diagonal structure against the unstructured one")
  }
  # the precision readings: conditional variances below the variances, and
  # machines A and C nearly uncorrelated given B while correlated marginally
  Om <- claims$Om
  S_pr <- solve(Om)
  pcor <- -Om[1, 3] / sqrt(Om[1, 1] * Om[3, 3])
  if (!near(S_pr, S_un, 1e-3) || !all(1 / diag(Om) < diag(S_pr)) ||
      !near(pcor, claims$pc_ac, 1e-12) || abs(pcor) > 0.1 ||
      S_pr[1, 3] / sqrt(S_pr[1, 1] * S_pr[3, 3]) < 0.5) {
    fail("the conditional variances and partial correlations")
  }
  S_ar <- param_value(ar1(5), eta = hy(fits$ar))
  vt <- glmmTMB::VarCorr(ref$ar)$cond$id
  if (!near(ml(fits$ar), logLik(ref$ar), 1e-6) ||
      !near(c(S_ar[1, 1], S_ar[1, 2] / S_ar[1, 1]), c(vt[1, 1], vt[1, 2] / vt[1, 1]), 1e-3) ||
      !near(exp(coef(fits$ar)$sigma[[1]]), sigma(ref$ar), 1e-3)) {
    fail("the AR(1) structure against glmmTMB()")
  }
  if (!(abs(sqrt(S_ar[1, 1]) - 1.5) < 0.4 && abs(S_ar[1, 2] / S_ar[1, 1] - 0.7) < 0.15)) {
    fail("the recovery of the AR(1) truth")
  }
  for (k in c("un", "cs", "dg", "pr", "ar")) {
    if (state(fits[[k]]) != "converged") fail(paste("the certificate of", k))
  }

  # --- nested and crossed groupings ------------------------------------------
  vn <- nlme::VarCorr(ref$ne)
  sd_ne <- as.numeric(vn[c(2, 4), "StdDev"])
  if (!near(hy(fits$ne), sd_ne, 1e-4) || !near(ml(fits$ne), logLik(ref$ne), 1e-6) ||
      state(fits$ne) != "converged") {
    fail("the nested groupings against lme()")
  }
  tn <- hy(fits$ne)
  if (!near(c(tn[1]^2 + tn[2]^2, tn[1]^2 / (tn[1]^2 + tn[2]^2)),
            c(S_cs[1, 1], S_cs[1, 2] / S_cs[1, 1]), 1e-3) ||
      !near(ml(fits$ne), ml(fits$cs), 1e-6)) {
    fail("the nested model as compound symmetry")
  }
  vx <- lme4::VarCorr(ref$cr)
  sd_cr <- c(attr(vx$plate, "stddev"), attr(vx$sample, "stddev"))
  if (!near(hy(fits$cr), sd_cr, 1e-4) || !near(ml(fits$cr), logLik(ref$cr), 1e-6) ||
      !near(exp(coef(fits$cr)$sigma[[1]]), sigma(ref$cr), 1e-4) ||
      state(fits$cr) != "converged") {
    fail("the crossed groupings against lmer()")
  }
  if (!(hy(fits$cr)[2] > hy(fits$cr)[1] &&
        hy(fits$cr)[1] > exp(coef(fits$cr)$sigma[[1]]))) {
    fail("the order of the crossed standard deviations")
  }

  # --- random effects on sigma -----------------------------------------------
  vd <- glmmTMB::VarCorr(ref$ds)
  sd_ds <- c(attr(vd$cond$Subject, "stddev"), attr(vd$disp$Subject, "stddev"))
  if (!near(hy(fits$ds), sd_ds, 1e-3) || !near(ml(fits$ds), logLik(ref$ds), 1e-5) ||
      state(fits$ds) != "converged") {
    fail("the random effect on sigma against glmmTMB()")
  }
  if (!(coef(fits$ds)$mu[["I(age - 11)"]] < coef(fits$ri)$mu[["I(age - 11)"]])) {
    fail("the slope in age under a random effect on sigma")
  }
  tb <- summary(fits$sh)@classes[[1]]$table
  rc <- tb[grepl("^cor", tb$name), ]
  if (!(nrow(rc) == 1L && rc$lower < 0 && rc$upper > 0 &&
        pchisq(claims$lr_sh, df = 1, lower.tail = FALSE) > 0.5 &&
        state(fits$sh) == "converged")) {
    fail("the shared covariance of mu and sigma")
  }

  # --- random effects with heavy tails ---------------------------------------
  b <- data$b_true
  rg <- coef(fits$hg)$mu; rg <- rg[grepl("random", names(rg))]
  rt2 <- coef(fits$ht)$mu; rt2 <- rt2[grepl("random", names(rt2))]
  e_g <- sqrt(mean((rg - b)^2)); e_t <- sqrt(mean((rt2 - b)^2))
  ht <- hy(fits$ht)
  if (!(hy(fits$hg) > 1.5 && ht[2] < 3 && ht[1] < 0.5 && e_t < e_g / 2 &&
        AIC(fits$ht) < AIC(fits$hg))) {
    fail("the Student t prior on the contaminated groups")
  }
  if (max(abs(coef(fits$hl)$mu - coef(fits$hg)$mu)) > 1e-4) fail("the t prior at its gaussian limit")
  re_edf <- function(f) sum(f@edf$edf[grepl("random", f@edf$term, fixed = TRUE)])
  if (!(re_edf(fits$ht) < re_edf(fits$hg))) fail("the degrees of freedom under the t prior")

  # --- other distributions on the real line ----------------------------------
  rmse <- function(f) {
    re <- coef(f)$mu[grepl("random", names(coef(f)$mu))]
    sqrt(mean((re - b)^2))
  }
  r_lo <- rmse(fits$lo); r_ca <- rmse(fits$ca); r_ph <- rmse(fits$ph)
  r_st <- rmse(fits$st)
  hst <- hy(fits$st)
  re_st <- coef(fits$st)$mu[grepl("random", names(coef(fits$st)$mu))]
  re_t <- coef(fits$ht)$mu[grepl("random", names(coef(fits$ht)$mu))]
  aics <- sapply(list(fits$hg, fits$lo, fits$ht, fits$ca, fits$ph, fits$st), AIC)
  rms <- c(e_g, r_lo, e_t, r_ca, r_ph, r_st)
  if (!(r_lo < e_g && r_lo > e_t && abs(AIC(fits$ca) - AIC(fits$ht)) < 0.1 &&
        abs(r_ca - e_t) < 0.01 && hy(fits$ph)[2] < 0.05 && r_ph < 0.15 &&
        hst[2] > 0 && which.min(aics) == 6L && which.min(rms) == 6L &&
        mean(re_st) > 0 && coef(fits$st)$mu[[1]] < coef(fits$ht)$mu[[1]] &&
        abs((coef(fits$st)$mu[[1]] + mean(re_st)) -
              (coef(fits$ht)$mu[[1]] + mean(re_t))) < 0.01)) {
    fail("the comparison of the distributions on the real line")
  }
  for (k in c("lo", "ca", "ph", "st")) {
    if (state(fits[[k]]) != "converged") fail(paste("the certificate of", k))
  }

  # --- the exact marginal likelihood of the t prior --------------------------
  op <- claims$op
  pt <- claims$par_tml
  l_exact <- -op$value
  l_lap <- claims$exact_loglik(pt)
  far <- c(3L, 9L, 16L)
  me <- claims$mode_b(op$par)
  mlap <- claims$mode_b(pt)
  if (!(exp(op$par[5]) < 1 && exp(op$par[4]) < exp(pt[4]) &&
        l_exact - l_lap > 1 && max(abs(op$par[1:3] - pt[1:3])) < 0.02 &&
        claims$rmse_exact > claims$rmse_laplace &&
        max(abs(me[far] - b[far])) < 0.3 && max(abs(mlap[far] - b[far])) < 0.3 &&
        sum(abs(me[-far]) < 0.01) == 16L && state(fits$tml) == "converged")) {
    fail("the comparison with the exact marginal likelihood")
  }

  # --- a multivariate t distribution -----------------------------------------
  hm <- hy(fits$mt)
  S_mt <- param_value(ar1(5), eta = hm[1:2])
  S_mg <- param_value(ar1(5), eta = hy(fits$mg))
  sd_t <- sqrt(hm[3] / (hm[3] - 2) * S_mt[1, 1])
  bt <- as.vector(t(data$b_mvt))
  rm_mt <- function(f) {
    re <- coef(f)$mu[grepl("random", names(coef(f)$mu))]
    sqrt(mean((re - bt)^2))
  }
  if (!(state(fits$mt) == "converged" && state(fits$mg) == "converged" &&
        AIC(fits$mt) < AIC(fits$mg) && rm_mt(fits$mt) < rm_mt(fits$mg) &&
        abs(sd_t - sqrt(S_mg[1, 1])) < 0.1 &&
        abs(S_mt[1, 2] / S_mt[1, 1] - S_mg[1, 2] / S_mg[1, 1]) < 0.05 &&
        sd_t < sqrt(3) && sqrt(S_mg[1, 1]) < sqrt(3) &&
        sqrt(mean(data$b_mvt^2)) < sqrt(3) &&
        hm[3] > 3 && all(claims$nu_seeds > 3))) {
    fail("the multivariate t distribution with an AR(1) structure")
  }

  # --- predictions -----------------------------------------------------------
  old <- data.frame(Days = c(0, 5, 9), Subject = "308")
  new <- data.frame(Days = c(0, 5, 9), Subject = "new")
  pc <- predict(fits$sc, what = "mu", newdata = old)
  if (!near(pc, stats::predict(ref$sc, newdata = old), 1e-5)) fail("the conditional prediction")
  pz <- claims$p_zero
  pm <- predict(fits$sc, what = "mu", newdata = new, se = TRUE, random = "marginal")
  if (!near(pz$fit, stats::predict(ref$sc, newdata = new, re.form = NA), 1e-6) ||
      !near(pm$fit, pz$fit, 1e-12) || !near(pm$se, pz$se, 1e-12)) {
    fail("the zero and marginal predictions")
  }
  Z <- cbind(1, new$Days)
  se_g <- sqrt(pz$se^2 + rowSums((Z %*% S_sc) * Z))
  v_ls <- vcov(fits$sc)["sigma:(Intercept)", "sigma:(Intercept)"]
  se_p <- sqrt(se_g^2 + exp(2 * coef(fits$sc)$sigma[[1]] + 2 * v_ls))
  if (!near(claims$p_group$se, se_g, 1e-6) || !near(claims$p_pred$se, se_p, 1e-6)) {
    fail("the group and prediction standard errors")
  }
  if (!all(diff(se_g) > 0)) fail("the group interval growing with the days")
  kap <- hy(fits$ds)[2]
  if (!near(claims$s_marg$fit, claims$s_zero$fit * exp(kap^2 / 2), 1e-6) ||
      !(claims$s_marg$fit > claims$s_zero$fit)) {
    fail("the marginal prediction of sigma")
  }

  # --- a new group under a distribution that is not gaussian ------------------
  hp <- claims$heavy_pred
  nw <- data.frame(x = 0.5, g = "new")
  par_at <- function(f, what) predict(f, what = what, newdata = nw, random = "zero")
  ends_g <- par_at(fits$hg, "mu") + c(-1, 1) * stats::qnorm(0.975) *
    sqrt(par_at(fits$hg, "sigma")^2 + hy(fits$hg)^2)
  if (!near(hp["gaussian", c("lower", "upper")], ends_g, 1e-6)) {
    fail("the gaussian prediction interval of a new group")
  }
  ht <- hy(fits$ht)
  mt <- par_at(fits$ht, "mu")
  st <- par_at(fits$ht, "sigma")
  Gt <- function(y) stats::integrate(function(w) {
    stats::pnorm((y - mt) / sqrt(st^2 + ht[1]^2 / w)) *
      stats::dgamma(w, shape = ht[2] / 2, rate = ht[2] / 2)
  }, lower = 0, upper = Inf, rel.tol = 1e-10)$value
  if (abs(Gt(hp["Student t", "lower"]) - 0.025) > 1e-6 ||
      abs(Gt(hp["Student t", "upper"]) - 0.975) > 1e-6 ||
      !is.na(hp["Student t", "se"]) || !(ht[2] < 2)) {
    fail("the Student t prediction interval of a new group")
  }
  # the logistic ends: G = 0.025 and 0.975 by integrate() over b, and the
  # standard error sqrt(sigma^2 + pi^2 s^2 / 3)
  ml <- par_at(fits$lo, "mu")
  sl <- par_at(fits$lo, "sigma")
  Gl <- function(y) stats::integrate(function(b)
    stats::pnorm((y - ml - b) / sl) * stats::dlogis(b, 0, hy(fits$lo)),
    -Inf, Inf, rel.tol = 1e-10)$value
  if (abs(Gl(hp["logistic", "lower"]) - 0.025) > 1e-5 ||
      abs(Gl(hp["logistic", "upper"]) - 0.975) > 1e-5 ||
      abs(hp["logistic", "se"] - sqrt(sl^2 + pi^2 * hy(fits$lo)^2 / 3)) >
        0.01 * hp["logistic", "se"]) {
    fail("the logistic prediction interval of a new group")
  }
  if (!(hp["gaussian", "width"] > hp["Student t", "width"])) {
    fail("the gaussian interval being the widest")
  }
  # the group interval of a new group's mean: G = 0.025 and 0.975 at the ends
  # of the t row, with the estimation error inside the scale mixture, and the
  # logistic row's standard error sqrt(se0^2 + pi^2 s^2 / 3)
  gh <- claims$group_heavy
  ct <- predict(fits$ht, what = "mu", newdata = nw, random = "zero", se = TRUE)
  Gg <- function(e) stats::integrate(function(v) {
    w <- exp(v)
    stats::pnorm((e - ct$fit) / sqrt(ct$se^2 + ht[1]^2 / w)) *
      stats::dgamma(w, shape = ht[2] / 2, rate = ht[2] / 2) * w
  }, -60, 15, rel.tol = 1e-11, subdivisions = 2000L)$value
  cl <- predict(fits$lo, what = "mu", newdata = nw, random = "zero", se = TRUE)
  if (abs(Gg(gh["Student t", "lower"]) - 0.025) > 1e-6 ||
      abs(Gg(gh["Student t", "upper"]) - 0.975) > 1e-6 ||
      !is.na(gh["Student t", "se"]) ||
      !near(gh["logistic", "se"], sqrt(cl$se^2 + pi^2 * hy(fits$lo)^2 / 3), 1e-8)) {
    fail("the group interval of a new group under the t and logistic priors")
  }

  invisible(TRUE)
}
