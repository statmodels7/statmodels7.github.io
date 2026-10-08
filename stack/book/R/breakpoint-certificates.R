# Certificates for Part III, chapter 11 (break-points).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
# segmented is always called with ::, never attached: it exports a seg() of
# its own.
#
# The claims, section by section:
#
#   seg      the GAG fit is segmented's: position, coefficients, the
#            position's standard error, and the RSS recomputed by lm.fit at
#            the position; the score process of the test is recomputed by
#            hand for a gaussian mean, with Davies' bound, and rejects; with
#            two break-points our RSS is not larger, the second position
#            agrees and the first sits on an observed age; the Poisson fit on
#            `down` is segmented's on glm (logLik, psi);
#   jump     the Nile step has stepmented's RSS, sits between 1898 and 1899,
#            its delta has lm()'s standard error with the indicator written
#            as a column, and sigma is lm()'s REML estimate; every bootstrap
#            position lies between two years and the share printed is
#            recomputed; the bootstrap test gives 1/100 with its supremum in
#            1898-1899; the coal step is stepmented's; the Gamma jseg is near
#            the simulation; the BIC is smallest for jseg() on the simulated
#            data;
#   locate   the default fit and the fit from (1880, 1950) without restarts
#            have segmented's RSS and positions; the table of the number of
#            break-points is recomputed by lm.fit at the fitted positions, the
#            BIC chooses four and the AIC five, selgmented chooses fewer
#            because its four-break-point fit has a larger RSS than ours;
#   by       the positions and the RSS of by = ~ 0 + group are segreg's; the
#            likelihood ratio is recomputed and the two models differ by
#            four degrees of freedom;
#   random   the random break-point fit certifies converged and agrees with
#            segmented.lme on the fixed effects (1e-3), the standard
#            deviations (one per cent) and the subject positions (0.005);
#            random effects on psi and gamma1: the shared variance is
#            rejected; the t prior estimates heavy tails and a smaller scale
#            than the gaussian;
#   smooth   the smoothed Nile fits give the position a finite standard
#            error inside 1898-1899, levels within 10 of the sharp step's
#            with larger standard errors, and a log-likelihood within one
#            unit of it; the width table: h = 0.1 does not converge, the
#            standard error grows with h from 0.5 on and the log-likelihood
#            falls; on the panel the standard deviation of the positions
#            grows with h and the change of level is too large at h = 2;
#   marginal the marginal fit certifies converged on the differenced
#            criterion, and agrees with the smoothed random step on the
#            position (0.3) and the change of level (0.2), the smoothed
#            standard deviation being the larger; its fitted values are the
#            posterior means, between the two levels of each subject;
#   sigma    the jseg in sigma's equation is at 14.7, its AIC is the
#            smallest, and it certifies converged; the test of the seg in
#            sigma's equation rejects.

assert_breakpoint_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part III chapter 11: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  # an empty or mismatched comparison is a failure, never a vacuous pass
  rel <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (!length(b) || length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
    max(abs(a - b)) / max(1, max(abs(b)))
  }
  rss_lm <- function(X, y) sum(stats::lm.fit(X, y)$residuals^2)
  state <- function(f) statmod_certificate(f)$state

  # --- seg -------------------------------------------------------------------
  g <- data$GAGurine
  cf <- coef(fits$gag)$mu
  sp <- ref$s_gag$psi
  cs <- stats::coef(ref$s_gag)
  if (rel(cf[["seg.psi1"]], sp[1L, "Est."]) > 1e-5 ||
      rel(cf[c("(Intercept)", "seg.beta", "seg.gamma1")],
          cs[c("(Intercept)", "Age", "U1.Age")]) > 1e-5) {
    fail("the GAG break-point against segmented")
  }
  if (rel(confint(fits$gag)["mu:seg.psi1", "se"], sp[1L, "St.Err"]) > 1e-3) {
    fail("the standard error of the GAG break-point")
  }
  X1 <- cbind(1, g$Age, pmax(g$Age - cf[["seg.psi1"]], 0))
  if (rel(claims$rss_gag, rss_lm(X1, g$lGAG)) > 1e-8 ||
      rel(claims$rss_gag, sum(stats::residuals(ref$s_gag)^2)) > 1e-6) {
    fail("the GAG residual sum of squares")
  }
  # the score process of the change of slope, written out for a gaussian
  # mean: least squares on the line, sigma at its maximum likelihood value
  tg <- claims$test_gag
  X0 <- cbind(1, g$Age)
  e0 <- stats::lm.fit(X0, g$lGAG)$residuals
  S <- vapply(tg$process$psi, function(q) {
    z <- pmax(g$Age - q, 0)
    sum(z * e0)^2 / (mean(e0^2) * sum(stats::lm.fit(X0, z)$residuals^2))
  }, 1)
  M <- max(S)
  pb <- stats::pchisq(M, 1, lower.tail = FALSE) +
    sum(abs(diff(sqrt(S)))) * exp(-M / 2) / sqrt(2 * pi)
  if (rel(tg$process$score, S) > 1e-6 || rel(tg$p.value, min(1, pb)) > 1e-6 ||
      !(tg$p.value < 1e-10)) {
    fail("the test for a change of slope in the GAG data")
  }
  cf2 <- coef(fits$gag2)$mu
  if (claims$rss_gag2 > sum(stats::residuals(ref$s_gag2)^2) + 1e-10 ||
      rel(cf2[["seg.psi2"]], ref$s_gag2$psi[2L, "Est."]) > 1e-4 ||
      !any(abs(g$Age - cf2[["seg.psi1"]]) < 1e-6)) {
    fail("the GAG fit with two break-points")
  }
  cfd <- coef(fits$down)$mu
  if (rel(as.numeric(logLik(fits$down)), as.numeric(stats::logLik(ref$s_down))) > 1e-6 ||
      rel(cfd[["seg.psi1"]], ref$s_down$psi[1L, "Est."]) > 1e-4) {
    fail("the Poisson break-point against segmented")
  }

  # --- jump ------------------------------------------------------------------
  nile <- data$nile
  cfn <- coef(fits$nile)$mu
  psi_n <- cfn[["jump.psi1"]]
  if (rel(claims$rss_nile, sum(stats::residuals(ref$s_nile)^2)) > 1e-8 ||
      !(psi_n > 1898 && psi_n < 1899)) {
    fail("the Nile step against stepmented")
  }
  ln <- stats::lm(flow ~ I(year > psi_n), data = nile)
  if (rel(confint(fits$nile)["mu:jump.delta1", "se"],
          sqrt(diag(stats::vcov(ln)))[2L]) > 1e-6 ||
      rel(exp(coef(fits$nile)$sigma[[1L]]), summary(ln)$sigma) > 1e-6) {
    fail("the conditional standard error and the REML sigma of the Nile step")
  }
  ps <- claims$psi_star
  if (length(ps) != 100L || any(!is.finite(ps)) ||
      any(abs(ps - round(ps)) < 1e-8) ||
      rel(claims$share, mean(floor(ps) == 1898)) > 0) {
    fail("the bootstrap of the Nile step")
  }
  tn <- claims$test_nile
  if (!isTRUE(all.equal(tn$p.value, 1 / 100)) ||
      !(tn$estimate > 1898 && tn$estimate < 1899)) {
    fail("the bootstrap test of the Nile step")
  }
  pc <- coef(fits$coal)$mu[["jump.psi1"]]
  if (rel(as.numeric(logLik(fits$coal)), as.numeric(stats::logLik(ref$s_coal))) > 1e-6 ||
      !(pc > 1891 && pc < 1892) || !identical(state(fits$coal), "converged")) {
    fail("the Poisson step on the coal disasters against stepmented")
  }
  cg <- coef(fits$gj)$mu
  if (abs(cg[["jseg.psi1"]] - 4) > 0.3 || abs(cg[["jseg.delta1"]] - 0.8) > 0.2 ||
      abs(cg[["jseg.gamma1"]] + 0.15) > 0.05 || !(BIC(fits$gj) < BIC(fits$gs)) ||
      !identical(state(fits$gj), "converged")) {
    fail("the Gamma jseg")
  }
  if (!(BIC(fits$jseg) < BIC(fits$seg0) && BIC(fits$jseg) < BIC(fits$jump0))) {
    fail("the BIC ordering of jseg, seg and jump")
  }

  # --- locate ----------------------------------------------------------------
  ga <- data$globTempAnom
  rss_g <- function(f) sum((ga$Anomaly - fitted(f, what = "mu"))^2)
  s_glob <- sum(stats::residuals(ref$s_glob)^2)
  for (f in list(fits$glob, fits$glob_0)) {
    if (abs(rss_g(f) - s_glob) > 1e-5 ||
        rel(coef(f)$mu[c("seg.psi1", "seg.psi2")], ref$s_glob$psi[, "Est."]) > 1e-4) {
      fail("the two starts on globTempAnom against segmented")
    }
  }
  tab <- claims$tab_npsi
  for (k in 0:5) {
    f <- fits$glob_k[[k + 1L]]
    X <- cbind(1, ga$Year)
    if (k) {
      q <- coef(f)$mu[paste0("seg.psi", seq_len(k))]
      X <- cbind(X, sapply(q, function(v) pmax(ga$Year - v, 0)))
    }
    if (rel(tab$rss[k + 1L], rss_lm(X, ga$Anomaly)) > 1e-8 ||
        abs(tab$edf[k + 1L] - (3 + 2 * k)) > 1e-8 ||
        !identical(state(f), "converged")) {
      fail(sprintf("the fit with %d break-points on globTempAnom", k))
    }
  }
  if (claims$k_bic != 4L || claims$k_aic != 5L ||
      !(nrow(ref$s_sel$psi) < claims$k_bic) ||
      !(sum(stats::residuals(ref$s_glob4)^2) > tab$rss[5L] + 0.1)) {
    fail("the choice of the number of break-points")
  }

  # --- by --------------------------------------------------------------------
  cb <- coef(fits$by)$mu
  rp <- ref$s_plant$psi[, "Est."]
  if (rel(cb[c("seg.psi1.groupRKV", "seg.psi1.groupRKW", "seg.psi1.groupRWC")],
          rp) > 1e-3 ||
      rel(sum((data$plant$y - fitted(fits$by, what = "mu"))^2),
          sum(stats::residuals(ref$s_plant)^2)) > 1e-6) {
    fail("the per-plant break-points against segreg")
  }
  lr <- 2 * (as.numeric(logLik(fits$by)) - as.numeric(logLik(fits$psi)))
  if (rel(claims$lr, lr) > 1e-10 ||
      abs(sum(fits$by@edf$edf) - sum(fits$psi@edf$edf) - 4) > 1e-6) {
    fail("the likelihood ratio between the two plant models")
  }

  # --- random ----------------------------------------------------------------
  if (!identical(state(fits$pan), "converged")) fail("the random break-point fit")
  sl <- ref$s_pan
  cp <- coef(fits$pan)$mu
  hp <- hyper(fits$pan)
  sd_lme <- as.numeric(nlme::VarCorr(sl$lme.fit)[, "StdDev"])
  fx <- nlme::fixef(sl$lme.fit)
  if (rel(c(cp[["seg.psi1.(Intercept)"]], cp[["seg.gamma1"]]),
          fx[c("G0", "U")]) > 1e-3 ||
      abs(hp$estimate[grepl("seg", hp$term)] / sd_lme[2L] - 1) > 0.01 ||
      abs(hp$estimate[!grepl("seg", hp$term)] / sd_lme[1L] - 1) > 0.01 ||
      max(abs(claims$psi_hat - sl$psi.i)) > 0.005) {
    fail("the random break-point against segmented.lme")
  }
  hs <- hyper(fits$pg_id)$estimate[grepl("seg", hyper(fits$pg_id)$term)]
  lr_g <- 2 * (as.numeric(logLik(fits$pg, type = "marginal")) -
                 as.numeric(logLik(fits$pg_id, type = "marginal")))
  if (!identical(state(fits$pg), "converged") ||
      !identical(state(fits$pg_id), "converged") ||
      length(hs) != 2L || abs(hs[1] - hs[2]) > 1e-12 ||
      rel(claims$lr_g, lr_g) > 1e-12 ||
      stats::pchisq(lr_g, df = 1, lower.tail = FALSE) > 0.05) {
    fail("the random effects on the position and on the change of slope")
  }
  # the t prior: heavy tails estimated, a scale well below the gaussian
  # standard deviation, and subject positions about as good
  ho <- hyper(fits$og)
  ht <- hyper(fits$ot)
  sc_g <- ho$estimate[ho$name == "sigma" & grepl("seg", ho$term)]
  sc_t <- ht$estimate[ht$name == "sigma" & grepl("seg", ht$term)]
  psi_of <- function(f) {
    cf <- coef(f)$mu
    cf[["seg.psi1.(Intercept)"]] + cf[paste0("seg.psi1.random.", 1:20)]
  }
  rm_g <- sqrt(mean((psi_of(fits$og) - claims$psi_o)^2))
  rm_t <- sqrt(mean((psi_of(fits$ot) - claims$psi_o)^2))
  if (!(claims$nu_ot < 4) || !(sc_t < 0.75 * sc_g) ||
      abs(rm_t / rm_g - 1) > 0.25) {
    fail("the t prior on the random positions")
  }

  # --- smooth ----------------------------------------------------------------
  ci_s <- confint(fits$nile)
  for (f in list(fits$np, fits$nq, fits$nh)) {
    ci_f <- confint(f)
    se <- ci_f["mu:jump.psi1", "se"]
    pe <- ci_f["mu:jump.psi1", "estimate"]
    lv <- c("mu:(Intercept)", "mu:jump.delta1")
    if (!is.finite(se) || se <= 0 || pe < 1898 || pe > 1899 ||
        any(abs(ci_f[lv, "estimate"] - ci_s[lv, "estimate"]) > 10) ||
        any(ci_f[lv, "se"] <= ci_s[lv, "se"]) ||
        abs(as.numeric(logLik(f)) - as.numeric(logLik(fits$nile))) >= 1) {
      fail("the smoothed steps against the sharp one")
    }
  }
  th <- as.data.frame(claims$tab_h)
  big <- th$h >= 0.5
  if (th$converged[th$h == 0.1] != 0 || any(th$converged[th$h >= 0.25] != 1) ||
      any(diff(th$se_psi[big]) <= 0) || any(diff(th$logLik[th$h >= 0.25]) >= 0) ||
      !all(th$se_psi[big] < 5)) {
    fail("the table of widths on the Nile")
  }
  tr <- as.data.frame(claims$tab_hr)
  if (any(diff(tr$sd_position) <= 0) || !(tr$delta[tr$h == 2] > 1.6)) {
    fail("the widths on the panel")
  }

  # --- marginal --------------------------------------------------------------
  ce <- statmod_certificate(fits$marg)
  if (!identical(ce$state, "converged") ||
      !identical(ce$curvature, "differenced criterion")) {
    fail("the certificate of the marginal break-point")
  }
  cm <- coef(fits$marg)$mu
  cr <- coef(fits$srand)$mu
  hr <- hyper(fits$srand)
  sd_r <- hr$estimate[grepl("jump", hr$term)]
  if (abs(cm[["jump.psi1.mean"]] - cr[["jump.psi1.(Intercept)"]]) > 0.3 ||
      abs(cm[["jump.delta1"]] - cr[["jump.delta1"]]) > 0.2 ||
      length(sd_r) != 1L || !(sd_r > cm[["jump.psi1.sd"]])) {
    fail("the marginal step against the smoothed random step")
  }
  lat <- statmod_latent(fits$marg)
  if (nrow(lat) != 20L || any(!is.finite(lat$mean)) || any(!is.finite(lat$sd)) ||
      any(lat$sd <= 0) || any(lat$mean < 0 | lat$mean > 10)) {
    fail("the posterior positions of the marginal break-point")
  }
  # the fitted values rise from one level to the other within each subject:
  # the increments over time lie between 0 and the change of level
  pj <- data$panel_j
  fv <- fitted(fits$marg, what = "mu")
  inc <- unlist(lapply(split(seq_along(fv), pj$id), function(ix) {
    diff(fv[ix[order(pj$t[ix])]])
  }))
  if (any(inc < -1e-8) || any(inc > cm[["jump.delta1"]] + 1e-8)) {
    fail("the fitted values of the marginal break-point")
  }

  # --- sigma -----------------------------------------------------------------
  if (abs(coef(fits$sjs)$sigma[["jseg.psi1"]] - 14.7) > 1e-8 ||
      !(AIC(fits$sjs) < AIC(fits$ss) && AIC(fits$sjs) < AIC(fits$sj)) ||
      !identical(state(fits$sjs), "converged")) {
    fail("the break-point in the equation of sigma")
  }
  ts <- statmod_breakpoint_test(fits$ss)
  if (!(ts$p.value < 1e-3)) fail("the test of the change of slope of sigma")
  invisible(TRUE)
}
