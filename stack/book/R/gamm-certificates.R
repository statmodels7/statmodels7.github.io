# Certificates for Part II, chapter 9 (generalized additive mixed models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   intercept  the random intercept is mgcv's s(Chick, bs = "re"), standard
#              deviation sigma / sqrt(sp) and fitted values; the random slope
#              is glmmTMB's (both standard deviations, the correlation, which
#              is strongly negative, and sigma); the slope's likelihood ratio
#              is large; the REML constant between the two packages is
#              log|a|, a the ratio of the columns of the smooth's free part;
#   band       the Bayesian and frequentist bands of the typical curve are
#              mgcv's Vp and Ve; the frequentist band is narrower everywhere,
#              the unconditional one at most 10% wider than the Bayesian, and
#              from day 7 the Bayesian standard error is within 10% of the
#              standard deviation of the mean deviation of 50 chicks; the
#              contribution of the smooth is mgcv's type = "terms";
#   by         the four smoothing parameters do not differ significantly,
#              diet 3's is the smallest, the uncorrelated version is mgcv's
#              with id = 1 (fitted values and edf);
#   chick      the REML criterion with smooth deviations is the exact gaussian
#              restricted likelihood; the deviations are significant, spend
#              100 to 160 edf and cut sigma to between a quarter and 0.4 of
#              its value; the curve of diet 3 moves by less than 5 g and its
#              band widens; chicks 33 and 11 stop growing;
#   sigma      the REML criterion of both sigma models is the Laplace
#              approximation written out with the analytic Hessian; the random
#              effects in sigma are significant; the chicks' standard
#              deviation in the mean is below 1 g with sigma modelled and above
#              20 g with sigma constant; glmmTMB's standard
#              deviations agree within 3% and its fitted values differ by more
#              than 0.01 and less than 1 g; under ML each package's maximum is
#              the exact gaussian ML with its own matrices, and the two differ;
#              the group interval of sigma is exp(eta0 +- 1.96 s);
#   counts     theta is mgcv's and glmmTMB's, the predictor is mgcv's; the
#              calls peak near 22:40 and dip near midnight and 3:00, and the
#              deprived nestlings call more, by a constant factor, at
#              every hour; with random effects on theta the standard
#              deviations are glmmTMB's and the likelihood ratio is the same in
#              both packages;
#   id         the shared parameter is mgcv's id = 1; the free parameters are
#              preferred and the same day's curve is the smoothest; the two
#              models' lag curves differ by less than 0.005 on the log scale
#              on 98% of the days, most at the ends of lag 0; the hottest day
#              is in July 1995; the
#              grouping of diets 1, 2 and 4 is not rejected.

assert_gamm_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 9: ", what, " no longer agrees with the package.",
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
  tmb_sd <- function(m, comp, grp) {
    vc <- glmmTMB::VarCorr(m)[[comp]][[grp]]
    as.numeric(attr(vc, "stddev"))
  }
  CW <- data$ChickWeight

  # --- a smooth and a random intercept ---------------------------------------
  g <- ref$g_ri
  if (!near(hy(fits$ri)[2], sqrt(g$sig2 / g$sp[2]), 1e-4) ||
      max(abs(fitted(fits$ri) - fitted(g))) > 1e-3) {
    fail("the random intercept against mgcv")
  }
  # the random slope: covariance from the chart, against glmmTMB
  h <- hyper(fits$rs)
  S <- parameters7::param_value(parameters7::log_cholesky(2),
                                eta = h$estimate[h$name != "lambda"])
  t <- ref$tmb_rs
  vt <- glmmTMB::VarCorr(t)$cond$Chick
  cor_rs <- S[1, 2] / sqrt(S[1, 1] * S[2, 2])
  if (!near(sqrt(diag(S)), attr(vt, "stddev"), 1e-3) ||
      !near(cor_rs, attr(vt, "correlation")[1, 2], 1e-3) ||
      !near(exp(coef(fits$rs)$sigma[[1]]), glmmTMB::sigma(t), 1e-4) ||
      cor_rs > -0.8) {
    fail("the random slope against glmmTMB")
  }
  # the variance of b0 + b1 t, minimized on a grid rather than by the formula
  tt <- seq(from = 0, to = 21, by = 0.001)
  vv <- S[1, 1] + 2 * tt * S[1, 2] + tt^2 * S[2, 2]
  if (!near(claims$t_min, tt[which.min(vv)], 1e-3) || claims$t_min < 1 ||
      claims$t_min > 6) {
    fail("the time at which the chicks are closest")
  }
  if (!(claims$lr_rs > 100 &&
        near(claims$lr_rs, 2 * (ml(fits$rs) - ml(fits$ri)), 1e-12))) {
    fail("the likelihood ratio of the random slope")
  }
  # the REML constant is the log of the ratio of the free columns, computed
  # here by a projection rather than by lm()
  xs <- as.numeric(model.matrix(fits$rs)[, "s(Time).lin"])
  xt <- as.numeric(glmmTMB::getME(t, "X")[, "s(Time)1"])
  a <- sum(xs * xt) / sum(xt * xt)
  if (max(abs(xs - a * xt)) > 1e-8 * max(abs(xs)) ||
      !near(log(abs(a)), as.numeric(logLik(t)) - ml(fits$rs), 1e-6) ||
      !near(claims$a_lin, a, 1e-8)) {
    fail("the REML constant of the smooth")
  }

  # --- the band of the typical curve -------------------------------------------
  # the Bayesian and frequentist bands of the random-intercept model are mgcv's
  # Vp and Ve, computed here from mgcv's own lpmatrix
  tg <- seq(from = 0, to = 21, length.out = 101)
  ndt <- data.frame(Time = tg, Chick = factor("1", levels = levels(CW$Chick)))
  Xg <- stats::predict(g, newdata = ndt, type = "lpmatrix")
  Xg[, grepl("Chick", colnames(Xg))] <- 0
  sem <- function(V) sqrt(rowSums((Xg %*% V) * Xg))
  ses <- function(f, ty) {
    predict(f, what = "mu", newdata = ndt, random = "zero", se = TRUE, type = ty)$se
  }
  if (max(abs(ses(fits$ri, "bayesian") / sem(g$Vp) - 1)) > 1e-4 ||
      max(abs(ses(fits$ri, "frequentist") / sem(g$Ve) - 1)) > 1e-4) {
    fail("the bands of the typical curve against mgcv")
  }
  bb <- claims$bands
  sb <- bb$bayesian$se; sf <- bb$frequentist$se; su <- bb$unconditional$se
  sdm <- sqrt(vapply(tg, function(t) drop(c(1, t) %*% S %*% c(1, t)), 0) / 50)
  if (!(all(sf < sb) && all(su >= sb * (1 - 1e-6)) && max(su / sb) < 1.1 &&
        min(sf) > 0.8 && max(sf) < 1.6 && sb[101] > 4 * sb[1] &&
        all(abs(sdm[tg >= 7] / sb[tg >= 7] - 1) < 0.1) &&
        near(sb, ses(fits$rs, "bayesian"), 1e-12))) {
    fail("the three bands of the typical curve")
  }
  # the contribution of the smooth is mgcv's type = "terms", computed here from
  # the coefficients and vcov() rather than by part_band()
  Xr <- as.matrix(model.matrix(fits$ri))
  ks <- startsWith(colnames(Xr), "s(Time).")
  Vs <- as.matrix(vcov(fits$ri)[paste0("mu:", colnames(Xr)[ks]),
                                paste0("mu:", colnames(Xr)[ks])])
  ft <- drop(Xr[, ks] %*% fits$ri@coefficients$mu[ks])
  st <- sqrt(rowSums((Xr[, ks] %*% Vs) * Xr[, ks]))
  tt_ <- stats::predict(g, type = "terms", se.fit = TRUE)
  if (max(abs(ft - tt_$fit[, "s(Time)"])) > 1e-3 ||
      max(abs(st - tt_$se.fit[, "s(Time)"])) > 1e-4 ||
      !near(claims$part_ri$fit, ft, 1e-10) || !near(claims$part_ri$se, st, 1e-10)) {
    fail("the contribution of the smooth against mgcv")
  }

  # --- a smooth for each group -------------------------------------------------
  if (!(stats::pchisq(claims$lr_by, df = 3, lower.tail = FALSE) > 0.05 &&
        near(claims$lr_by, 2 * (ml(fits$lv) - ml(fits$sh)), 1e-12))) {
    fail("the test of the four smoothing parameters")
  }
  lam <- hy(fits$lv)[1:4]
  q3 <- lam[3] / lam[1:2]
  if (which.min(lam) != 3L || any(q3 < 0.15) || any(q3 > 0.4) ||
      !(lam[4] > lam[3] && lam[4] < min(lam[1:2]))) {
    fail("the four smoothing parameters")
  }
  edf_mu <- function(f) sum(f@edf$edf[f@edf$parameter == "mu"])
  if (max(abs(fitted(fits$shu) - fitted(ref$g_shu))) > 1e-4 ||
      !near(edf_mu(fits$shu), sum(ref$g_shu$edf), 1e-4)) {
    fail("the uncorrelated slopes with a factor by against mgcv")
  }
  # the curves of the diets: diet 3 heaviest and diet 1 lightest at day 21
  w21 <- vapply(1:4, function(d) {
    predict(fits$sh, what = "mu", random = "zero",
            newdata = data.frame(Time = 21, Diet = factor(d, levels = levels(CW$Diet)),
                                 Chick = factor("1", levels = levels(CW$Chick))))
  }, 0)
  if (which.max(w21) != 3L || which.min(w21) != 1L) fail("the order of the diets")

  # --- a curve for each chick ------------------------------------------------------
  # the REML criterion of both models is the exact gaussian restricted
  # likelihood with the covariance written out
  y <- CW$weight
  if (!near(exact_reml(fits$gs, y), ml(fits$gs), 1e-8) ||
      !near(exact_reml(fits$sh, y), ml(fits$sh), 1e-8)) {
    fail("the REML criterion of the model with smooth deviations")
  }
  e_gs <- fits$gs@edf$edf[4]
  sg_ratio <- exp(coef(fits$gs)$sigma[[1]] - coef(fits$sh)$sigma[[1]])
  if (!(near(claims$lr_gs, 2 * (ml(fits$gs) - ml(fits$sh)), 1e-12) &&
        claims$lr_gs > 500 && e_gs > 100 && e_gs < 160 &&
        sg_ratio < 0.4 && sg_ratio > 0.25)) {
    fail("the test of the smooth deviations")
  }
  # the free line of each diet multiplies the standardized time on that
  # diet's rows; the slopes are ordered diet 3 largest, diet 1 smallest
  Xgs <- as.matrix(model.matrix(fits$gs))
  zt <- (CW$Time - mean(CW$Time)) / stats::sd(CW$Time)
  lin_nm <- paste0("s(Time).", 1:4, ".lin")
  for (d in 1:4) {
    r <- CW$Diet == d
    if (max(abs(Xgs[r, lin_nm[d]] - zt[r])) > 1e-8 || any(Xgs[!r, lin_nm[d]] != 0)) {
      fail("the free line of the diet curves")
    }
  }
  sl <- coef(fits$gs)$mu[lin_nm]
  if (which.max(sl) != 3L || which.min(sl) != 1L) fail("the slopes of the diets")
  s_rows <- summary(fits$gs)@tables$mu
  s_rows <- unlist(lapply(s_rows, function(b) b$table$name[b$table$role == "coefficient"]))
  if (!all(lin_nm %in% s_rows)) fail("the free lines in the summary")
  d3 <- claims$diet3
  if (max(abs(d3$straight$fit - d3$smooth$fit)) > 5 ||
      mean(d3$smooth$se > d3$straight$se) < 0.75) {
    fail("the curve of diet 3 with smooth deviations")
  }
  # chicks 33 and 11 stop growing; chick 33's straight curve misses its birth
  # weight by about 20 g
  for (ck in c("33", "11")) {
    nd <- data.frame(Time = c(18, 21), Diet = CW$Diet[CW$Chick == ck][1],
                     Chick = factor(ck, levels = levels(CW$Chick)))
    p <- predict(fits$gs, what = "mu", newdata = nd)
    if (p[2] >= p[1]) fail("the chicks that stop growing")
  }
  n33 <- data.frame(Time = 0, Diet = factor(3, levels = levels(CW$Diet)),
                    Chick = factor("33", levels = levels(CW$Chick)))
  miss <- predict(fits$sh, what = "mu", newdata = n33) -
    CW$weight[CW$Chick == "33" & CW$Time == 0]
  if (miss < 15 || miss > 25) fail("the weight at birth of chick 33")

  # --- smooths and random effects in sigma --------------------------------------
  # the REML criterion is the Laplace approximation written out here with the
  # analytic Hessian of the gaussian location-scale log-likelihood, at the
  # fitted coefficients, sigma's free coefficients held
  y <- CW$weight
  for (f in list(fits$sg0, fits$sg)) {
    if (!near(laplace_ls(f, y), ml(f), 1e-9)) fail("the REML criterion of the sigma model")
  }
  hs <- hy(fits$sg)
  if (!(near(claims$lr_sg, 2 * (ml(fits$sg) - ml(fits$sg0)), 1e-12) &&
        0.5 * stats::pchisq(claims$lr_sg, df = 1, lower.tail = FALSE) < 1e-6 &&
        hs[2] < 1 && hy(fits$sg0)[2] < 1 && hy(fits$ri)[2] > 20)) {
    fail("the test of the random effects in sigma")
  }
  ts <- ref$tmb_sg
  d_fit <- max(abs(fitted(fits$sg) - fitted(ts)))
  sd_t <- c(tmb_sd(ts, "cond", "Chick"), tmb_sd(ts, "disp", "Chick"))
  if (length(sd_t) != 2L || any(abs(hs[c(2, 4)] - sd_t) / sd_t >= 0.03) ||
      d_fit > 1 || d_fit < 0.01) {
    fail("the random effects in sigma against glmmTMB")
  }
  # the explanation: under ML the two packages' smooths are different models,
  # each package's maximum being the exact gaussian ML with its own matrices
  gm <- ref$g_ml
  Xs7 <- as.matrix(model.matrix(fits$ml))
  ex7 <- exact_ml(y, Xs7[, 1:2], Xs7[, -(1:2)])
  Xg <- stats::predict(gm, type = "lpmatrix")
  sm <- gm$smooth[[1]]
  ii <- sm$first.para:sm$last.para
  e <- eigen(sm$S[[1]], symmetric = TRUE)
  r <- sum(e$values > max(e$values) * 1e-10)
  Xu <- Xg[, ii] %*% e$vectors
  exg <- exact_ml(y, cbind(Xg[, 1], Xu[, -(1:r), drop = FALSE]),
                  sweep(Xu[, 1:r], 2, sqrt(e$values[1:r]), "/"))
  if (!near(ex7, ml(fits$ml), 1e-6) || !near(exg, -gm$gcv.ubre, 1e-6) ||
      abs(ex7 - exg) < 0.1) {
    fail("the two maximum likelihoods of the smooth")
  }
  nd <- data.frame(Time = c(2, 10, 20),
                   Diet = factor(1, levels = levels(CW$Diet)), Chick = "new")
  gs <- predict(fits$sg, what = "sigma", newdata = nd, random = "zero",
                interval = "group")
  cs <- predict(fits$sg, what = "link:sigma", newdata = nd, random = "zero",
                se = TRUE)
  s <- sqrt(cs$se^2 + hs[4]^2)
  if (!near(c(gs$lower, gs$upper),
            exp(c(cs$fit - stats::qnorm(0.975) * s, cs$fit + stats::qnorm(0.975) * s)),
            1e-8)) {
    fail("the group interval of sigma")
  }

  # --- a model for counts --------------------------------------------------------
  Ow <- data$Owls
  th <- exp(coef(fits$ow)$theta[[1]])
  if (!near(th, ref$g_ow$family$getTheta(TRUE), 1e-4) ||
      !near(th, glmmTMB::sigma(ref$tmb_ow), 1e-4) ||
      max(abs(predict(fits$ow, what = "link:mu") -
                predict(ref$g_ow, type = "link"))) > 1e-3 ||
      !near(hy(fits$ow)[2], tmb_sd(ref$tmb_ow, "cond", "Nest"), 1e-3)) {
    fail("the owls' smooth against mgcv and glmmTMB")
  }
  ag <- seq(from = min(Ow$ArrivalTime), to = max(Ow$ArrivalTime), length.out = 200)
  curve <- function(food) {
    nd <- data.frame(ArrivalTime = ag,
                     FoodTreatment = factor(food, levels(Ow$FoodTreatment)),
                     SexParent = factor("Male", levels(Ow$SexParent)),
                     logBroodSize = 0, Nest = "new")
    predict(fits$ow, what = "mu", newdata = nd, random = "zero")
  }
  cd <- curve("Deprived")
  cs2 <- curve("Satiated")
  # the peak near 22:40, the lowest points near midnight and near 3:00 (local
  # minima of the curve), a constant ratio between the two treatments
  peak <- ag[which.max(cd)]
  mins <- ag[which(diff(sign(diff(cd))) > 0) + 1]
  rng <- range(Ow$ArrivalTime)
  ratio <- cd / cs2
  if (!(peak > 22.4 && peak < 23 && length(mins) >= 2 &&
        any(abs(mins - 24) < 0.4) && any(abs(mins - 27) < 0.4) &&
        all(cd > cs2) && diff(range(ratio)) < 1e-8 &&
        abs(rng[1] - (21 + 40 / 60)) < 0.1 && abs(rng[2] - (29 + 15 / 60)) < 0.1)) {
    fail("the shape of the owls' curves")
  }
  t2 <- ref$tmb_ow2
  if (!near(hy(fits$ow2)[2:3], c(tmb_sd(t2, "cond", "Nest"), tmb_sd(t2, "disp", "Nest")),
            1e-3) ||
      !near(coef(fits$ow2)$theta[1:2], glmmTMB::fixef(t2)$disp, 1e-3) ||
      !near(claims$lr_ow, claims$lr_ow_tmb, 1e-4) ||
      !near(claims$lr_ow, 2 * (ml(fits$ow2) - ml(fits$ow)), 1e-12)) {
    fail("the random effects on theta against glmmTMB")
  }

  # --- shared smoothing parameters -----------------------------------------------
  hid <- hy(fits$id)
  if (!(length(unique(hid[2:4])) == 1L &&
        max(abs(predict(fits$id, what = "link:mu") -
                  predict(ref$g_id, type = "link"))) < 1e-4)) {
    fail("the shared parameter of the lags against mgcv")
  }
  hf <- hy(fits$free)
  if (!(stats::pchisq(claims$lr_id, df = 2, lower.tail = FALSE) < 0.01 &&
        which.max(hf[2:4]) == 1L &&
        identical(fits$free@edf$term[3], "s(temp, bspline_smooth(k = 8))") &&
        all(fits$free@edf$edf[3] < fits$free@edf$edf[4:5] - 2) &&
        near(claims$lr_id, 2 * (ml(fits$free) - ml(fits$id)), 1e-12))) {
    fail("the test of the shared parameter of the lags")
  }
  # the curves of the lags: the two models agree at lags 1 and 2 and differ
  # only at the ends of lag 0; the rate ratios quoted in the text
  ch <- data$chicagoNMMAPS
  part <- function(f, lab) {
    X <- as.matrix(model.matrix(f))
    k <- startsWith(colnames(X), paste0(lab, "."))
    drop(X[, k] %*% f@coefficients$mu[k])
  }
  dd <- vapply(c("s(temp)", "s(temp1)", "s(temp2)"), function(l) {
    a <- abs(part(fits$id, l) - part(fits$free, l))
    c(max = max(a), q98 = as.numeric(stats::quantile(a, 0.98)))
  }, c(max = 0, q98 = 0))
  hot <- which.max(ch$temp)
  r1 <- exp(part(fits$free, "s(temp1)"))[which.max(ch$temp1)]
  r2 <- exp(part(fits$free, "s(temp2)"))[which.max(ch$temp2)]
  rc <- exp(part(fits$free, "s(temp2)"))[which.min(abs(ch$temp2 + 15))]
  r0c <- exp(part(fits$free, "s(temp)"))[which.min(abs(ch$temp + 15))]
  if (!(all(dd["q98", ] < 0.005) && dd["max", 1] > 1.5 * max(dd["max", 2:3]) &&
        format(ch$date[hot], "%Y-%m") == "1995-07" &&
        ch$temp1[which.max(ch$temp1)] == ch$temp[hot] &&
        near(c(r1, r2, rc), claims$lag_ratios, 1e-10) &&
        r1 > 1.3 && r2 > r1 && rc > 1.05 && r0c < 1)) {
    fail("the curves of the lags")
  }

  hl <- hy(fits$lid)
  if (!(hl[1] == hl[2] && hl[1] == hl[4] && hl[3] != hl[1] &&
        stats::pchisq(claims$lr_lid, df = 2, lower.tail = FALSE) > 0.05 &&
        near(claims$lr_lid, 2 * (ml(fits$lv) - ml(fits$lid)), 1e-12))) {
    fail("the grouping of the diets")
  }

  if (!all(vapply(fits, state, "") %in% c("converged", "boundary"))) {
    fail("the certificates of chapter 9")
  }
  invisible(TRUE)
}


# The REML criterion of a gaussian location-scale fit, by Laplace's method
# written out: the log-likelihood and the log prior of the penalized
# coefficients at the fitted coefficients, minus half the log-determinant of
# the analytic negative Hessian over the integrated coefficients (every
# coefficient of mu, and the penalized ones of sigma). A smooth's wiggly
# coordinates carry the penalty lambda * I, a random effect 1 / sd^2.
laplace_ls <- function(fit, y) {
  Xm <- as.matrix(model.matrix(fit, what = "mu"))
  Xs <- as.matrix(model.matrix(fit, what = "sigma"))
  b <- fit@coefficients$mu
  g <- fit@coefficients$sigma
  h <- hyper(fit)
  est <- h$estimate
  pen_of <- function(nm, par) {
    lam <- numeric(length(nm))
    z <- grepl(".z", nm, fixed = TRUE)
    re <- grepl("random", nm, fixed = TRUE)
    if (any(z)) lam[z] <- est[h$parameter == par & h$name == "lambda"]
    if (any(re)) lam[re] <- 1 / est[h$parameter == par & h$name == "sigma"]^2
    lam
  }
  lm_ <- pen_of(colnames(Xm), "mu")
  ls_ <- pen_of(colnames(Xs), "sigma")
  ints <- ls_ > 0
  r <- y - drop(Xm %*% b)
  es <- drop(Xs %*% g)
  s2 <- exp(2 * es)
  ll <- sum(-es - 0.5 * log(2 * pi) - r^2 / (2 * s2))
  pr <- function(cf, lam) {
    k <- lam > 0
    sum(-lam[k] / 2 * cf[k]^2 + 0.5 * log(lam[k]) - 0.5 * log(2 * pi))
  }
  Xi <- Xs[, ints, drop = FALSE]
  Hmm <- crossprod(Xm, Xm / s2) + diag(lm_)
  Hms <- crossprod(Xm, Xi * (2 * r / s2))
  Hss <- crossprod(Xi, Xi * (2 * r^2 / s2)) + diag(ls_[ints], sum(ints))
  H <- rbind(cbind(Hmm, Hms), cbind(t(Hms), Hss))
  ll + pr(b, lm_) + pr(g, ls_) - 0.5 * as.numeric(determinant(H)$modulus) +
    ncol(H) / 2 * log(2 * pi)
}

# The exact gaussian restricted log-likelihood of a ChickWeight model whose
# penalized parts are the wiggly coordinates of the diet curves (lambda_1 I),
# the random intercept and slope of each chick (Sigma_b) and, if present, the
# smooth deviations of the chicks (lambda_2 I), with sigma at its estimate and
# a flat prior on the coefficients without a penalty.
exact_reml <- function(fit, y) {
  X <- as.matrix(model.matrix(fit))
  nm <- colnames(X)
  h <- hyper(fit)
  est <- h$estimate
  lam <- est[h$name == "lambda"]
  S <- parameters7::param_value(parameters7::log_cholesky(2),
                                eta = est[h$name != "lambda"])
  sz <- startsWith(nm, "s(Time).") & grepl(".z", nm, fixed = TRUE)
  dv <- startsWith(nm, "dev.")
  ri <- startsWith(nm, "random.")
  fx <- !(sz | dv | ri)
  V <- exp(2 * coef(fit)$sigma[[1]]) * diag(length(y)) +
    tcrossprod(X[, sz]) / lam[1]
  if (any(dv)) V <- V + tcrossprod(X[, dv]) / lam[2]
  Zr <- X[, ri]
  chick <- sub("random.", "", colnames(Zr), fixed = TRUE)
  chick <- sub(".(Intercept)", "", sub(".Time", "", chick, fixed = TRUE), fixed = TRUE)
  for (cc in unique(chick)) {
    Zc <- Zr[, chick == cc]
    V <- V + Zc %*% S %*% t(Zc)
  }
  R <- chol(V)
  W <- backsolve(R, X[, fx], transpose = TRUE)
  u <- backsolve(R, y, transpose = TRUE)
  rr <- u - W %*% qr.coef(qr(W), u)
  -0.5 * ((length(y) - sum(fx)) * log(2 * pi) + 2 * sum(log(diag(R))) +
            as.numeric(determinant(crossprod(W))$modulus) + sum(rr^2))
}

# The exact gaussian log-likelihood of y ~ N(Xf beta, sigma^2 I + Z Z' / lambda),
# maximized over beta, sigma and lambda.
exact_ml <- function(y, Xf, Z) {
  n <- length(y)
  one <- function(p) {
    V <- exp(2 * p[2]) * diag(n) + tcrossprod(Z) / exp(p[1])
    R <- chol(V)
    W <- backsolve(R, Xf, transpose = TRUE)
    u <- backsolve(R, y, transpose = TRUE)
    rr <- u - W %*% qr.coef(qr(W), u)
    -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(R))) + sum(rr^2))
  }
  -stats::optim(c(log(0.04), log(38)), function(p) -one(p),
                control = list(reltol = 1e-12))$value
}
