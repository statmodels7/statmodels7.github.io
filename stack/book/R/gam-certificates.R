# Certificates for Part II, chapter 6 (generalized additive models).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   bspline    the B-spline smooth is mgcv's bs = "bs" (edf, sigma, fitted
#              values) at the orders 1, 2 and 3; the edf grow from k = 10 to
#              20 by more than from 20 to 40; the order-3 smooth has the free
#              coefficients lin and poly2; the order-1 smooth has the most
#              edf; the shrinkage removes the smooth of a covariate without
#              effect where the default keeps its straight line;
#   families   the P-spline is mgcv's bs = "ps" and its AIC is close to the
#              B-spline's; the cyclic and Fourier smooths take the same value
#              at 0 and 366 where the B-spline smooth differs at the two ends
#              by more than a degree; the cyclic smooth is mgcv's "cc" with
#              k = 11 and equally spaced knots, with the same number of
#              coefficients; the Fourier smooth has an AIC within 2 of the
#              cyclic one and fewer edf; the B-spline smooth has a smaller AIC
#              than the Legendre one; the adaptive smooth is mgcv's bs = "ad",
#              its smoothing parameters span orders of magnitude, and it has
#              fewer edf and a smaller AIC than the P-spline;
#   operators  the null spaces printed; the Fourier smooth has the free
#              coefficients sin1 and cos1, the two-harmonic one sin2 and cos2
#              too, the second-derivative one none, and the harmonic operator
#              has the fewest edf and the smallest AIC of the three; the
#              B-spline smooth over the days has a free sine and cosine, the
#              second derivative at k = 20 is nearly unpenalized and loses to
#              the operator, and at k = 60 it has the smallest AIC of the
#              four; the exponential operator on wtloss gives 2 edf and the
#              nls curve, with a smaller AIC than the second derivative;
#              the co2 operator has order 4 and a
#              null space of dimension 4, and the fit with it has a much
#              smaller sigma and AIC than the fit on the second derivative;
#   by         the shared diet curves have one smoothing parameter and are
#              mgcv's with id = 1; the per-level ones have four; diet 3 is the
#              heaviest at day 21; the varying coefficient is mgcv's, carries
#              the constant, and its fitted coefficient of z is close to the
#              true one;
#   te         the anisotropic tensor product has two smoothing parameters,
#              the additive model a much larger AIC; the design of te() spans
#              the same functions as mgcv's and the fits differ;
#   sigma      the smooth sigma lowers the AIC, is larger in January than in
#              July, and agrees with gaulss to 0.1 degrees;
#   counts     the Poisson smooth is mgcv's; the Pearson ratio is above one
#              and the negative binomial has the smaller AIC;
#   lasso      the lasso keeps some and drops the others of the 18 penalized
#              coefficients, its edf are that count plus two, lambda is
#              chosen by BIC, its AIC is smaller than the roughness smooth's,
#              and the certificate is unknown.

assert_gam_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 6: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  edf <- function(f) sum(f@edf$edf[f@edf$parameter == "mu"])
  gap <- function(f, g) max(abs(fitted(f) - stats::fitted(g)))
  same <- function(f, g, tol_edf = 1e-3, tol_fit = 1e-2) {
    abs(edf(f) - sum(g$edf)) / sum(g$edf) < tol_edf && gap(f, g) < tol_fit
  }

  # --- B-spline smooths -----------------------------------------------------
  if (!same(fits$b20, ref$b20) ||
      abs(exp(coef(fits$b20)$sigma[[1]]) - sqrt(ref$b20$sig2)) > 1e-4) {
    fail("the B-spline smooth against mgcv")
  }
  if (!same(fits$o1, ref$o1) || !same(fits$o3, ref$o3)) {
    fail("the orders 1 and 3 against mgcv")
  }
  e <- c(edf(fits$b10), edf(fits$b20), edf(fits$b40))
  if (!(e[2] - e[1] > 1 && abs(e[3] - e[2]) < 0.5)) fail("the edf against k")
  if (!all(c("s(lstat).lin", "s(lstat).poly2") %in% names(coef(fits$o3)$mu))) {
    fail("the free coefficients at order 3")
  }
  if (!(edf(fits$o1) > edf(fits$b20) && edf(fits$b20) > edf(fits$o3))) {
    fail("the order of the edf of the three penalties")
  }
  ek <- fits$sk@edf$edf; es <- fits$ss@edf$edf
  if (!(abs(ek[4] - 1) < 1e-3 && es[4] < 0.01 && all(abs(ek[2:3] - es[2:3]) < 1) &&
        AIC(fits$ss) < AIC(fits$sk))) {
    fail("the shrinkage of the indus smooth")
  }

  # --- other smoothers ------------------------------------------------------
  if (!same(fits$ps, ref$ps)) fail("the P-spline against mgcv")
  if (abs(AIC(fits$ps) - AIC(fits$b20)) > 2) fail("the P-spline against the B-spline")
  ends0 <- data.frame(doy = c(0, 366))
  for (k in c("tc", "tf")) {
    p <- predict(fits[[k]], what = "mu", newdata = ends0)
    if (abs(p[1] - p[2]) > 1e-8) fail(paste("the periodicity of", k))
  }
  pb <- predict(fits$tb, what = "mu", newdata = data.frame(doy = c(1, 366)))
  if (!(abs(pb[1] - pb[2]) > 1)) fail("the two ends of the B-spline smooth")
  if (!same(fits$tc, ref$cc, tol_edf = 1e-5, tol_fit = 1e-4) ||
      length(coef(fits$tc)$mu) != length(stats::coef(ref$cc))) {
    fail("the cyclic smooth against mgcv's cc with k = 11")
  }
  if (!(abs(AIC(fits$tf) - AIC(fits$tc)) < 2 && edf(fits$tf) < edf(fits$tc))) {
    fail("the Fourier smooth against the cyclic one")
  }
  if (!(AIC(fits$b20) < AIC(fits$lg))) fail("the Legendre smooth's AIC")
  if (!same(fits$ad, ref$ad, tol_fit = 0.05)) fail("the adaptive smooth against mgcv")
  la <- hyper(fits$ad)$estimate
  if (!(max(la) / min(la) > 1e6)) fail("the spread of the adaptive smoothing parameters")
  if (!(edf(fits$ad) < edf(fits$p40) && AIC(fits$ad) < AIC(fits$p40))) {
    fail("the adaptive smooth against the P-spline")
  }

  # --- operators ------------------------------------------------------------
  if (!all(c("s(doy).sin1", "s(doy).cos1") %in% names(coef(fits$tf)$mu))) {
    fail("the free coefficients of the Fourier smooth")
  }
  if (!all(c("s(doy).sin2", "s(doy).cos2") %in% names(coef(fits$tf3)$mu)) ||
      any(grepl("sin|cos", names(coef(fits$tf2)$mu)))) {
    fail("the free coefficients of the Fourier penalties")
  }
  if (!(AIC(fits$tf) < AIC(fits$tf2) && AIC(fits$tf) < AIC(fits$tf3) &&
        edf(fits$tf) < edf(fits$tf2) && edf(fits$tf) < edf(fits$tf3))) {
    fail("the harmonic operator against the other Fourier penalties")
  }
  if (!all(c("s(day).sin1", "s(day).cos1") %in% names(coef(fits$h20)$mu))) {
    fail("the free sine and cosine of the B-spline smooth over the days")
  }
  a4 <- c(AIC(fits$d20), AIC(fits$h20), AIC(fits$d60), AIC(fits$h60))
  if (!(edf(fits$d20) > 19 && a4[2] < a4[1] - 50 && which.min(a4) == 3L)) {
    fail("the yearly cycle on the B-spline basis")
  }
  if (!(abs(edf(fits$we) - 2) < 1e-3 && gap(fits$we, ref$nls) < 1e-4 &&
        edf(fits$w2) > 4 && AIC(fits$we) < AIC(fits$w2))) {
    fail("the exponential operator on the wtloss data")
  }
  op <- claims$op_co2
  if (operator_order(op) != 4L) fail("the order of the co2 operator")
  s1 <- exp(coef(fits$c1)$sigma[[1]]); s2 <- exp(coef(fits$c2)$sigma[[1]])
  s3 <- exp(coef(fits$c3)$sigma[[1]])
  if (!(s2 < s1 / 2 && AIC(fits$c2) < AIC(fits$c1) - 500 &&
        abs(s3 - s1) < 0.1 * s1)) {
    fail("the co2 fit with the operator")
  }

  # --- by -------------------------------------------------------------------
  if (nrow(hyper(fits$cs)) != 1L || nrow(hyper(fits$cl)) != 4L) {
    fail("the smoothing parameters of the diet curves")
  }
  if (!same(fits$cs0, ref$cs0)) fail("the shared diet curves against mgcv")
  cw <- data$cw
  d21 <- vapply(1:4, function(d) {
    predict(fits$cs, what = "mu",
            newdata = data.frame(Time = 21, Diet = factor(d, levels = levels(cw$Diet))))
  }, numeric(1))
  if (which.max(d21) != 3L) fail("the heaviest diet at day 21")
  if (!same(fits$vc, ref$vc, tol_edf = 1e-5, tol_fit = 1e-5)) {
    fail("the varying coefficient against mgcv")
  }
  if (!("s(x):z.const" %in% names(coef(fits$vc)$mu))) fail("the constant of the varying coefficient")
  xv <- claims$xv
  if (max(abs(claims$beta_vc - (1 + 2 * xv^2))) > 0.4) fail("the fitted coefficient of z")

  # --- tensor products ------------------------------------------------------
  if (nrow(hyper(fits$te)) != 2L || nrow(hyper(fits$ti)) != 1L) {
    fail("the smoothing parameters of the tensor products")
  }
  if (!(AIC(fits$ta) > AIC(fits$te) + 20)) fail("the additive model against the tensor product")
  Xo <- model.matrix(fits$te, what = "mu")
  Xm <- stats::predict(ref$te, type = "lpmatrix")
  res <- Xm - Xo %*% qr.solve(Xo, Xm)
  if (max(abs(res)) > 1e-8) fail("the span of te() against mgcv's")
  if (!(gap(fits$te, ref$te) > 1)) fail("the difference of the two tensor products")

  # --- sigma ----------------------------------------------------------------
  days <- data.frame(doy = c(15, 105, 196, 288))
  sg <- predict(fits$ts, what = "sigma", newdata = days)
  sm <- 1 / stats::predict(ref$ts, newdata = days, type = "response")[, 2]
  if (!(AIC(fits$ts) < AIC(fits$tc) && sg[1] > sg[3] && max(abs(sg - sm)) < 0.1)) {
    fail("the smooth standard deviation of the temperature")
  }

  # --- counts ---------------------------------------------------------------
  if (!same(fits$dp, ref$dp, tol_fit = 0.1)) fail("the Poisson smooth against mgcv")
  mu <- fitted(fits$dp); dd <- data$chic
  pr <- sum((dd$death - mu)^2 / mu) / (nrow(dd) - edf(fits$dp))
  if (!(pr > 1 && AIC(fits$dn) < AIC(fits$dp))) fail("the overdispersion of the deaths")
  th <- exp(coef(fits$dn)$theta[[1]])
  if (abs(th - ref$dn$family$getTheta(TRUE)) / th > 1e-3 || !same(fits$dn, ref$dn, tol_fit = 0.1)) {
    fail("the negative binomial smooth against mgcv")
  }

  # --- lasso ----------------------------------------------------------------
  cl <- coef(fits$ls)$mu
  nz <- sum(cl[-(1:2)] != 0)
  if (!(nz > 0 && nz < 18)) fail("the coefficients kept by the lasso")
  if (abs(edf(fits$ls) - (nz + 2)) > 1e-6) fail("the edf of the lasso")
  if (hyper(fits$ls)$source != "bic") fail("the criterion of the lasso")
  if (!(AIC(fits$ls) < AIC(fits$r20))) fail("the AIC of the lasso")
  if (statmod_certificate(fits$ls)$state != "unknown") fail("the certificate of the lasso")

  invisible(TRUE)
}
