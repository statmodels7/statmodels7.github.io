# Certificates for Part III, chapter 13 (time series).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   gas       the Nile log-likelihood is a filter written out here, started at
#             zero, and no point near the estimates is higher (so sigma is
#             maximized too); the interval for xi1 is the rhobit interval of
#             the partial autocorrelation mapped back, inside (-1, 1); the
#             ARMA reading is close to arima(); with q = 2 the same, and the
#             likelihood ratio is recomputed and not significant;
#             at d = 1 the same likelihood, the loading divided by sigma^2;
#   scale     the DAX log-likelihood is the log-scale filter written out, at a
#             maximum; rugarch's log-likelihood is a GARCH recursion written
#             out at its estimates; GARCH is the larger with the gaussian
#             distribution, and the two standard deviations correlate;
#             the spelling with omega free and no intercept gives the same
#             likelihood at omega = beta0 (1 - beta);
#   families  the Student t filter written out with its bounded score; the t
#             model beats the gaussian by AIC and the t GARCH by
#             log-likelihood; after August 1991 the gaussian standard
#             deviation rises more and stays high longer; the Poisson and
#             negative binomial filters written out, the Poisson at a
#             maximum and both above tscount; the exponential filter at a
#             maximum and on ACDm's estimates through the omega mapping; the
#             gamma and beta filters written out, both below their static
#             models by AIC, the law with a negative coefficient;
#   scaling   the Poisson filters scaled by mu^-1/2 and mu^-1 written out,
#             d = 1 at a maximum and the largest; the GARCH, INGARCH and
#             linear ACD models are the scaled filters on the identity link,
#             written out both ways for GARCH, at a maximum, and within 5e-3
#             of rugarch, tscount and ACDm through the omega mapping;
#   panels    the Ovary log-likelihoods are one filter per mare written out,
#             the shared one at a maximum and close to gls() with ARMA(1, 1)
#             errors; the shared model is the model with tau at zero (its
#             marginal likelihood at tau = 1e-5 is the shared maximum); the
#             likelihood ratio is recomputed and not significant under the
#             boundary mixture; the forecast of mare 1 is her recursion
#             continued without the score; the random intercept goes to the
#             shared model and has the smaller marginal likelihood;
#   forecast  the continuation is the recursion without the score, decaying
#             to the intercept; its standard error comes with a warning and
#             is the delta method on the hand-written forecast, rising and
#             then settling near the intercept's;
#   regime    every log-likelihood is a forward recursion written out from
#             the stationary distribution, with the transition matrix built
#             from the log-ratios by hand; the probabilities the summary
#             reports are that matrix; the first row of the geyser chain sits
#             at the edge and carries no standard error; statmod_latent() is a
#             forward-backward pass written out, and depmixS4's forward pass
#             agrees with it; the criteria the text reads are recomputed; the
#             predicted scale of the DAX chain is the exponential of the
#             predictor averaged over the states.

assert_timeseries_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part III chapter 13: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  # an empty or mismatched comparison is a failure, never a vacuous pass
  rel <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (!length(b) || length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
    max(abs(a - b)) / max(1, max(abs(b)))
  }
  ll <- function(f) as.numeric(logLik(f))
  # no point near p is higher than the value at p: the fit is at a maximum
  # of the hand-written log-likelihood, over every parameter
  at_max <- function(nll, p, tol = 1e-5) {
    o <- stats::optim(p, nll, method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
    nll(p) - o$value < tol
  }

  # --- gas on the mean: Nile --------------------------------------------------
  y <- data$nile$flow
  n <- length(y)
  # f_t = a s_{t-1} + b1 f_{t-1} + b2 f_{t-2}, s = (y - mu)/s2, f_1 = f_0 = 0
  nile_ll <- function(c0, a, b1, b2, s) {
    f <- numeric(n + 1L)            # f[t + 1] is f_t, f[1] is f_0
    out <- 0
    for (t in seq_len(n)) {
      if (t > 1L) {
        sc <- (y[t - 1L] - c0 - f[t]) / s^2
        f[t + 1L] <- a * sc + b1 * f[t] + b2 * f[t - 1L]
      }
      out <- out + stats::dnorm(y[t], c0 + f[t + 1L], s, log = TRUE)
    }
    out
  }
  c1 <- coef(fits$nile)$mu
  s1 <- exp(coef(fits$nile)$sigma[["(Intercept)"]])
  if (abs(nile_ll(c1[["(Intercept)"]], c1[["gas.kappa1"]], c1[["gas.xi1"]], 0, s1) -
          ll(fits$nile)) > 1e-8) {
    fail("the Nile log-likelihood")
  }
  nll1 <- function(p) -nile_ll(p[1], exp(p[2]), tanh(p[3]), 0, exp(p[4]))
  if (!at_max(nll1, c(c1[["(Intercept)"]], log(c1[["gas.kappa1"]]),
                      atanh(c1[["gas.xi1"]]), log(s1)))) {
    fail("the maximum of the Nile log-likelihood")
  }
  # the interval of xi1: the rhobit interval of pacf1, mapped back
  V <- vcov(fits$nile, readable = FALSE)
  zp <- fits$nile@structural[[1]]$unconstrained[["pacf1"]]
  sz <- sqrt(V["mu:gas.pacf1", "mu:gas.pacf1"])
  q <- stats::qnorm(0.975)
  ci_b <- confint(fits$nile)["mu:gas.xi1", c("lower", "upper")]
  # the rhobit chart is atanh: the coordinate is atanh(xi1)
  if (abs(zp - atanh(c1[["gas.xi1"]])) > 1e-10 ||
      rel(unlist(ci_b), tanh(zp + c(-1, 1) * q * sz)) > 1e-10 ||
      !(ci_b$lower > -1 && ci_b$upper < 1)) {
    fail("the interval of the persistence")
  }
  # the spelling without an intercept: omega free, the filter started at
  # omega/(1 - beta); at omega = beta0 (1 - beta) it is the same likelihood
  om_ll <- function(om, a, b, s) {
    f <- om / (1 - b)
    out <- 0
    for (t in seq_len(n)) {
      if (t > 1L) f <- om + a * (y[t - 1L] - fp) / s^2 + b * fp
      fp <- f
      out <- out + stats::dnorm(y[t], f, s, log = TRUE)
    }
    out
  }
  if (abs(om_ll(c1[["(Intercept)"]] * (1 - c1[["gas.xi1"]]), c1[["gas.kappa1"]],
                c1[["gas.xi1"]], s1) - ll(fits$nile)) > 1e-8) {
    fail("the two spellings of the level")
  }
  # d = 1 on the gaussian mean: the information is constant, so the same
  # likelihood, with the loading divided by sigma^2
  c1d1 <- coef(fits$nile_1)$mu
  if (abs(ll(fits$nile_1) - ll(fits$nile)) > 1e-6 ||
      rel(c1d1[["gas.kappa1"]], c1[["gas.kappa1"]] / s1^2) > 1e-4 ||
      rel(c1d1[["gas.xi1"]], c1[["gas.xi1"]]) > 1e-4) {
    fail("the Nile model at d = 1")
  }
  # the ARMA reading
  s2 <- s1^2
  ar <- stats::coef(ref$arima_nile)
  if (abs(c1[["gas.xi1"]] - ar[["ar1"]]) > 0.05 ||
      abs(c1[["gas.kappa1"]] / s2 - c1[["gas.xi1"]] - ar[["ma1"]]) > 0.05 ||
      abs(ll(fits$nile) - as.numeric(stats::logLik(ref$arima_nile))) > 0.5) {
    fail("the ARMA(1, 1) reading of the Nile model")
  }
  # q = 2
  c2 <- coef(fits$nile2)$mu
  s2b <- exp(coef(fits$nile2)$sigma[["(Intercept)"]])
  if (abs(nile_ll(c2[["(Intercept)"]], c2[["gas.kappa1"]], c2[["gas.xi1"]],
                  c2[["gas.xi2"]], s2b) - ll(fits$nile2)) > 1e-8) {
    fail("the Nile log-likelihood with two lags")
  }
  # the stationary chart of q = 2: b1 = r1 (1 - r2), b2 = r2
  r2 <- c2[["gas.xi2"]]
  r1 <- c2[["gas.xi1"]] / (1 - r2)
  nll2 <- function(p) {
    a1 <- tanh(p[3]); a2 <- tanh(p[4])
    -nile_ll(p[1], exp(p[2]), a1 * (1 - a2), a2, exp(p[5]))
  }
  if (!at_max(nll2, c(c2[["(Intercept)"]], log(c2[["gas.kappa1"]]),
                      atanh(r1), atanh(r2), log(s2b))) ||
      rel(claims$lr_nile, 2 * (ll(fits$nile2) - ll(fits$nile))) > 1e-12 ||
      # with two lags the interval of xi1 is the symmetric one on the
      # coefficient's own scale, and here it leaves (-1, 1)
      rel(unlist(confint(fits$nile2)["mu:gas.xi1", c("lower", "upper")]),
          c2[["gas.xi1"]] + c(-1, 1) * q *
            confint(fits$nile2)["mu:gas.xi1", "se"]) > 1e-10 ||
      !(confint(fits$nile2)["mu:gas.xi1", "upper"] > 1) ||
      !(stats::pchisq(claims$lr_nile, 1, lower.tail = FALSE) > 0.05)) {
    fail("the comparison of one and two lags")
  }

  # --- gas on the scale: DAX ---------------------------------------------------
  r <- data$dax$r
  m <- length(r)
  # log sigma_t = g0 + f_t, f_t = a s_{t-1} + b f_{t-1}, f_1 = 0
  scale_ll <- function(mu, g0, a, b, score, dens) {
    f <- 0
    out <- 0
    for (t in seq_len(m)) {
      if (t > 1L) f <- a * score((r[t - 1L] - mu) / sp) + b * f
      sp <- exp(g0 + f)
      out <- out + dens((r[t] - mu) / sp) - log(sp)
    }
    out
  }
  gs <- function(z) z^2 - 1
  gd <- function(z) stats::dnorm(z, log = TRUE)
  cd <- coef(fits$dax)
  pd <- c(cd$mu[["(Intercept)"]], cd$sigma[["(Intercept)"]],
          cd$sigma[["gas.kappa1"]], cd$sigma[["gas.xi1"]])
  if (abs(scale_ll(pd[1], pd[2], pd[3], pd[4], gs, gd) - ll(fits$dax)) > 1e-7) {
    fail("the DAX log-likelihood")
  }
  nlld <- function(p) -scale_ll(p[1], p[2], exp(p[3]), tanh(p[4]), gs, gd)
  if (!at_max(nlld, c(pd[1], pd[2], log(pd[3]), atanh(pd[4])))) {
    fail("the maximum of the DAX log-likelihood")
  }
  # rugarch's GARCH(1, 1), started at the mean squared residual
  garch_ll <- function(cf, dens) {
    e <- r - cf[["mu"]]
    h <- mean(e^2)
    out <- 0
    for (t in seq_len(m)) {
      if (t > 1L) h <- cf[["omega"]] + cf[["alpha1"]] * e[t - 1L]^2 +
        cf[["beta1"]] * h
      out <- out + dens(e[t], sqrt(h))
    }
    out
  }
  gn <- rugarch::coef(ref$garch_n)
  if (abs(garch_ll(gn, function(e, s) stats::dnorm(e, 0, s, log = TRUE)) -
          rugarch::likelihood(ref$garch_n)) > 1e-3 ||
      !(rugarch::likelihood(ref$garch_n) > ll(fits$dax)) ||
      !(stats::cor(predict(fits$dax, what = "sigma"),
                   as.numeric(rugarch::sigma(ref$garch_n))) > 0.8)) {
    fail("the comparison with the gaussian GARCH model")
  }

  # --- families: Student t and Poisson ------------------------------------------
  ct <- coef(fits$dax_t)
  nu <- exp(ct$nu[["(Intercept)"]])
  ts_ <- function(z) (nu + 1) * z^2 / (nu + z^2) - 1
  td <- function(z) stats::dt(z, nu, log = TRUE)
  pt_ <- c(ct$mu[["(Intercept)"]], ct$sigma[["(Intercept)"]],
           ct$sigma[["gas.kappa1"]], ct$sigma[["gas.xi1"]])
  # rugarch's std is the t standardized to unit variance
  gt <- rugarch::coef(ref$garch_t)
  sh <- gt[["shape"]]
  std_d <- function(e, s) {
    k <- sqrt(sh / (sh - 2))
    stats::dt(e / s * k, sh, log = TRUE) + log(k) - log(s)
  }
  if (abs(scale_ll(pt_[1], pt_[2], pt_[3], pt_[4], ts_, td) - ll(fits$dax_t)) > 1e-7 ||
      max(abs(ts_(c(1e3, -1e4)) - nu)) > 1e-3 ||
      !(AIC(fits$dax_t) < AIC(fits$dax)) ||
      abs(garch_ll(gt, std_d) - rugarch::likelihood(ref$garch_t)) > 1e-3 ||
      !(ll(fits$dax_t) > rugarch::likelihood(ref$garch_t))) {
    fail("the Student t score-driven model")
  }
  w <- data$disc$n
  k <- length(w)
  pois_ll <- function(c0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_len(k)) {
      if (t > 1L) f <- a * (w[t - 1L] - exp(c0 + fp)) + b * f
      fp <- f
      out <- out + stats::dpois(w[t], exp(c0 + f), log = TRUE)
    }
    out
  }
  cp <- coef(fits$disc)$mu
  if (abs(pois_ll(cp[["(Intercept)"]], cp[["gas.kappa1"]], cp[["gas.xi1"]]) -
          ll(fits$disc)) > 1e-8 ||
      !at_max(function(p) -pois_ll(p[1], exp(p[2]), tanh(p[3])),
              c(cp[["(Intercept)"]], log(cp[["gas.kappa1"]]),
                atanh(cp[["gas.xi1"]]))) ||
      !(AIC(fits$disc) < AIC(fits$disc0) && AIC(fits$disc) < AIC(fits$disc_nb))) {
    fail("the Poisson score-driven model")
  }
  # the gaussian and t standard deviations after August 1991: the gaussian
  # peak is more than twice the t peak, and thirty days later the gaussian is
  # still above its median while the t is below its own
  sd_g <- predict(fits$dax, what = "sigma")
  sd_t <- predict(fits$dax_t, what = "sigma") * sqrt(nu / (nu - 2))
  j <- which.max(sd_g[1:200])
  if (!(max(sd_g[1:200]) > 2 * max(sd_t[1:200])) ||
      !(sd_g[j + 30L] > stats::median(sd_g) && sd_t[j + 30L] < stats::median(sd_t))) {
    fail("the two standard deviations after August 1991")
  }
  # the negative binomial filter: score theta (y - mu)/(theta + mu)
  cn <- coef(fits$disc_gnb)
  th <- exp(cn$theta[["(Intercept)"]])
  nb_ll <- function(c0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_len(k)) {
      if (t > 1L) {
        m0 <- exp(c0 + fp)
        f <- a * th * (w[t - 1L] - m0) / (th + m0) + b * f
      }
      fp <- f
      out <- out + stats::dnbinom(w[t], size = th, mu = exp(c0 + f), log = TRUE)
    }
    out
  }
  cnm <- cn$mu
  if (abs(nb_ll(cnm[["(Intercept)"]], cnm[["gas.kappa1"]], cnm[["gas.xi1"]]) -
          ll(fits$disc_gnb)) > 1e-8 ||
      !(AIC(fits$disc_gnb) < AIC(fits$disc))) {
    fail("the negative binomial score-driven model")
  }
  # tscount: its log-likelihoods are the sums of the densities at its fitted
  # means, and ours are the larger
  tp <- ref$tsc_p
  tn <- ref$tsc_nb
  mt <- as.numeric(stats::fitted(tp))
  if (abs(sum(stats::dpois(w, mt, log = TRUE)) - as.numeric(stats::logLik(tp))) > 1e-6 ||
      !(ll(fits$disc) > as.numeric(stats::logLik(tp))) ||
      !(ll(fits$disc_gnb) > as.numeric(stats::logLik(tn)))) {
    fail("the comparison with tscount")
  }
  # the logarithmic ACD: the exponential filter written out, at a maximum,
  # its omega mapped onto ACDm's, the estimates within 1e-3
  xd <- data$dur$x
  nd_ <- length(xd)
  ca <- coef(fits$acd)$mu
  acd_ll <- function(c0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_len(nd_)) {
      if (t > 1L) f <- a * (xd[t - 1L] / exp(c0 + fp) - 1) + b * f
      fp <- f
      out <- out + stats::dexp(xd[t], rate = exp(-(c0 + f)), log = TRUE)
    }
    out
  }
  ap <- ref$acd$mPar
  om_acd <- ca[["(Intercept)"]] * (1 - ca[["gas.xi1"]]) - ca[["gas.kappa1"]]
  if (abs(acd_ll(ca[["(Intercept)"]], ca[["gas.kappa1"]], ca[["gas.xi1"]]) -
          ll(fits$acd)) > 1e-8 ||
      !at_max(function(p) -acd_ll(p[1], exp(p[2]), tanh(p[3])),
              c(ca[["(Intercept)"]], log(ca[["gas.kappa1"]]),
                atanh(ca[["gas.xi1"]]))) ||
      max(abs(c(om_acd, ca[["gas.kappa1"]], ca[["gas.xi1"]]) - ap)) > 1e-3 ||
      abs(ll(fits$acd) - ref$acd$goodnessOfFit["LogLikelihood", "value"]) > 0.1) {
    fail("the logarithmic ACD model against ACDm")
  }

  # --- scaling the score ------------------------------------------------------
  # the Poisson filter on the log link driven by u = (y - mu)/mu^d, written
  # out at d = 1/2 and d = 1; d = 1 at a maximum and the largest of the three
  yd <- data$disc$n
  nn <- length(yd)
  pois_d_ll <- function(c0, a, b, d) {
    f <- 0
    out <- 0
    for (t in seq_len(nn)) {
      if (t > 1L) {
        m <- exp(c0 + fp)
        f <- a * (yd[t - 1L] - m) / m^d + b * f
      }
      fp <- f
      out <- out + stats::dpois(yd[t], exp(c0 + f), log = TRUE)
    }
    out
  }
  ch <- coef(fits$disc_h)$mu
  c1d <- coef(fits$disc_1)$mu
  if (abs(pois_d_ll(ch[["(Intercept)"]], ch[["gas.kappa1"]], ch[["gas.xi1"]], 0.5) -
          ll(fits$disc_h)) > 1e-8 ||
      abs(pois_d_ll(c1d[["(Intercept)"]], c1d[["gas.kappa1"]], c1d[["gas.xi1"]], 1) -
          ll(fits$disc_1)) > 1e-8 ||
      !at_max(function(p) -pois_d_ll(p[1], exp(p[2]), tanh(p[3]), 1),
              c(c1d[["(Intercept)"]], log(c1d[["gas.kappa1"]]),
                atanh(c1d[["gas.xi1"]]))) ||
      !(ll(fits$disc_1) > ll(fits$disc_h) && ll(fits$disc_1) > ll(fits$disc))) {
    fail("the scaled Poisson filters")
  }
  # GARCH(1, 1): the variance filter driven by (r - mu)^2 - sigma2, written
  # out; the same likelihood as the GARCH recursion at omega = g0 (1 - xi),
  # kappa, xi - kappa, started at g0; at a maximum; within 5e-3 of rugarch's
  # estimates and 0.1 of its log-likelihood; the three coefficients positive
  rr <- data$dax$r
  nr <- length(rr)
  gas_v_ll <- function(mu, g0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_len(nr)) {
      if (t > 1L) f <- a * ((rr[t - 1L] - mu)^2 - (g0 + fp)) + b * f
      fp <- f
      v <- g0 + f
      if (!(v > 0)) return(-1e10)
      out <- out + stats::dnorm(rr[t], mu, sqrt(v), log = TRUE)
    }
    out
  }
  garch_map_ll <- function(mu, om, ag, bg, v1) {
    v <- v1
    out <- 0
    for (t in seq_len(nr)) {
      if (t > 1L) v <- om + ag * (rr[t - 1L] - mu)^2 + bg * v
      out <- out + stats::dnorm(rr[t], mu, sqrt(v), log = TRUE)
    }
    out
  }
  mg <- coef(fits$garch)$mu[["(Intercept)"]]
  cg <- coef(fits$garch)$sigma2
  g_map <- c(mg, cg[["(Intercept)"]] * (1 - cg[["gas.xi1"]]), cg[["gas.kappa1"]],
             cg[["gas.xi1"]] - cg[["gas.kappa1"]])
  if (abs(gas_v_ll(mg, cg[["(Intercept)"]], cg[["gas.kappa1"]], cg[["gas.xi1"]]) -
          ll(fits$garch)) > 1e-8 ||
      abs(garch_map_ll(g_map[1], g_map[2], g_map[3], g_map[4], cg[["(Intercept)"]]) -
          ll(fits$garch)) > 1e-8 ||
      !at_max(function(p) -gas_v_ll(p[1], p[2], exp(p[3]), tanh(p[4])),
              c(mg, cg[["(Intercept)"]], log(cg[["gas.kappa1"]]),
                atanh(cg[["gas.xi1"]]))) ||
      max(abs(g_map - rugarch::coef(ref$garch_n))) > 5e-3 ||
      abs(ll(fits$garch) - rugarch::likelihood(ref$garch_n)) > 0.1 ||
      !all(g_map[2:4] > 0)) {
    fail("the GARCH model as a scaled score-driven model")
  }
  # INGARCH(1, 1) and ACD(1, 1): the filters driven by y - mu and x - psi on
  # the identity link, written out, at a maximum, within 5e-3 of tscount's
  # and ACDm's estimates through the same mapping
  lin_ll <- function(z, dens) function(g0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_along(z)) {
      if (t > 1L) f <- a * (z[t - 1L] - (g0 + fp)) + b * f
      fp <- f
      m <- g0 + f
      if (!(m > 0)) return(-1e10)
      out <- out + dens(z[t], m)
    }
    out
  }
  ing_ll <- lin_ll(yd, function(y, m) stats::dpois(y, m, log = TRUE))
  acl_ll <- lin_ll(data$dur$x, function(x, m) stats::dexp(x, rate = 1 / m, log = TRUE))
  ci <- coef(fits$ingarch)$mu
  cl <- coef(fits$acd_lin)$mu
  lin_map <- function(cc) c(cc[["(Intercept)"]] * (1 - cc[["gas.xi1"]]),
                            cc[["gas.kappa1"]], cc[["gas.xi1"]] - cc[["gas.kappa1"]])
  lin_ok <- function(fun, cc, fit) {
    abs(fun(cc[["(Intercept)"]], cc[["gas.kappa1"]], cc[["gas.xi1"]]) - ll(fit)) < 1e-8 &&
      at_max(function(p) -fun(p[1], exp(p[2]), tanh(p[3])),
             c(cc[["(Intercept)"]], log(cc[["gas.kappa1"]]), atanh(cc[["gas.xi1"]])))
  }
  if (!lin_ok(ing_ll, ci, fits$ingarch) ||
      max(abs(lin_map(ci) - stats::coef(ref$tsc_id))) > 5e-3 ||
      abs(ll(fits$ingarch) - as.numeric(stats::logLik(ref$tsc_id))) > 1e-3) {
    fail("the INGARCH model against tscount")
  }
  if (!lin_ok(acl_ll, cl, fits$acd_lin) ||
      max(abs(lin_map(cl) - ref$acd_lin$mPar)) > 5e-3 ||
      abs(ll(fits$acd_lin) - ref$acd_lin$goodnessOfFit["LogLikelihood", "value"]) > 0.1) {
    fail("the linear ACD model against ACDm")
  }
  # the gamma filter on lh: shape 1/phi, mean mu, score (y/mu - 1)/phi
  yl <- data$lhd$y
  cgm <- coef(fits$lh)
  phl <- exp(cgm$phi[["(Intercept)"]])
  gm_ll <- function(c0, a, b) {
    f <- 0
    out <- 0
    for (t in seq_along(yl)) {
      if (t > 1L) f <- a * (yl[t - 1L] / exp(c0 + fp) - 1) / phl + b * f
      fp <- f
      out <- out + stats::dgamma(yl[t], shape = 1 / phl,
                                 scale = exp(c0 + f) * phl, log = TRUE)
    }
    out
  }
  if (abs(gm_ll(cgm$mu[["(Intercept)"]], cgm$mu[["gas.kappa1"]],
                cgm$mu[["gas.xi1"]]) - ll(fits$lh)) > 1e-8 ||
      !(AIC(fits$lh) < AIC(fits$lh0))) {
    fail("the gamma score-driven model")
  }
  # the beta filter on the front-seat share: score with respect to logit mu
  bt <- data$belts
  cb <- coef(fits$belt)
  phb <- exp(cb$phi[["(Intercept)"]])
  b_ll <- function(c0, cl, a, b) {
    f <- 0
    out <- 0
    for (t in seq_len(nrow(bt))) {
      if (t > 1L) {
        m0 <- stats::plogis(c0 + cl * bt$law[t - 1L] + fp)
        sc <- phb * m0 * (1 - m0) *
          (stats::qlogis(bt$share[t - 1L]) - (digamma(m0 * phb) - digamma((1 - m0) * phb)))
        f <- a * sc + b * f
      }
      fp <- f
      m <- stats::plogis(c0 + cl * bt$law[t] + f)
      out <- out + stats::dbeta(bt$share[t], m * phb, (1 - m) * phb, log = TRUE)
    }
    out
  }
  if (abs(b_ll(cb$mu[["(Intercept)"]], cb$mu[["law"]], cb$mu[["gas.kappa1"]],
               cb$mu[["gas.xi1"]]) - ll(fits$belt)) > 1e-8 ||
      !(cb$mu[["law"]] < 0) || !(AIC(fits$belt) < AIC(fits$belt0))) {
    fail("the beta score-driven model")
  }

  # --- panels: Ovary -----------------------------------------------------------
  ov <- data$ov
  # one filter per mare, each from f = w_i / (1 - b); mu = xb + f
  panel_ll <- function(xb, w, a, b, s) {
    out <- 0
    fl <- list()
    for (mm in levels(ov$Mare)) {
      rr <- which(ov$Mare == mm)
      rr <- rr[order(ov$Time[rr])]
      wi <- w[[mm]]
      f <- wi / (1 - b)
      for (t in seq_along(rr)) {
        if (t > 1L) f <- wi + a * (ov$follicles[rr[t - 1L]] - xb[rr[t - 1L]] - fp) / s^2 +
          b * fp
        fp <- f
        out <- out + stats::dnorm(ov$follicles[rr[t]], xb[rr[t]] + f, s, log = TRUE)
      }
      fl[[mm]] <- list(rows = rr, f = f, wi = wi)
    }
    list(ll = out, last = fl)
  }
  zero_w <- stats::setNames(as.list(rep(0, nlevels(ov$Mare))), levels(ov$Mare))
  co <- coef(fits$ov)$mu
  so <- exp(coef(fits$ov)$sigma[["(Intercept)"]])
  sn <- sin(2 * pi * ov$Time); cs <- cos(2 * pi * ov$Time)
  kS <- "sin(2 * pi * Time)"; kC <- "cos(2 * pi * Time)"
  xbo <- function(p) p[1] + p[2] * sn + p[3] * cs
  if (abs(panel_ll(xbo(co[c("(Intercept)", kS, kC)]), zero_w, co[["gas.kappa1"]],
                   co[["gas.xi1"]], so)$ll - ll(fits$ov)) > 1e-8 ||
      !at_max(function(p) -panel_ll(xbo(p[1:3]), zero_w, exp(p[4]), tanh(p[5]),
                                    exp(p[6]))$ll,
              c(co[["(Intercept)"]], co[[kS]], co[[kC]], log(co[["gas.kappa1"]]),
                atanh(co[["gas.xi1"]]), log(so)))) {
    fail("the shared panel filter")
  }
  ag <- stats::coef(ref$gls_ov$modelStruct$corStruct, unconstrained = FALSE)
  if (abs(co[["gas.xi1"]] - ag[["Phi1"]]) > 0.05 ||
      abs(co[["gas.kappa1"]] / so^2 - co[["gas.xi1"]] - ag[["Theta1"]]) > 0.05 ||
      abs(ll(fits$ov) - as.numeric(stats::logLik(ref$gls_ov))) > 2) {
    fail("the ARMA reading of the panel against gls()")
  }
  cm <- coef(fits$ov_om)$mu
  sm <- exp(coef(fits$ov_om)$sigma[["(Intercept)"]])
  w_om <- stats::setNames(lapply(levels(ov$Mare), function(mm)
    cm[["gas.omega.(Intercept)"]] + cm[[paste0("gas.omega.random.", mm)]]),
    levels(ov$Mare))
  xb_om <- cm[[kS]] * sn + cm[[kC]] * cs
  pom <- panel_ll(xb_om, w_om, cm[["gas.kappa1"]], cm[["gas.xi1"]], sm)
  if (abs(pom$ll - ll(fits$ov_om)) > 1e-8) fail("the panel with a level for each mare")
  # the shared model is the one with tau at zero: at a tiny tau the marginal
  # likelihood is the shared maximum
  f0 <- statmod(follicles ~ 0 + sin(2 * pi * Time) + cos(2 * pi * Time) +
                  gas(p = 1, q = 1, omega ~ 1 + random(~ 1 | Mare, hyper = c(sigma = 1e-5)),
                      by = Mare, time = Time),
                distrib = gaussian1_distrib(), data = ov)
  ns <- asNamespace("statmodels7")
  # integrated over the directions ml() integrates, the penalized ones alone
  d0 <- ns$statmod_design(f0@spec)
  m0 <- ns$statmod_marginal(f0@spec, d0, f0@coefficients, f0@hyper, ml(),
                            basis = ns$integrated_basis(f0@spec, d0, "ml"))
  lm1 <- as.numeric(logLik(fits$ov_ml, type = "marginal"))
  if (abs(m0$value - ll(fits$ov)) > 1e-3 ||
      rel(claims$lr_ov, 2 * (lm1 - ll(fits$ov))) > 1e-12 ||
      !(claims$lr_ov > 0) ||
      !(0.5 * stats::pchisq(claims$lr_ov, 1, lower.tail = FALSE) > 0.05)) {
    fail("the likelihood ratio for a level for each mare")
  }
  # the forecast of mare 1: the recursion without the score, from her state
  l1 <- pom$last[["1"]]
  r1 <- l1$rows
  nt <- claims$new_1
  xb1 <- cm[[kS]] * sin(2 * pi * nt$Time) + cm[[kC]] * cos(2 * pi * nt$Time)
  fh <- numeric(nrow(nt))
  fh[1] <- l1$wi + cm[["gas.kappa1"]] *
    (ov$follicles[r1[length(r1)]] - xb_om[r1[length(r1)]] - l1$f) / sm^2 +
    cm[["gas.xi1"]] * l1$f
  for (h in seq_along(fh)[-1]) fh[h] <- l1$wi + cm[["gas.xi1"]] * fh[h - 1L]
  if (rel(claims$fc_ov, xb1 + fh) > 1e-10) fail("the forecast of mare 1")
  # the random intercept goes to the shared model, the level of each mare has
  # the larger maximized marginal likelihood, and lev_ml is tau/(1 - beta)
  li <- as.numeric(logLik(fits$ov_int, type = "marginal"))
  if (!(claims$tau_int < 1e-2) ||
      !(lm1 > li) ||
      rel(claims$lev_ml, hyper(fits$ov_ml)$estimate /
            (1 - coef(fits$ov_ml)$mu[["gas.xi1"]])) > 1e-12 ||
      abs(li - ll(fits$ov)) > 0.1) {
    fail("the random intercept against the level of each mare")
  }

  # --- forecast ---------------------------------------------------------------
  b0 <- c1[["(Intercept)"]]
  f <- 0
  for (t in seq_len(n)) {
    if (t > 1L) f <- c1[["gas.kappa1"]] * (y[t - 1L] - b0 - fprev) / s2 +
      c1[["gas.xi1"]] * f
    fprev <- f
  }
  fn1 <- c1[["gas.kappa1"]] * (y[n] - b0 - f) / s2 + c1[["gas.xi1"]] * f
  hand_fc <- b0 + fn1 * c1[["gas.xi1"]]^(seq_along(claims$fc_nile) - 1L)
  if (rel(claims$fc_nile, hand_fc) > 1e-10) fail("the forecast of the Nile flow")
  # the standard error of the forecast: a warning, and the delta method on the
  # hand-written forecast differentiated numerically in the coordinates the
  # covariance is written in (intercept, log kappa, atanh xi, log sigma)
  wn <- NULL
  pse <- withCallingHandlers(
    predict(fits$nile, what = "mu", newdata = data.frame(year = 1971:1990), se = TRUE),
    warning = function(w) {
      wn <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })
  fc_of <- function(p) {
    b0 <- p[1]; a <- exp(p[2]); b <- tanh(p[3]); s2_ <- exp(2 * p[4])
    f <- 0
    for (t in seq_len(n)) {
      if (t > 1L) f <- a * (y[t - 1L] - b0 - fp) / s2_ + b * f
      fp <- f
    }
    b0 + (a * (y[n] - b0 - f) / s2_ + b * f) * b^(0:19)
  }
  Vn <- vcov(fits$nile, readable = FALSE)
  kn <- c("mu:(Intercept)", "mu:gas.kappa1", "mu:gas.pacf1", "sigma:(Intercept)")
  if (!all(kn %in% rownames(Vn))) fail("the coordinates of the Nile covariance")
  zn <- fits$nile@structural[[1]]$unconstrained
  p0 <- c(c1[["(Intercept)"]], zn[["kappa1"]], zn[["pacf1"]], log(s1))
  Jn <- numDeriv::jacobian(fc_of, p0)
  se_ref <- sqrt(rowSums((Jn[, seq_along(kn)] %*% Vn[kn, kn]) * Jn[, seq_along(kn)]))
  se_int <- confint(fits$nile)["mu:(Intercept)", "se"]
  if (is.null(wn) || !grepl("parameters alone", wn, fixed = TRUE) ||
      rel(pse$fit, hand_fc) > 1e-10 ||
      rel(pse$se, se_ref) > 1e-5 ||
      rel(claims$se_fc, pse$se) > 1e-12 ||
      !(pse$se[1] < max(pse$se)) ||
      abs(pse$se[20] - se_int) / se_int > 0.15) {
    fail("the standard error of the Nile forecast")
  }

  # --- regime -------------------------------------------------------------------
  # the forward recursion from the stationary distribution, and the smoothed
  # probabilities by a backward pass
  hmm <- function(D, P) {
    nn <- nrow(D); kk <- ncol(D)
    ev <- eigen(t(P))
    d0 <- Re(ev$vectors[, which.min(abs(ev$values - 1))])
    d0 <- d0 / sum(d0)
    al <- matrix(0, nn, kk)
    a <- d0 * D[1, ]
    lsum <- log(sum(a)); al[1, ] <- a / sum(a)
    for (t in seq_len(nn)[-1]) {
      a <- as.vector(al[t - 1L, ] %*% P) * D[t, ]
      lsum <- lsum + log(sum(a)); al[t, ] <- a / sum(a)
    }
    be <- matrix(1, nn, kk)
    for (t in rev(seq_len(nn - 1L))) {
      b <- as.vector(P %*% (D[t + 1L, ] * be[t + 1L, ]))
      be[t, ] <- b / sum(b)
    }
    g <- al * be
    list(ll = lsum, gamma = g / rowSums(g))
  }
  # the transition matrix from additive log-ratios, row by row
  alr_P <- function(z, kk) {
    t(vapply(seq_len(kk), function(i) {
      e <- c(exp(z[sprintf("alr%d.%d", i, seq_len(kk - 1L))]), 1)
      e / sum(e)
    }, numeric(kk)))
  }
  levels_of <- function(z, kk) {
    cumsum(c(0, exp(z[sprintf("gap%d", seq.int(2L, kk))])))
  }
  check_regime <- function(fit, yv, mu_static, s_static, kk, where, eq = "mu") {
    z <- fit@structural[[1]]$unconstrained
    P <- alr_P(z, kk)
    lv <- levels_of(z, kk)
    D <- vapply(seq_len(kk), function(j) {
      if (eq == "mu") stats::dnorm(yv, mu_static + lv[j], s_static)
      else stats::dnorm(yv, mu_static, exp(s_static + lv[j]))
    }, numeric(length(yv)))
    h <- hmm(D, P)
    pn <- sprintf("%s:regime.p%d.%d", eq, rep(seq_len(kk), each = kk),
                  rep(seq_len(kk), kk))
    ci <- confint(fit)
    if (abs(h$ll - ll(fit)) > 1e-7 ||
        rel(ci[pn, "estimate"], as.numeric(t(P))) > 1e-10) {
      fail(where)
    }
    list(h = h, ci = ci, pn = pn, P = P, lv = lv)
  }
  gy <- data$geyser$waiting
  cg <- coef(fits$gey)
  rg <- check_regime(fits$gey, gy, cg$mu[["(Intercept)"]],
                     exp(cg$sigma[["(Intercept)"]]), 2L, "the two-state geyser chain")
  # the first row at the edge, without a standard error; the second with one
  if (!(abs(fits$gey@structural[[1]]$unconstrained[["alr1.1"]]) > 8) ||
      !all(is.na(rg$ci[rg$pn[1:2], c("se", "lower", "upper")])) ||
      !all(is.finite(rg$ci[rg$pn[3:4], "se"])) ||
      !(rg$P[1, 2] > 0.999) ||
      rel(as.matrix(claims$st_gey), rg$h$gamma) > 1e-8 ||
      rel(predict(fits$gey, what = "mu"),
          cg$mu[["(Intercept)"]] + as.vector(rg$h$gamma %*% rg$lv)) > 1e-8) {
    fail("the edge and the states of the geyser chain")
  }
  # depmixS4's forward recursion at our estimates, its gap from the states,
  # and its own fit with a standard deviation for each state
  if (abs(ref$fb_gey$logLike - ll(fits$gey)) > 1e-6 ||
      rel(claims$gap_fb, max(abs(rg$h$gamma - ref$fb_gey$gamma))) > 1e-6 ||
      !(claims$gap_fb < 1e-10) ||
      abs(claims$bic_hmm - (-2 * depmixS4::forwardbackward(ref$hmm_fit)$logLike +
                              log(length(gy)) * depmixS4::freepars(ref$hmm_fit))) > 1e-8 ||
      !(claims$bic_hmm < BIC(fits$gey))) {
    fail("the comparison with depmixS4")
  }
  cg3 <- coef(fits$gey3)
  check_regime(fits$gey3, gy, cg3$mu[["(Intercept)"]],
               exp(cg3$sigma[["(Intercept)"]]), 3L, "the three-state geyser chain")
  if (!(BIC(fits$gey3) < BIC(fits$gey)) ||
      abs(BIC(fits$gey3) - (-2 * ll(fits$gey3) + log(299) * 10)) > 1e-6) {
    fail("the choice of three states")
  }
  gd2 <- data$gey
  cgd <- coef(fits$gey_d)
  check_regime(fits$gey_d, gd2$waiting,
               cgd$mu[["(Intercept)"]] + cgd$mu[["dlag"]] * gd2$dlag,
               exp(cgd$sigma[["(Intercept)"]]), 2L, "the geyser chain with a covariate")
  if (!(BIC(fits$gey_d) < BIC(fits$gey_d0)) ||
      !(cgd$mu[["regime.gap2"]] < cg$mu[["regime.gap2"]])) {
    fail("the covariate beside the chain")
  }
  cdr <- coef(fits$dax_r)
  rd <- check_regime(fits$dax_r, r, cdr$mu[["(Intercept)"]],
                     cdr$sigma[["(Intercept)"]], 2L, "the chain on the DAX scale",
                     eq = "sigma")
  if (rel(as.matrix(claims$st_dax), rd$h$gamma) > 1e-8 ||
      !(AIC(fits$dax_r) < AIC(fits$dax)) ||
      !(rd$P[1, 1] > 0.9 && rd$P[2, 2] > 0.9) ||
      rel(predict(fits$dax_r, what = "sigma"),
          exp(cdr$sigma[["(Intercept)"]] + as.vector(rd$h$gamma %*% rd$lv))) > 1e-8) {
    fail("the two volatility states")
  }
  invisible(TRUE)
}
