# Certificates for the statmodels7 chapter.
#
# Same contract as the other gates: everything the chapter asserts is
# re-derived here from a route the layer does not itself take, and the render
# stops if any of them disagrees. Nothing in this file is visible to the
# reader.
#
# Five claims are checked.
#
#   1. the penalized mode of a fit whose objective is quadratic agrees with
#      the closed form. With the scale held the objective is
#      (1/2)|y - Z b|^2 + (lambda/2)|b_R|^2, so the mode is
#      (Z'Z + S)^{-1} Z'y, and Z and S are transcribed by hand rather than
#      read off the design -- comparing the fit with the fit's own assembly
#      would prove nothing about the assembly. The scoring loop reaches the
#      same point through an augmented QR, which shares no arithmetic with a
#      solve of the normal equations;
#   2. the exact gradient of the marginal criterion agrees with a central
#      difference of that criterion, the mode refitted at every probe, read
#      AWAY from the optimum where the gradient is not nearly zero;
#   3. the effective degrees of freedom agree with tr[(H+S)^-1 H] assembled
#      here out of the family's expected information and a penalty Hessian
#      taken straight from penalties7;
#   4. a structural level and a linear intercept in the same equation are
#      exactly confounded, so holding the level costs nothing: the two
#      parametrizations reach the same maximum and the held fit's intercept
#      is the free fit's stationary level omega/(1 - sum b);
#   5. prediction reapplies a term's blueprint rather than rebuilding it,
#      pinned by the identity that predicting on rows the model was fitted to
#      returns the fitted values there. The subset matters: a rebuilt basis
#      agrees with the fitted values on the WHOLE data whatever its knots
#      are, so a check on all the rows passes for a defect it exists to
#      catch.
#
# Injection-checked when written, each caught and the gate clean again with
# the injection removed. The agreements asserted, against the tolerance each
# is held to:
#
#   1  the penalized mode                6.661e-16   (1e-8)
#      its stationarity residual         4.263e-14
#   2  the criterion's exact gradient    7.877e-07   (1e-4, relative)
#   3  the effective degrees of freedom  1.776e-15   (1e-8)
#   4  the two parametrizations' maxima  0.000e+00   (1e-3)
#      the intercept against omega/(1-b) 3.573e-12   (1e-2, relative)
#      the held level                    0.000e+00
#   5  prediction on fitted rows         0.000e+00   (1e-8)
#
# and the injections:
#
#   1a the closed form computed at 2*lambda                     2.294e-02
#   1b the penalty applied to the intercept as well             3.785e-02
#   2a the gradient read on the parameter scale, not the free   1.000e+00
#   2b the gradient read at the optimum instead of at lambda0   1.000e+00
#   2c the mode taken at the wrong lambda                       2.264e-02
#   3a the observed information in place of the expected        3.338e-04
#   3b the penalty omitted, S = 0                               9.702e-01
#   3c lambda doubled in the penalty Hessian                    4.230e-01
#   4a the stationary level read as omega alone                 4.650e-01
#   4b the stationary level read as omega/(1+b)                 5.812e-01
#   4c the FREE fit's omega compared against zero               3.098e-01
#   5a the prediction read at rows shifted by one               2.066e-01
#   5b the smooth's block rebuilt on the subset, not reapplied  2.895e+00
#
# 5b is the defect check 5 exists for, measured on the block itself rather
# than on the fit. 3a is not a defect but a different quantity, and it is
# what says the count is read on the information the fit inverted.
#
# The gate makes eight statmod() calls and takes 3.4 s from cold, 0.9 s with
# the caches warm. That was worth checking rather than assuming: the warm
# reading is faster than the four-check version's cold one, which looks like
# a check that stopped running, and counting the fits is what said otherwise.

assert_statmod_ok <- function() {

  fail <- function(what) {
    stop("statmodels7 chapter: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }

  # --- 1. a penalized mode that is available in closed form ----------------

  set.seed(20260828)
  n <- 80
  dd <- data.frame(x = stats::runif(n, -2, 2))
  dd$R <- matrix(stats::rnorm(n * 3), n, 3)
  dd$y <- 1 + 0.8 * dd$x +
    as.numeric(dd$R %*% c(0.6, -0.4, 0)) + stats::rnorm(n, sd = 0.7)

  lam <- 2.5
  # The scale is held, so the equation is one weighted least squares problem
  # with unit weights and the penalized objective is quadratic in beta.
  fam <- distributions7::fixed(distributions7::gaussian1_distrib(), sigma = 1)
  fit <- statmod(y ~ x + ridge(R, lambda = lam), fam, dd)

  # Transcribed from the printed model, not read off the fit: an unpenalized
  # intercept and slope, then the three penalized columns.
  Z <- cbind(1, dd$x, dd$R)
  S <- diag(c(0, 0, rep(lam, 3)))
  ref <- as.numeric(solve(crossprod(Z) + S, crossprod(Z, dd$y)))
  got <- unname(coef(fit)$mu)

  if (length(got) != length(ref) ||
      max(abs(got - ref)) > 1e-8 * max(1, max(abs(ref)))) {
    fail("the penalized mode against its closed form")
  }

  # The same statement read as a stationarity condition, which needs no
  # matrix inverse: the gradient of the penalized objective vanishes there.
  if (max(abs(crossprod(Z, dd$y - Z %*% got) - S %*% got)) > 1e-6) {
    fail("the stationarity of the reported mode")
  }

  # --- 2. the exact gradient of the marginal criterion, and the edf --------

  set.seed(20260830)
  nb <- 250
  dc <- data.frame(z = sort(stats::runif(nb, -3, 3)))
  dc$y <- sin(1.4 * dc$z) + stats::rnorm(nb, sd = 0.4)
  gau <- distributions7::gaussian1_distrib()

  fb <- statmod(y ~ s(z, k = 8), gau, dc)
  sp <- fb@spec
  de <- statmod_design(sp)
  key <- "s(z, k = 8)"
  idx <- statmodels7:::outer_hyper_index(
    sp, statmodels7:::statmod_blocks(sp, de))

  # The criterion at any smoothing parameter: the mode there, then the
  # Laplace value read at it. The base point is deliberately away from the
  # optimum -- at the optimum the gradient is zero and a relative comparison
  # would be between two nearly-zero numbers, which is the trap the fitting
  # gate of chapter 3 records.
  hy_at <- function(lam) {
    h <- fb@hyper
    h[["mu"]][[key]][["lambda"]] <- lam
    h
  }
  mode_at <- function(lam) {
    statmod(stats::as.formula(
      sprintf("y ~ s(z, k = 8, lambda = %.17g)", lam)), gau, dc)@coefficients
  }
  crit_at <- function(lam) {
    v <- statmodels7:::statmod_marginal(sp, de, mode_at(lam), hy_at(lam),
                                        reml())
    if (is.null(v)) NA_real_ else v$value
  }

  lam0 <- 0.5
  g_exact <- statmodels7:::statmod_marginal_grad(
    sp, de, mode_at(lam0), hy_at(lam0), reml(), idx)
  h <- 3e-3
  e0 <- log(lam0)
  g_fd <- (crit_at(exp(e0 + h)) - crit_at(exp(e0 - h))) / (2 * h)
  if (!is.finite(g_exact) || !is.finite(g_fd) ||
      abs(g_fd - g_exact) / max(1, abs(g_exact)) > 1e-4) {
    fail("the exact gradient of the marginal criterion")
  }

  # The effective degrees of freedom, against the trace assembled here from
  # the family's information and a penalty Hessian built straight from
  # penalties7. The information is the EXPECTED one, which is what the fit
  # inverted; reading the observed one instead moves the count by 3.3e-04,
  # which is a different quantity rather than a worse one.
  H <- -hessian(fb, expected = TRUE)
  cols <- de$mu$blocks[[key]]
  ptot <- sum(vapply(de, function(e) e$npar, integer(1)))
  S <- matrix(0, ptot, ptot)
  S[cols, cols] <- as.matrix(penalties7::penalty_hessian(
    modelterms7::term_penalty(sp@terms$mu[[key]]),
    fb@coefficients$mu[cols],
    list(lambda = fb@hyper$mu[[key]][["lambda"]])))
  if (abs(sum(diag(solve(H + S, H))) - sum(fb@edf$edf)) > 1e-8) {
    fail("the effective degrees of freedom")
  }

  # --- 3. a structural level and an intercept are the same model -----------

  # The chapter claims the two are exactly confounded and that holding the
  # term's own level costs nothing. The check is end to end and uses only
  # public calls: the two parametrizations must reach the same maximum, and
  # the intercept the held fit reports must BE the free fit's stationary
  # level omega/(1 - sum b), which is a different arithmetic route to the
  # same number.
  set.seed(20260831)
  ng <- 400
  om <- 0.3; a1 <- 0.4; b1 <- 0.6
  lev <- numeric(ng); yg <- numeric(ng); lev[1] <- om / (1 - b1)
  for (t in seq_len(ng)) {
    yg[t] <- stats::rnorm(1, lev[t], 1)
    if (t < ng) lev[t + 1] <- om + a1 * (yg[t] - lev[t]) + b1 * lev[t]
  }
  dg <- data.frame(y = yg)

  held <- statmod(y ~ gas(p = 1, q = 1), gau, dg)
  free <- statmod(y ~ gas(p = 1, q = 1) - 1, gau, dg)

  if (abs(as.numeric(logLik(held)) - as.numeric(logLik(free))) > 1e-3) {
    fail("the two parametrizations of a structural level")
  }
  pr <- free@structural[[1]]$parameter
  stationary <- pr[["omega"]] / (1 - pr[["pacf1"]])
  if (abs(unname(coef(held)$mu[["(Intercept)"]]) - stationary) >
      1e-2 * max(1, abs(stationary))) {
    fail("the held intercept against the free fit's stationary level")
  }
  # and the term's own level really is held at zero where the intercept won
  if (abs(held@structural[[1]]$parameter[["omega"]]) > 1e-12) {
    fail("the level held at zero beside an intercept")
  }

  # --- 4. prediction reapplies the blueprint -------------------------------

  set.seed(20260829)
  m <- 200
  db <- data.frame(z = sort(stats::runif(m, -3, 3)))
  db$y <- sin(1.4 * db$z) + stats::rnorm(m, sd = 0.4)

  # The smoothing parameter is held. What is being checked is the blueprint,
  # not the search over it, and holding it keeps the check independent of
  # where a criterion happens to stop.
  sf <- statmod(y ~ s(z, k = 8, lambda = 2),
                distributions7::gaussian1_distrib(), db)

  # A subset over which a rebuilt basis would put its knots somewhere else.
  rows <- which(abs(db$z) < 0.5)
  if (length(rows) < 10) fail("the subset the prediction check needs")
  pr <- predict(sf, "mu", newdata = db[rows, , drop = FALSE])
  if (max(abs(pr - sf@fitted$mu[rows])) > 1e-8) {
    fail("prediction on rows the model was fitted to")
  }

  invisible(TRUE)
}
