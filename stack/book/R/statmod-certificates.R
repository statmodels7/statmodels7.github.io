# Certificates for the statmodels7 chapter.
#
# Same contract as the other gates: everything the chapter asserts is
# re-derived here from a route the layer does not itself take, and the render
# stops if any of them disagrees. Nothing in this file is visible to the
# reader.
#
# Two claims are checked, one per section written so far.
#
#   1. the penalized mode of a fit whose objective is quadratic agrees with
#      the closed form. With the scale held the objective is
#      (1/2)|y - Z b|^2 + (lambda/2)|b_R|^2, so the mode is
#      (Z'Z + S)^{-1} Z'y, and Z and S are transcribed by hand rather than
#      read off the design -- comparing the fit with the fit's own assembly
#      would prove nothing about the assembly. The scoring loop reaches the
#      same point through an augmented QR, which shares no arithmetic with a
#      solve of the normal equations;
#   2. prediction reapplies a term's blueprint rather than rebuilding it,
#      pinned by the identity that predicting on rows the model was fitted to
#      returns the fitted values there. The subset matters: a rebuilt basis
#      agrees with the fitted values on the WHOLE data whatever its knots
#      are, so a check on all the rows passes for a defect it exists to
#      catch.
#
# Injection-checked when written, against a tolerance of 1e-8, each caught
# and the gate clean again with the injection removed:
#
#   the closed form computed at 2*lambda                        2.294e-02
#   the penalty applied to the intercept as well as the block   3.785e-02
#   the prediction read at rows shifted by one                  2.066e-01
#   the smooth's block rebuilt on the subset, not reapplied     2.895e+00
#
# against the agreements the gate asserts, which are 6.661e-16 on the mode
# and exactly 0 on the prediction. The last line is the defect check 2
# exists for, measured on the block itself rather than on the fit.

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

  # --- 2. prediction reapplies the blueprint -------------------------------

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
