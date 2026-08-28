# Certificates for the modelterms7 chapter.
#
# Each claim the chapter prints is checked against a route the package does
# not use: blocks against hand-built references, the penalized terms' penalty
# against the penalties7 constructor called directly, the random block and
# its Kronecker-replicated precision against assembled matrices, and the edf
# trace against the eigenvalue form, and the two break-point constructions
# against the printed formulas transcribed by hand. Injection-checked when
# written: a block construction transposed in the reference, a 5% inflation
# of the edf, a penalty compared at a different theta, the sign of the
# break-point column, the jump identity 2.5% wrong, a break-point 2% wrong
# and psi read off the wrong combination were each caught.
#
# Extended 2026-08-28 with four more, and their measurements. The
# agreements: standardization leaves the design block untouched and the
# standardized penalty at beta equals the plain one at the scaled beta to
# 0.000e+00; the three states of a hyperparameter argument -- free, held,
# a written-out grid -- are distinct; a score-driven term's simulated path
# is reproduced to 0.000e+00 by running the filter on the response it drew,
# which is the identity that says it drew AS the recursion ran rather than
# after it; and a smoothed step agrees with the sharp one to 0.000e+00
# twenty to forty widths out, with the smoothed construction reporting a
# Jacobian block where the sharp discontinuous one does not.
#
# Their injections, each caught: the standardized penalty compared at the
# unscaled beta 1.070e+00 and at inverted spreads 4.484e+00; the filter
# re-run with alpha1 off by 0.01, 3.347e-02; the smoothed step read a
# quarter of a width from the corner instead of twenty, 4.013e-01.

assert_terms_ok <- function() {
  set.seed(20260808)
  n <- 30
  dd <- data.frame(x1 = rnorm(n), x2 = runif(n),
                   g = factor(rep(letters[1:5], 6)))
  dd$L <- matrix(rnorm(n * 3), n, 3)

  fail <- function(what) {
    stop("modelterms7 chapter: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }

  # 1. the collected parametric block is the GLM design matrix
  out <- interpret_formula(x1 ~ x2 + log(x2) + lasso(L) + random(~ x1 | g), dd)
  lp <- term_build(out$terms$linpar, dd)
  ref <- stats::model.matrix(~ x2 + log(x2), dd)
  if (!isTRUE(all.equal(unname(term_matrix(lp)), unname(ref),
                        tolerance = 1e-12, check.attributes = FALSE))) {
    fail("the collected linpar block")
  }

  # 2. the lasso term's penalty is the penalties7 object, value for value
  bl <- term_build(out$terms[["lasso(L)"]], dd)
  pen <- term_penalty(bl)
  refp <- penalties7::lasso_penalty(n_coef = 3)
  beta <- c(0.7, -0.2, 1.1)
  if (!identical(penalties7::penalty_value(pen, beta, list(lambda = 1.3)),
                 penalties7::penalty_value(refp, beta, list(lambda = 1.3)))) {
    fail("the lasso term's penalty")
  }
  if (term_smooth(bl) || !term_smooth(lp)) {
    fail("the smoothness split")
  }

  # 3. the random block against a hand-built interaction, and its precision
  #    against the assembled Kronecker product
  br <- term_build(out$terms[["random(~x1 | g)"]], dd)
  m <- nlevels(dd$g)
  Zref <- matrix(0, n, 2 * m)
  for (i in seq_len(n)) {
    j <- as.integer(dd$g[i])
    Zref[i, (j - 1) * 2 + 1] <- 1
    Zref[i, (j - 1) * 2 + 2] <- dd$x1[i]
  }
  # the block is built sparse -- a row belongs to one group, so the density is
  # 1/m -- and what is asserted here is which column each observation lands
  # in, which is a statement about the entries and not about the container
  Zgot <- as.matrix(term_matrix(br))
  if (!isTRUE(all.equal(unname(Zgot), Zref, tolerance = 1e-12,
                        check.attributes = FALSE))) {
    fail("the random-effect block")
  }
  st <- parameters7::log_cholesky(2)
  eta <- c(0.2, -0.1, 0.3)
  full <- parameters7::param_value(
    parameters7::kron_identity(st, m), eta)
  if (!isTRUE(all.equal(unname(unclass(full)),
                        kronecker(diag(m),
                                  unclass(parameters7::param_value(st, eta))),
                        tolerance = 1e-12, check.attributes = FALSE))) {
    fail("the Kronecker replication of the per-group precision")
  }

  # 4. the printed edf trace against the eigenvalue form, computed apart
  brg <- term_build(ridge(~ x1 + x2), dd)
  H <- crossprod(term_matrix(brg))
  # the ridge's hyperparameter IS the precision, so it enters the trace as it
  # stands rather than as one over a square
  lambda <- 1 / 1.4^2
  got <- edf(brg, coef = c(0.5, -1), hessian = H,
             theta = list(lambda = lambda))
  d <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (abs(got - sum(d / (d + lambda))) > 1e-10) {
    fail("the edf trace")
  }
  if (edf(term_build(lasso(~x1), dd), coef = 0) != 0 ||
      edf(lp) != ncol(term_matrix(lp))) {
    fail("the edf counting rules")
  }

  # 5. the two break-point constructions, transcribed by hand from the
  #    printed formulas rather than read off the package
  sx <- data.frame(x = seq(0.5, 9.5, length.out = 41))
  bsg <- term_build(seg(x, psi = 4.7), sx)
  cfs <- c(0.4, 1.7, 5.3)
  Xs <- term_matrix(term_refresh(bsg, cfs))
  jac <- cbind(sx$x, pmax(sx$x - cfs[3], 0),
               -cfs[2] * (sx$x > cfs[3]))
  if (!isTRUE(all.equal(unname(Xs), jac, tolerance = 1e-12,
                        check.attributes = FALSE))) {
    fail("the segmented Jacobian")
  }

  # the identity 1(x>psi) = 1/2 + (x-psi)/(2|x-psi|), away from the band
  bjp <- term_build(jump(x, psi = 4), sx)
  kap <- 2.5
  psi <- 4
  cfj <- c(kap, -kap * psi)
  Xj <- term_matrix(term_refresh(bjp, cfj))
  away <- abs(sx$x - psi) > 0.02 * diff(range(sx$x))
  step <- 0.5 + (sx$x - psi) / (2 * abs(sx$x - psi))
  if (!isTRUE(all.equal(as.numeric(Xj %*% cfj)[away],
                        (kap * step)[away], tolerance = 1e-10))) {
    fail("the jump identity")
  }
  if (abs(seg_psi(bjp, cfj) - psi) > 1e-12) {
    fail("the break-point read off two coefficients")
  }

  # 6. the blueprint: a droplevels subset reproduces the rows
  idx <- which(dd$g != "e")
  sub <- droplevels(dd[idx, , drop = FALSE])
  if (!isTRUE(all.equal(unname(term_predict(lp, sub)),
                        unname(term_matrix(lp)[idx, , drop = FALSE]),
                        tolerance = 1e-12, check.attributes = FALSE))) {
    fail("blueprint prediction on a subset")
  }

  # 7. standardization is a map on the PENALTY and not on the design
  bp <- term_build(lasso(~ x1 + x2), dd)
  bs <- term_build(lasso(~ x1 + x2, standardize = TRUE), dd)
  if (!isTRUE(all.equal(unname(as.matrix(term_matrix(bp))),
                        unname(as.matrix(term_matrix(bs))),
                        tolerance = 1e-12, check.attributes = FALSE))) {
    fail("standardization left the design block alone")
  }
  sdev <- bs@blueprint$standardize
  bb <- c(0.7, -1.3)
  th <- list(lambda = 1.1)
  # the standardized penalty at beta is the plain one at the scaled beta,
  # built here straight from penalties7 rather than read off the term
  plain <- penalties7::lasso_penalty(n_coef = 2)
  if (abs(penalties7::penalty_value(term_penalty(bs), bb, th) -
          penalties7::penalty_value(plain, sdev * bb, th)) > 1e-10) {
    fail("the standardized penalty against the plain one at scaled coefficients")
  }

  # 8. the three states of a hyperparameter argument are distinct
  if (length(term_hyper(lasso(~x1))) != 0L ||
      term_hyper(lasso(~x1, lambda = 3))[[1]]$lambda != 3 ||
      length(term_values(lasso(~x1, lambda = c(0.1, 1, 10)))[[1]]$lambda) != 3L ||
      length(term_hyper(lasso(~x1, lambda = c(0.1, 1, 10)))) != 0L) {
    fail("the three states of a hyperparameter argument")
  }

  # 9. a score-driven term draws AS the recursion runs, so filtering the
  #    response it produced reproduces the predictor it produced. The two
  #    routes share no arithmetic beyond the recursion itself: one draws
  #    forward, the other reads a response it is given.
  dt <- data.frame(t = 1:200)
  tg <- term_build(gas(p = 1, q = 1, time = t), dt)
  psi <- list(omega = 0.4, alpha1 = 0.3, pacf1 = 0.6)
  set.seed(77)
  sim <- term_simulate(tg, psi, rep(0, 200),
                       draw = function(e, i) stats::rnorm(1, e, 1))
  back <- term_filter(tg, eta = rep(0, 200), y = sim$y,
                      score = function(e, i) sim$y[i] - e,
                      curvature = function(e, i) -1, psi = psi)
  if (max(abs(back$eta - sim$eta)) > 1e-10) {
    fail("the filter's simulated path against the filter run on its own draw")
  }

  # 10. the smoothed construction is a Jacobian block, and away from the
  #     break-point its step agrees with the sharp one
  sm <- penalties7::smooth_probit()
  bsm <- term_build(seg(x, psi = 4.7, smoothed = sm), sx)
  if (!term_jacobian_block(bsm) || term_jacobian_block(term_build(jump(x), sx))) {
    fail("which constructions report a Jacobian block")
  }
  # the width the term resolved from the covariate's spacing, which is what
  # it reports and what the smoothing bias lives on
  hw <- bsm@blueprint$smooth$width
  uu <- c(-40, -20, 20, 40) * hw
  sharp <- as.numeric(uu > 0)
  smoothstep <- (1 + penalties7::smoother_deriv(sm, uu, hw, 1L)) / 2
  if (max(abs(smoothstep - sharp)) > 1e-6) {
    fail("the smoothed step against the sharp one away from the corner")
  }

  invisible(TRUE)
}
