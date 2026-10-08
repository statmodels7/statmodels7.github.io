# Certificates for Part I, chapter 1 (a first model, the motorcycle data).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   summary      the smooth has k - 1 columns summing to zero over the
#                observed times, its first column `lin` is the time
#                standardized to mean zero and sd one, and the interval of
#                the smoothing parameter is symmetric on the log scale and
#                positive; hyper() reports the lambda the summary prints;
#   coef         vcov() is 21 by 21, its block between mu and sigma is
#                exactly zero and so is the block of the expected information,
#                the information of the mean is X' diag(1/sigma^2) X; the
#                confint() rows are the estimate +- qnorm(0.975) se with se the
#                square root of the diagonal of vcov();
#   residuals    fitted() is predict(mu); the quantile residuals are
#                qnorm(pnorm(y, mu, sigma)), taken from stats::pnorm; the
#                residuals of the first model outside +-2 all lie between 15
#                and 40 ms; those of the second model have sd within
#                [0.7, 1.2] in each of the three phases;
#   predict      the interval of the mean is fit +- 1.96 se; the interval of
#                sigma is symmetric on the log scale and positive, and in the
#                first model constant;
#   sigma        the linear part of log sigma is positive with p < 0.001;
#   compare      the log-likelihood is the sum of stats::dnorm with mu and
#                sigma read back through predict(); each term's edf is the
#                trace of its block of (H + S)^{-1} H with S assembled here
#                from penalties7; cAIC and cBIC are -2 l + 2 tau and
#                -2 l + log(n) tau; the second model has the lower cAIC and
#                cBIC; the constant sigma is sqrt(RSS/(n - tau_mu)) and larger
#                than sqrt(RSS/n);
#   certificate  both fits are `converged`, with decrement and mode error
#                below 1e-6 and nothing at a boundary;
#   simulate     the truth rstatmod() returns is predict() of the second
#                model, with no hyperparameter reported; simulate() returns
#                three columns of the right length.

assert_model_ok <- function(f0, f1, data, claims) {

  fail <- function(what) {
    stop("Part I chapter 1: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  z <- stats::qnorm(0.975)

  # --- summary: the constrained smooth and its linear column --------------
  sp <- f1@spec
  de <- statmod_design(sp)
  Xmu <- as.matrix(de$mu$X)
  key <- names(de$mu$blocks)[vapply(de$mu$blocks, length, integer(1)) > 1L]
  kk <- as.integer(sub(".*k = ([0-9]+).*", "\\1", key))
  Xs <- Xmu[, de$mu$blocks[[key]], drop = FALSE]
  if (ncol(Xs) != kk - 1L) fail("the number of coefficients of the smooth")
  if (max(abs(colSums(Xs))) > 1e-10 * nrow(Xs) * max(abs(Xs))) {
    fail("the sum-to-zero constraint of the smooth")
  }
  lin <- Xs[, 1]
  if (!endsWith(colnames(Xmu)[de$mu$blocks[[key]][1]], ".lin") ||
      abs(mean(lin)) > 1e-10 || abs(stats::sd(lin) - 1) > 1e-10 ||
      max(abs(stats::resid(stats::lm(lin ~ data$times)))) > 1e-10) {
    fail("the standardized linear column")
  }
  tb <- summary(f0)@tables$mu[[2]]$table
  lam <- tb[tb$name == "lambda", ]
  if (nrow(lam) != 1L || !(lam$lower > 0) ||
      abs((log(lam$lower) + log(lam$upper)) / 2 - log(lam$estimate)) > 1e-8) {
    fail("the log-scale interval of the smoothing parameter")
  }
  h0 <- hyper(f0)
  if (nrow(h0) != 1L || rel(h0$estimate, lam$estimate) > 1e-10 ||
      h0$source != "reml" || isTRUE(h0$held)) {
    fail("hyper()")
  }

  # --- coef: vcov, the zero block, confint --------------------------------
  V <- stats::vcov(f0)
  if (!identical(dim(V), c(21L, 21L))) fail("the size of vcov()")
  imu <- grep("^mu:", rownames(V)); isg <- grep("^sigma:", rownames(V))
  if (any(V[imu, isg] != 0)) fail("the zero block of vcov()")
  H <- -hessian(f1, expected = TRUE)
  offs <- cumsum(c(0L, vapply(de, function(e) e$npar, integer(1))))
  i_mu <- offs[1] + seq_len(de$mu$npar)
  i_sg <- offs[2] + seq_len(de$sigma$npar)
  if (max(abs(H[i_mu, i_sg])) > 1e-12 * max(abs(H))) {
    fail("the zero block between the mean and the spread")
  }
  sg1 <- predict(f1, what = "sigma", data)
  if (rel(H[i_mu, i_mu], crossprod(Xmu, Xmu / sg1^2)) > 1e-8) {
    fail("the weights of the information of the mean")
  }
  ci0 <- claims$ci0
  if (rel(ci0$se, sqrt(diag(V[claims$keep, claims$keep]))) > 1e-10 ||
      rel(ci0$lower, ci0$estimate - z * ci0$se) > 1e-10 ||
      rel(ci0$upper, ci0$estimate + z * ci0$se) > 1e-10) {
    fail("confint()")
  }

  # --- residuals ------------------------------------------------------------
  mu0 <- predict(f0, what = "mu", data); sg0 <- predict(f0, what = "sigma", data)
  if (rel(stats::fitted(f0), mu0) > 1e-12) fail("fitted()")
  r0 <- stats::qnorm(stats::pnorm(data$accel, mu0, sg0))
  if (rel(stats::residuals(f0), r0) > 1e-8) fail("the quantile residuals")
  middle <- data$times >= 15 & data$times < 40
  if (!all(middle[abs(r0) > 2])) fail("the residuals outside +-2")
  r1 <- stats::residuals(f1)
  for (ph in list(data$times < 15, middle, data$times >= 40)) {
    s <- stats::sd(r1[ph])
    if (s < 0.7 || s > 1.2) fail("the residuals of the second model")
  }

  # --- predict ---------------------------------------------------------------
  p0 <- claims$p0
  if (rel(c(p0$lower, p0$upper), c(p0$fit - z * p0$se, p0$fit + z * p0$se)) >
      1e-10) fail("the interval of the mean")
  for (f in list(f0, f1)) {
    ps <- predict(f, what = "sigma", claims$new, se = TRUE)
    if (any(ps$lower <= 0) ||
        max(abs((log(ps$lower) + log(ps$upper)) / 2 - log(ps$fit))) > 1e-8) {
      fail("the log-scale interval of sigma")
    }
  }
  if (diff(range(predict(f0, what = "sigma", claims$grid))) > 1e-10) {
    fail("the constant sigma of the first model")
  }

  # --- sigma: the linear part of log sigma ---------------------------------
  ts <- summary(f1)@tables$sigma[[2]]$table
  ln <- ts[endsWith(ts$name, ".lin"), ]
  if (nrow(ln) != 1L || !(ln$estimate > 0) || !(ln$p_value < 1e-3)) {
    fail("the linear part of log sigma")
  }

  # --- compare ---------------------------------------------------------------
  for (f in list(f0, f1)) {
    ll <- sum(stats::dnorm(data$accel, predict(f, what = "mu", data),
                           predict(f, what = "sigma", data), log = TRUE))
    if (rel(ll, as.numeric(logLik(f))) > 1e-10) fail("the log-likelihood")
    e <- sum(f@edf$edf)
    l <- as.numeric(logLik(f))
    if (rel(-2 * l + 2 * e, AIC(f)) > 1e-10) fail("cAIC")
    if (rel(-2 * l + log(nrow(data)) * e, BIC(f)) > 1e-10) fail("cBIC")
  }
  if (!(AIC(f1) < AIC(f0)) || !(BIC(f1) < BIC(f0))) {
    fail("the comparison of the two models")
  }
  ptot <- nrow(H)
  S <- matrix(0, ptot, ptot)
  shares <- list()
  for (k in seq_along(de)) {
    par <- names(de)[k]
    for (bk in names(de[[k]]$blocks)) {
      tm <- sp@terms[[par]][[bk]]
      pen <- modelterms7::term_penalty(tm)
      cols <- offs[k] + de[[k]]$blocks[[bk]]
      if (!is.null(pen)) {
        S[cols, cols] <- as.matrix(penalties7::penalty_hessian(
          pen, f1@coefficients[[par]][de[[k]]$blocks[[bk]]],
          list(lambda = f1@hyper[[par]][[bk]][["lambda"]])))
      }
      shares[[paste(par, bk)]] <- cols
    }
  }
  Fm <- solve(H + S, H)
  got <- vapply(shares, function(cols) sum(diag(Fm)[cols]), numeric(1))
  ref <- stats::setNames(f1@edf$edf, paste(f1@edf$parameter, f1@edf$term))
  ref <- ref[names(got)]
  if (anyNA(ref) || rel(got, ref) > 1e-8) fail("the effective degrees of freedom")
  rss <- sum((data$accel - mu0)^2)
  tmu <- sum(f0@edf$edf[f0@edf$parameter == "mu"])
  s0 <- unique(round(sg0, 12))
  if (length(s0) != 1L || rel(s0, sqrt(rss / (nrow(data) - tmu))) > 1e-5 ||
      !(s0 > sqrt(rss / nrow(data)))) fail("the constant scale")

  # --- certificate -------------------------------------------------------------
  for (f in list(f0, f1)) {
    ce <- statmod_certificate(f)
    if (!identical(ce$state, "converged") || !(ce$decrement < 1e-6) ||
        !(ce$mode_error < 1e-6) || length(ce$boundary) ||
        !identical(ce$curvature, "analytic")) {
      fail("the certificate")
    }
  }

  # --- simulate ----------------------------------------------------------------
  sim <- claims$sim
  if (length(sim$data) != 3L || nrow(sim$hyper) != 0L) fail("rstatmod()")
  for (k in 1:3) {
    if (rel(sim$theta[[k]]$sigma, sg1) > 1e-10 ||
        rel(sim$theta[[k]]$mu, predict(f1, what = "mu", data)) > 1e-10) {
      fail("the truth rstatmod() returns")
    }
  }
  if (!identical(dim(claims$sims), c(nrow(data), 3L))) fail("simulate()")

  invisible(TRUE)
}
