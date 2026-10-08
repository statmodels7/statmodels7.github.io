# Certificates for Part V, chapter 16, sections 16.4-16.5 (writing an optimizer).
#
# Same contract as every gate of the book: what the chapter states is re-derived
# by a route the reader's code does not take, and the render stops on any
# disagreement. Nothing in this file is visible to the reader.
#
# The claims:
#
#   subproblem  the step tr_step() returns satisfies the first-order conditions
#               of the trust-region subproblem, written here as matrix
#               identities with no eigendecomposition; its model value is not
#               beaten by a random search over the ball (the global minimum);
#   battery     check_optimizer() passes every check of the trust-region method;
#   gamma       fit_distrib() with the trust region reaches the logLik of the
#               default method;
#   inner       every inner run of the statmod() fit meets the rule, and the
#               logLik equals the default inner method's;
#   outer       the three outer methods reach the same logLik and hyperparameters;
#               the trust region takes fewer iterations than the unlimited
#               Newton method; every outer run converged;
#   radius      the iteration counts for the three starting radii are within a
#               factor of two of one another.

assert_own_optim_ok <- function(fits, claims, data, ref) {

  fail <- function(what) {
    stop("Chapter 'Writing a model term' (sections 16.4-16.5): ", what,
         " no longer agrees with the package.", call. = FALSE)
  }
  ll <- function(f) as.numeric(stats::logLik(f))

  # --- subproblem ----------------------------------------------------------------
  set.seed(16)
  n_bad <- 0L
  worst_model_gap <- -Inf
  for (i in 1:60) {
    M <- matrix(stats::rnorm(4), 2)
    H <- (M + t(M)) * 1.5
    g <- stats::rnorm(2)
    D <- stats::runif(1, 0.05, 3)
    p <- tr_step(g, H, D)
    mu <- attr(p, "mu")
    res_stat <- max(abs(drop((H + mu * diag(2)) %*% p) + g))
    res_comp <- abs(mu * (D - sqrt(sum(p^2))))
    res_feas <- max(0, sqrt(sum(p^2)) - D)
    psd <- !inherits(try(chol(H + (mu + 1e-10) * diag(2)), silent = TRUE), "try-error")
    if (res_stat > 1e-7 || res_comp > 1e-7 || res_feas > 1e-9 || mu < -1e-12 || !psd) {
      n_bad <- n_bad + 1L
    }
    m_tr <- sum(g * p) + 0.5 * sum(p * drop(H %*% p))
    # random directions inside the ball: none may lower the model
    for (j in 1:400) {
      q <- stats::rnorm(2)
      q <- q / sqrt(sum(q^2)) * D * stats::runif(1)^0.5
      m_q <- sum(g * q) + 0.5 * sum(q * drop(H %*% q))
      worst_model_gap <- max(worst_model_gap, m_tr - m_q)
    }
  }
  if (n_bad > 0L) fail("the first-order conditions of the trust-region step")
  if (worst_model_gap > 1e-9) fail("the trust-region step against a random search on the ball")

  # --- battery -------------------------------------------------------------------
  if (!all(unlist(claims$chk_tr$checks))) fail("check_optimizer() on the trust-region method")

  # --- gamma ---------------------------------------------------------------------
  if (abs(fits$gamma_tr@loglik - fits$gamma_default@loglik) > 1e-6 ||
      abs(fits$gamma_tr@loglik - (-1065.306)) > 1e-3) {
    fail("fit_distrib() gamma with the trust region")
  }

  # --- inner ---------------------------------------------------------------------
  if (length(claims$inner_log$runs) == 0L) fail("the inner runs were not kept")
  not_conv <- sum(vapply(claims$inner_log$runs, function(r) !r$converged, TRUE))
  if (not_conv > 0L) fail("the inner runs of the trust region")
  f_default_inner <- statmod(accel ~ s(times) | sigma ~ s(times),
                             distrib = gaussian1_distrib(), data = data$mcycle)
  if (abs(ll(fits$inner) - ll(f_default_inner)) > 1e-5) {
    fail("the inner fit against the default inner method")
  }

  # --- outer ---------------------------------------------------------------------
  ot <- claims$outer_table
  if (max(abs(ot$log_likelihood - ot$log_likelihood[1])) > 1e-5) {
    fail("the logLik of the three outer methods")
  }
  hy <- lapply(fits$outer, function(f) hyper(f)$estimate)
  if (max(abs(unlist(hy[-1]) - unlist(hy[[1]]))) > 1e-4) {
    fail("the hyperparameters of the three outer methods")
  }
  for (f in fits$outer) {
    if (statmod_certificate(f)$state != "converged") fail("an outer run")
  }
  if (!(ot$iterations[1] < ot$iterations[2])) {
    fail("the iteration count of the trust region against the unlimited Newton method")
  }

  # --- radius --------------------------------------------------------------------
  r0 <- claims$r0_counts
  if (any(!is.finite(r0)) || max(r0) > 2 * min(r0)) {
    fail("the iteration counts for the three starting radii")
  }
  invisible(TRUE)
}
