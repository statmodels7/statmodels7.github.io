# Consistency gate for the fitting chapter.
#
# Same contract as the other three: everything Chapter 5 asserts about a fit is
# re-derived here from quantities the chapter does not itself use, and the
# render stops if any of them disagrees. Nothing in this file is visible to the
# reader.
#
# Four claims are checked, one per section of the chapter:
#
#   1. the reported optimum really is one -- the link-scale score vanishes there;
#   2. the delta method holds -- V_theta = diag(h') V_eta diag(h'), with h'
#      recomputed from linkfunctions7 rather than taken from the fit;
#   3. every reported interval lies inside its parameter's domain, and its ends
#      are in increasing order even for a decreasing link;
#   4. the three optimisers agree about where the maximum is.

# The link-scale score of the whole sample at a given theta, summed over
# observations. Deliberately assembled here from distrib_gradient(scale =
# "link") rather than from the package's internal fit_score(), so that the two
# routes to the same quantity are independent.
.fit_score_at <- function(d, y, theta) {
  g <- distrib_gradient(d, y, theta, scale = "link")
  vapply(g, sum, numeric(1))
}

.fit_theta_at <- function(d, eta) {
  ps <- d@params
  stats::setNames(
    lapply(seq_along(ps), function(i) {
      linkfunctions7::linkinv(d@link_params[[ps[i]]], eta[i])
    }), ps)
}

# One fitted model, checked in every way the chapter claims.
.certify_one_fit <- function(d, theta_true, n = 400, seed = 1, tol = 1e-4) {
  set.seed(seed)
  y <- distrib_rng(d, n, theta_true)
  fit <- fit_distrib(d, y)
  ps <- d@params
  out <- character()

  if (!fit@converged) {
    return(sprintf("%s: fit did not converge", d@distrib_name))
  }

  # (1) the score vanishes at the reported optimum. Scaled by n, since the
  # summed score grows with the sample and an absolute threshold would mean
  # different things for different n.
  U <- .fit_score_at(d, y, .fit_theta_at(d, fit@eta))
  if (max(abs(U)) / n > tol) {
    out <- c(out, sprintf("%s: score at optimum is %.3e per observation",
                          d@distrib_name, max(abs(U)) / n))
  }

  # (2) the delta method, with h' taken straight from linkfunctions7
  J <- vapply(seq_along(ps), function(i) {
    linkfunctions7::dlinkinv(d@link_params[[ps[i]]], fit@eta[i])
  }, numeric(1))
  V_expect <- diag(J, length(ps)) %*% fit@vcov_eta %*% diag(J, length(ps))
  err <- max(abs(V_expect - fit@vcov) / pmax(abs(fit@vcov), 1e-8))
  if (is.finite(err) && err > 1e-8) {
    out <- c(out, sprintf("%s: delta method disagrees, relative error %.3e",
                          d@distrib_name, err))
  }

  # (3) intervals lie inside the domain and are ordered. Read through confint(),
  # which is what the chapter displays, rather than through the stored slot: a
  # check on a quantity the reader never sees would not cover the one shown.
  ci <- stats::confint(fit)
  ci_eta <- stats::confint(fit, scale = "link")

  for (i in seq_along(ps)) {
    b <- d@params_bounds[[ps[i]]]
    lo <- ci[i, 1]
    hi <- ci[i, 2]
    if (!is.finite(lo) || !is.finite(hi)) next
    if (lo <= b[1] || hi >= b[2]) {
      out <- c(out, sprintf("%s: interval for %s escapes (%s, %s)",
                            d@distrib_name, ps[i], format(b[1]), format(b[2])))
    }
    if (lo > hi) {
      out <- c(out, sprintf("%s: interval for %s is inverted", d@distrib_name, ps[i]))
    }

    # the link-scale interval is symmetric about the estimate, and the
    # parameter-scale one is its image under g^{-1}. Mapped here independently
    # of the package's own mapping, and sorted, since a link may decrease.
    if (!all(is.finite(ci_eta[i, ]))) next
    half <- c(fit@eta[i] - ci_eta[i, 1], ci_eta[i, 2] - fit@eta[i])
    if (abs(half[1] - half[2]) > 1e-10 * max(1, abs(fit@eta[i]))) {
      out <- c(out, sprintf("%s: link-scale interval for %s is not symmetric",
                            d@distrib_name, ps[i]))
    }
    mapped <- sort(linkfunctions7::linkinv(d@link_params[[ps[i]]], ci_eta[i, ]))
    if (max(abs(mapped - ci[i, ]) / pmax(abs(ci[i, ]), 1e-8)) > 1e-10) {
      out <- c(out, sprintf("%s: interval for %s is not the image of the link one",
                            d@distrib_name, ps[i]))
    }
  }

  # (4) the three optimisers find the same maximum. Compared on the
  # log-likelihood rather than on the estimates: a flat direction can move the
  # estimate without moving the objective, and it is the objective that the
  # methods are all maximising.
  ll <- c(fisher = fit@loglik)
  for (m in c("newton", "bfgs")) {
    alt <- tryCatch(fit_distrib(d, y, method = m), error = function(e) NULL)
    if (!is.null(alt) && alt@converged) ll[m] <- alt@loglik
  }
  if (length(ll) > 1 && diff(range(ll)) > 1e-4 * (abs(ll[[1]]) + 1)) {
    out <- c(out, sprintf("%s: methods disagree on the maximum (%s)",
                          d@distrib_name,
                          paste(sprintf("%s=%.6f", names(ll), ll), collapse = ", ")))
  }

  out
}

# The gate itself. A spread of families rather than one: a two-parameter
# continuous one, a one-parameter discrete one, a bounded one, and the Laplace,
# whose location is exactly the non-regular case the chapter appeals to.
assert_fit_ok <- function() {
  cases <- list(
    list(d = gaussian_distrib(), theta = list(mu = 1, sigma = 2)),
    list(d = gamma_distrib(),    theta = list(mu = 4, sigma2 = 6)),
    list(d = poisson_distrib(),  theta = list(mu = 3)),
    list(d = beta_distrib(),     theta = list(mu = 0.4, phi = 8)),
    list(d = laplace_distrib(),  theta = list(mu = 0, b = 1.5))
  )

  problems <- character()
  for (k in seq_along(cases)) {
    cse <- cases[[k]]
    res <- tryCatch(
      .certify_one_fit(cse$d, cse$theta, seed = 100 + k),
      error = function(e) sprintf("%s: error -- %s",
                                  cse$d@distrib_name, conditionMessage(e))
    )
    problems <- c(problems, res)
  }

  problems <- c(problems, .certify_starting_values())

  if (length(problems)) {
    stop("Chapter 5 consistency gate failed:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}


# The section on starting values makes two claims: that a family which knows
# an estimator returns it, and that starting there is what lets the fit
# converge. Both are checked against arithmetic done here, and the second is
# checked the only way it can be -- by starting somewhere else and watching
# the same fit fail.
.certify_starting_values <- function() {
  out <- character()
  d <- mvgaussian_distrib(4)
  y <- as.matrix(datasets::iris[, 1:4])

  st <- distrib_start(d, y)[[1]]
  if (max(abs(as.numeric(mv_location(d, st)) - colMeans(y))) > 1e-10) {
    out <- c(out, "the proposed start is not the sample mean")
  }
  s_hat <- crossprod(sweep(y, 2L, colMeans(y))) / nrow(y)
  if (max(abs(unname(mv_sigma(d, st)) - unname(s_hat))) > 1e-8) {
    out <- c(out, "the proposed start is not the sample second moment")
  }

  # The exact maximum, in closed form, and the fit that must reach it.
  ll_hat <- sum(distrib_pdf(d, y, st, log = TRUE))
  f <- fit_distrib(d, y)
  if (!f@converged) {
    out <- c(out, "the four-dimensional gaussian fit did not converge")
  }
  if (abs(as.numeric(logLik(f)) - ll_hat) > 1e-8) {
    out <- c(out, sprintf(
      "the fit reached %s where the closed-form maximum is %s",
      format(as.numeric(logLik(f)), digits = 10), format(ll_hat, digits = 10)
    ))
  }
  if (f@iterations > 3) {
    out <- c(out, sprintf(
      "the fit took %d iterations from its own start, not the one it claims",
      f@iterations
    ))
  }

  # And the counterfactual: the same everything, started at the origin, does
  # not get there. Without this the claim would be untested -- a fit that
  # converges quickly from a good start proves nothing about the start unless
  # a bad one is shown to fail.
  zero <- as.list(stats::setNames(rep(0, d@n_params), d@params))
  fz <- fit_distrib(d, y, start = zero, maxit = 500)
  if (isTRUE(fz@converged) &&
      abs(as.numeric(logLik(fz)) - ll_hat) < 1e-8) {
    out <- c(out, paste0(
      "the fit started at the origin now reaches the maximum too, so the ",
      "section's claim about starting values no longer has a counterexample"
    ))
  }
  out
}
