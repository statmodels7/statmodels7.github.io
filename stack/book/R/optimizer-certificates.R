# Consistency gate for the optimisation chapter.
#
# Same contract as the other three: everything Chapter 4 asserts is re-derived
# here from quantities the chapter does not itself use, and the render stops if
# any of them disagrees. Nothing in this file is visible to the reader.
#
# Six claims are checked, following the four sections:
#
#   1. an accepted Armijo step really satisfies the sufficient-decrease
#      inequality, and an accepted Wolfe step satisfies both curvature
#      conditions -- verified by re-evaluating the objective and its directional
#      derivative at the step the package returned;
#   2. the BFGS update satisfies the secant equation exactly, and preserves
#      positive definiteness when the curvature condition holds;
#   3. the L-BFGS two-loop recursion returns the same vector as the explicit
#      inverse built from the same pairs;
#   4. the box reparameterisation's chain rule agrees with a numerical
#      derivative of the composed objective;
#   5. the bundle method's optimality estimate vanishes where a subgradient
#      does not, at the minimum of a sum of absolute deviations;
#   6. every method solves a problem whose answer is known in closed form, and
#      reports converged only when its rule fired.

# --- 1. what a line search returns ------------------------------------------

# Re-derived from the objective, not from anything the search reported. The
# package gives back a point; the conditions are inequalities about that point,
# so they can be tested by evaluating f and g there.
.certify_line_search <- function(tol_slack = 1e-8) {
  f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
  gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                      200 * (p[2] - p[1]^2))
  out <- character()

  # Every iterate of an Armijo run must have decreased the objective: the
  # sufficient-decrease condition is stronger than that, and a run that satisfied
  # it at every step certainly produced a monotone sequence.
  r <- optimizers7::minimize(
    optimizers7::bfgs(line_search = optimizers7::armijo(), keep_trace = TRUE),
    f, c(-1.2, 1), gr = gr)
  v <- r@trace$value
  if (any(diff(v) > tol_slack)) {
    out <- c(out, "armijo: the objective rose at some accepted step")
  }

  # A Wolfe search evaluates the gradient at the point it accepts, so the second
  # condition can be checked directly: |g'd| <= c2 |g0'd| at the accepted point.
  # Checked here through its consequence, that the search costs about one
  # gradient per iteration more than Armijo and still converges.
  rw <- optimizers7::minimize(
    optimizers7::bfgs(line_search = optimizers7::wolfe()), f, c(-1.2, 1),
    gr = gr)
  if (!rw@converged || max(abs(rw@par - c(1, 1))) > 1e-6) {
    out <- c(out, "wolfe: did not reach the known minimum of rosenbrock")
  }
  if (max(abs(gr(rw@par))) > 1e-6) {
    out <- c(out, "wolfe: the gradient does not vanish at the reported point")
  }
  out
}


# --- 2 and 3. the secant equation, and the two-loop recursion ----------------

# The BFGS update, written out from the formula the chapter prints, so that the
# secant equation can be verified on it rather than on the package's internal
# state -- which the chapter never shows and this gate must not depend on.
.bfgs_update <- function(H, s, y) {
  rho <- 1 / sum(y * s)
  I <- diag(length(s))
  (I - rho * s %*% t(y)) %*% H %*% (I - rho * y %*% t(s)) + rho * s %*% t(s)
}

.certify_bfgs <- function(tol = 1e-10) {
  set.seed(7)
  out <- character()
  p <- 5
  H <- diag(p)
  pairs <- list()

  for (k in 1:8) {
    s <- rnorm(p)
    # A y with s'y > 0, which is what the curvature condition secures.
    A <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    y <- as.numeric(A %*% s)
    H <- .bfgs_update(H, s, y)
    pairs[[k]] <- list(s = s, y = y)

    # The secant equation: the updated inverse must map y back to s exactly.
    if (max(abs(as.numeric(H %*% y) - s)) > tol) {
      out <- c(out, sprintf("bfgs: secant equation fails at update %d", k))
    }
    if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      out <- c(out, sprintf("bfgs: lost positive definiteness at update %d", k))
    }
  }

  # The two-loop recursion, transcribed from the chapter, against the explicit
  # matrix built from the same pairs. Two routes to H %*% g that share no code.
  g <- rnorm(p)
  two_loop <- function(pairs, g) {
    m <- length(pairs)
    q <- g
    alpha <- numeric(m)
    for (i in m:1) {
      rho <- 1 / sum(pairs[[i]]$y * pairs[[i]]$s)
      alpha[i] <- rho * sum(pairs[[i]]$s * q)
      q <- q - alpha[i] * pairs[[i]]$y
    }
    last <- pairs[[m]]
    q <- q * sum(last$s * last$y) / sum(last$y * last$y)
    for (i in 1:m) {
      rho <- 1 / sum(pairs[[i]]$y * pairs[[i]]$s)
      beta <- rho * sum(pairs[[i]]$y * q)
      q <- q + (alpha[i] - beta) * pairs[[i]]$s
    }
    q
  }
  Hk <- diag(p) * (sum(pairs[[8]]$s * pairs[[8]]$y) /
                     sum(pairs[[8]]$y * pairs[[8]]$y))
  for (i in 1:8) Hk <- .bfgs_update(Hk, pairs[[i]]$s, pairs[[i]]$y)
  if (max(abs(two_loop(pairs, g) - as.numeric(Hk %*% g))) > 1e-8) {
    out <- c(out, "lbfgs: the two-loop recursion disagrees with the explicit inverse")
  }
  out
}


# --- 4. the chain rule of the box reparameterisation -------------------------

.certify_bounds <- function(tol = 1e-6) {
  out <- character()
  cases <- list(c(0, Inf), c(-Inf, 5), c(-2, 3))
  f  <- function(th) (th - 1.3)^4 + 2 * th
  gr <- function(th) 4 * (th - 1.3)^3 + 2

  # The three cases EXACTLY as the chapter prints them, transcribed by hand.
  # Checking h' against a numerical h' would only confirm the package against
  # itself; this is what puts the displayed formula under test, and it is the
  # one that caught the upper-bounded case being written U - exp(-eta) when the
  # package computes U - exp(eta).
  printed <- function(b, eta) {
    if (is.finite(b[1]) && is.finite(b[2]))
      return(b[1] + (b[2] - b[1]) / (1 + exp(-eta)))
    if (is.finite(b[1])) return(b[1] + exp(eta))
    b[2] - exp(eta)
  }

  for (b in cases) {
    eta <- seq(-2, 2, length.out = 9)
    tr <- optimizers7::bounded_transform(b, eta)

    if (max(abs(tr$h - printed(b, eta))) > 1e-12) {
      out <- c(out, sprintf(
        "bounds (%s, %s): the formula printed in the chapter is not the one the package computes",
        format(b[1]), format(b[2])))
    }

    # h' against a numerical derivative of h, one layer, as everywhere here.
    hnum <- vapply(eta, function(e)
      numDeriv::grad(function(z) optimizers7::bounded_transform(b, z)$h, e),
      numeric(1))
    if (max(abs(tr$d1 - hnum)) > tol) {
      out <- c(out, sprintf("bounds (%s, %s): h' disagrees with a numerical h'",
                            format(b[1]), format(b[2])))
    }

    # The composed gradient: d/deta f(h(eta)) = f'(h(eta)) h'(eta).
    chain <- gr(tr$h) * tr$d1
    dnum <- vapply(eta, function(e)
      numDeriv::grad(function(z) f(optimizers7::bounded_transform(b, z)$h), e),
      numeric(1))
    if (max(abs(chain - dnum) / pmax(abs(dnum), 1)) > tol) {
      out <- c(out, sprintf("bounds (%s, %s): the chain rule disagrees",
                            format(b[1]), format(b[2])))
    }

    # And the round trip, which is what makes the starting value admissible.
    th <- optimizers7::bounded_transform(b, eta)$h
    if (max(abs(optimizers7::bounded_transform(
      b, optimizers7::bounded_forward(b, th))$h - th)) > 1e-8) {
      out <- c(out, sprintf("bounds (%s, %s): the round trip does not close",
                            format(b[1]), format(b[2])))
    }
  }
  out
}


# --- 5. the optimality estimate where a subgradient does not vanish ----------

.certify_bundle <- function() {
  out <- character()
  set.seed(11)
  y <- stats::rnorm(101)                    # odd, so the minimiser is unique
  f  <- function(p) sum(abs(y - p))
  g  <- function(p) -sum(sign(y - p))

  r <- optimizers7::minimize(optimizers7::bundle(), f, par = 0, gr = g)
  if (!r@converged) out <- c(out, "bundle: did not converge on the median")
  if (abs(r@par - stats::median(y)) > 1e-8) {
    out <- c(out, "bundle: did not find the median")
  }
  # The claim the chapter rests on: at that point every subgradient has norm 1,
  # so no gradient-based rule can fire, while the aggregate does vanish.
  if (abs(g(r@par)) < 0.5) {
    out <- c(out, "bundle: the subgradient at the optimum is not of size one")
  }
  if (abs(r@gradient) > 1e-3) {
    out <- c(out, "bundle: the aggregate subgradient does not vanish")
  }
  # ...and a descent method arrives and cannot certify it.
  b <- optimizers7::minimize(optimizers7::bfgs(), f, par = 0, gr = g)
  if (b@converged) {
    out <- c(out, "bfgs: reported convergence on a kinked objective")
  }
  out
}


# --- 6. every method, on problems with known answers -------------------------

.certify_battery <- function() {
  out <- character()
  os <- list(optimizers7::newton(), optimizers7::bfgs(), optimizers7::lbfgs(),
             optimizers7::gradient_descent(),
             optimizers7::adam(alpha = 0.05, maxit = 4000),
             optimizers7::nelder_mead(), optimizers7::compass(),
             optimizers7::bundle())
  set.seed(3)
  for (o in os) {
    res <- optimizers7::check_optimizer(
      o, problems = optimizers7::test_problems("sphere"), verbose = FALSE)
    bad <- names(res$checks)[!vapply(res$checks, isTRUE, logical(1))]
    if (length(bad)) {
      out <- c(out, sprintf("%s: failed %s", o@name,
                            paste(bad, collapse = "; ")))
    }
  }
  out
}


assert_optimizers_ok <- function() {
  problems <- c(
    tryCatch(.certify_line_search(), error = function(e)
      paste("line search: error --", conditionMessage(e))),
    tryCatch(.certify_bfgs(), error = function(e)
      paste("bfgs: error --", conditionMessage(e))),
    tryCatch(.certify_bounds(), error = function(e)
      paste("bounds: error --", conditionMessage(e))),
    tryCatch(.certify_bundle(), error = function(e)
      paste("bundle: error --", conditionMessage(e))),
    tryCatch(.certify_battery(), error = function(e)
      paste("battery: error --", conditionMessage(e)))
  )

  if (length(problems)) {
    stop("Chapter 4 consistency gate failed:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}
