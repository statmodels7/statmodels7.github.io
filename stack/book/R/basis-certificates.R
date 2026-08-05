# Consistency gate for the basis chapter.
#
# Same contract as the other four: everything Chapter 5 prints is re-derived
# here from a route the chapter does not itself take, and the render stops if
# any of them disagrees. Nothing in this file is visible to the reader.
#
# The discipline that matters is the one link-formulas.R established: a gate
# that compares the package with itself proves nothing. Every displayed formula
# below is transcribed BY HAND from the chapter -- the Cox-de Boor recurrence,
# the B-spline derivative recurrence, the Fourier phase shift, the Legendre
# polynomials written out as explicit polynomials, the closed-form Gram
# diagonals, the Khatri-Rao column ordering and the polyadic contraction -- and
# compared against what the package computes. Where a claim is structural
# rather than a formula (orthonormality, the anchored integral, the three
# Demmler-Reinsch properties, separability of the Gram matrix), it is checked
# against an independent quadrature or an independent numerical derivative,
# never against the package's own route to the same number.

# --- shared helpers ---------------------------------------------------------

# A grid strictly inside the interval and away from the knots, so that neither
# the half-open convention of the recurrence nor a one-sided finite-difference
# stencil is what is being measured.
.inner <- function(lo, hi, n = 37) {
  seq(lo + 0.017 * (hi - lo), hi - 0.013 * (hi - lo), length.out = n)
}

.maxdiff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

# One numerical derivative, never nested, of a matrix-valued function.
.num_deriv_matrix <- function(f, x, h = 1e-5) {
  (f(x + h) - f(x - h)) / (2 * h)
}

# The full knot vector of a package B-spline basis, assembled from what the
# object records: the boundary knots repeated degree + 1 times each.
.full_knots <- function(b) {
  p <- b@basis_params
  c(rep(b@lower, p$degree + 1L), p$knots, rep(b@upper, p$degree + 1L))
}


# --- 1. Section 5.1: the expansion, the anchor, the stencil -----------------

# Equation (fd-weights), transcribed: sum_j w_j s_j^i must be 0 for i != d and
# d! for i = d. This puts the printed system under test rather than the solver.
.certify_fd_weights <- function() {
  out <- character()
  for (d in 1:4) {
    off <- numericals7::fd_offsets(d)
    for (nm in c("central", "forward", "backward")) {
      s <- off[[nm]]
      w <- numericals7::fd_weights(s, d)
      for (i in seq_along(s) - 1L) {
        want <- if (i == d) factorial(d) else 0
        got <- sum(w * s^i)
        if (abs(got - want) > 1e-6 * max(1, abs(want))) {
          out <- c(out, sprintf(
            "fd weights (order %d, %s): sum w s^%d is %s, the chapter requires %s",
            d, nm, i, format(got), format(want)
          ))
        }
      }
    }
  }
  out
}

.certify_anchor <- function() {
  out <- character()
  bases <- list(
    bspline = basis7::bspline_basis(lower = -1, upper = 2, dimension = 7),
    fourier = basis7::fourier_basis(lower = -1, upper = 3, dimension = 5),
    legendre = basis7::poly_basis(lower = 0, upper = 2, dimension = 5)
  )
  bases$orthonorm <- basis7::orthonorm_basis(bases$bspline)
  bases$constrained <- basis7::constrain_basis(
    bases$bspline, colSums(basis7::basis_eval(bases$bspline, .inner(-1, 2)))
  )

  for (nm in names(bases)) {
    b <- bases[[nm]]
    v <- basis7::basis_int(b, b@lower)
    # Equation (basis-anchor): exactly zero, not nearly.
    if (max(abs(v)) != 0) {
      out <- c(out, sprintf(
        "%s: the integral at the lower endpoint is %s, not exactly zero",
        nm, format(max(abs(v)))
      ))
    }
    # ...and it really is the integral: differentiating it returns the basis.
    x <- .inner(b@lower, b@upper, 11)
    num <- .num_deriv_matrix(function(z) basis7::basis_int(b, z), x)
    if (.maxdiff(num, basis7::basis_eval(b, x)) > 1e-6) {
      out <- c(out, sprintf("%s: the integral does not differentiate back", nm))
    }
  }
  out
}

# The Gaussian bumps of the chapter, whose numerical derivative is compared
# with a closed form written out here and known to no part of the package.
.certify_fallback_accuracy <- function() {
  out <- character()
  Bumps <- S7::new_class("Bumps", parent = basis7::basis, package = NULL)
  # A local name for the generic: S7::method<- mutates it in place, so the
  # method is registered globally all the same, and `pkg::generic` on the left
  # of a replacement is not something R can assign to.
  gen_eval <- basis7::basis_eval
  S7::method(gen_eval, Bumps) <- function(basis, x, ...) {
    p <- basis@basis_params
    out <- exp(-0.5 * outer(x, p$centers, "-")^2 / p$width^2)
    colnames(out) <- basis7::basis_colnames(basis)
    out
  }
  bumps <- Bumps(
    basis_name = "bumps", dimension = 5L, lower = 0, upper = 1,
    basis_params = list(centers = seq(0.1, 0.9, length.out = 5), width = 0.15)
  )

  flags <- basis7::basis_is_numerical(bumps)
  if (!all(flags)) {
    out <- c(out, "the bumps basis is not reported as fully numerical")
  }

  x <- .inner(0, 1, 15)
  z <- outer(x, bumps@basis_params$centers, "-") / bumps@basis_params$width
  ev <- basis7::basis_eval(bumps, x)
  exact <- list(
    -z / bumps@basis_params$width * ev,
    (z^2 - 1) / bumps@basis_params$width^2 * ev,
    (-z^3 + 3 * z) / bumps@basis_params$width^3 * ev
  )
  # Relative to the size of the quantity itself: a third derivative of a bump of
  # width 0.15 is of order 500, so an absolute tolerance would be a statement
  # about the width rather than about the stencil.
  tol <- c(1e-8, 1e-6, 1e-4)
  for (d in 1:3) {
    scale <- max(1, max(abs(exact[[d]])))
    err <- .maxdiff(basis7::basis_deriv(bumps, x, order = d), exact[[d]]) / scale
    if (!is.finite(err) || err > tol[d]) {
      out <- c(out, sprintf(
        "the numerical derivative of order %d is off by %s in relative terms",
        d, format(err)
      ))
    }
  }
  out
}

assert_expansions_ok <- function() {
  .stop_on(c(
    tryCatch(.certify_fd_weights(), error = function(e)
      paste("fd weights: error --", conditionMessage(e))),
    tryCatch(.certify_anchor(), error = function(e)
      paste("anchor: error --", conditionMessage(e))),
    tryCatch(.certify_fallback_accuracy(), error = function(e)
      paste("fallbacks: error --", conditionMessage(e)))
  ), "Section 5.1")
}


# --- 2. Section 5.2: the three families -------------------------------------

# Equations (bspline-zero) and (cox-de-boor), transcribed literally. Nothing
# here comes from splines2, which is what the package calls.
.cox_de_boor <- function(x, knots, degree) {
  nb <- length(knots) - 1L
  cur <- matrix(0, length(x), nb)
  for (i in seq_len(nb)) {
    cur[, i] <- as.numeric(x >= knots[i] & x < knots[i + 1L])
  }
  if (degree > 0L) {
    for (r in seq_len(degree)) {
      nb <- nb - 1L
      nxt <- matrix(0, length(x), nb)
      for (i in seq_len(nb)) {
        d1 <- knots[i + r] - knots[i]
        d2 <- knots[i + r + 1L] - knots[i + 1L]
        t1 <- if (d1 > 0) (x - knots[i]) / d1 * cur[, i] else 0
        t2 <- if (d2 > 0) (knots[i + r + 1L] - x) / d2 * cur[, i + 1L] else 0
        nxt[, i] <- t1 + t2
      }
      cur <- nxt
    }
  }
  cur
}

# Equation (bspline-derivative), applied d times. Its argument is the degree,
# and the number of functions of that degree follows from the knot vector.
.cox_deriv <- function(x, knots, degree, d) {
  if (d == 0L) return(.cox_de_boor(x, knots, degree))
  k <- length(knots) - degree - 1L
  if (degree == 0L) return(matrix(0, length(x), k))
  lower <- .cox_deriv(x, knots, degree - 1L, d - 1L)
  out <- matrix(0, length(x), k)
  for (i in seq_len(k)) {
    d1 <- knots[i + degree] - knots[i]
    d2 <- knots[i + degree + 1L] - knots[i + 1L]
    t1 <- if (d1 > 0) lower[, i] / d1 else 0
    t2 <- if (d2 > 0) lower[, i + 1L] / d2 else 0
    out[, i] <- degree * (t1 - t2)
  }
  out
}

.certify_bspline <- function() {
  out <- character()
  for (cfg in list(
    list(lo = 0, hi = 1, k = 6L, m = 3L),
    list(lo = -2, hi = 3, k = 9L, m = 2L),
    list(lo = 0, hi = 1, k = 5L, m = 4L)
  )) {
    b <- basis7::bspline_basis(cfg$lo, cfg$hi, cfg$k, cfg$m)
    kn <- .full_knots(b)
    x <- .inner(cfg$lo, cfg$hi)

    err <- .maxdiff(.cox_de_boor(x, kn, cfg$m), basis7::basis_eval(b, x))
    if (err > 1e-12) {
      out <- c(out, sprintf(
        "bspline(K=%d, m=%d): the recurrence printed in the chapter differs from the package by %s",
        cfg$k, cfg$m, format(err)
      ))
    }

    for (d in seq_len(cfg$m)) {
      err <- .maxdiff(
        .cox_deriv(x, kn, cfg$m, d), basis7::basis_deriv(b, x, order = d)
      )
      if (err > 1e-10) {
        out <- c(out, sprintf(
          "bspline(K=%d, m=%d): the derivative recurrence at order %d differs by %s",
          cfg$k, cfg$m, d, format(err)
        ))
      }
    }

    # Above the degree the derivative is the zero function, not a refusal.
    if (any(basis7::basis_deriv(b, x, order = cfg$m + 1L) != 0)) {
      out <- c(out, sprintf(
        "bspline(K=%d, m=%d): the derivative above the degree is not zero",
        cfg$k, cfg$m
      ))
    }

    # Partition of unity.
    if (max(abs(rowSums(basis7::basis_eval(b, x)) - 1)) > 1e-12) {
      out <- c(out, sprintf(
        "bspline(K=%d, m=%d): the functions do not sum to one", cfg$k, cfg$m
      ))
    }

    # Local support: B_i vanishes outside [t_i, t_{i+m+1}].
    ev <- basis7::basis_eval(b, x)
    for (i in seq_len(cfg$k)) {
      outside <- x < kn[i] | x > kn[i + cfg$m + 1L]
      if (any(outside) && max(abs(ev[outside, i])) > 1e-12) {
        out <- c(out, sprintf(
          "bspline(K=%d, m=%d): function %d is non-zero outside its support",
          cfg$k, cfg$m, i
        ))
      }
    }
  }
  out
}

# Equation (fourier-shift), transcribed. Order -1 gives the antiderivative, and
# the chapter says the value at the lower endpoint is subtracted from it.
.fourier_hand <- function(b, x, k) {
  p <- b@basis_params
  z <- 2 * pi * (x - b@lower) / p$omega
  cols <- lapply(seq_len(p$n_pairs), function(j) {
    s <- (2 * pi * j / p$omega)^k
    cbind(sin(j * z + k * pi / 2) * s, cos(j * z + k * pi / 2) * s)
  })
  const <- if (k == 0L) rep(1, length(x)) else rep(0, length(x))
  cbind(const, do.call(cbind, cols))
}

.certify_fourier <- function() {
  out <- character()
  for (cfg in list(
    list(lo = 0, hi = 1, k = 5L), list(lo = -1, hi = 3, k = 7L)
  )) {
    b <- basis7::fourier_basis(cfg$lo, cfg$hi, cfg$k)
    x <- .inner(cfg$lo, cfg$hi)

    for (d in 0:4) {
      err <- .maxdiff(
        .fourier_hand(b, x, d),
        if (d == 0L) basis7::basis_eval(b, x) else
          basis7::basis_deriv(b, x, order = d)
      )
      if (err > 1e-10) {
        out <- c(out, sprintf(
          "fourier(K=%d): the shift identity at order %d differs by %s",
          cfg$k, d, format(err)
        ))
      }
    }

    # The antiderivative, anchored as the chapter describes: order -1 minus its
    # value at the lower endpoint, with x - lower for the constant column.
    anti <- .fourier_hand(b, x, -1L)
    at_lo <- .fourier_hand(b, b@lower, -1L)
    hand <- sweep(anti, 2L, as.numeric(at_lo), "-")
    hand[, 1L] <- x - b@lower
    err <- .maxdiff(hand, basis7::basis_int(b, x))
    if (err > 1e-10) {
      out <- c(out, sprintf(
        "fourier(K=%d): the anchored antiderivative differs by %s",
        cfg$k, format(err)
      ))
    }
  }
  out
}

# The first five Legendre polynomials written out, so that the recurrence the
# package runs is compared against something that is not a recurrence.
.legendre_hand <- function(t, d) {
  p0 <- list(
    function(t) rep(1, length(t)),
    function(t) t,
    function(t) (3 * t^2 - 1) / 2,
    function(t) (5 * t^3 - 3 * t) / 2,
    function(t) (35 * t^4 - 30 * t^2 + 3) / 8
  )
  p1 <- list(
    function(t) rep(0, length(t)),
    function(t) rep(1, length(t)),
    function(t) 3 * t,
    function(t) (15 * t^2 - 3) / 2,
    function(t) (35 * t^3 - 15 * t) / 2
  )
  p2 <- list(
    function(t) rep(0, length(t)),
    function(t) rep(0, length(t)),
    function(t) rep(3, length(t)),
    function(t) 15 * t,
    function(t) (105 * t^2 - 15) / 2
  )
  fs <- list(p0, p1, p2)[[d + 1L]]
  do.call(cbind, lapply(fs, function(f) f(t)))
}

.certify_legendre <- function() {
  out <- character()
  for (cfg in list(list(lo = 0, hi = 1), list(lo = -2, hi = 4))) {
    b <- basis7::poly_basis(cfg$lo, cfg$hi, 5L)
    x <- .inner(cfg$lo, cfg$hi)
    half <- (cfg$hi - cfg$lo) / 2
    t <- (x - cfg$lo) / half - 1

    for (d in 0:2) {
      # The chain rule for the shift contributes (2/(u - l))^d.
      hand <- .legendre_hand(t, d) / half^d
      got <- if (d == 0L) basis7::basis_eval(b, x) else
        basis7::basis_deriv(b, x, order = d)
      err <- .maxdiff(hand, got)
      if (err > 1e-10) {
        out <- c(out, sprintf(
          "legendre on (%s, %s): the polynomials at order %d differ by %s",
          format(cfg$lo), format(cfg$hi), d, format(err)
        ))
      }
    }

    # Equation (legendre-integral), transcribed with P_5 written out.
    p <- cbind(.legendre_hand(t, 0L), (63 * t^5 - 70 * t^3 + 15 * t) / 8)
    hand <- matrix(0, length(x), 5L)
    hand[, 1L] <- x - cfg$lo
    for (n in 1:4) hand[, n + 1L] <- (p[, n + 2L] - p[, n]) / (2 * n + 1) * half
    err <- .maxdiff(hand, basis7::basis_int(b, x))
    if (err > 1e-10) {
      out <- c(out, sprintf(
        "legendre on (%s, %s): the integral identity differs by %s",
        format(cfg$lo), format(cfg$hi), format(err)
      ))
    }
  }
  out
}

assert_families_ok <- function() {
  .stop_on(c(
    tryCatch(.certify_bspline(), error = function(e)
      paste("bspline: error --", conditionMessage(e))),
    tryCatch(.certify_fourier(), error = function(e)
      paste("fourier: error --", conditionMessage(e))),
    tryCatch(.certify_legendre(), error = function(e)
      paste("legendre: error --", conditionMessage(e)))
  ), "Section 5.2")
}


# --- 3. Section 5.3: inner products -----------------------------------------

# Every entry of a Gram matrix by adaptive quadrature, which shares no code
# with the exact rules the package uses.
.gram_by_integrate <- function(b, order = 0L, weight = NULL) {
  k <- b@dimension
  g <- matrix(0, k, k)
  fn <- function(a, c) function(t) {
    m <- if (order == 0L) basis7::basis_eval(b, t) else
      basis7::basis_deriv(b, t, order = order)
    w <- if (is.null(weight)) 1 else weight(t)
    m[, a] * m[, c] * w
  }
  for (a in seq_len(k)) {
    for (c in a:k) {
      v <- stats::integrate(fn(a, c), b@lower, b@upper,
        rel.tol = 1e-11, subdivisions = 500L
      )$value
      g[a, c] <- v
      g[c, a] <- v
    }
  }
  g
}

.certify_gram_exact <- function() {
  out <- character()
  cases <- list(
    list(b = basis7::bspline_basis(0, 1, 6L, 3L), d = 0:2, tol = 1e-9),
    list(b = basis7::bspline_basis(-1, 2, 7L, 2L), d = 0:2, tol = 1e-8),
    list(b = basis7::fourier_basis(0, 1, 5L), d = 0:1, tol = 1e-8),
    list(b = basis7::poly_basis(0, 2, 5L), d = 0:2, tol = 1e-8)
  )
  for (cs in cases) {
    for (d in cs$d) {
      err <- .maxdiff(
        basis7::basis_gram(cs$b, order = d), .gram_by_integrate(cs$b, d)
      )
      if (err > cs$tol) {
        out <- c(out, sprintf(
          "%s: the order-%d Gram matrix differs from an independent quadrature by %s",
          cs$b@basis_name, d, format(err)
        ))
      }
    }
  }
  out
}

.certify_gram_forms <- function() {
  out <- character()

  # Equation (fourier-gram), transcribed by hand.
  for (cfg in list(list(lo = 0, hi = 1), list(lo = -1, hi = 3))) {
    b <- basis7::fourier_basis(cfg$lo, cfg$hi, 7L)
    omega <- cfg$hi - cfg$lo
    for (d in 0:2) {
      j <- seq_len(3L)
      diag_hand <- c(
        if (d == 0L) omega else 0,
        rep(omega / 2 * (2 * pi * j / omega)^(2 * d), each = 2L)
      )
      hand <- diag(diag_hand, nrow = 7L)
      err <- .maxdiff(hand, basis7::basis_gram(b, order = d))
      if (err > 1e-8 * max(1, max(abs(diag_hand)))) {
        out <- c(out, sprintf(
          "fourier on (%s, %s): the closed-form Gram diagonal at order %d differs by %s",
          format(cfg$lo), format(cfg$hi), d, format(err)
        ))
      }
    }
  }

  # The Legendre diagonal (u - l)/(2n + 1) at order zero.
  for (cfg in list(list(lo = 0, hi = 1), list(lo = -2, hi = 4))) {
    b <- basis7::poly_basis(cfg$lo, cfg$hi, 6L)
    hand <- diag((cfg$hi - cfg$lo) / (2 * (seq_len(6L) - 1L) + 1), nrow = 6L)
    err <- .maxdiff(hand, basis7::basis_gram(b))
    if (err > 1e-10) {
      out <- c(out, sprintf(
        "legendre on (%s, %s): the closed-form Gram diagonal differs by %s",
        format(cfg$lo), format(cfg$hi), format(err)
      ))
    }
  }

  # The band structure: zero whenever |a - b| > m.
  for (m in 1:3) {
    b <- basis7::bspline_basis(dimension = 9L, degree = m)
    g <- basis7::basis_gram(b, order = 0L)
    idx <- abs(outer(seq_len(9L), seq_len(9L), "-")) > m
    if (max(abs(g[idx])) > 1e-14) {
      out <- c(out, sprintf(
        "bspline(m=%d): the Gram matrix is not banded at width %d", m, m
      ))
    }
  }
  out
}

.certify_gram_measures <- function() {
  out <- character()
  b <- basis7::bspline_basis(dimension = 6L, degree = 3L)

  # Equation (gram-quadratic-form), against an independent integration of the
  # squared derivative of the fitted function.
  set.seed(11)
  beta <- stats::rnorm(6)
  for (d in 0:2) {
    qf <- drop(t(beta) %*% basis7::basis_gram(b, order = d) %*% beta)
    fn <- function(t) {
      m <- if (d == 0L) basis7::basis_eval(b, t) else
        basis7::basis_deriv(b, t, order = d)
      as.numeric(m %*% beta)^2
    }
    iv <- stats::integrate(fn, 0, 1, rel.tol = 1e-11,
      subdivisions = 500L)$value
    if (abs(qf - iv) > 1e-7 * max(1, abs(iv))) {
      out <- c(out, sprintf(
        "the quadratic form at order %d is %s against an integral of %s",
        d, format(qf), format(iv)
      ))
    }
  }

  # Equation (gram-empirical), transcribed: B'B / n.
  set.seed(3)
  xs <- stats::runif(500)
  for (d in 0:1) {
    m <- if (d == 0L) basis7::basis_eval(b, xs) else
      basis7::basis_deriv(b, xs, order = d)
    hand <- crossprod(m) / length(xs)
    if (.maxdiff(hand, basis7::basis_gram(b, order = d, at = xs)) > 1e-12) {
      out <- c(out, sprintf(
        "the empirical Gram matrix at order %d is not B'B/n", d
      ))
    }
  }

  # A weighted Lebesgue measure, against adaptive quadrature of the same
  # integrand with the weight inside it.
  w <- function(t) stats::dbeta(t, 2, 2)
  err <- .maxdiff(
    basis7::basis_gram(b, weight = w), .gram_by_integrate(b, 0L, weight = w)
  )
  if (err > 1e-8) {
    out <- c(out, sprintf("the weighted Gram matrix differs by %s", format(err)))
  }
  out
}

assert_gram_ok <- function() {
  .stop_on(c(
    tryCatch(.certify_gram_exact(), error = function(e)
      paste("gram exactness: error --", conditionMessage(e))),
    tryCatch(.certify_gram_forms(), error = function(e)
      paste("gram closed forms: error --", conditionMessage(e))),
    tryCatch(.certify_gram_measures(), error = function(e)
      paste("gram measures: error --", conditionMessage(e)))
  ), "Section 5.3")
}


# --- 4. Section 5.4: linear transformations ---------------------------------

.certify_transform_algebra <- function() {
  out <- character()
  b <- basis7::bspline_basis(0, 1, 6L, 3L)
  tr <- list(
    orthonorm = basis7::orthonorm_basis(b),
    constrained = basis7::constrain_basis(
      b, colSums(basis7::basis_eval(b, seq(0, 1, length.out = 200)))
    )
  )
  set.seed(5)
  tr$dr <- basis7::dr_basis(b, sort(stats::runif(200)))

  for (nm in names(tr)) {
    tb <- tr[[nm]]
    x <- .inner(0, 1, 13)

    # The transformed Gram matrix, by adaptive quadrature of the transformed
    # functions themselves rather than by the congruence the package applies.
    err <- .maxdiff(basis7::basis_gram(tb), .gram_by_integrate(tb, 0L))
    if (err > 1e-8) {
      out <- c(out, sprintf(
        "%s: the congruence disagrees with a direct integration by %s",
        nm, format(err)
      ))
    }

    # Derivatives and integrals transform by the same matrix, checked against
    # one numerical derivative of the transformed evaluation.
    num <- .num_deriv_matrix(function(z) basis7::basis_eval(tb, z), x)
    if (.maxdiff(num, basis7::basis_deriv(tb, x, order = 1L)) >
      1e-5 * max(1, max(abs(num)))) {
      out <- c(out, sprintf("%s: the derivative is not that of the basis", nm))
    }
    if (max(abs(basis7::basis_int(tb, tb@lower))) != 0) {
      out <- c(out, sprintf("%s: the anchor did not survive the transform", nm))
    }
  }

  # Equation (orthonorm): the Gram matrix of an orthonormalized basis is the
  # identity, verified by the independent quadrature above.
  err <- .maxdiff(.gram_by_integrate(tr$orthonorm, 0L), diag(6L))
  if (err > 1e-8) {
    out <- c(out, sprintf(
      "orthonormalization: the integrated inner products are off the identity by %s",
      format(err)
    ))
  }

  # Transforms compose by multiplication: the product is one object, and its
  # matrix is the product of the two.
  o2 <- basis7::orthonorm_basis(tr$constrained)
  if (!S7::S7_inherits(o2@parent_basis, basis7::BsplineBasis)) {
    out <- c(out, "a chain of transforms nested instead of composing")
  }
  if (.maxdiff(
    basis7::basis_eval(o2, .inner(0, 1, 9)),
    basis7::basis_eval(b, .inner(0, 1, 9)) %*% o2@transform
  ) > 1e-10) {
    out <- c(out, "the composed transform is not the product of the two")
  }
  out
}

.certify_constraint <- function() {
  out <- character()
  b <- basis7::bspline_basis(dimension = 8L, degree = 3L)
  grid <- seq(0, 1, length.out = 200)

  for (cm in list(
    matrix(colSums(basis7::basis_eval(b, grid)), nrow = 1L),
    rbind(colSums(basis7::basis_eval(b, grid)),
          basis7::basis_eval(b, 0.5))
  )) {
    r <- qr(cm)$rank
    cs <- basis7::constrain_basis(b, cm)
    if (cs@dimension != b@dimension - r) {
      out <- c(out, sprintf(
        "a constraint of rank %d reduced the dimension to %d, not %d",
        r, cs@dimension, b@dimension - r
      ))
    }
    # The defining property: C T = 0, so every admissible coefficient vector
    # satisfies the constraint exactly.
    if (max(abs(cm %*% cs@transform)) > 1e-10) {
      out <- c(out, "the constrained basis does not satisfy C T = 0")
    }
  }
  out
}

.certify_dr <- function() {
  out <- character()
  set.seed(9)
  cases <- list(
    list(b = basis7::bspline_basis(dimension = 10L), x = sort(stats::runif(300))),
    # A design that has lost rank, which is what the section claims survives.
    list(b = basis7::bspline_basis(dimension = 14L),
         x = sort(c(seq(0, 0.25, length.out = 120),
                    seq(0.78, 1, length.out = 120))))
  )
  for (cs in cases) {
    b <- cs$b
    x <- cs$x
    p <- basis7::basis_gram(b, order = 2L)
    d <- basis7::dr_basis(b, x)
    tm <- d@transform
    z <- basis7::basis_eval(d, x)

    # Equation (dr-properties), all three, computed here from Z and P.
    zz <- crossprod(z)
    off <- max(abs(zz[upper.tri(zz)])) / max(diag(zz))
    if (off > 1e-8) {
      out <- c(out, sprintf(
        "dr(K=%d): Z'Z is not diagonal (relative off-diagonal %s)",
        b@dimension, format(off)
      ))
    }
    tpt <- crossprod(tm, p %*% tm)
    scal <- mean(diag(tpt))
    if (.maxdiff(tpt, diag(scal, nrow(tpt))) > 1e-8 * max(1, abs(scal))) {
      out <- c(out, sprintf(
        "dr(K=%d): T'PT is not a multiple of the identity", b@dimension
      ))
    }
    if (max(abs(crossprod(cbind(1, x), z))) > 1e-8) {
      out <- c(out, sprintf(
        "dr(K=%d): the basis is not empirically orthogonal to 1 and x",
        b@dimension
      ))
    }

    # The consequence the section states: the penalized normal equations are
    # diagonal, so each coefficient is shrunk on its own.
    a <- zz + 0.37 * tpt
    if (max(abs(a[upper.tri(a)])) / max(diag(a)) > 1e-8) {
      out <- c(out, sprintf(
        "dr(K=%d): the penalized normal equations are not diagonal", b@dimension
      ))
    }
    if (b@dimension - d@dimension != 2L) {
      out <- c(out, sprintf(
        "dr(K=%d): the dimension fell by %d, not by the rank of the constraint",
        b@dimension, b@dimension - d@dimension
      ))
    }
  }

  # The rank deficiency of the second case is real, so the claim is not vacuous.
  if (qr(basis7::basis_eval(cases[[2]]$b, cases[[2]]$x))$rank >=
    cases[[2]]$b@dimension) {
    out <- c(out, "the rank-deficient example is not rank deficient")
  }
  out
}

assert_transforms_ok <- function() {
  .stop_on(c(
    tryCatch(.certify_transform_algebra(), error = function(e)
      paste("transform algebra: error --", conditionMessage(e))),
    tryCatch(.certify_constraint(), error = function(e)
      paste("constraint: error --", conditionMessage(e))),
    tryCatch(.certify_dr(), error = function(e)
      paste("demmler-reinsch: error --", conditionMessage(e)))
  ), "Section 5.4")
}


# --- 5. Section 5.5: products of bases --------------------------------------

# Equation (khatri-rao), transcribed literally: column (a - 1) K2 + c is the
# product of column a of the first and column c of the second. This is what
# fixes the ordering, and pairing a coefficient with the wrong function is a
# mistake that returns a perfectly finite number.
.hand_khatri <- function(parts) {
  Reduce(function(A, B) {
    pa <- ncol(A)
    pb <- ncol(B)
    out <- matrix(0, nrow(A), pa * pb)
    for (a in seq_len(pa)) {
      for (c in seq_len(pb)) {
        out[, (a - 1L) * pb + c] <- A[, a] * B[, c]
      }
    }
    out
  }, parts)
}

# The same ordering, applied to a coefficient array: the column index of
# (a_1, ..., a_D) built by the mixed-radix fold with the LAST index fastest.
.hand_flatten <- function(cf, dims) {
  idx <- as.matrix(expand.grid(lapply(dims, seq_len)))
  out <- numeric(prod(dims))
  for (r in seq_len(nrow(idx))) {
    a <- idx[r, ]
    col <- 0L
    for (j in seq_along(dims)) col <- col * dims[j] + (a[j] - 1L)
    out[col + 1L] <- cf[matrix(a, nrow = 1L)]
  }
  out
}

.certify_tensor_algebra <- function() {
  out <- character()
  m1 <- basis7::bspline_basis(0, 1, 4L, 3L)
  m2 <- basis7::fourier_basis(-1, 2, 3L)
  m3 <- basis7::poly_basis(0, 2, 3L)
  tb <- basis7::tensor_basis(m1, m2, m3)

  set.seed(21)
  n <- 12L
  x <- cbind(stats::runif(n, 0.05, 0.95), stats::runif(n, -0.9, 1.9),
    stats::runif(n, 0.05, 1.95))

  # Equation (khatri-rao).
  parts <- list(basis7::basis_eval(m1, x[, 1]), basis7::basis_eval(m2, x[, 2]),
    basis7::basis_eval(m3, x[, 3]))
  if (.maxdiff(.hand_khatri(parts), basis7::basis_eval(tb, x)) > 1e-12) {
    out <- c(out, "the tensor design matrix is not the printed Khatri-Rao product")
  }

  # Equation (tensor-derivative): each marginal to its own order.
  for (ord in list(c(1L, 0L, 0L), c(0L, 2L, 1L), c(1L, 1L, 1L))) {
    hand <- .hand_khatri(list(
      basis7::basis_deriv(m1, x[, 1], order = ord[1]),
      basis7::basis_deriv(m2, x[, 2], order = ord[2]),
      basis7::basis_deriv(m3, x[, 3], order = ord[3])
    ))
    if (.maxdiff(hand, basis7::basis_deriv(tb, x, order = ord)) > 1e-10) {
      out <- c(out, sprintf(
        "the partial derivative at order (%s) does not separate",
        paste(ord, collapse = ", ")
      ))
    }
  }

  # A single non-zero order is refused; zero is not.
  refused <- tryCatch({
    basis7::basis_deriv(tb, x, order = 1L)
    FALSE
  }, error = function(e) TRUE)
  if (!refused) out <- c(out, "a single non-zero order was accepted")
  if (.maxdiff(basis7::basis_deriv(tb, x, order = 0L),
    basis7::basis_eval(tb, x)) != 0) {
    out <- c(out, "order zero is not the evaluation")
  }

  # The anchored integral over the box, at the lower corner.
  lo <- matrix(c(m1@lower, m2@lower, m3@lower), nrow = 1L)
  if (max(abs(basis7::basis_int(tb, lo))) != 0) {
    out <- c(out, "the tensor integral is not zero at the lower corner")
  }
  out
}

# Equation (tensor-gram): checked against a two-dimensional adaptive
# integration of the actual product functions, which never forms a Kronecker
# product of anything.
.certify_tensor_gram <- function() {
  out <- character()
  m1 <- basis7::bspline_basis(0, 1, 4L, 3L)
  m2 <- basis7::fourier_basis(0, 1, 3L)
  tb <- basis7::tensor_basis(m1, m2)
  g <- basis7::basis_gram(tb)

  # The pairs are drawn from the entries the Kronecker structure predicts to be
  # non-zero, and from the entries it predicts to vanish, in both cases against
  # the same two-dimensional integral. Drawing at random would mostly pick the
  # zeros, since the Fourier marginal is orthogonal, and a reference scaled by a
  # constant agrees with zero however wrong it is.
  set.seed(31)
  big <- which(abs(g) > 1e-3, arr.ind = TRUE)
  small <- which(abs(g) < 1e-12, arr.ind = TRUE)
  pairs <- rbind(
    big[sample.int(nrow(big), 6L), , drop = FALSE],
    small[sample.int(nrow(small), 3L), , drop = FALSE]
  )
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1]
    c <- pairs[r, 2]
    inner <- function(s) {
      vapply(s, function(si) {
        stats::integrate(function(t) {
          m <- basis7::basis_eval(tb, cbind(rep(si, length(t)), t))
          m[, a] * m[, c]
        }, 0, 1, rel.tol = 1e-10)$value
      }, numeric(1))
    }
    v <- stats::integrate(inner, 0, 1, rel.tol = 1e-10)$value
    if (abs(v - g[a, c]) > 1e-7 * max(1, abs(v))) {
      out <- c(out, sprintf(
        "the tensor Gram entry (%d, %d) is %s against an integral of %s",
        a, c, format(g[a, c]), format(v)
      ))
    }
  }
  out
}

.certify_contraction <- function() {
  out <- character()
  m1 <- basis7::bspline_basis(0, 1, 5L, 3L)
  m2 <- basis7::bspline_basis(0, 1, 4L, 2L)
  m3 <- basis7::fourier_basis(0, 1, 3L)
  tb <- basis7::tensor_basis(m1, m2, m3)
  dims <- c(5L, 4L, 3L)

  set.seed(41)
  n <- 20L
  x <- matrix(stats::runif(n * 3, 0.02, 0.98), n, 3L)
  cf <- array(stats::rnorm(prod(dims)), dim = dims)

  # The array route, against the design matrix times the coefficients flattened
  # by the hand-written index map rather than by the package's aperm().
  hand <- drop(basis7::basis_eval(tb, x) %*% .hand_flatten(cf, dims))
  if (.maxdiff(basis7::basis_contract(tb, x, cf), hand) > 1e-10) {
    out <- c(out, "the array contraction pairs coefficients with the wrong functions")
  }
  # ...and with a block size below n, so the blocking is exercised.
  if (.maxdiff(basis7::basis_contract(tb, x, cf, block = 7L), hand) > 1e-10) {
    out <- c(out, "the blocked array contraction disagrees with the unblocked one")
  }

  # Equation (cp-contraction), transcribed as the sum over f of the product
  # over j of the marginal evaluations times a column of the factor matrix.
  fr <- 3L
  gam <- lapply(dims, function(k) matrix(stats::rnorm(k * fr), k, fr))
  ev <- list(basis7::basis_eval(m1, x[, 1]), basis7::basis_eval(m2, x[, 2]),
    basis7::basis_eval(m3, x[, 3]))
  hand_cp <- numeric(n)
  for (f in seq_len(fr)) {
    term <- rep(1, n)
    for (j in 1:3) term <- term * as.numeric(ev[[j]] %*% gam[[j]][, f])
    hand_cp <- hand_cp + term
  }
  if (.maxdiff(basis7::basis_contract(tb, x, gam), hand_cp) > 1e-10) {
    out <- c(out, "the polyadic contraction differs from the printed formula")
  }

  # The two shapes describe the same function: the array assembled from the
  # factor matrices must contract to the same values.
  full <- array(0, dim = dims)
  for (f in seq_len(fr)) {
    full <- full + outer(outer(gam[[1]][, f], gam[[2]][, f]), gam[[3]][, f])
  }
  if (.maxdiff(basis7::basis_contract(tb, x, full), hand_cp) > 1e-10) {
    out <- c(out, "the two coefficient shapes do not describe the same function")
  }

  # And the cost claim is at least structurally true: a product whose design
  # matrix cannot be formed is contracted all the same.
  big <- basis7::tensor_basis(rep(list(basis7::bspline_basis(dimension = 10L)), 6L))
  if (big@dimension != 1e6) {
    out <- c(out, "the six-variable product does not have a million functions")
  }
  gb <- replicate(6L, matrix(stats::rnorm(20L), 10L, 2L), simplify = FALSE)
  xb <- matrix(stats::runif(50 * 6, 0.02, 0.98), 50L, 6L)
  vb <- basis7::basis_contract(big, xb, gb)
  if (length(vb) != 50L || !all(is.finite(vb))) {
    out <- c(out, "the six-variable contraction did not produce finite values")
  }
  out
}

assert_tensor_ok <- function() {
  .stop_on(c(
    tryCatch(.certify_tensor_algebra(), error = function(e)
      paste("tensor algebra: error --", conditionMessage(e))),
    tryCatch(.certify_tensor_gram(), error = function(e)
      paste("tensor gram: error --", conditionMessage(e))),
    tryCatch(.certify_contraction(), error = function(e)
      paste("contraction: error --", conditionMessage(e)))
  ), "Section 5.5")
}


# --- the whole chapter ------------------------------------------------------

.stop_on <- function(problems, where) {
  if (length(problems)) {
    stop(where, " consistency gate failed:\n  ",
      paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_bases_ok <- function() {
  assert_expansions_ok()
  assert_families_ok()
  assert_gram_ok()
  assert_transforms_ok()
  assert_tensor_ok()
  invisible(TRUE)
}
