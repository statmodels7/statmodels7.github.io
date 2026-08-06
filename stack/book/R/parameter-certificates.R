# Consistency gate for chapter 6, the parameters7 package.
#
# Same contract as the other gates. Every formula the chapter prints is
# transcribed BY HAND here and compared against what the package computes;
# every structural claim is checked against a route the package does not take.
# Nothing in this file is visible to the reader.
#
# The transcriptions are the log-Cholesky map and its two derivative formulas,
# the linearity of its log-determinant, the diagonal family's derivatives, the
# scaled family's whole set, the pseudo-determinant as a sum of logs of
# non-zero eigenvalues, the stationary point lambda = r / (b'Pb), and the
# identity d log|M| = tr(M^{-1} dM). Where a claim is about behavior rather
# than a formula -- that a rank-deficient solve is refused, that the eigenvalue
# count is not scale invariant while the stacked rank is -- it is exercised
# directly.

.cs_maxdiff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

# One numerical derivative, never nested.
.cs_grad_matrix <- function(f, eta, k, h = 1e-5) {
  up <- dn <- eta
  up[k] <- eta[k] + h
  dn[k] <- eta[k] - h
  (f(up) - f(dn)) / (2 * h)
}

.cs_eta3 <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)

# The lower-triangular factor, assembled here from the printed description of
# the free vector: the logs of the diagonal first, then the below-diagonal
# entries column by column.
.cs_factor <- function(eta, p) {
  l <- matrix(0, p, p)
  diag(l) <- exp(eta[seq_len(p)])
  k <- p
  for (j in seq_len(p - 1L)) {
    for (i in seq(j + 1L, p)) {
      k <- k + 1L
      l[i, j] <- eta[k]
    }
  }
  l
}

# The single-entry derivative of the factor, equation (logchol-derivs).
.cs_dfactor <- function(eta, p, k) {
  lk <- matrix(0, p, p)
  if (k <= p) {
    lk[k, k] <- exp(eta[k])
  } else {
    idx <- p
    for (j in seq_len(p - 1L)) {
      for (i in seq(j + 1L, p)) {
        idx <- idx + 1L
        if (idx == k) lk[i, j] <- 1
      }
    }
  }
  lk
}


# --- 1. The log-Cholesky map ------------------------------------------------

.certify_logchol_map <- function() {
  unlist(lapply(list(
    list(p = 3L, eta = .cs_eta3),
    # p = 4 as well, and not for coverage. At p = 3 the below-diagonal entries
    # are (2,1), (3,1), (3,2) whether they are read down the columns or across
    # the rows, so the ordering @eq-logchol-free states is not under test
    # there. At p = 4 the two orderings differ.
    list(p = 4L, eta = c(0.2, -0.1, 0.4, -0.3,
                         0.5, -0.4, 0.2, 0.3, -0.6, 0.1))
  ), .certify_logchol_map_at))
}

.certify_logchol_map_at <- function(case) {
  out <- character()
  p <- case$p
  s <- parameters7::log_cholesky(p)
  eta <- case$eta
  l <- .cs_factor(eta, p)

  # M = L L', with L transcribed above rather than taken from the package
  err <- .cs_maxdiff(parameters7::param_value(s, eta), tcrossprod(l))
  if (err > 1e-12) {
    out <- c(out, sprintf(
      "the printed log-Cholesky map differs from the package by %s",
      format(err, digits = 3)
    ))
  }

  # the free vector really is the one the chapter describes: the logs of the
  # diagonal first, then the below-diagonal entries column by column
  want_names <- c(
    paste0("log_L", seq_len(p)),
    unlist(lapply(seq_len(p - 1L), function(j)
      paste0("L", seq(j + 1L, p), ".", j)))
  )
  if (!identical(s@free_names, want_names)) {
    out <- c(out, sprintf(
      "the free names at p = %d are not the ordering the chapter prints", p
    ))
  }

  # equation (logchol-derivs), first order
  for (k in seq_len(s@n_free)) {
    lk <- .cs_dfactor(eta, p, k)
    want <- tcrossprod(lk, l) + t(tcrossprod(lk, l))
    got <- parameters7::param_d1(s, eta)[[k]]
    if (.cs_maxdiff(want, got) > 1e-12) {
      out <- c(out, sprintf(
        "the printed first derivative in direction %d disagrees with the package", k
      ))
    }
    # and against a numerical derivative of the map, which shares no code
    num <- .cs_grad_matrix(function(e) parameters7::param_value(s, e), eta, k)
    if (.cs_maxdiff(want, num) > 1e-6) {
      out <- c(out, sprintf(
        "the printed first derivative in direction %d disagrees with finite differences", k
      ))
    }
  }

  # equation (logchol-derivs), second order
  idx <- parameters7::param_tuple_indices(s)
  d2 <- parameters7::param_d2(s, eta)
  for (i in seq_along(idx)) {
    k <- idx[[i]][1L]
    m <- idx[[i]][2L]
    lk <- .cs_dfactor(eta, p, k)
    lm <- .cs_dfactor(eta, p, m)
    want <- tcrossprod(lk, lm) + t(tcrossprod(lk, lm))
    if (k == m && k <= p) {
      b <- tcrossprod(lk, l)
      want <- want + b + t(b)
    }
    if (.cs_maxdiff(want, d2[[i]]) > 1e-12) {
      out <- c(out, sprintf(
        "the printed second derivative (%d, %d) disagrees with the package", k, m
      ))
    }
  }

  # equation (logchol-logdet): linear in the free vector, gradient 2 and 0
  for (e in list(eta, rev(eta) * 1.7, rep(0, s@n_free))) {
    want_ld <- 2 * sum(e[seq_len(p)])
    if (abs(parameters7::param_logdet(s, e) - want_ld) > 1e-10) {
      out <- c(out, "the log-Cholesky log-determinant is not 2 sum log L_ii")
    }
    want_g <- c(rep(2, p), rep(0, s@n_free - p))
    if (.cs_maxdiff(parameters7::param_dlogdet(s, e), want_g) > 1e-12) {
      out <- c(out, "the log-determinant gradient is not the constant the chapter prints")
    }
    if (max(abs(parameters7::param_d2logdet(s, e))) > 1e-12) {
      out <- c(out, "the log-determinant Hessian does not vanish")
    }
  }

  # the inverse map closes exactly
  if (.cs_maxdiff(parameters7::param_free(s, parameters7::param_value(s, eta)),
                  eta) > 1e-12) {
    out <- c(out, "the log-Cholesky round trip does not close")
  }
  out
}


# --- 2. Equation (parameter-dlogdet), the trace identity --------------------

.certify_dlogdet_identity <- function() {
  out <- character()
  cases <- list(
    list(s = parameters7::log_cholesky(3), eta = .cs_eta3),
    list(s = parameters7::diagonal_matrix(3), eta = c(0.2, -0.4, 0.1)),
    list(s = parameters7::diagonal_matrix(2, link = linkfunctions7::sqrt_link()),
         eta = c(1.3, 2.1))
  )
  for (cs in cases) {
    s <- cs$s
    eta <- cs$eta
    m <- parameters7::param_value(s, eta)
    a <- parameters7::param_d1(s, eta)
    want <- vapply(a, function(ak) sum(diag(solve(m, ak))), numeric(1))
    got <- parameters7::param_dlogdet(s, eta)
    if (.cs_maxdiff(want, got) > 1e-8) {
      out <- c(out, sprintf(
        "tr(M^{-1} dM) and the log-determinant gradient of '%s' differ by %s",
        s@param_name, format(.cs_maxdiff(want, got), digits = 3)
      ))
    }
  }
  out
}


# --- 3. The diagonal family, equation (diag-struct) -------------------------

.certify_diag_family <- function() {
  out <- character()
  for (lnk in list(linkfunctions7::log_link(), linkfunctions7::sqrt_link())) {
    s <- parameters7::diagonal_matrix(3, link = lnk)
    eta <- c(0.3, -0.7, 1.1)
    want <- diag(linkfunctions7::linkinv(lnk, eta))
    if (.cs_maxdiff(parameters7::param_value(s, eta), want) > 1e-12) {
      out <- c(out, sprintf(
        "the diagonal map with the %s link is not diag(h(eta))", lnk@link_name
      ))
    }
    a <- parameters7::param_d1(s, eta)
    for (k in 1:3) {
      wk <- matrix(0, 3, 3)
      wk[k, k] <- linkfunctions7::dlinkinv(lnk, eta[k])
      if (.cs_maxdiff(a[[k]], wk) > 1e-10) {
        out <- c(out, sprintf(
          "the diagonal first derivative in direction %d is not h'(eta_k) e_k e_k'", k
        ))
      }
    }
    idx <- parameters7::param_tuple_indices(s)
    d2 <- parameters7::param_d2(s, eta)
    for (i in seq_along(idx)) {
      k <- idx[[i]][1L]
      l <- idx[[i]][2L]
      wk <- matrix(0, 3, 3)
      if (k == l) wk[k, k] <- linkfunctions7::d2linkinv(lnk, eta[k])
      if (.cs_maxdiff(d2[[i]], wk) > 1e-10) {
        out <- c(out, sprintf(
          "the diagonal second derivative (%d, %d) is not what the chapter prints", k, l
        ))
      }
    }
    if (abs(parameters7::param_logdet(s, eta) -
            sum(log(linkfunctions7::linkinv(lnk, eta)))) > 1e-10) {
      out <- c(out, "the diagonal log-determinant is not the sum of logs")
    }
  }
  out
}


# --- 4. The scaled family, equations (scaled-struct) and (scaled-logdet) ----

.certify_scaled_family <- function() {
  out <- character()
  dm <- diff(diag(6), differences = 2)
  p_mat <- crossprod(dm)
  s <- parameters7::scaled_matrix(p_mat)

  # the pseudo-determinant, written out as the sum of the logs of the non-zero
  # eigenvalues of the assembled matrix
  ev <- eigen(p_mat, symmetric = TRUE, only.values = TRUE)$values
  keep <- ev > 1e-10 * max(ev)
  if (sum(keep) != s@rank) {
    out <- c(out, "the declared rank and the eigenvalue count of P disagree")
  }
  for (eta in c(-3, 0, 2.5)) {
    lambda <- exp(eta)
    want_ld <- sum(log(lambda * ev[keep]))
    if (abs(parameters7::param_logdet(s, eta) - want_ld) > 1e-8) {
      out <- c(out, sprintf(
        "the scaled pseudo-determinant at eta = %s is not sum log(lambda nu_j)",
        format(eta)
      ))
    }
    # equation (scaled-logdet): the derivative is the rank, at any scale
    if (abs(parameters7::param_dlogdet(s, eta) - s@rank) > 1e-10) {
      out <- c(out, "the derivative of the scaled pseudo-determinant is not the rank")
    }
    if (abs(parameters7::param_d2logdet(s, eta)) > 1e-12) {
      out <- c(out, "the second derivative of the scaled pseudo-determinant is not zero")
    }
    # equation (scaled-struct): dM = M and d2M = M under the log link
    m <- parameters7::param_value(s, eta)
    if (.cs_maxdiff(parameters7::param_d1(s, eta)[[1L]], m) > 1e-10 ||
        .cs_maxdiff(parameters7::param_d2(s, eta)[[1L]], m) > 1e-10) {
      out <- c(out, "the scaled derivatives are not the matrix itself")
    }
    if (.cs_maxdiff(m, lambda * p_mat) > 1e-10) {
      out <- c(out, "the scaled map is not lambda P")
    }
  }

  # equation (scaled-lambda): the stationary point, against a numerical
  # minimization of the printed objective
  set.seed(3)
  beta <- as.numeric(stats::rnorm(6))
  quad <- drop(t(beta) %*% p_mat %*% beta)
  obj <- function(lambda) 0.5 * lambda * quad - 0.5 * s@rank * log(lambda)
  num <- stats::optimize(obj, c(1e-8, 1e8))$minimum
  if (abs(s@rank / quad - num) / (s@rank / quad) > 1e-4) {
    out <- c(out, sprintf(
      "lambda = r / (b'Pb) is %s while the numerical minimizer is %s",
      format(s@rank / quad, digits = 6), format(num, digits = 6)
    ))
  }

  # a fully known matrix has no free value at all
  if (parameters7::scaled_matrix(diag(3), link = NULL)@n_free != 0L) {
    out <- c(out, "a scaled structure with no link still has a free value")
  }
  out
}


# --- 5. Rank, the null space, and the refused solve ------------------------

.certify_rank_claims <- function() {
  out <- character()
  d2 <- function(k) crossprod(diff(diag(k), differences = 2))
  p1 <- kronecker(d2(4), diag(8))
  p2 <- kronecker(diag(4), d2(8))

  count_ev <- function(m, tol = 1e-10) {
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    sum(ev > tol * max(ev))
  }

  # The chapter's claim has two halves: the eigenvalue count of the assembled
  # sum DEGRADES with the ratio, and the stacked rank does not. Both are
  # checked, because the second alone would not distinguish a stacked route
  # that is right from an eigenvalue route that is also right.
  ns <- parameters7::param_null_basis(list(p1, p2))
  if (ns$rank != 28L) {
    out <- c(out, sprintf(
      "the stacked rank of the tensor penalty is %d, the chapter says 28", ns$rank
    ))
  }
  if (count_ev(p1 + p2) != 28L) {
    out <- c(out, "the eigenvalue count already disagrees at equal weights")
  }
  if (count_ev(p1 + 1e10 * p2) >= 28L) {
    out <- c(out, "the eigenvalue count does not degrade at a ratio of 1e10")
  }
  # and the null residual does not move, which is what the chapter shows
  for (f in c(1, 1e8, 1e12)) {
    m <- p1 + f * p2
    if (max(abs(m %*% ns$null_basis)) / max(abs(m)) > 1e-12) {
      out <- c(out, sprintf(
        "the stored null basis leaves a residual at a ratio of %s", format(f)
      ))
    }
  }

  # the null space of a second-difference penalty is the constants and the
  # lines, checked by the penalty itself rather than by the stored basis
  p1 <- d2(8)
  pen <- parameters7::scaled_matrix(p1)
  nb <- pen@null_basis
  if (ncol(nb) != 2L) {
    out <- c(out, "the second-difference penalty does not have a two-dimensional null space")
  }
  line <- 1 + 2 * seq_len(8)
  if (drop(t(line) %*% p1 %*% line) > 1e-18) {
    out <- c(out, "a straight line is penalized by the second-difference penalty")
  }
  if (max(abs(line - nb %*% crossprod(nb, line))) > 1e-9) {
    out <- c(out, "a straight line is not in the stored null basis")
  }
  if (drop(t(seq_len(8)^2) %*% p1 %*% seq_len(8)^2) < 1)  {
    out <- c(out, "a quadratic is not penalized by the second-difference penalty")
  }

  # the solve is refused rather than returning a pseudo-inverse
  refused <- tryCatch({
    parameters7::param_solve(pen, 0)
    FALSE
  }, error = function(e) TRUE)
  if (!refused) {
    out <- c(out, "a rank-deficient structure returned a solve instead of refusing")
  }

  # and a full-rank one solves, agreeing with base R
  s <- parameters7::log_cholesky(3)
  m <- parameters7::param_value(s, .cs_eta3)
  if (.cs_maxdiff(parameters7::param_solve(s, .cs_eta3), solve(m)) > 1e-9) {
    out <- c(out, "the full-rank solve disagrees with base R")
  }

  # a rank-deficient structure is not a density, and the distribution says so
  refused2 <- tryCatch({
    mvgaussian_distrib(8, sigma = pen)
    FALSE
  }, error = function(e) TRUE)
  if (!refused2) {
    out <- c(out, "a rank-deficient structure was accepted as a covariance")
  }
  out
}


# --- 5b. What a free name denotes -------------------------------------------

# Section 6.1 states a convention and names three coordinates as instances of
# it. The three are transcribed here rather than derived from the package's
# own tagging helper, which would compare the package with itself; and the
# property behind the convention -- that a labeled coordinate is genuinely
# unrestricted -- is exercised separately, since a spelling that happened to
# be right for the wrong reason would pass the first check alone.
.certify_free_names <- function() {
  out <- character()

  printed <- list(
    list(what = "log_scale", got = parameters7::scalar_matrix(3)@free_names),
    list(what = "z_rho", got = parameters7::ar1(5)@free_names[2L]),
    list(what = "L2.1", got = parameters7::log_cholesky(2)@free_names[3L])
  )
  for (p in printed) {
    if (!identical(p$got, p$what)) {
      out <- c(out, sprintf(
        "the section prints the coordinate %s and the package calls it %s",
        p$what, paste(p$got, collapse = ", ")
      ))
    }
  }

  # the tag names the transformation, not the bound: compound symmetry is
  # definite only above -1/(p-1), and no part of that number reaches the label
  cs <- parameters7::compound_symmetry(4)@free_names
  if (!identical(cs, c("log_scale", "logit_rho"))) {
    out <- c(out, sprintf(
      "compound symmetry names its coordinates %s", paste(cs, collapse = ", ")
    ))
  }

  # a labeled coordinate ranges over the whole line
  fams <- list(parameters7::scalar_matrix(3), parameters7::ar1(5),
               parameters7::compound_symmetry(4),
               parameters7::autoregressive(6, order = 2),
               parameters7::log_cholesky(3))
  for (s in fams) {
    for (k in seq_len(s@n_free)) {
      for (v in c(-12, 12)) {
        eta <- numeric(s@n_free)
        eta[k] <- v
        m <- parameters7::param_value(s, eta)
        ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
        if (!all(is.finite(m)) || min(ev) < -1e-13 * max(ev)) {
          out <- c(out, sprintf(
            "%s leaves its set at %s = %g", s@param_name, s@free_names[k], v
          ))
        }
      }
    }
  }

  # and the consumer really does flatten it into identity links, which is what
  # makes a label promising a bounded quantity a misreport rather than a
  # cosmetic complaint
  d <- distributions7::mvgaussian_distrib(3, sigma = parameters7::ar1(3))
  tags <- vapply(d@link_params, function(l) l@link_name, character(1))
  if (!all(tags == "identity")) {
    out <- c(out, "the multivariate parameters do not all carry identity links")
  }
  out
}


# --- 6. The gates ----------------------------------------------------------

assert_parameter_ok <- function() {
  problems <- c(.certify_logchol_map(), .certify_dlogdet_identity(),
                .certify_free_names())
  if (length(problems)) {
    stop("Section 6.1 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_parameter_families_ok <- function() {
  problems <- c(.certify_logchol_map(), .certify_diag_family(),
                .certify_scaled_family(), .certify_new_families(),
                .certify_phase2_families(),
                .certify_arp())
  if (length(problems)) {
    stop("Section 6.2 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_parameter_rank_ok <- function() {
  problems <- c(.certify_scaled_family(), .certify_rank_claims())
  if (length(problems)) {
    stop("Section 6.3 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}


# The claims of sections 6.2.4 to 6.2.6: the simplex identity the chapter
# prints, the two exactnesses of the matrix logarithm, the row-locality of a
# transition matrix, and the closed higher orders against one product stencil
# on the map -- a route none of the closed forms share.
.certify_new_families <- function() {
  out <- character()

  s <- parameters7::simplex(4)
  eta <- c(0.5, -0.2, 0.8)
  v <- parameters7::param_value(s, eta)
  if (abs(sum(v) - 1) > 1e-12 || any(v <= 0)) {
    out <- c(out, "the simplex value printed by 6.2.4 is not on the simplex")
  }
  for (o in 1:4) {
    d <- switch(o, parameters7::param_d1, parameters7::param_d2,
                parameters7::param_d3, parameters7::param_d4)(s, eta)
    if (max(abs(vapply(d, sum, numeric(1)))) > 1e-10) {
      out <- c(out, sprintf(
        "a simplex derivative of order %d does not sum to zero", o
      ))
    }
  }

  m <- parameters7::matrix_log(3)
  em <- c(0.2, -0.3, 0.4, 0.5, -0.1, 0.3)
  if (abs(parameters7::param_logdet(m, em) - sum(em[1:3])) > 1e-12) {
    out <- c(out, "the matrix logarithm's log-determinant is not the trace")
  }
  if (max(abs(parameters7::param_solve(m, em) -
              solve(parameters7::param_value(m, em)))) > 1e-10) {
    out <- c(out, "the matrix logarithm's inverse is not exp(-S)")
  }
  a3 <- parameters7::param_d3(m, em)
  n3 <- parameters7::numerical_d3(m, em)
  worst <- max(mapply(function(a, b) {
    max(abs(a - b)) / max(1, max(abs(b)))
  }, a3, n3))
  if (worst > 1e-4) {
    out <- c(out, sprintf(
      "matrix_log third derivatives disagree with one stencil by %.1e", worst
    ))
  }

  tm <- parameters7::transition_matrix(3)
  set.seed(5)
  et <- stats::rnorm(tm@n_free, sd = 0.6)
  d2 <- parameters7::param_d2(tm, et)
  idx <- parameters7::param_tuple_indices(tm, 2L)
  row_of <- rep(seq_len(3), each = 2)
  for (i in seq_along(idx)) {
    rows <- row_of[idx[[i]]]
    if (rows[1L] != rows[2L] && max(abs(d2[[i]])) != 0) {
      out <- c(out, "a cross-row transition component is not exactly zero")
      break
    }
  }
  out
}


# The claims of sections 6.2.4 and 6.2.5. Each is checked against a route the
# chapter does not use: the correlation's unit diagonal against the value
# itself, its log-determinant against the spectrum, the two economical
# families' eigenvalue formulas against an eigendecomposition, the definiteness
# bound of compound symmetry against the matrix at a correlation just inside
# and just outside it, and the AR(1) precision's bandwidth against a general
# inverse.
.certify_phase2_families <- function() {
  out <- character()

  s <- parameters7::correlation_matrix(4)
  eta <- c(0.4, -0.2, 0.6, 0.1, -0.5, 0.3)
  r <- parameters7::param_value(s, eta)
  if (max(abs(diag(r) - 1)) > 1e-12) {
    out <- c(out, "the correlation matrix printed by 6.2.4 has no unit diagonal")
  }
  ev <- eigen(r, symmetric = TRUE, only.values = TRUE)$values
  if (abs(parameters7::param_logdet(s, eta) - sum(log(ev))) > 1e-9) {
    out <- c(out, "the correlation's log-determinant disagrees with its spectrum")
  }

  # the two distinct eigenvalues compound symmetry claims
  p <- 5
  cs <- parameters7::compound_symmetry(p)
  ecs <- c(log(1.7), 0.6)
  m <- parameters7::param_value(cs, ecs)
  v <- m[1L, 1L]
  rho <- m[2L, 1L] / v
  want <- sort(c(v * (1 + (p - 1) * rho), rep(v * (1 - rho), p - 1)))
  if (max(abs(sort(eigen(m, symmetric = TRUE, only.values = TRUE)$values) -
              want)) > 1e-9) {
    out <- c(out, "compound symmetry's eigenvalues are not the two 6.2.5 states")
  }

  # the bound -1/(p-1): just inside it the matrix is definite, just outside
  # it is not, so the interval the link maps onto is the right one
  pat <- function(rho) (1 - rho) * diag(p) + rho * matrix(1, p, p)
  inside <- min(eigen(pat(-1 / (p - 1) + 1e-6), symmetric = TRUE,
                      only.values = TRUE)$values)
  outside <- min(eigen(pat(-1 / (p - 1) - 1e-3), symmetric = TRUE,
                       only.values = TRUE)$values)
  if (!(inside > 0 && outside < 0)) {
    out <- c(out, "the definiteness bound -1/(p-1) is not where 6.2.5 says")
  }
  bounds <- cs@param_params$link_rho@link_bounds
  if (max(abs(unname(bounds) - c(-1 / (p - 1), 1))) > 1e-12) {
    out <- c(out, "the correlation link of compound symmetry uses another bound")
  }

  # AR(1): the determinant formula and the tridiagonal precision
  a <- parameters7::ar1(5)
  ea <- c(log(2), atanh(0.6))
  ma <- parameters7::param_value(a, ea)
  if (abs(parameters7::param_logdet(a, ea) -
          (5 * log(2) + 4 * log(1 - 0.6^2))) > 1e-9) {
    out <- c(out, "the AR(1) log-determinant is not p log s + (p-1) log(1-r^2)")
  }
  pr <- parameters7::param_solve(a, ea)
  far <- abs(row(pr) - col(pr)) > 1
  if (max(abs(pr[far])) != 0 ||
      max(abs(pr - solve(ma))) > 1e-10) {
    out <- c(out, "the AR(1) precision is not the tridiagonal inverse")
  }

  # and the compound-symmetric inverse is compound symmetric, which is the
  # closure under sides the section claims
  pc <- parameters7::param_solve(cs, ecs)
  if (diff(range(diag(pc))) > 1e-10 ||
      diff(range(pc[lower.tri(pc)])) > 1e-10) {
    out <- c(out, "the compound-symmetric inverse is not compound symmetric")
  }
  out
}


# The claims of section 6.2.6, each against a route the chapter does not use:
# the stationary region against the roots of the polynomial the chapter
# writes, the log-determinant against a spectrum, the band against a general
# inverse, and the order-one case against the independently written ar1().
.certify_arp <- function() {
  out <- character()
  s <- parameters7::autoregressive(7, order = 2L)
  eta <- c(log(2), atanh(0.7), atanh(-0.3))

  phi <- parameters7:::ar_taylor(s, eta)$phi[, 1L]
  if (min(Mod(polyroot(c(1, -phi)))) <= 1) {
    out <- c(out, "the chart of 6.2.6 left the stationary region")
  }
  # and the triangle the section names is the region at order two
  if (!(abs(phi[2L]) < 1 && phi[2L] + phi[1L] < 1 && phi[2L] - phi[1L] < 1)) {
    out <- c(out, "the order-two coefficients are outside the stated triangle")
  }

  m <- parameters7::param_value(s, eta)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  if (abs(parameters7::param_logdet(s, eta) - sum(log(ev))) > 1e-9) {
    out <- c(out, "the AR(q) log-determinant disagrees with its spectrum")
  }

  pr <- parameters7::param_solve(s, eta)
  if (max(abs(pr[abs(row(pr) - col(pr)) > 2L])) != 0 ||
      max(abs(pr - solve(m))) > 1e-9) {
    out <- c(out, "the AR(q) precision is not the banded exact inverse")
  }

  # the two implementations of the order-one case
  s1 <- parameters7::autoregressive(6, order = 1L)
  a1 <- parameters7::ar1(6)
  e1 <- c(log(2), atanh(0.6))
  gap <- max(abs(parameters7::param_value(s1, e1) -
                 parameters7::param_value(a1, e1)))
  if (gap > 1e-12) {
    out <- c(out, sprintf(
      "autoregressive(order = 1) and ar1() disagree by %.1e", gap
    ))
  }
  out
}
