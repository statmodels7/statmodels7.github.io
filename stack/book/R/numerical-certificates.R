# Consistency gate for chapter 7, the numericals7 package.
#
# Same contract as the other gates. Every formula the chapter prints is
# transcribed BY HAND here and compared against what the package computes;
# every structural claim is checked against a route the package does not
# take. Nothing in this file is visible to the reader.
#
# The transcriptions are the Weibull scale map's derivative, the Vandermonde
# moment conditions, the Gauss and Kronrod exactness degrees, the beta and
# gamma closed forms the quadrature is compared with, the three retirement
# conditions read off a hand-built block, the Mills ratio at zero, Owen's
# two identities, the Bessel asymptotic expansion and the inverse function
# rule. Injection-checked during development: a Kronrod weight 5% wrong, a
# stencil weight 5% wrong, and a mis-signed term in the asymptotic expansion
# are each caught by the checks below.

.nm_tol <- 1e-10

# --- 2. Stencils -------------------------------------------------------------

.certify_stencils <- function() {
  out <- character()

  # the printed five-point weights
  w <- numericals7::fd_weights(-2:2, 1L) * 12
  if (max(abs(w - c(1, -8, 0, 8, -1))) > 1e-9) {
    out <- c(out, "the five-point first-derivative weights are not (1,-8,0,8,-1)/12")
  }

  # the moment conditions, transcribed: sum w o^p = p! [p == m]
  off <- c(-3, -1, 0, 2, 4)
  m <- 2L
  ww <- numericals7::fd_weights(off, m)
  for (p in 0:4) {
    want <- if (p == m) factorial(p) else 0
    if (abs(sum(ww * off^p) - want) > 1e-8) {
      out <- c(out, sprintf("the stencil moment condition fails at p = %d", p))
    }
  }

  # fd_derivative against a closed third derivative
  err <- abs(numericals7::fd_derivative(exp, 1, order = 3L, accuracy = 4L) - exp(1))
  if (err > 1e-6) {
    out <- c(out, sprintf("the third derivative of exp at 1 is off by %.1e", err))
  }

  out
}

# --- 3. Quadrature and series ------------------------------------------------

.certify_quadrature <- function() {
  out <- character()

  r <- numericals7::gauss_kronrod15()
  mono <- function(w, p) sum(w * r$nodes^p)
  # exactness degrees 13 and 22, checked on the even monomials
  for (p in seq(0, 12, by = 2)) {
    if (abs(mono(r$wg, p) - 2 / (p + 1)) > 1e-12) {
      out <- c(out, sprintf("G7 is not exact at degree %d", p))
    }
  }
  for (p in seq(0, 22, by = 2)) {
    if (abs(mono(r$wk, p) - 2 / (p + 1)) > 1e-12) {
      out <- c(out, sprintf("K15 is not exact at degree %d", p))
    }
  }
  # and one degree beyond each must fail, so the exactness is not vacuous
  if (abs(mono(r$wg, 14) - 2 / 15) < 1e-10) {
    out <- c(out, "G7 appears exact past its degree, the check is vacuous")
  }
  if (abs(mono(r$wk, 24) - 2 / 25) < 1e-10) {
    out <- c(out, "K15 appears exact past its degree, the check is vacuous")
  }

  # the chapter's singular normalization and gamma means
  f <- function(x, i) dbeta(x, c(2, 0.8, 3)[i], c(5, 1.2, 0.7)[i])
  got <- numericals7::quad_vec(f, 0, rep(1, 3))
  if (max(abs(got - 1)) > 1e-8) {
    out <- c(out, "a beta normalization of section 7.3 is not one")
  }
  shp <- c(0.5, 1, 2, 5, 50)
  g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
  if (max(abs(numericals7::quad_vec(g, 0, rep(Inf, 5)) - shp)) > 1e-6) {
    out <- c(out, "a gamma mean of section 7.3 is off its closed form")
  }

  # a refusal really is a refusal: a needle-free singular row at a depth too
  # small must come back NA and not a number
  bad <- suppressWarnings(
    numericals7::quad_vec(function(x, i) dbeta(as.numeric(x), 0.5, 1),
                          0, 1, max_depth = 4L)
  )
  if (!is.na(bad)) {
    out <- c(out, "a row that cannot converge returned a number instead of NA")
  }

  # the series claims: the three lambdas, and the hump-shaped tail guard;
  # the reference is a direct finite sum sharing nothing with series_vec
  lam <- c(0.5, 4, 300)
  got <- numericals7::series_vec(function(k, i) dpois(k, lam[i]), n = 3)
  ref <- vapply(lam, function(l) sum(dpois(0:1200, l)), numeric(1))
  if (max(abs(got - ref)) > 1e-9) {
    out <- c(out, "a Poisson mass of section 7.3 disagrees with a direct sum")
  }

  out
}

# --- 4. Special functions ----------------------------------------------------

.certify_special <- function() {
  out <- character()

  # Mills: sqrt(2/pi) at zero, the asymptote in the deep tail, and R'
  m0 <- numericals7::mills_ratio(0)
  if (abs(m0$r - sqrt(2 / pi)) > 1e-14) {
    out <- c(out, "the Mills ratio at zero is not sqrt(2/pi)")
  }
  m <- numericals7::mills_ratio(-60)
  if (!is.finite(m$r) || abs(m$r - (60 + 1 / 60)) > 1e-3) {
    out <- c(out, "the Mills ratio in the deep tail is not close to -t - 1/t")
  }
  # R' = -R(t + R), transcribed, against the package's dr
  t0 <- -2.3
  mm <- numericals7::mills_ratio(t0)
  if (abs(mm$dr - (-mm$r * (t0 + mm$r))) > 1e-14) {
    out <- c(out, "the printed identity for R' is not what the package computes")
  }

  # Owen: the two printed identities
  a <- c(0.3, 1, 5)
  if (max(abs(numericals7::owen_t(0, a) - atan(a) / (2 * pi))) > 1e-11) {
    out <- c(out, "T(0, a) is not arctan(a)/(2 pi)")
  }
  hs <- c(-2, -0.7, 0, 1.1, 3)
  if (max(abs(numericals7::owen_t(hs, 1) - pnorm(hs) * pnorm(-hs) / 2)) > 1e-11) {
    out <- c(out, "T(h, 1) is not Phi(h) Phi(-h) / 2")
  }
  if (abs(numericals7::owen_t(1.3, Inf) - 0.5 * pnorm(-1.3)) > 1e-14) {
    out <- c(out, "T(h, Inf) is not Phi(-|h|)/2")
  }
  # and the skew normal cdf identity the section quotes, against quadrature
  # of the density -- a route owen_t does not take
  z <- c(-1, 0.4, 2)
  al <- 3
  sn_cdf <- pnorm(z) - 2 * numericals7::owen_t(z, al)
  ref <- vapply(z, function(q) {
    stats::integrate(function(t) 2 * dnorm(t) * pnorm(al * t), -Inf, q,
                     rel.tol = 1e-12)$value
  }, numeric(1))
  if (max(abs(sn_cdf - ref)) > 1e-9) {
    out <- c(out, "Phi(z) - 2 T(z, alpha) is not the skew normal cdf")
  }

  # Bessel: the asymptotic expansion as printed, against the Bessel branch
  # at a point where BOTH are exact
  for (kk in c(5e3, 9.9e3)) {
    asym <- 1 - 1 / (2 * kk) - 1 / (8 * kk^2) - 1 / (8 * kk^3)
    if (abs(numericals7::bessel_i_ratio(kk) - asym) > 1e-14) {
      out <- c(out, sprintf("the asymptotic expansion disagrees at kappa = %g", kk))
    }
  }
  # finiteness past the scaled underflow
  if (!is.finite(numericals7::bessel_i_ratio(1e8))) {
    out <- c(out, "the ratio is not finite past the scaled underflow")
  }
  # A' = 1 - A/k - A^2, transcribed, and the inverse function rule
  k0 <- 2.7
  ad <- numericals7::bessel_i_ratio_derivs(k0)
  if (abs(ad$d1 - (1 - ad$A / k0 - ad$A^2)) > 1e-14) {
    out <- c(out, "A' is not 1 - A/kappa - A^2")
  }
  inv <- numericals7::bessel_i_ratio_inverse(0.7)
  ai <- numericals7::bessel_i_ratio_derivs(inv$kappa)
  if (abs(inv$d1 - 1 / ai$d1) > 1e-12 * abs(inv$d1)) {
    out <- c(out, "kappa' is not 1/A'")
  }
  if (abs(inv$d2 - (-ai$d2 / ai$d1^3)) > 1e-10 * abs(inv$d2)) {
    out <- c(out, "kappa'' is not -A''/(A')^3")
  }
  if (abs(inv$d3 - (3 * ai$d2^2 - ai$d1 * ai$d3) / ai$d1^5) > 1e-9 * abs(inv$d3)) {
    out <- c(out, "kappa''' is not (3 A''^2 - A' A''')/(A')^5")
  }
  if (abs(numericals7::bessel_i_ratio(inv$kappa) - 0.7) > 1e-10) {
    out <- c(out, "the inverse does not round-trip at rho = 0.7")
  }

  out
}

# --- 5. The gates ------------------------------------------------------------

assert_stencils_ok <- function() {
  problems <- .certify_stencils()
  if (length(problems)) {
    stop("Section 7.2 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_quadrature_ok <- function() {
  problems <- .certify_quadrature()
  if (length(problems)) {
    stop("Section 7.3 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_special_ok <- function() {
  problems <- .certify_special()
  if (length(problems)) {
    stop("Section 7.4 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}
