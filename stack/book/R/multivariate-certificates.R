# Consistency gate for section 3.5, the multivariate families.
#
# Same contract as the other gates: every formula the section prints is
# transcribed BY HAND below and compared against what the package computes,
# and every structural claim is checked against a route the package does not
# take. Nothing here is visible to the reader.
#
# What "by hand" means precisely. The matrix itself comes from the structure,
# since building a log-Cholesky factor is chapter 6's subject and chapter 6's
# gate; everything the section claims ON TOP of it -- the density, the score,
# the Hessian, the trace form of the information, the weight of the t and its
# score in nu -- is written out here from the printed equations and from
# nothing else. The scores and Hessians are ALSO compared against numerical
# derivatives of the hand-written log-density, which shares no code with the
# package at all.

.mv_maxdiff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

# One numerical derivative, never nested (the rule of check_link()).
.mv_grad <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(k) {
    up <- dn <- x
    up[k] <- x[k] + h
    dn[k] <- x[k] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

# A second derivative from ONE stencil, for the same reason.
.mv_hess1 <- function(f, x, k, l, h = 1e-4) {
  hk <- h * max(1, abs(x[k]))
  hl <- h * max(1, abs(x[l]))
  if (k == l) {
    up <- dn <- x
    up[k] <- x[k] + hk
    dn[k] <- x[k] - hk
    return((f(up) - 2 * f(x) + f(dn)) / hk^2)
  }
  pp <- pm <- mp <- mm <- x
  pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
  pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
  mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
  mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
  (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * hk * hl)
}

# The fixture every check below runs on: a two-dimensional family with a
# log-Cholesky scale, a handful of observations including one far out, and the
# ingredients the printed equations are written in.
.mv_fixture <- function(nu = NULL) {
  s <- covstructs7::log_cholesky(2)
  eta <- c(0.1, -0.2, 0.5)
  mu <- c(0.4, -0.3)
  sigma <- covstructs7::struct_matrix(s, eta)
  y <- rbind(c(0, 0), c(1, -1), c(-0.5, 0.8), c(3.2, 2.4))
  list(
    s = s, eta = eta, mu = mu, sigma = sigma, y = y, p = 2L, nu = nu,
    a = covstructs7::struct_dmatrix(s, eta),
    a2 = covstructs7::struct_d2matrix(s, eta)
  )
}

# --- 1. The gaussian density, equation (mvgaussian-density) ------------------

# Written from the printed equation, with the determinant and the solve taken
# from base R rather than from the structure.
.mv_dens_gauss <- function(y, mu, sigma) {
  p <- length(mu)
  apply(y, 1, function(row) {
    r <- row - mu
    -0.5 * (p * log(2 * pi) + log(det(sigma)) +
      drop(t(r) %*% solve(sigma) %*% r))
  })
}

.certify_mvgauss_density <- function() {
  f <- .mv_fixture()
  d <- mvgaussian_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta), d@params))
  err <- .mv_maxdiff(
    distrib_pdf(d, f$y, th, log = TRUE),
    .mv_dens_gauss(f$y, f$mu, f$sigma)
  )
  if (err > 1e-12) {
    return(sprintf(
      "the gaussian density of the chapter and the package differ by %s",
      format(err, digits = 3)
    ))
  }
  character()
}

# --- 2. The score, equation (mvgaussian-score) ------------------------------

.mv_score_gauss <- function(y, mu, sigma, a, dld) {
  si <- solve(sigma)
  r <- sweep(y, 2L, mu)
  w <- r %*% si
  mean_part <- lapply(seq_along(mu), function(j) w[, j])
  mat_part <- lapply(seq_along(a), function(k) {
    -0.5 * dld[[k]] + 0.5 * rowSums((w %*% a[[k]]) * w)
  })
  c(mean_part, mat_part)
}

.certify_mvgauss_score <- function() {
  f <- .mv_fixture()
  d <- mvgaussian_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta), d@params))
  dld <- covstructs7::struct_dlogdet(f$s, f$eta)

  book <- .mv_score_gauss(f$y, f$mu, f$sigma, f$a, dld)
  pkg <- distrib_gradient(d, f$y, th)
  out <- character()
  err <- max(vapply(seq_along(book), function(i)
    .mv_maxdiff(book[[i]], pkg[[i]]), numeric(1)))
  if (err > 1e-10) {
    out <- c(out, sprintf(
      "the printed gaussian score and the package differ by %s",
      format(err, digits = 3)
    ))
  }

  # and the same score against a numerical derivative of the HAND-WRITTEN
  # log-likelihood, which shares nothing with the package
  ll <- function(v) {
    sum(.mv_dens_gauss(f$y, v[1:2], covstructs7::struct_matrix(f$s, v[3:5])))
  }
  num <- .mv_grad(ll, c(f$mu, f$eta))
  err2 <- .mv_maxdiff(vapply(book, sum, numeric(1)), num)
  if (err2 > 1e-5) {
    out <- c(out, sprintf(
      "the printed gaussian score disagrees with finite differences by %s",
      format(err2, digits = 3)
    ))
  }
  out
}

# --- 3. The Hessian, equation (mvgaussian-hessian) --------------------------

.certify_mvgauss_hessian <- function() {
  f <- .mv_fixture()
  d <- mvgaussian_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta), d@params))
  si <- solve(f$sigma)
  r <- sweep(f$y, 2L, f$mu)
  w <- r %*% si
  d2ld <- covstructs7::struct_d2logdet(f$s, f$eta)
  pairs <- covstructs7::struct_pair_indices(f$s)
  pnames <- covstructs7::struct_pair_names(f$s)

  # The three blocks, transcribed from the printed equation.
  book <- list()
  for (a in 1:2) for (b in a:2) {
    book[[paste0("mu", a, "_mu", b)]] <- rep(-si[a, b], nrow(f$y))
  }
  for (a in 1:2) for (k in seq_along(f$a)) {
    book[[paste0("mu", a, "|", k)]] <- -((w %*% t(si %*% f$a[[k]]))[, a])
  }
  for (i in seq_along(pairs)) {
    k <- pairs[[i]][1L]
    l <- pairs[[i]][2L]
    book[[paste0("eta", k, "_eta", l)]] <-
      -0.5 * d2ld[[i]] + 0.5 * rowSums((w %*% f$a2[[i]]) * w) -
        rowSums((w %*% f$a[[l]] %*% si %*% f$a[[k]]) * w)
  }

  pkg <- distrib_hessian(d, f$y, th)
  # The DISTRIBUTION's names, not the structure's: a multivariate family
  # prefixes its structure's free values with the matrix they describe, and
  # a lookup built from the bare names would silently miss every component.
  free <- setdiff(d@params, paste0("mu", seq_len(f$p)))
  key <- function(nm) {
    # the package keys a component by parameter names; the transcription above
    # keys it by position, so the two are matched here explicitly
    nm
  }
  out <- character()
  lookup <- c(
    setNames(c("mu1_mu1", "mu1_mu2", "mu2_mu2"),
             c("mu1_mu1", "mu1_mu2", "mu2_mu2")),
    setNames(paste0("mu", rep(1:2, each = length(free)), "_",
                    rep(free, 2)),
             paste0("mu", rep(1:2, each = length(free)), "|",
                    rep(seq_along(free), 2))),
    setNames(vapply(pairs, function(p) paste0(free[p[1]], "_", free[p[2]]),
                    character(1)),
             vapply(pairs, function(p) paste0("eta", p[1], "_eta", p[2]),
                    character(1)))
  )
  for (nm in names(book)) {
    target <- lookup[[nm]]
    if (is.null(target) || is.null(pkg[[target]])) {
      out <- c(out, sprintf("Hessian component '%s' is not in the package", nm))
      next
    }
    err <- .mv_maxdiff(book[[nm]], pkg[[target]])
    if (err > 1e-9) {
      out <- c(out, sprintf(
        "the printed gaussian Hessian component %s differs by %s",
        target, format(err, digits = 3)
      ))
    }
  }

  # and against one mixed stencil on the hand-written log-likelihood
  ll <- function(v) {
    sum(.mv_dens_gauss(f$y, v[1:2], covstructs7::struct_matrix(f$s, v[3:5])))
  }
  v0 <- c(f$mu, f$eta)
  prs <- distributions7:::hess_pairs(d@params)
  for (nm in names(pkg)) {
    k <- match(prs[[nm]], d@params)
    num <- .mv_hess1(ll, v0, k[1], k[2])
    err <- abs(sum(pkg[[nm]]) - num) / max(1, abs(num))
    if (err > 1e-5) {
      out <- c(out, sprintf(
        "Hessian component %s disagrees with finite differences by %s",
        nm, format(err, digits = 3)
      ))
    }
  }
  out
}

# --- 4. The information, equation (mvgaussian-information) ------------------

.certify_mvgauss_information <- function() {
  f <- .mv_fixture()
  d <- mvgaussian_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta), d@params))
  si <- solve(f$sigma)
  # The DISTRIBUTION's names, not the structure's: a multivariate family
  # prefixes its structure's free values with the matrix they describe, and
  # a lookup built from the bare names would silently miss every component.
  free <- setdiff(d@params, paste0("mu", seq_len(f$p)))
  pkg <- distrib_expected_hessian(d, f$y[1, , drop = FALSE], th)

  out <- character()
  # the mean block and the exact zero of the mixed block
  for (a in 1:2) for (b in a:2) {
    nm <- paste0("mu", a, "_mu", b)
    if (abs(pkg[[nm]][1] + si[a, b]) > 1e-12) {
      out <- c(out, sprintf("the information block %s is not -Sigma^{-1}", nm))
    }
  }
  for (a in 1:2) for (fr in free) {
    nm <- if (!is.null(pkg[[paste0("mu", a, "_", fr)]])) {
      paste0("mu", a, "_", fr)
    } else {
      paste0(fr, "_mu", a)
    }
    if (abs(pkg[[nm]][1]) > 1e-14) {
      out <- c(out, sprintf("the mixed information block %s is not zero", nm))
    }
  }
  # the trace form, transcribed
  pairs <- covstructs7::struct_pair_indices(f$s)
  for (i in seq_along(pairs)) {
    k <- pairs[[i]][1L]
    l <- pairs[[i]][2L]
    want <- -0.5 * sum(diag(si %*% f$a[[k]] %*% si %*% f$a[[l]]))
    nm <- paste0(free[k], "_", free[l])
    err <- abs(pkg[[nm]][1] - want) / max(1, abs(want))
    if (err > 1e-10) {
      out <- c(out, sprintf(
        "the printed trace form for %s differs by %s", nm, format(err, digits = 3)
      ))
    }
  }

  # and the whole thing against a Monte Carlo average of the OBSERVED Hessian,
  # which is the second Bartlett identity taken as a fact about the family
  set.seed(4021)
  big <- distrib_rng(d, 2e5, th)
  hb <- distrib_hessian(d, big, th)
  # The tolerance is the sampling error of the average itself plus a small
  # relative allowance, rather than a flat number: a block that is exactly zero
  # is compared against noise, and a flat relative tolerance either passes it
  # for free or fails it for being noisy.
  for (nm in names(pkg)) {
    se <- stats::sd(hb[[nm]]) / sqrt(length(hb[[nm]]))
    if (abs(pkg[[nm]][1] - mean(hb[[nm]])) > 4 * se + 0.01 * abs(pkg[[nm]][1])) {
      out <- c(out, sprintf(
        "the information component %s is %s while the sample average is %s",
        nm, format(pkg[[nm]][1], digits = 4), format(mean(hb[[nm]]), digits = 4)
      ))
    }
  }
  out
}

# --- 5. The t: density, weight, score ---------------------------------------

.mv_dens_t <- function(y, mu, sigma, nu) {
  p <- length(mu)
  apply(y, 1, function(row) {
    r <- row - mu
    q <- drop(t(r) %*% solve(sigma) %*% r)
    lgamma((nu + p) / 2) - lgamma(nu / 2) - (p / 2) * log(nu * pi) -
      0.5 * log(det(sigma)) - ((nu + p) / 2) * log1p(q / nu)
  })
}

.certify_mvt_density <- function() {
  f <- .mv_fixture(nu = 4)
  d <- mvstudent_t_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta, 4), d@params))
  out <- character()
  err <- .mv_maxdiff(
    distrib_pdf(d, f$y, th, log = TRUE),
    .mv_dens_t(f$y, f$mu, f$sigma, 4)
  )
  if (err > 1e-12) {
    out <- c(out, sprintf(
      "the t density of the chapter and the package differ by %s",
      format(err, digits = 3)
    ))
  }

  # equation (mvt-mixture): the density is the chi-squared mixture of gaussians,
  # integrated here rather than taken from any formula above
  row <- c(1, -1)
  mix <- stats::integrate(function(v) {
    vapply(v, function(vi) {
      sc <- (4 / vi) * f$sigma
      exp(.mv_dens_gauss(matrix(row, 1L), f$mu, sc)) * stats::dchisq(vi, 4)
    }, numeric(1))
  }, 0, Inf, rel.tol = 1e-10)$value
  err2 <- abs(mix - exp(.mv_dens_t(matrix(row, 1L), f$mu, f$sigma, 4)))
  if (err2 > 1e-8) {
    out <- c(out, sprintf(
      "the scale-mixture representation is off by %s", format(err2, digits = 3)
    ))
  }
  out
}

.certify_mvt_score <- function() {
  f <- .mv_fixture(nu = 4)
  nu <- 4
  p <- 2
  d <- mvstudent_t_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta, nu), d@params))

  si <- solve(f$sigma)
  r <- sweep(f$y, 2L, f$mu)
  w <- r %*% si
  q <- rowSums(r * w)
  cw <- (nu + p) / (nu + q)                      # equation (mvt-weight)
  dld <- covstructs7::struct_dlogdet(f$s, f$eta)

  book <- c(
    lapply(1:2, function(j) cw * w[, j]),        # equation (mvt-score)
    lapply(seq_along(f$a), function(k) {
      -0.5 * dld[[k]] + 0.5 * cw * rowSums((w %*% f$a[[k]]) * w)
    }),
    list(0.5 * (                                 # equation (mvt-score-nu)
      digamma((nu + p) / 2) - digamma(nu / 2) - p / nu -
        log1p(q / nu) + (nu + p) * q / (nu * (nu + q))
    ))
  )

  pkg <- distrib_gradient(d, f$y, th)
  out <- character()
  err <- max(vapply(seq_along(book), function(i)
    .mv_maxdiff(book[[i]], pkg[[i]]), numeric(1)))
  if (err > 1e-10) {
    out <- c(out, sprintf(
      "the printed t score and the package differ by %s", format(err, digits = 3)
    ))
  }

  ll <- function(v) {
    sum(.mv_dens_t(f$y, v[1:2], covstructs7::struct_matrix(f$s, v[3:5]), v[6]))
  }
  num <- .mv_grad(ll, c(f$mu, f$eta, nu))
  err2 <- .mv_maxdiff(vapply(book, sum, numeric(1)), num)
  if (err2 > 1e-5) {
    out <- c(out, sprintf(
      "the printed t score disagrees with finite differences by %s",
      format(err2, digits = 3)
    ))
  }

  # the link scale is the chain rule and nothing else: d/d(log nu) = nu * d/dnu
  gl <- distrib_gradient(d, f$y, th, scale = "link")
  if (.mv_maxdiff(gl[["nu"]], nu * pkg[["nu"]]) > 1e-10) {
    out <- c(out, "the link-scale score in nu is not nu times the parameter one")
  }
  out
}

# --- 6. Moments, marginals and the precision form ---------------------------

.certify_mv_moments <- function() {
  f <- .mv_fixture()
  d <- mvstudent_t_distrib(2)
  out <- character()

  th <- as.list(stats::setNames(c(f$mu, f$eta, 6), d@params))
  if (.mv_maxdiff(variance(d, th), f$sigma * 6 / 4) > 1e-12) {
    out <- c(out, "the t covariance is not nu Sigma / (nu - 2)")
  }
  if (.mv_maxdiff(mv_sigma(d, th), f$sigma) > 1e-12) {
    out <- c(out, "mv_sigma() does not return the scale matrix")
  }
  th2 <- th; th2$nu <- 1.5
  if (!all(is.infinite(variance(d, th2)))) {
    out <- c(out, "the covariance is finite at nu = 1.5")
  }
  if (!is.finite(distrib_pdf(d, c(0, 0), th2))) {
    out <- c(out, "the density is not finite at nu = 1.5")
  }
  th3 <- th; th3$nu <- 0.5
  if (!all(is.nan(mean(d, th3)))) {
    out <- c(out, "the mean is not NaN at nu = 0.5")
  }

  # the marginal: same family, block of the matrix, and for the t the same nu
  d3 <- mvstudent_t_distrib(3)
  s3 <- covstructs7::log_cholesky(3)
  th3f <- as.list(stats::setNames(
    c(0, 1, -1, 0.1, -0.1, 0.2, 0.5, -0.2, 0.3, 7), d3@params
  ))
  m <- mv_marginal(d3, th3f, c(1L, 3L))
  blk <- mv_sigma(d3, th3f)[c(1, 3), c(1, 3)]
  if (.mv_maxdiff(mv_sigma(m$distrib, m$theta), blk) > 1e-12) {
    out <- c(out, "the t marginal is not the block of the scale matrix")
  }
  if (abs(m$theta$nu - 7) > 1e-12) {
    out <- c(out, "the t marginal changed the degrees of freedom")
  }

  g3 <- mvgaussian_distrib(3)
  thg <- as.list(stats::setNames(
    c(0, 1, -1, 0.1, -0.1, 0.2, 0.5, -0.2, 0.3), g3@params
  ))
  mg <- mv_marginal(g3, thg, c(1L, 2L))
  if (.mv_maxdiff(mv_sigma(mg$distrib, mg$theta),
                  mv_sigma(g3, thg)[1:2, 1:2]) > 1e-12) {
    out <- c(out, "the gaussian marginal is not the block of the covariance")
  }

  # the two parametrisations of the gaussian describe the same law
  ds <- mvgaussian_distrib(2)
  do <- mvgaussian_distrib(2, struct_omega = covstructs7::log_cholesky(2))
  ths <- as.list(stats::setNames(c(f$mu, f$eta), ds@params))
  etao <- covstructs7::struct_free(do@struct, solve(mv_sigma(ds, ths)))
  tho <- as.list(stats::setNames(c(f$mu, unname(etao)), do@params))
  if (.mv_maxdiff(distrib_pdf(do, f$y, tho, log = TRUE),
                  distrib_pdf(ds, f$y, ths, log = TRUE)) > 1e-10) {
    out <- c(out, "the precision form is not the same law as the covariance form")
  }
  out
}

# --- 7. What a fit reports: the delta method and its Jacobian ---------------

.certify_mv_reporting <- function() {
  f <- .mv_fixture()
  d <- mvgaussian_distrib(2)
  th <- as.list(stats::setNames(c(f$mu, f$eta), d@params))
  out <- character()

  # The decomposition itself, written out here.
  s <- unname(f$sigma)
  sd_book <- sqrt(diag(s))
  rho_book <- s[1, 2] / (sd_book[1] * sd_book[2])
  der <- mv_derived(d, th)
  if (.mv_maxdiff(der$value, c(sd_book, rho_book)) > 1e-12) {
    out <- c(out, "the printed decomposition Sigma = D R D disagrees with the package")
  }

  # Equation (mv-sdcor-jacobian), transcribed, against the package's Jacobian.
  a <- covstructs7::struct_dmatrix(f$s, f$eta)
  jb <- matrix(0, 3, length(th))
  for (l in seq_along(a)) {
    al <- unname(a[[l]])
    col <- 2L + l
    jb[1, col] <- al[1, 1] / (2 * sd_book[1])
    jb[2, col] <- al[2, 2] / (2 * sd_book[2])
    jb[3, col] <- al[1, 2] / (sd_book[1] * sd_book[2]) -
      (rho_book / 2) * (al[1, 1] / s[1, 1] + al[2, 2] / s[2, 2])
  }
  if (.mv_maxdiff(der$jacobian, jb) > 1e-12) {
    out <- c(out, "the printed Jacobian of (sd, cor) disagrees with the package")
  }
  # and against a numerical derivative of the decomposition, which shares no
  # code with either
  num <- numDeriv::jacobian(function(v) {
    m <- mv_sigma(d, as.list(stats::setNames(v, d@params)))
    sdv <- sqrt(diag(m))
    c(sdv, m[1, 2] / (sdv[1] * sdv[2]))
  }, unlist(th))
  if (.mv_maxdiff(der$jacobian, num) > 1e-6) {
    out <- c(out, "the Jacobian of (sd, cor) disagrees with finite differences")
  }

  # Equation (mv-delta): the standard errors are J V J', and the intervals are
  # built on log and on Fisher's z, so neither can leave its own set.
  set.seed(5107)
  y <- distrib_rng(d, 600, th)
  fit <- fit_distrib(d, y)
  tab <- mv_summary(fit)
  dh <- mv_derived(d, as.list(coef(fit)))
  se <- sqrt(diag(dh$jacobian %*% vcov(fit) %*% t(dh$jacobian)))
  if (.mv_maxdiff(tab[["Std. Error"]], se) > 1e-12) {
    out <- c(out, "the reported standard errors are not J V J'")
  }
  if (any(tab[1:2, 3] <= 0)) {
    out <- c(out, "a standard deviation's interval reaches zero")
  }
  if (abs(tab[3, 3]) >= 1 || abs(tab[3, 4]) >= 1) {
    out <- c(out, "a correlation's interval leaves (-1, 1)")
  }
  # the interval really is the image of a symmetric one on the stated scale
  z <- stats::qnorm(0.975)
  want_lo <- exp(log(tab[1, 1]) - z * tab[1, 2] / tab[1, 1])
  if (abs(tab[1, 3] - want_lo) > 1e-10) {
    out <- c(out, "the standard deviation's interval is not built on the log scale")
  }
  want_lo_r <- tanh(atanh(tab[3, 1]) - z * tab[3, 2] / (1 - tab[3, 1]^2))
  if (abs(tab[3, 3] - want_lo_r) > 1e-10) {
    out <- c(out, "the correlation's interval is not built on Fisher's z")
  }

  # The precision form: same law, so the same standard deviations and
  # correlations, plus the readings that are its own.
  do <- mvgaussian_distrib(3, struct_omega = covstructs7::log_cholesky(3))
  ds <- mvgaussian_distrib(3)
  set.seed(5108)
  ths <- generate_random_theta(ds)
  sig3 <- unname(mv_sigma(ds, ths))
  om3 <- solve(sig3)
  tho <- as.list(stats::setNames(
    c(unlist(ths)[1:3], unname(covstructs7::struct_free(do@struct, om3))),
    do@params
  ))
  vs <- mv_derived(ds, ths)$value
  vo <- mv_derived(do, tho)$value
  if (.mv_maxdiff(vo[names(vs)], vs) > 1e-8) {
    out <- c(out, "the two parametrisations report different standard deviations")
  }
  # Equation (mv-conditional-variance), and its regression reading: the ratio
  # of the conditional to the marginal variance is 1 - R^2 for the regression
  # of that coordinate on all the others, which is computed here by least
  # squares and nowhere in the package.
  if (.mv_maxdiff(vo[paste0("cvar_v", 1:3)], 1 / diag(om3)) > 1e-8) {
    out <- c(out, "the conditional variances are not 1/Omega_jj")
  }
  r2 <- vapply(1:3, function(j) {
    # R^2 of Y_j on the others, from the covariance alone
    s_jj <- sig3[j, j]
    s_jo <- sig3[j, -j, drop = FALSE]
    s_oo <- sig3[-j, -j, drop = FALSE]
    drop(s_jo %*% solve(s_oo, t(s_jo))) / s_jj
  }, numeric(1))
  if (.mv_maxdiff(unname(vo[paste0("cvar_v", 1:3)]) / diag(sig3), 1 - r2) > 1e-8) {
    out <- c(out, "the conditional variance is not (1 - R^2) times the marginal one")
  }
  # Equation (mv-partial-correlation)
  if (abs(vo[["pcor_v1_v2"]] +
          om3[1, 2] / sqrt(om3[1, 1] * om3[2, 2])) > 1e-8) {
    out <- c(out, "the partial correlation is not the negated normalised precision")
  }

  # The two parametrisations describe the same law, so a fit of each must
  # reach the same maximised likelihood and the same information criteria.
  # A disagreement means one of them has not converged, which is exactly the
  # failure this invariance is worth asserting against.
  set.seed(5109)
  y3 <- distrib_rng(ds, 400, ths)
  fs_ <- fit_distrib(ds, y3)
  fo_ <- fit_distrib(do, y3)
  if (abs(as.numeric(logLik(fs_)) - as.numeric(logLik(fo_))) > 1e-6 ||
      abs(fs_@aic - fo_@aic) > 1e-6 || abs(fs_@bic - fo_@bic) > 1e-6) {
    out <- c(out, sprintf(
      "the covariance and precision fits disagree: logLik %s against %s",
      format(as.numeric(logLik(fs_)), digits = 10),
      format(as.numeric(logLik(fo_)), digits = 10)
    ))
  }
  if (!fs_@converged || !fo_@converged) {
    out <- c(out, "a four-parameter gaussian fit did not converge")
  }

  # The t reports the SCALE standard deviations, and says so in the name.
  dt <- mvstudent_t_distrib(2)
  tht <- as.list(stats::setNames(c(f$mu, f$eta, 6), dt@params))
  vt <- mv_derived(dt, tht)$value
  if (!all(c("scale_sd_v1", "scale_sd_v2", "cor_v1_v2") %in% names(vt))) {
    out <- c(out, "the t does not name its diagonal quantities as scale quantities")
  }
  if (.mv_maxdiff(vt[1:2], sqrt(diag(s))) > 1e-12) {
    out <- c(out, "the t's scale standard deviations are not the scale matrix's")
  }
  # its correlations are the response's, the covariance being a multiple
  if (abs(vt[["cor_v1_v2"]] - rho_book) > 1e-12) {
    out <- c(out, "the t's correlation differs from its scale matrix's")
  }
  out
}


# --- 8. The gate -----------------------------------------------------------

assert_multivariate_ok <- function() {
  problems <- c(
    .certify_mvgauss_density(),
    .certify_mvgauss_score(),
    .certify_mvgauss_hessian(),
    .certify_mvgauss_information(),
    .certify_mvt_density(),
    .certify_mvt_score(),
    .certify_mv_moments(),
    .certify_mv_reporting()
  )
  if (length(problems)) {
    stop("Section 3.5 disagrees with the packages:\n  ",
         paste(problems, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}
