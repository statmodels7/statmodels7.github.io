# Certificates for Chapter 4.
#
# The wrappers are the part of the toolkit that cannot be checked against a
# reference implementation, because they compose: the Hessian of a zero-inflated
# distribution is assembled at run time from whatever parent it was handed. So
# every check here is either structural (an identity that must hold between two
# separately implemented quantities) or an independent re-derivation written in
# this file.

# ---------------------------------------------------------------------------
# The twelve transformers
# ---------------------------------------------------------------------------

TRANSFORMERS <- list(
  list(id = "log", title = "Logarithm", ctor = "log_transform()",
       obj = function() log_transform(),
       g = "\\log x", ginv = "e^{y}", jac = "e^{y}",
       domain = "x > 0", decreasing = FALSE),
  list(id = "exp", title = "Exponential", ctor = "exp_transform()",
       obj = function() exp_transform(),
       g = "e^{x}", ginv = "\\log y", jac = "1/y",
       domain = "x \\in \\mathbb{R}", decreasing = FALSE),
  list(id = "sqrt", title = "Square root", ctor = "sqrt_transform()",
       obj = function() sqrt_transform(),
       g = "\\sqrt{x}", ginv = "y^{2}", jac = "2y",
       domain = "x \\geq 0", decreasing = FALSE),
  list(id = "inverse", title = "Reciprocal", ctor = "inverse_transform()",
       obj = function() inverse_transform(),
       g = "1/x", ginv = "1/y", jac = "1/y^{2}",
       domain = "0 \\notin \\operatorname{supp}", decreasing = TRUE),
  list(id = "power", title = "Power", ctor = "power_transform(p = 2)",
       obj = function() power_transform(2),
       g = "x^{p}", ginv = "\\operatorname{sign}(y)\\,\\lvert y\\rvert^{1/p}",
       jac = "\\dfrac{\\lvert y\\rvert^{1/p - 1}}{\\lvert p\\rvert}",
       domain = "p\\text{-dependent}", decreasing = FALSE),
  list(id = "bc", title = "Box-Cox", ctor = "bc_transform(lambda = 0.5)",
       obj = function() bc_transform(0.5),
       g = "\\dfrac{x^{\\lambda}-1}{\\lambda}", ginv = "(\\lambda y + 1)^{1/\\lambda}",
       jac = "(\\lambda y + 1)^{1/\\lambda - 1}",
       domain = "x \\geq 0", decreasing = FALSE),
  list(id = "yj", title = "Yeo-Johnson", ctor = "yj_transform(lambda = 0.5)",
       obj = function() yj_transform(0.5),
       g = "\\begin{cases}\\frac{(x+1)^{\\lambda}-1}{\\lambda} & x \\geq 0\\\\ -\\frac{(1-x)^{2-\\lambda}-1}{2-\\lambda} & x < 0\\end{cases}",
       ginv = "\\text{piecewise inverse}", jac = "(\\lvert x\\rvert+1)^{\\operatorname{sign}(x)(\\lambda-1)}",
       domain = "x \\in \\mathbb{R}", decreasing = FALSE),
  list(id = "softplus", title = "Softplus", ctor = "softplus_transform(a = 1)",
       obj = function() softplus_transform(1),
       g = "\\tfrac{1}{a}\\log\\!\\left(e^{ax}-1\\right)",
       ginv = "\\tfrac{1}{a}\\log\\!\\left(1+e^{ay}\\right)",
       jac = "\\sigma(ay)", domain = "x \\geq 0", decreasing = FALSE),
  list(id = "asinh", title = "Inverse hyperbolic sine", ctor = "asinh_transform()",
       obj = function() asinh_transform(),
       g = "\\operatorname{arsinh} x", ginv = "\\sinh y", jac = "\\cosh y",
       domain = "x \\in \\mathbb{R}", decreasing = FALSE),
  list(id = "logit", title = "Logit", ctor = "logit_transform()",
       obj = function() logit_transform(),
       g = "\\log\\dfrac{x}{1-x}", ginv = "\\sigma(y)", jac = "\\sigma(y)(1-\\sigma(y))",
       domain = "x \\in (0,1)", decreasing = FALSE),
  list(id = "expit", title = "Expit", ctor = "expit_transform()",
       obj = function() expit_transform(),
       g = "\\sigma(x)", ginv = "\\log\\dfrac{y}{1-y}", jac = "\\dfrac{1}{y(1-y)}",
       domain = "x \\in \\mathbb{R}", decreasing = FALSE),
  list(id = "affine", title = "Affine", ctor = "affine_transform(loc = 1, scale = 2)",
       obj = function() affine_transform(1, 2),
       g = "c + s\\,x", ginv = "\\dfrac{y-c}{s}", jac = "1/\\lvert s\\rvert",
       domain = "x \\in \\mathbb{R}", decreasing = FALSE)
)

# A parent whose support is compatible with each transformer, and a grid in x.
transformer_parent <- function(id) {
  switch(id,
    log = , sqrt = , inverse = , power = , bc = , softplus =
      list(d = gamma_distrib(), th = list(mu = 2, sigma2 = 1), x = seq(0.4, 5, length.out = 30)),
    logit =
      list(d = beta_distrib(), th = list(mu = 0.4, phi = 6), x = seq(0.06, 0.94, length.out = 30)),
    list(d = gaussian_distrib(), th = list(mu = 0.5, sigma = 1),
         x = seq(-2.2, 2.6, length.out = 30))
  )
}

# Certificate: the transformer's own contract.
#   round trip   g^{-1}(g(x)) = x
#   jacobian     |J(y)| against a numerical derivative of g^{-1}
#   density      f_Y(y) against f_X(g^{-1}(y)) |J(y)|, assembled here
#   score        the transformed score against the parent's at g^{-1}(y)
transformer_certificate <- function() {
  rows <- lapply(TRANSFORMERS, function(tr) {
    tf <- tr$obj()
    pp <- transformer_parent(tr$id)
    d <- pp$d; th <- pp$th; x <- pp$x

    e_round <- tryCatch(rel_err(tf@trans_inv(tf@trans_fun(x)), x), error = function(e) NA_real_)

    y <- tf@trans_fun(x)
    y <- sort(y)
    e_jac <- tryCatch({
      num <- num_grad1(function(v) tf@trans_inv(v), y)
      rel_err(tf@trans_abs_jac(y, log = FALSE), abs(num))
    }, error = function(e) NA_real_)

    td <- tryCatch(transformation(d, tf), error = function(e) NULL)
    if (is.null(td)) {
      e_dens <- NA_real_; e_score <- NA_real_; e_mass <- NA_real_
    } else {
      # f_Y(y) = f_X(g^{-1}(y)) |J(y)|, assembled independently
      e_dens <- tryCatch({
        book <- distrib_pdf(d, tf@trans_inv(y), th, log = TRUE) + tf@trans_abs_jac(y, log = TRUE)
        rel_err(book, distrib_pdf(td, y, th, log = TRUE))
      }, error = function(e) NA_real_)

      e_score <- tryCatch({
        a <- distrib_gradient(td, y, th)
        b <- distrib_gradient(d, tf@trans_inv(y), th)
        max(vapply(names(a), function(k) rel_err(a[[k]], b[[k]]), numeric(1)))
      }, error = function(e) NA_real_)

      e_mass <- tryCatch({
        bb <- td@bounds
        m <- suppressWarnings(distrib_quantile(td, 0.5, th))
        kn <- sort(unique(c(bb[1], m, bb[2])))
        tot <- sum(vapply(seq_len(length(kn) - 1L), function(k)
          stats::integrate(function(t) distrib_pdf(td, t, th), kn[k], kn[k + 1],
                           rel.tol = 1e-9)$value, numeric(1)))
        abs(tot - 1)
      }, error = function(e) NA_real_)
    }

    data.frame(transformer = tr$title,
               round_trip = e_round, jacobian = e_jac, density = e_dens,
               score = e_score, mass = e_mass, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

format_transformer_certificate <- function(tab) {
  worst <- pmax(tab$round_trip, tab$jacobian, tab$density, tab$score, na.rm = TRUE)
  data.frame(
    transformer = tab$transformer,
    `round trip` = fmt_err(tab$round_trip),
    `Jacobian` = fmt_err(tab$jacobian),
    `density` = fmt_err(tab$density),
    `score` = fmt_err(tab$score),
    `total mass` = fmt_err(tab$mass),
    verdict = ifelse(is.na(worst) | worst > 1e-6, "DISAGREE", "agree"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Certificates for the zero wrappers
# ---------------------------------------------------------------------------

zero_cases <- function() {
  list(
    zip = list(title = "Zero-inflated Poisson",
               d = zero_inflated(poisson_distrib()), th = list(mu = 3, zi = 0.25)),
    zinb = list(title = "Zero-inflated NB2",
                d = zero_inflated(negbin_distrib()), th = list(mu = 4, theta = 1.5, zi = 0.3)),
    zap = list(title = "Hurdle Poisson",
               d = zero_adjusted(poisson_distrib()), th = list(mu = 3, za = 0.4)),
    zanb = list(title = "Hurdle NB2",
                d = zero_adjusted(negbin_distrib()), th = list(mu = 4, theta = 1.5, za = 0.35)),
    zaga = list(title = "Zero-adjusted gamma",
                d = zero_adjusted(gamma_distrib()), th = list(mu = 3, sigma2 = 2, za = 0.3)),
    zabe = list(title = "Zero-adjusted beta",
                d = zero_adjusted(beta_distrib()), th = list(mu = 0.4, phi = 6, za = 0.25))
  )
}

# Independent re-implementations of the derived pmf/pdf, written from the
# formulas printed in the chapter.
book_zi_logpdf <- function(parent, y, th_parent, zi) {
  f0 <- distrib_pdf(parent, 0, th_parent)
  log(zi * (y == 0) + (1 - zi) * distrib_pdf(parent, y, th_parent))
}
book_za_discrete_logpdf <- function(parent, y, th_parent, za) {
  f0 <- distrib_pdf(parent, 0, th_parent)
  ifelse(y == 0, log(za),
         log(1 - za) + distrib_pdf(parent, y, th_parent, log = TRUE) - log1p(-f0))
}
book_za_continuous_logpdf <- function(parent, y, th_parent, za) {
  ifelse(y == 0, log(za), log(1 - za) + distrib_pdf(parent, y, th_parent, log = TRUE))
}

# Certificate 1: the derived density, against the implementation.
zero_density_certificate <- function() {
  cs <- zero_cases()
  rows <- lapply(names(cs), function(nm) {
    cse <- cs[[nm]]
    d <- cse$d; th <- cse$th
    parent <- d@parent_distrib
    p_th <- th[seq_len(d@n_params - 1L)]
    mix <- th[[d@n_params]]
    is_disc <- S7::S7_inherits(d, discrete_distrib)
    y <- if (is_disc) 0:25 else c(0, seq(0.02, 12, length.out = 40))
    if (nm == "zabe") y <- c(0, seq(0.02, 0.98, length.out = 40))

    book <- if (grepl("^zi", nm)) book_zi_logpdf(parent, y, p_th, mix)
            else if (is_disc) book_za_discrete_logpdf(parent, y, p_th, mix)
            else book_za_continuous_logpdf(parent, y, p_th, mix)

    e_dens <- rel_err(book, distrib_pdf(d, y, th, log = TRUE))

    # total mass: sum over the support, or integral plus the atom
    tot <- if (is_disc) {
      numerical_series(function(k) distrib_pdf(d, k, th), d@bounds[1], d@bounds[2])
    } else {
      at <- distrib_atoms(d, th)
      bb <- d@bounds
      kn <- sort(unique(c(bb[1], suppressWarnings(distrib_quantile(d, 0.75, th)), bb[2])))
      sum(at$p) + sum(vapply(seq_len(length(kn) - 1L), function(k)
        stats::integrate(function(t) distrib_pdf(d, t, th), kn[k], kn[k + 1],
                         rel.tol = 1e-10)$value, numeric(1)))
    }

    data.frame(model = cse$title, density = e_dens, mass = abs(tot - 1),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Certificate 2: derivatives. The score and Hessian against finite differences of
# the log-density, and the expected Hessian against the expectation of the
# observed one -- the latter being a structural identity, since the two are
# implemented by different formulas.
zero_derivative_certificate <- function(seed = 7) {
  set.seed(seed)
  cs <- zero_cases()
  rows <- lapply(names(cs), function(nm) {
    cse <- cs[[nm]]
    d <- cse$d; th <- cse$th
    y <- distrib_rng(d, 40, th)
    if (!any(y == 0)) y[1] <- 0            # make sure both branches are exercised

    fd_g <- function(h) {
      out <- lapply(seq_along(th), function(j) {
        tp <- tm <- th
        hh <- h * max(1, abs(th[[j]]))
        tp[[j]] <- th[[j]] + hh; tm[[j]] <- th[[j]] - hh
        (distrib_pdf(d, y, tp, log = TRUE) - distrib_pdf(d, y, tm, log = TRUE)) / (2 * hh)
      })
      names(out) <- names(th)
      out
    }
    ag <- distrib_gradient(d, y, th)
    ng <- fd_g(1e-5)
    e_grad <- max(vapply(names(th), function(k) rel_err(ag[[k]], ng[[k]]), numeric(1)))

    ah <- distrib_hessian(d, y, th)
    e_hess <- max(vapply(names(ah), function(k) {
      pr <- book_hess_pair(names(th), k)
      i <- match(pr[1], names(th)); j <- match(pr[2], names(th))
      hi <- 1e-4 * max(1, abs(th[[i]])); hj <- 1e-4 * max(1, abs(th[[j]]))
      shift <- function(a, b) { t2 <- th; t2[[i]] <- t2[[i]] + a * hi; t2[[j]] <- t2[[j]] + b * hj; t2 }
      num <- (distrib_pdf(d, y, shift(1, 1), log = TRUE) -
              distrib_pdf(d, y, shift(1, -1), log = TRUE) -
              distrib_pdf(d, y, shift(-1, 1), log = TRUE) +
              distrib_pdf(d, y, shift(-1, -1), log = TRUE)) / (4 * hi * hj)
      rel_err(ah[[k]], num)
    }, numeric(1)))

    ae <- distrib_expected_hessian(d, 0, th)
    e_exp <- max(vapply(names(ae), function(k) {
      num <- expectation(d, function(y, theta) distrib_hessian(d, y, theta)[[k]], th)
      rel_err(ae[[k]][1], num)
    }, numeric(1)))

    data.frame(model = cse$title, score = e_grad, hessian = e_hess,
               expected = e_exp, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Certificate 3: the closed-form moments derived in the chapter.
zero_moment_certificate <- function() {
  mu <- 3; zi <- 0.25
  zip <- zero_inflated(poisson_distrib()); th_zip <- list(mu = mu, zi = zi)
  za <- 0.4
  zap <- zero_adjusted(poisson_distrib()); th_zap <- list(mu = mu, za = za)
  zag <- zero_adjusted(gamma_distrib()); th_zag <- list(mu = 3, sigma2 = 2, za = 0.3)
  f0 <- dpois(0, mu)

  data.frame(
    model = c("ZI Poisson", "ZI Poisson", "Hurdle Poisson", "ZA gamma", "ZA gamma"),
    quantity = c("$\\mathbb{E}[Y] = (1-\\zeta)\\mu$",
                 "$\\operatorname{Var}(Y) = (1-\\zeta)\\mu + \\zeta(1-\\zeta)\\mu^{2}$",
                 "$\\mathbb{E}[Y] = (1-\\pi)\\mu/(1-f_0)$",
                 "$\\mathbb{E}[Y] = (1-\\pi)\\mu$",
                 "$\\operatorname{Var}(Y) = (1-\\pi)\\sigma^{2} + \\pi(1-\\pi)\\mu^{2}$"),
    error = c(
      abs(mean(zip, th_zip) - (1 - zi) * mu),
      abs(variance(zip, th_zip) - ((1 - zi) * mu + zi * (1 - zi) * mu^2)),
      abs(mean(zap, th_zap) - (1 - za) * mu / (1 - f0)),
      abs(mean(zag, th_zag) - 0.7 * 3),
      abs(variance(zag, th_zag) - (0.7 * 2 + 0.3 * 0.7 * 9))
    ),
    stringsAsFactors = FALSE
  )
}

# Certificate 4: check_distrib() on the wrappers.
zero_check_certificate <- function(n = 40, nsim = 5e4, orders = 1:2, seed = 20260726) {
  set.seed(seed)
  cs <- zero_cases()
  rows <- lapply(names(cs), function(nm) {
    cse <- cs[[nm]]
    out <- check_distrib(cse$d, theta = cse$th, n = n, nsim = nsim,
                         orders = orders, verbose = FALSE)
    failed <- out$check[out$status != "OK"]
    data.frame(model = cse$title, checks = nrow(out), passed = sum(out$status == "OK"),
               failures = if (length(failed)) paste(failed, collapse = "; ") else "none",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# Certificates for truncation
# ---------------------------------------------------------------------------

truncation_cases <- function() {
  list(
    ztp  = list(title = "Zero-truncated Poisson",
                d = truncated(poisson_distrib(), lower = 1), th = list(mu = 2.5)),
    tpu  = list(title = "Poisson, truncated above at 8",
                d = truncated(poisson_distrib(), upper = 8), th = list(mu = 3)),
    tnb  = list(title = "NB2 on [1, 15]",
                d = truncated(negbin_distrib(), 1, 15), th = list(mu = 4, theta = 1.5)),
    tn   = list(title = "Gaussian on [-1, 2]",
                d = truncated(gaussian_distrib(), -1, 2), th = list(mu = 0.5, sigma = 1.5)),
    tg   = list(title = "Gamma on [0.5, 8]",
                d = truncated(gamma_distrib(), 0.5, 8), th = list(mu = 3, sigma2 = 2)),
    tzag = list(title = "Zero-adjusted gamma, capped at 5",
                d = truncated(zero_adjusted(gamma_distrib()), upper = 5),
                th = list(mu = 3, sigma2 = 2, za = 0.3))
  )
}

# Two truncated families have textbook closed forms. They are the only external
# reference available for this wrapper, and they are used.
truncated_closed_form_certificate <- function() {
  mu <- 2
  ztp <- truncated(poisson_distrib(), lower = 1)
  lo <- -1; up <- 2; m <- 0.5; s <- 1.5
  tn <- truncated(gaussian_distrib(), lower = lo, upper = up)
  a <- (lo - m) / s; b <- (up - m) / s; Z <- pnorm(b) - pnorm(a)

  data.frame(
    model = c("Zero-truncated Poisson", "Zero-truncated Poisson",
              "Truncated Gaussian", "Truncated Gaussian", "Truncated Gaussian"),
    quantity = c("$P(Y=y) = f(y)/(1-f_0)$",
                 "$\\mathbb{E}[Y] = \\mu/(1-e^{-\\mu})$",
                 "$f_T(y) = \\varphi_{\\mu,\\sigma}(y)/Z$",
                 "$\\mathbb{E}[Y] = \\mu + \\sigma\\frac{\\varphi(a)-\\varphi(b)}{Z}$",
                 "$\\operatorname{Var}(Y)$ (standard formula)"),
    error = c(
      rel_err(distrib_pdf(ztp, 1:8, list(mu = mu)), dpois(1:8, mu) / (1 - dpois(0, mu))),
      abs(mean(ztp, list(mu = mu)) - mu / (1 - exp(-mu))),
      rel_err(distrib_pdf(tn, seq(-0.9, 1.9, length.out = 20), list(mu = m, sigma = s)),
              dnorm(seq(-0.9, 1.9, length.out = 20), m, s) / Z),
      abs(mean(tn, list(mu = m, sigma = s)) - (m + s * (dnorm(a) - dnorm(b)) / Z)),
      abs(variance(tn, list(mu = m, sigma = s)) -
            s^2 * (1 + (a * dnorm(a) - b * dnorm(b)) / Z - ((dnorm(a) - dnorm(b)) / Z)^2))
    ),
    stringsAsFactors = FALSE
  )
}

# Derivatives: score and Hessian against finite differences of the log-density;
# E[score] = 0 (the first Bartlett identity, which the recentring makes automatic);
# and the closed-form expected Hessian -Cov_T(s) against the expectation of the
# observed Hessian.
truncated_derivative_certificate <- function(seed = 19) {
  set.seed(seed)
  cs <- truncation_cases()
  rows <- lapply(names(cs), function(nm) {
    cse <- cs[[nm]]
    d <- cse$d; th <- cse$th
    y <- distrib_rng(d, 25, th)

    ag <- distrib_gradient(d, y, th)
    ng <- lapply(seq_along(th), function(j) {
      tp <- tm <- th; hh <- 1e-5 * max(1, abs(th[[j]]))
      tp[[j]] <- th[[j]] + hh; tm[[j]] <- th[[j]] - hh
      (distrib_pdf(d, y, tp, log = TRUE) - distrib_pdf(d, y, tm, log = TRUE)) / (2 * hh)
    })
    names(ng) <- names(th)
    e_grad <- max(vapply(names(th), function(k) rel_err(ag[[k]], ng[[k]]), numeric(1)))

    ah <- distrib_hessian(d, y, th)
    e_hess <- max(vapply(names(ah), function(k) {
      pr <- book_hess_pair(names(th), k)
      i <- match(pr[1], names(th)); j <- match(pr[2], names(th))
      hi <- 1e-4 * max(1, abs(th[[i]])); hj <- 1e-4 * max(1, abs(th[[j]]))
      sh <- function(p, q) { t2 <- th; t2[[i]] <- t2[[i]] + p * hi; t2[[j]] <- t2[[j]] + q * hj; t2 }
      num <- (distrib_pdf(d, y, sh(1, 1), log = TRUE) - distrib_pdf(d, y, sh(1, -1), log = TRUE) -
              distrib_pdf(d, y, sh(-1, 1), log = TRUE) + distrib_pdf(d, y, sh(-1, -1), log = TRUE)) /
             (4 * hi * hj)
      rel_err(ah[[k]], num)
    }, numeric(1)))

    e_score0 <- max(vapply(names(th), function(k) {
      abs(expectation(d, function(y, theta) distrib_gradient(d, y, theta)[[k]], th))
    }, numeric(1)))

    ae <- distrib_expected_hessian(d, 0, th)
    e_exp <- max(vapply(names(ae), function(k) {
      num <- expectation(d, function(y, theta) distrib_hessian(d, y, theta)[[k]], th)
      rel_err(ae[[k]][1], num)
    }, numeric(1)))

    data.frame(model = cse$title, score = e_grad, hessian = e_hess,
               score_mean = e_score0, expected = e_exp, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

truncated_check_certificate <- function(n = 40, nsim = 3e4, orders = 1:2, seed = 20260726) {
  set.seed(seed)
  cs <- truncation_cases()
  rows <- lapply(names(cs), function(nm) {
    cse <- cs[[nm]]
    out <- check_distrib(cse$d, theta = cse$th, n = n, nsim = nsim,
                         orders = orders, verbose = FALSE)
    failed <- out$check[out$status != "OK"]
    data.frame(model = cse$title, checks = nrow(out), passed = sum(out$status == "OK"),
               failures = if (length(failed)) paste(failed, collapse = "; ") else "none",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


# Silent consistency gate for Chapter 4. Stops the render if any derived
# formula -- density, score, Hessian, expected Hessian, moment, closed form --
# disagrees with the implementation, or if check_distrib() fails on any wrapper.
assert_transformations_ok <- function() {
  tc <- transformer_certificate()
  worst <- pmax(tc$round_trip, tc$jacobian, tc$density, tc$score, na.rm = TRUE)
  if (any(is.na(worst) | worst > 1e-6)) {
    stop("A transformer disagrees with its printed contract: ",
         paste(tc$transformer[is.na(worst) | worst > 1e-6], collapse = ", "), call. = FALSE)
  }
  zd <- zero_density_certificate()
  if (any(zd$density > 1e-10 | zd$mass > 1e-8)) {
    stop("A zero-wrapper density disagrees with its derivation.", call. = FALSE)
  }
  zdc <- zero_derivative_certificate()
  if (any(zdc$score > 1e-6 | zdc$hessian > 1e-4 | zdc$expected > 1e-6)) {
    stop("A zero-wrapper derivative disagrees with its derivation.", call. = FALSE)
  }
  zm <- zero_moment_certificate()
  if (any(zm$error > 1e-6)) stop("A zero-wrapper moment disagrees.", call. = FALSE)
  tcf <- truncated_closed_form_certificate()
  if (any(tcf$error > 1e-6)) stop("Truncation disagrees with a closed form.", call. = FALSE)
  tdc <- truncated_derivative_certificate()
  if (any(tdc$score > 1e-6 | tdc$hessian > 1e-4 |
          tdc$score_mean > 1e-6 | tdc$expected > 1e-6)) {
    stop("A truncated derivative disagrees with its derivation.", call. = FALSE)
  }
  zc <- zero_check_certificate()
  tcc <- truncated_check_certificate()
  if (any(zc$checks != zc$passed) || any(tcc$checks != tcc$passed)) {
    stop("check_distrib() failed on a wrapper model.", call. = FALSE)
  }

  # Third and fourth derivatives: the partition machinery against finite
  # differences, and the closed form for the pure zi component at zero.
  hi <- list(
    list(d = transformation(gaussian_distrib(), exp_transform()),
         th = list(mu = 0.5, sigma = 1.1), y = c(0.6, 1.4, 3.0)),
    list(d = zero_inflated(negbin_distrib()),
         th = list(mu = 4, theta = 1.5, zi = 0.3), y = c(0, 2, 6)),
    list(d = zero_adjusted(poisson_distrib()),
         th = list(mu = 3, za = 0.4), y = c(0, 1, 5)),
    list(d = zero_adjusted(gamma_distrib()),
         th = list(mu = 3, sigma2 = 2, za = 0.3), y = c(0, 1.2, 5)),
    list(d = truncated(gaussian_distrib(), -1, 3),
         th = list(mu = 0.5, sigma = 1.5), y = c(-0.5, 0.3, 2.4)),
    list(d = truncated(poisson_distrib(), lower = 1), th = list(mu = 2.5), y = c(1, 3, 6))
  )
  for (cs in hi) {
    for (k in 3:4) {
      a <- if (k == 3) distrib_deriv3(cs$d, cs$y, cs$th) else distrib_deriv4(cs$d, cs$y, cs$th)
      n <- if (k == 3) numerical_deriv3(cs$d, cs$y, cs$th) else numerical_deriv4(cs$d, cs$y, cs$th)
      e <- max(vapply(names(n), function(j) rel_err(a[[j]], n[[j]]), numeric(1)))
      if (!is.finite(e) || e > 1e-4) {
        stop(sprintf("order-%d derivatives of '%s' disagree with finite differences (%.2e).",
                     k, cs$d@distrib_name, e), call. = FALSE)
      }
    }
  }
  zd <- zero_inflated(poisson_distrib())
  f0 <- dpois(0, 3); L0 <- 0.25 + 0.75 * f0; r <- (1 - f0) / L0
  if (abs(distrib_deriv3(zd, 0, list(mu = 3, zi = 0.25))[["zi_zi_zi"]] - 2 * r^3) > 1e-10) {
    stop("The closed form for the pure zi third derivative no longer matches.", call. = FALSE)
  }
  invisible(TRUE)
}
