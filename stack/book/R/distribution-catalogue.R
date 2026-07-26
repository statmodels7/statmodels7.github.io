# The distribution catalogue, and the certificates for Chapter 3.
#
# As in R/link-formulas.R, each record keeps the printed formula and an
# independent R transcription of it side by side. The transcription is written
# from the density as displayed in the book, using only lgamma/lbeta/besselK and
# arithmetic -- never the `d*` function the package itself calls -- so that the
# comparison actually tests the parametrisation rather than restating it.

DISTRIBS <- list()

DISTRIBS$gaussian <- list(
  title = "Gaussian",
  ctor = "gaussian_distrib()",
  obj = function() gaussian_distrib(),
  theta = list(mu = 1.5, sigma = 2),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{\\sigma\\sqrt{2\\pi}}\\exp\\!\\left\\{-\\frac{(y-\\mu)^{2}}{2\\sigma^{2}}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}",
  ld = function(y, th) {
    mu <- th$mu; s <- th$sigma
    -log(s) - 0.5 * log(2 * pi) - (y - mu)^2 / (2 * s^2)
  },
  grid = function(th) seq(-4, 7, length.out = 40),
  note = "The reference case: $\\sigma$ is the standard deviation, as in
  `dnorm(y, mu, sigma)`."
)

DISTRIBS$cauchy <- list(
  title = "Cauchy",
  ctor = "cauchy_distrib()",
  obj = function() cauchy_distrib(),
  theta = list(mu = 0.5, sigma = 1.4),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{\\pi\\sigma\\left[1 + \\left(\\frac{y-\\mu}{\\sigma}\\right)^{2}\\right]}",
  moments = "no moments exist",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    -log(pi) - log(th$sigma) - log1p(z^2)
  },
  grid = function(th) seq(-12, 13, length.out = 40),
  note = "Neither the mean nor the variance exists, which is why
  `check_distrib()` validates the generator against the **cdf** rather than
  against the first two moments: a moment-based check would have nothing to
  compare here, and the cdf is both always available and the stronger claim."
)

DISTRIBS$logistic <- list(
  title = "Logistic",
  ctor = "logistic_distrib()",
  obj = function() logistic_distrib(),
  theta = list(mu = 0.5, sigma = 1.4),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{e^{-z}}{\\sigma\\left(1+e^{-z}\\right)^{2}}, \\qquad z = \\frac{y-\\mu}{\\sigma}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\pi^{2}\\sigma^{2}/3",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    -z - log(th$sigma) - 2 * log1p(exp(-z))
  },
  grid = function(th) seq(-9, 10, length.out = 40),
  note = "Writing $t = \\sigma(z) = (1+e^{-z})^{-1}$ and $u = 1-t$, the
  log-density is $-\\log\\sigma + \\log t + \\log u$ and every observed derivative
  is a polynomial in $t$: with $g_1 = 1-2t$, $g_2 = -2tu$, $g_3 = g_2 g_1$,
  $g_4 = g_2(1-6tu)$, the recursions
  $\\partial_\\mu[A/\\sigma^{k}] = -A'/\\sigma^{k+1}$ and
  $\\partial_\\sigma[A/\\sigma^{k}] = -(zA' + kA)/\\sigma^{k+1}$ generate all of
  them, overflow-free out to $|z| = 4000$. The **expected** derivatives were
  rederived here because symbolic algebra systems return, for this family,
  expressions containing polylogarithms, the imaginary unit and a dependence on
  $\\mu$ --- which is impossible for a location-scale family and is a useful
  reminder that a computer-algebra answer is not automatically a correct one.
  Seven of the nine expected third-order components are known exactly
  ($\\mathbb{E}[\\ell_{\\mu\\mu\\sigma}] = 1/(2\\sigma^{3})$,
  $\\mathbb{E}[\\ell_{\\sigma\\sigma\\sigma}] = (\\pi^{2}+2)/(2\\sigma^{3})$,
  $\\mathbb{E}[\\ell_{\\mu^{4}}] = 1/(15\\sigma^{4})$, the rest zero by symmetry);
  the two that remain require $\\int w^{k}\\operatorname{sech}^{4}w\\tanh^{2}w\\,dw$
  for $k = 2, 4$ and are **open**. In the meantime `approx = \"integrate\"`
  supplies them to within Monte Carlo noise."
)

DISTRIBS$student_t <- list(
  title = "Student t",
  ctor = "student_t_distrib()",
  obj = function() student_t_distrib(),
  theta = list(mu = 0.5, sigma = 1.3, nu = 6),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{\\Gamma\\!\\left(\\frac{\\nu+1}{2}\\right)}{\\sigma\\sqrt{\\nu\\pi}\\,\\Gamma\\!\\left(\\frac{\\nu}{2}\\right)}\\left(1+\\frac{z^{2}}{\\nu}\\right)^{-\\frac{\\nu+1}{2}}, \\qquad z = \\frac{y-\\mu}{\\sigma}",
  moments = "\\mathbb{E}[Y] = \\mu\\ (\\nu>1), \\qquad \\operatorname{Var}(Y) = \\sigma^{2}\\nu/(\\nu-2)\\ (\\nu>2)",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    lgamma((th$nu + 1) / 2) - lgamma(th$nu / 2) - 0.5 * log(th$nu * pi) -
      log(th$sigma) - ((th$nu + 1) / 2) * log1p(z^2 / th$nu)
  },
  grid = function(th) seq(-7, 8, length.out = 40),
  note = "A location-scale family built on the standard $t$: $\\sigma$ is a scale,
  **not** the standard deviation, which is $\\sigma\\sqrt{\\nu/(\\nu-2)}$ and
  exists only for $\\nu > 2$. The degrees of freedom $\\nu$ are continuous and
  estimable."
)

DISTRIBS$laplace <- list(
  title = "Laplace",
  ctor = "laplace_distrib()",
  obj = function() laplace_distrib(),
  theta = list(mu = 0.5, b = 1.4),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{2b}\\exp\\!\\left\\{-\\frac{|y-\\mu|}{b}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = 2b^{2}",
  ld = function(y, th) -log(2 * th$b) - abs(y - th$mu) / th$b,
  grid = function(th) seq(-8, 9, length.out = 40),
  note = "The one non-regular member of the catalogue: $\\ell$ has a kink at
  $y = \\mu$, the observed information in $\\mu$ is identically zero, and the
  Fisher information $1/b^{2}$ can only be recovered from the score. See
  @sec-nonregular. `params_smooth` records this, and it is the reason that
  property exists at all."
)

DISTRIBS$pseudohuber <- list(
  title = "Pseudo-Huber",
  ctor = "pseudohuber_distrib()",
  obj = function() pseudohuber_distrib(),
  theta = list(mu = 0.5, sigma = 1.2, nu = 3),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{\\exp(-D)}{2\\sigma\\sqrt{\\nu}\\,K_{1}(\\sqrt{\\nu})}, \\qquad D = \\sqrt{\\nu + \\left(\\frac{y-\\mu}{\\sigma}\\right)^{2}}",
  moments = "\\mathbb{E}[Y] = \\mu; variance has no elementary closed form",
  ld = function(y, th) {
    D <- sqrt(th$nu + ((y - th$mu) / th$sigma)^2)
    sn <- sqrt(th$nu)
    # log K_1(x) via the exponentially scaled Bessel function: K_1(x) = e^{-x} * scaled
    lognorm <- log(2) + log(th$sigma) + 0.5 * log(th$nu) +
      log(besselK(sn, 1, expon.scaled = TRUE)) - sn
    -D - lognorm
  },
  grid = function(th) seq(-8, 9, length.out = 40),
  note = "The normalising constant is exact, not approximate:
  $\\int_{-\\infty}^{\\infty} e^{-\\sqrt{a^{2}+z^{2}}}\\,dz = 2aK_{1}(a)$, so with
  $a = \\sqrt{\\nu}$ and $dy = \\sigma\\,dz$ the constant is
  $2\\sigma\\sqrt{\\nu}K_{1}(\\sqrt{\\nu})$. The Bessel terms are
  degree-homogeneous, which is why the **exponentially scaled**
  `besselK(x, nu, expon.scaled = TRUE)` may be used: the scaling cancels between
  numerator and denominator, so it is exact rather than approximate, and it
  avoids overflow out to $\\nu = 2000$. $\\nu \\to \\infty$ gives a Gaussian limit
  and small $\\nu$ a Laplace-like one, so the family interpolates between the two
  --- with, unlike the Laplace, a genuinely smooth log-density everywhere."
)

DISTRIBS$gamma <- list(
  title = "Gamma",
  ctor = "gamma_distrib()",
  obj = function() gamma_distrib(),
  theta = list(mu = 3, sigma2 = 2),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{\\beta^{\\alpha}}{\\Gamma(\\alpha)}y^{\\alpha-1}e^{-\\beta y}, \\qquad \\alpha = \\frac{\\mu^{2}}{\\sigma^{2}}, \\quad \\beta = \\frac{\\mu}{\\sigma^{2}}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}",
  ld = function(y, th) {
    a <- th$mu^2 / th$sigma2
    b <- th$mu / th$sigma2
    a * log(b) - lgamma(a) + (a - 1) * log(y) - b * y
  },
  grid = function(th) seq(0.05, 15, length.out = 40),
  note = "**Not** the base-R parametrisation. `dgamma()` takes a shape and a rate;
  this object takes the **mean and the variance** directly, which is what a
  modelling framework wants, since it is $\\mu$ that gets a linear predictor. The
  map is $\\alpha = \\mu^{2}/\\sigma^{2}$, $\\beta = \\mu/\\sigma^{2}$, whence
  $\\alpha/\\beta = \\mu$ and $\\alpha/\\beta^{2} = \\sigma^{2}$ as claimed."
)

DISTRIBS$invgauss <- list(
  title = "Inverse Gaussian",
  ctor = "invgauss_distrib()",
  obj = function() invgauss_distrib(),
  theta = list(mu = 2, phi = 0.7),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{\\sqrt{2\\pi\\phi y^{3}}}\\exp\\!\\left\\{-\\frac{(y-\\mu)^{2}}{2\\phi\\mu^{2}y}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\phi\\mu^{3}",
  ld = function(y, th) {
    -0.5 * log(2 * pi * th$phi * y^3) - (y - th$mu)^2 / (2 * th$phi * th$mu^2 * y)
  },
  grid = function(th) seq(0.05, 12, length.out = 40),
  note = "In the **dispersion** parametrisation, $\\phi = 1/\\lambda$ where
  $\\lambda$ is the shape used by some references. The variance function
  $\\phi\\mu^{3}$ is what makes the inverse-square link canonical here (see
  @sec-link-inverse-sq)."
)

DISTRIBS$lognormal <- list(
  title = "Lognormal",
  ctor = "lognormal_distrib()",
  obj = function() lognormal_distrib(),
  theta = list(mu = 0.5, sigma2 = 1.3),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{y\\sqrt{2\\pi\\sigma^{2}}}\\exp\\!\\left\\{-\\frac{(\\log y - \\mu)^{2}}{2\\sigma^{2}}\\right\\}",
  moments = "\\mathbb{E}[Y] = e^{\\mu + \\sigma^{2}/2}, \\qquad \\operatorname{Var}(Y) = \\left(e^{\\sigma^{2}}-1\\right)e^{2\\mu+\\sigma^{2}}",
  ld = function(y, th) {
    -log(y) - 0.5 * log(2 * pi * th$sigma2) - (log(y) - th$mu)^2 / (2 * th$sigma2)
  },
  grid = function(th) seq(0.02, 20, length.out = 40),
  note = "$\\mu$ and $\\sigma^{2}$ are the mean and **variance on the log scale**,
  not of $Y$ itself: `dlnorm()` takes `sdlog`, and this object takes its square.
  The same law can also be built as
  `transformation(gaussian_distrib(), exp_transform())`, by the change of
  variables of @sec-transformation, and the two constructions coincide."
)

DISTRIBS$beta <- list(
  title = "Beta",
  ctor = "beta_distrib()",
  obj = function() beta_distrib(),
  theta = list(mu = 0.4, phi = 6),
  support = "y \\in (0, 1)",
  pdf_latex = "f(y) = \\frac{y^{a-1}(1-y)^{b-1}}{B(a,b)}, \\qquad a = \\mu\\phi, \\quad b = (1-\\mu)\\phi",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\dfrac{\\mu(1-\\mu)}{1+\\phi}",
  ld = function(y, th) {
    a <- th$mu * th$phi
    b <- (1 - th$mu) * th$phi
    (a - 1) * log(y) + (b - 1) * log1p(-y) - (lgamma(a) + lgamma(b) - lgamma(a + b))
  },
  grid = function(th) seq(0.02, 0.98, length.out = 40),
  note = "The **mean-precision** parametrisation of beta regression, not the
  `shape1`/`shape2` one: $a = \\mu\\phi$, $b = (1-\\mu)\\phi$, so
  $a/(a+b) = \\mu$ and $a+b = \\phi$. Larger $\\phi$ means less dispersion."
)

DISTRIBS$bernoulli <- list(
  title = "Bernoulli",
  ctor = "bernoulli_distrib()",
  obj = function() bernoulli_distrib(),
  theta = list(mu = 0.35),
  support = "y \\in \\{0, 1\\}",
  pdf_latex = "P(Y=y) = \\mu^{y}(1-\\mu)^{1-y}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu(1-\\mu)",
  ld = function(y, th) y * log(th$mu) + (1 - y) * log1p(-th$mu),
  grid = function(th) 0:1,
  note = "Its support has two points, which is the whole content of
  @sec-identifiability: a Bernoulli has exactly one free probability, so neither
  zero-inflating nor zero-adjusting it produces an identified model."
)

DISTRIBS$binomial <- list(
  title = "Binomial",
  ctor = "binomial_distrib(size = 10)",
  obj = function() binomial_distrib(size = 10),
  theta = list(mu = 0.35),
  support = "y \\in \\{0, 1, \\dots, n\\}",
  pdf_latex = "P(Y=y) = \\binom{n}{y}\\mu^{y}(1-\\mu)^{n-y}",
  moments = "\\mathbb{E}[Y] = n\\mu, \\qquad \\operatorname{Var}(Y) = n\\mu(1-\\mu)",
  ld = function(y, th) {
    n <- 10
    lchoose(n, y) + y * log(th$mu) + (n - y) * log1p(-th$mu)
  },
  grid = function(th) 0:10,
  note = "The number of trials $n$ is a **property of the object**, not a
  parameter: it is known, not estimated, so it belongs with `size` rather than in
  `theta`. Shown at $n = 10$."
)

DISTRIBS$poisson <- list(
  title = "Poisson",
  ctor = "poisson_distrib()",
  obj = function() poisson_distrib(),
  theta = list(mu = 4),
  support = "y \\in \\{0, 1, 2, \\dots\\}",
  pdf_latex = "P(Y=y) = \\frac{e^{-\\mu}\\mu^{y}}{y!}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu",
  ld = function(y, th) -th$mu + y * log(th$mu) - lgamma(y + 1),
  grid = function(th) 0:20,
  note = "Equidispersed by construction, which is what makes it the natural parent
  for the zero-inflation and hurdle wrappers of @sec-transformations: excess
  zeros are one of the two standard ways real count data depart from it."
)

DISTRIBS$negbin <- list(
  title = "Negative binomial (NB2)",
  ctor = "negbin_distrib()",
  obj = function() negbin_distrib(),
  theta = list(mu = 4, theta = 1.7),
  support = "y \\in \\{0, 1, 2, \\dots\\}",
  pdf_latex = "P(Y=y) = \\frac{\\Gamma(y+\\theta)}{\\Gamma(\\theta)\\,y!}\\left(\\frac{\\theta}{\\theta+\\mu}\\right)^{\\theta}\\left(\\frac{\\mu}{\\theta+\\mu}\\right)^{y}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu + \\dfrac{\\mu^{2}}{\\theta}",
  ld = function(y, th) {
    lgamma(y + th$theta) - lgamma(th$theta) - lgamma(y + 1) +
      th$theta * (log(th$theta) - log(th$theta + th$mu)) +
      y * (log(th$mu) - log(th$theta + th$mu))
  },
  grid = function(th) 0:40,
  note = "The NB2 parametrisation: $\\theta$ is the size, and dispersion
  **decreases** as $\\theta$ grows, with the Poisson recovered as
  $\\theta \\to \\infty$. The excess kurtosis is
  $6/\\theta + \\theta/(\\mu(\\theta+\\mu))$. That expression is worth stating
  explicitly because the predecessor package this one replaces had
  $6/\\theta + (\\theta+\\mu)/(\\mu\\theta)$, which is wrong; it was found by
  numerical comparison, not by reading. Every formula inherited from that package
  was subsequently treated as unverified until checked."
)


# ---------------------------------------------------------------------------
# Catalogue rendering
# ---------------------------------------------------------------------------

# hess_pairs() is internal to distributions7; the book needs the same map, so it
# is recomputed here rather than reached for across the namespace boundary.
book_hess_pair <- function(params, nm) {
  hits <- which(vapply(params, function(p) startsWith(nm, paste0(p, "_")), logical(1)))
  first <- params[hits[which.max(nchar(params[hits]))]]
  c(first, substring(nm, nchar(first) + 2L))
}

fmt_bounds <- function(b, discrete) {
  f <- function(x) if (is.infinite(x)) (if (x > 0) "\\infty" else "-\\infty") else format(x)
  if (discrete) {
    if (is.infinite(b[2])) sprintf("$\\{%s, %s, \\dots\\}$", format(b[1]), format(b[1] + 1))
    else sprintf("$\\{%s, \\dots, %s\\}$", format(b[1]), format(b[2]))
  } else {
    sprintf("$(%s, %s)$", f(b[1]), f(b[2]))
  }
}

catalogue_table <- function() {
  rows <- lapply(names(DISTRIBS), function(id) {
    d <- DISTRIBS[[id]]$obj()
    is_disc <- S7::S7_inherits(d, discrete_distrib)
    data.frame(
      distribution = DISTRIBS[[id]]$title,
      kind = if (is_disc) "discrete" else "continuous",
      support = fmt_bounds(d@bounds, is_disc),
      parameters = paste(sprintf("`%s`", d@params), collapse = ", "),
      links = paste(vapply(d@params, function(p) d@link_params[[p]]@link_name, character(1)),
                    collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# The formula-vs-implementation certificate for one distribution.
certify_distrib_density <- function(rec) {
  d <- rec$obj()
  y <- rec$grid(rec$theta)
  book <- rec$ld(y, rec$theta)
  pkg <- distrib_pdf(d, y, rec$theta, log = TRUE)
  rel_err(book, pkg)
}

render_distrib_entry <- function(id, rec) {
  cat(sprintf("\n### %s {#sec-dist-%s}\n\n", rec$title, gsub("_", "-", id)))
  cat(sprintf("``` r\n%s\n```\n\n", rec$ctor))
  cat(sprintf("Support: $%s$.\n\n", rec$support))
  cat(sprintf("$$%s$$\n\n", rec$pdf_latex))
  cat(sprintf("$$%s$$\n\n", rec$moments))
  cat(gsub("\n\\s+", "\n", rec$note), "\n\n")
}

# Silent consistency gate for Chapter 3: every printed density against the
# implementation, the Bartlett identities, the link-scale machinery, and the
# package's own validator on all fourteen distributions. Stops the render on
# any disagreement; never rendered.
assert_distributions_ok <- function() {
  for (id in names(DISTRIBS)) {
    err <- certify_distrib_density(DISTRIBS[[id]])
    if (is.na(err) || err > 1e-10) {
      stop(sprintf("Distribution '%s': the printed density disagrees with distrib_pdf().", id),
           call. = FALSE)
    }
  }
  bt <- bartlett_certificate()
  if (any(bt$verdict != "holds")) {
    stop("A Bartlett identity fails: ", paste(bt$identity[bt$verdict != "holds"], collapse = "; "),
         call. = FALSE)
  }
  ls <- link_scale_certificate()
  if (any(ls$verdict != "agree")) {
    stop("Link-scale derivatives disagree with finite differences in eta.", call. = FALSE)
  }
  dc <- distrib_certificate()
  if (!dc$all_ok) {
    stop("check_distrib() failed on: ",
         paste(dc$table$distribution[dc$table$checks != dc$table$passed], collapse = ", "),
         call. = FALSE)
  }

  # The non-regular section makes three specific numerical claims; pin them.
  # (a) the shipped Laplace has a closed-form E[H], so approx is ignored and all
  #     strategies return -1/b^2;
  # (b) on a bare Laplace without closed forms, bartlett recovers -1/b^2 while
  #     integrate and mc, which both average the a.e.-zero observed l_mumu,
  #     return zero.
  lap <- laplace_distrib()
  th <- list(mu = 0, b = 2)
  shipped <- vapply(c("bartlett", "integrate", "mc"), function(a) {
    distrib_expected_hessian(lap, 0, th, approx = a)$mu_mu
  }, numeric(1))
  if (max(abs(shipped + 0.25)) > 1e-12) {
    stop("The shipped Laplace no longer ignores 'approx'; the non-regular section is stale.",
         call. = FALSE)
  }
  BareLap <- S7::new_class("BareLapGate", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, BareLap) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  S7::method(distrib_gradient, BareLap) <- function(distrib, y, theta,
                                                    scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]; b <- theta[[2]]
    list(mu = sign(r) / b, b = (abs(r) / b - 1) / b)
  }
  S7::method(distrib_hessian, BareLap) <- function(distrib, y, theta,
                                                   scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]; b <- theta[[2]]; n <- length(y)
    list(mu_mu = rep(0, n), b_b = (b - 2 * abs(r)) / b^3, mu_b = -sign(r) / b^2)
  }
  bare <- BareLap(
    distrib_name = "bare laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "location", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(mu = identity_link(), b = log_link()),
    params_smooth = c(mu = FALSE, b = TRUE)
  )
  set.seed(1)
  strat <- c(
    bartlett  = distrib_expected_hessian(bare, 0, th, approx = "bartlett")$mu_mu,
    integrate = distrib_expected_hessian(bare, 0, th, approx = "integrate")$mu_mu,
    mc        = distrib_expected_hessian(bare, 0, th, approx = "mc")$mu_mu
  )
  if (abs(strat[["bartlett"]] + 0.25) > 1e-6 ||
      abs(strat[["integrate"]]) > 1e-8 || abs(strat[["mc"]]) > 1e-8) {
    stop("The bare-Laplace strategy comparison no longer matches the prose of the ",
         "non-regular section: ", paste(names(strat), round(strat, 4), collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Certificate: the Bartlett identities
#
# Both sides computed here, from the observed derivatives, by routes that share
# nothing but those derivatives: the left-hand side is what the package returns
# for the expected derivative, the right-hand side is the partition sum of
# expectations of products, assembled in this file.
# ---------------------------------------------------------------------------

bartlett_certificate <- function(ids = c("gaussian", "gamma", "beta", "poisson", "negbin")) {
  rows <- list()
  for (id in ids) {
    rec <- DISTRIBS[[id]]
    d <- rec$obj()
    th <- rec$theta
    params <- d@params

    # k = 1: E[l_i] = 0
    e1 <- max(vapply(params, function(p) {
      abs(expectation(d, function(y, theta) distrib_gradient(d, y, theta)[[p]], th))
    }, numeric(1)))

    # k = 2: E[l_ij] = -E[l_i l_j]
    e2 <- max(vapply(hess_names(params), function(nm) {
      pr <- book_hess_pair(params, nm)
      lhs <- distrib_expected_hessian(d, 0, th)[[nm]][1]
      rhs <- -expectation(d, function(y, theta) {
        g <- distrib_gradient(d, y, theta); g[[pr[1]]] * g[[pr[2]]]
      }, th)
      abs(lhs - rhs) / max(1, abs(rhs))
    }, numeric(1)))

    # k = 3: E[l_ijk] = -(E[l_ij l_k] + E[l_ik l_j] + E[l_jk l_i] + E[l_i l_j l_k])
    e3 <- max(vapply(deriv_names(params, 3), function(nm) {
      idx <- strsplit(nm, "_", fixed = TRUE)[[1]]
      i <- idx[1]; j <- idx[2]; k <- idx[3]
      key <- function(a, b) {
        cand <- c(paste(a, b, sep = "_"), paste(b, a, sep = "_"))
        cand[cand %in% hess_names(params)][1]
      }
      lhs <- distrib_deriv3(d, 0, th, expected = TRUE, approx = "integrate")[[nm]][1]
      rhs <- -expectation(d, function(y, theta) {
        g <- distrib_gradient(d, y, theta)
        H <- distrib_hessian(d, y, theta)
        H[[key(i, j)]] * g[[k]] + H[[key(i, k)]] * g[[j]] + H[[key(j, k)]] * g[[i]] +
          g[[i]] * g[[j]] * g[[k]]
      }, th)
      abs(lhs - rhs) / max(1, abs(rhs))
    }, numeric(1)))

    rows[[length(rows) + 1L]] <- data.frame(
      distribution = rec$title,
      identity = c("$k=1$: $\\mathbb{E}[\\ell_i]=0$",
                   "$k=2$: information matrix equality",
                   "$k=3$: partition sum"),
      error = c(e1, e2, e3),
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  data.frame(
    distribution = tab$distribution,
    identity = tab$identity,
    `max deviation` = fmt_err(tab$error),
    verdict = ifelse(tab$error < 1e-4, "holds", "VIOLATED"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Certificate: link-scale derivatives
#
# Orders 1..4 of the diagonal link-scale derivative, each against ONE numerical
# differentiation of the analytical order below. Order 1 is compared against a
# numerical derivative of the log-density itself, so the chain is anchored on
# something the link-scale machinery had no part in producing.
# ---------------------------------------------------------------------------

link_scale_certificate <- function(ids = c("gaussian", "gamma", "beta", "negbin"),
                                   y = c(0.5, 1.2)) {
  rows <- list()
  for (id in ids) {
    rec <- DISTRIBS[[id]]
    d <- rec$obj()
    th <- rec$theta
    params <- d@params
    yy <- rec$grid(th)[c(5, 12)]

    eta0 <- vapply(seq_along(params), function(i) {
      linkfun(d@link_params[[params[i]]], th[[i]])
    }, numeric(1))
    theta_of <- function(e) {
      out <- lapply(seq_along(params), function(i) linkinv(d@link_params[[params[i]]], e[i]))
      names(out) <- params
      out
    }
    # analytic diagonal link-scale derivative of order k, in parameter i
    ana <- function(i, k, e) {
      th_e <- theta_of(e)
      nm <- paste(rep(params[i], k), collapse = "_")
      switch(k,
        distrib_gradient(d, yy, th_e, scale = "link")[[params[i]]],
        distrib_hessian(d, yy, th_e, scale = "link")[[nm]],
        distrib_deriv3(d, yy, th_e, scale = "link")[[nm]],
        distrib_deriv4(d, yy, th_e, scale = "link")[[nm]])
    }
    worst <- numeric(4)
    for (k in 1:4) {
      errs <- vapply(seq_along(params), function(i) {
        num <- vapply(seq_along(yy), function(m) {
          f <- function(v) {
            e <- eta0; e[i] <- v
            if (k == 1L) distrib_pdf(d, yy[m], theta_of(e), log = TRUE) else ana(i, k - 1L, e)[m]
          }
          numDeriv::grad(f, eta0[i])
        }, numeric(1))
        rel_err(ana(i, k, eta0), num)
      }, numeric(1))
      worst[k] <- max(errs)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      distribution = rec$title,
      order = paste0("$\\partial^{", 1:4, "}\\ell/\\partial\\eta^{", 1:4, "}$"),
      error = worst, stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  data.frame(
    distribution = tab$distribution,
    order = tab$order,
    `max relative error` = fmt_err(tab$error),
    verdict = ifelse(tab$error < 1e-4, "agree", "DISAGREE"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Certificate: check_distrib() on the whole catalogue
# ---------------------------------------------------------------------------

distrib_certificate <- function(n = 40, nsim = 5e4, orders = 1:4, seed = 20260726) {
  set.seed(seed)
  rows <- lapply(names(DISTRIBS), function(id) {
    rec <- DISTRIBS[[id]]
    out <- check_distrib(rec$obj(), theta = rec$theta, n = n, nsim = nsim,
                         orders = orders, verbose = FALSE)
    failed <- out$check[out$status != "OK"]
    data.frame(
      distribution = rec$title,
      checks = nrow(out),
      passed = sum(out$status == "OK"),
      failures = if (length(failed)) paste(failed, collapse = "; ") else "none",
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  list(table = tab,
       n_checks = sum(tab$checks),
       n_ok = sum(tab$passed),
       all_ok = all(tab$checks == tab$passed))
}
