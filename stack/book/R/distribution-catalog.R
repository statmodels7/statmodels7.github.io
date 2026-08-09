# The distribution catalog, and the certificates for Chapter 3.
#
# As in R/link-formulas.R, each record keeps the printed formula and an
# independent R transcription of it side by side. The transcription is written
# from the density as displayed in the book, using only lgamma/lbeta/besselK and
# arithmetic -- never the `d*` function the package itself calls -- so that the
# comparison actually tests the parametrization rather than restating it.

DISTRIBS <- list()

DISTRIBS$gaussian <- list(
  title = "Gaussian, mean and standard deviation",
  ctor = "gaussian1_distrib()",
  obj = function() gaussian1_distrib(),
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
  moments = "",
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
  title = "Student t, scale",
  ctor = "student_t1_distrib()",
  obj = function() student_t1_distrib(),
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
  theta = list(mu = 0.5, sigma = 1.4),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{2\\sigma}\\exp\\!\\left\\{-\\frac{|y-\\mu|}{\\sigma}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = 2\\sigma^{2}",
  ld = function(y, th) -log(2 * th$sigma) - abs(y - th$mu) / th$sigma,
  grid = function(th) seq(-8, 9, length.out = 40),
  note = "The one non-regular member of the catalog: $\\ell$ has a kink at
  $y = \\mu$, the observed information in $\\mu$ is identically zero, and the
  Fisher information $1/\\sigma^{2}$ can only be recovered from the score. See
  @sec-nonregular. `params_smooth` records this, and it is the reason that
  property exists at all. `laplace2` carries the same law in location and
  rate."
)

DISTRIBS$laplace2 <- list(
  title = "Laplace, location and rate",
  ctor = "laplace2_distrib()",
  obj = function() laplace2_distrib(),
  theta = list(mu = 0.5, lambda = 0.7),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{\\lambda}{2}\\exp\\!\\left\\{-\\lambda|y-\\mu|\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = 2/\\lambda^{2}",
  ld = function(y, th) log(th$lambda / 2) - th$lambda * abs(y - th$mu),
  grid = function(th) seq(-8, 9, length.out = 40),
  note = "`laplace` in the rate $\\lambda = 1/\\sigma$. This is the form a
  penalty consumes: at fixed $\\mu = 0$ the negative log-density is
  $\\lambda\\lvert y \\rvert$ up to a constant, so the lasso penalty is
  linear in $\\lambda$ and every derivative in $\\lambda$ beyond the first
  is free of the data. The kink at $y = \\mu$ and its handling are those of
  `laplace`."
)

DISTRIBS$enet <- list(
  title = "Elastic net",
  ctor = "enet_distrib()",
  obj = function() enet_distrib(),
  theta = list(mu = 0.3, lambda = 1.5, alpha = 0.6),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{Z}\\exp\\!\\left\\{-a|y-\\mu| - \\frac{c}{2}(y-\\mu)^{2}\\right\\}, \\qquad a = \\lambda\\alpha, \\quad c = \\lambda(1-\\alpha), \\quad Z = \\frac{2M(a/\\sqrt{c})}{\\sqrt{c}}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\frac{1 + xG}{c}, \\qquad x = \\frac{a}{\\sqrt{c}}, \\quad G = \\frac{\\mathrm{d}\\log M}{\\mathrm{d}x}",
  ld = function(y, th) {
    a <- th$lambda * th$alpha
    cc <- th$lambda * (1 - th$alpha)
    x <- a / sqrt(cc)
    # M is the Mills ratio Phi(-x)/phi(x), written here from pnorm and
    # dnorm rather than through the package's own stable route
    logM <- stats::pnorm(-x, log.p = TRUE) - stats::dnorm(x, log = TRUE)
    logZ <- log(2) - log(cc) / 2 + logM
    -a * abs(y - th$mu) - cc * (y - th$mu)^2 / 2 - logZ
  },
  grid = function(th) seq(-6, 7, length.out = 40),
  note = "The product of a Laplace and a Gaussian at the same location,
  normalized, and so the density whose negative logarithm is the
  elastic-net penalty. $M$ is the Mills ratio $\\Phi(-x)/\\varphi(x)$, which
  makes the constant closed and finite at both ends: $Z \\to 2/a$ as
  $\\alpha \\to 1$, where the family is `laplace2`, and
  $Z \\to \\sqrt{2\\pi/c}$ as $\\alpha \\to 0$, where it is the Gaussian. Both
  of those remain families of their own and $\\alpha$ is confined to the
  open interval. Every derivative in $a$ and $c$ is a polynomial in $x$
  and $G$, with $G' = 1 + xG - G^{2}$, and the chain to
  $(\\lambda, \\alpha)$ is bilinear, so the mixed second derivative picks up
  the map's own cross term. The kink at $y = \\mu$ is the Laplace's."
)

DISTRIBS$pseudohuber <- list(
  title = "Pseudo-Huber",
  ctor = "pseudohuber_distrib()",
  obj = function() pseudohuber_distrib(),
  theta = list(mu = 0.5, sigma = 1.2, nu = 3),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{\\exp(-D)}{2\\sigma\\sqrt{\\nu}\\,K_{1}(\\sqrt{\\nu})}, \\qquad D = \\sqrt{\\nu + \\left(\\frac{y-\\mu}{\\sigma}\\right)^{2}}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}\\sqrt{\\nu}\\,\\frac{K_{2}(\\sqrt{\\nu})}{K_{1}(\\sqrt{\\nu})}",
  ld = function(y, th) {
    D <- sqrt(th$nu + ((y - th$mu) / th$sigma)^2)
    sn <- sqrt(th$nu)
    # log K_1(x) via the exponentially scaled Bessel function: K_1(x) = e^{-x} * scaled
    lognorm <- log(2) + log(th$sigma) + 0.5 * log(th$nu) +
      log(besselK(sn, 1, expon.scaled = TRUE)) - sn
    -D - lognorm
  },
  grid = function(th) seq(-8, 9, length.out = 40),
  note = "The normalizing constant is exact, not approximate:
  $\\int_{-\\infty}^{\\infty} e^{-\\sqrt{a^{2}+z^{2}}}\\,dz = 2aK_{1}(a)$, so with
  $a = \\sqrt{\\nu}$ and $dy = \\sigma\\,dz$ the constant is
  $2\\sigma\\sqrt{\\nu}K_{1}(\\sqrt{\\nu})$. The Bessel terms are
  degree-homogeneous, which is why the **exponentially scaled**
  `besselK(x, nu, expon.scaled = TRUE)` may be used: the scaling cancels between
  numerator and denominator, so it is exact rather than approximate, and it
  avoids overflow out to $\\nu = 2000$. Differentiating the same identity in
  $a$ gives the variance printed above, computed through the same scaled
  ratio. $\\nu \\to \\infty$ gives a Gaussian limit
  and small $\\nu$ a Laplace-like one, so the family interpolates between the two
  --- with, unlike the Laplace, a genuinely smooth log-density everywhere."
)

DISTRIBS$gamma <- list(
  title = "Gamma, mean and variance",
  ctor = "gamma2_distrib()",
  obj = function() gamma2_distrib(),
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
  note = "**Not** the base-R parametrization. `dgamma()` takes a shape and a rate;
  this object takes the **mean and the variance** directly, which is what a
  modeling framework wants, since it is $\\mu$ that gets a linear predictor. The
  map is $\\alpha = \\mu^{2}/\\sigma^{2}$, $\\beta = \\mu/\\sigma^{2}$, whence
  $\\alpha/\\beta = \\mu$ and $\\alpha/\\beta^{2} = \\sigma^{2}$ as claimed."
)

DISTRIBS$invgauss <- list(
  title = "Inverse Gaussian, mean and dispersion",
  ctor = "invgauss1_distrib()",
  obj = function() invgauss1_distrib(),
  theta = list(mu = 2, phi = 0.7),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{\\sqrt{2\\pi\\phi y^{3}}}\\exp\\!\\left\\{-\\frac{(y-\\mu)^{2}}{2\\phi\\mu^{2}y}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\phi\\mu^{3}",
  ld = function(y, th) {
    -0.5 * log(2 * pi * th$phi * y^3) - (y - th$mu)^2 / (2 * th$phi * th$mu^2 * y)
  },
  grid = function(th) seq(0.05, 12, length.out = 40),
  note = "In the **dispersion** parametrization, $\\phi = 1/\\lambda$ where
  $\\lambda$ is the shape used by some references. The variance function
  $\\phi\\mu^{3}$ is what makes the inverse-square link canonical here (see
  @sec-link-inverse-sq)."
)

DISTRIBS$lognormal <- list(
  title = "Lognormal, log-scale parameters",
  ctor = "lognormal1_distrib()",
  obj = function() lognormal1_distrib(),
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
  `transformation(gaussian1_distrib(), exp_transform())`, by the change of
  variables of @sec-transformation, and the two constructions coincide."
)

DISTRIBS$beta <- list(
  title = "Beta, mean and precision",
  ctor = "beta1_distrib()",
  obj = function() beta1_distrib(),
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
  note = "The **mean-precision** parametrization of beta regression, not the
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
  title = "Negative binomial, NB2",
  ctor = "negbin2_distrib()",
  obj = function() negbin2_distrib(),
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
  note = "The NB2 parametrization: $\\theta$ is the size, and dispersion
  **decreases** as $\\theta$ grows, with the Poisson recovered as
  $\\theta \\to \\infty$. The excess kurtosis is
  $6/\\theta + \\theta/(\\mu(\\theta+\\mu))$. That expression is worth stating
  explicitly because the predecessor package this one replaces had
  $6/\\theta + (\\theta+\\mu)/(\\mu\\theta)$, which is wrong; it was found by
  numerical comparison, not by reading. Every formula inherited from that package
  was subsequently treated as unverified until checked."
)


# --- the remaining univariate families --------------------------------------

DISTRIBS$gaussian2 <- list(
  title = "Gaussian, mean and variance",
  ctor = "gaussian2_distrib()",
  obj = function() gaussian2_distrib(),
  theta = list(mu = 1.5, sigma2 = 4),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{\\sqrt{2\\pi\\sigma^{2}}}\\exp\\!\\left\\{-\\frac{(y-\\mu)^{2}}{2\\sigma^{2}}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}",
  ld = function(y, th) -0.5 * log(2 * pi * th$sigma2) - (y - th$mu)^2 / (2 * th$sigma2),
  grid = function(th) seq(-4, 7, length.out = 40),
  note = "The variance parametrization of the Gaussian: a linear predictor on
  the second parameter models $\\sigma^{2}$, where `gaussian1` models
  $\\sigma$."
)

DISTRIBS$gaussian3 <- list(
  title = "Gaussian, mean and precision",
  ctor = "gaussian3_distrib()",
  obj = function() gaussian3_distrib(),
  theta = list(mu = 1.5, tau = 0.25),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\sqrt{\\frac{\\tau}{2\\pi}}\\exp\\!\\left\\{-\\frac{\\tau(y-\\mu)^{2}}{2}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = 1/\\tau",
  ld = function(y, th) 0.5 * log(th$tau) - 0.5 * log(2 * pi) - th$tau * (y - th$mu)^2 / 2,
  grid = function(th) seq(-4, 7, length.out = 40),
  note = "The precision parametrization, the conjugate scale of Bayesian
  work."
)

DISTRIBS$gamma1 <- list(
  title = "Gamma, mean and dispersion",
  ctor = "gamma1_distrib()",
  obj = function() gamma1_distrib(),
  theta = list(mu = 3, phi = 0.4),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{y^{1/\\varphi - 1}e^{-y/(\\varphi\\mu)}}{(\\varphi\\mu)^{1/\\varphi}\\,\\Gamma(1/\\varphi)}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\varphi\\mu^{2}",
  ld = function(y, th) {
    k <- 1 / th$phi; r <- 1 / (th$phi * th$mu)
    k * log(r) + (k - 1) * log(y) - r * y - lgamma(k)
  },
  grid = function(th) seq(0.2, 12, length.out = 40),
  note = "The dispersion form: the squared coefficient of variation is
  $\\varphi$, constant in the mean, so shape and rate are
  $1/\\varphi$ and $1/(\\varphi\\mu)$. `gamma2` carries the variance
  instead."
)

DISTRIBS$negbin1 <- list(
  title = "Negative binomial, NB1",
  ctor = "negbin1_distrib()",
  obj = function() negbin1_distrib(),
  theta = list(mu = 4, theta = 1.5),
  support = "y \\in \\{0, 1, \\dots\\}",
  pdf_latex = "P(Y=y) = \\frac{\\Gamma(y + \\mu/\\theta)}{\\Gamma(\\mu/\\theta)\\,y!}\\left(\\frac{1}{1+\\theta}\\right)^{\\mu/\\theta}\\left(\\frac{\\theta}{1+\\theta}\\right)^{y}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu(1+\\theta)",
  ld = function(y, th) {
    r <- th$mu / th$theta
    lgamma(y + r) - lgamma(r) - lfactorial(y) - r * log1p(th$theta) +
      y * (log(th$theta) - log1p(th$theta))
  },
  grid = function(th) 0:15,
  note = "The NB1 of Cameron and Trivedi: the variance is linear in the mean,
  where `negbin2` makes it quadratic. The expected information is a series
  against the exact mass."
)

DISTRIBS$weibull1 <- list(
  title = "Weibull, scale and shape",
  ctor = "weibull1_distrib()",
  obj = function() weibull1_distrib(),
  theta = list(mu = 3, sigma = 1.7),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{\\sigma}{\\mu}\\left(\\frac{y}{\\mu}\\right)^{\\sigma-1}\\exp\\!\\left\\{-\\left(\\frac{y}{\\mu}\\right)^{\\sigma}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu\\,\\Gamma(1+1/\\sigma), \\qquad \\operatorname{Var}(Y) = \\mu^{2}\\left[\\Gamma(1+2/\\sigma) - \\Gamma(1+1/\\sigma)^{2}\\right]",
  ld = function(y, th) log(th$sigma) - log(th$mu) +
    (th$sigma - 1) * (log(y) - log(th$mu)) - (y / th$mu)^th$sigma,
  grid = function(th) seq(0.2, 9, length.out = 40),
  note = "The `WEI` of gamlss: $\\mu$ is the scale and not the mean. With
  $u = (Y/\\mu)^{\\sigma}$ standard exponential, every expectation the
  information needs is a derivative of $\\Gamma$ at 2, so the expected
  third and fourth orders are closed as well."
)

DISTRIBS$weibull3 <- list(
  title = "Weibull, mean and shape",
  ctor = "weibull3_distrib()",
  obj = function() weibull3_distrib(),
  theta = list(mean = 3, sigma = 1.7),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{\\sigma}{b}\\left(\\frac{y}{b}\\right)^{\\sigma-1}e^{-(y/b)^{\\sigma}}, \\qquad b = \\frac{\\mu}{\\Gamma(1+1/\\sigma)}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = b^{2}\\left[\\Gamma(1+2/\\sigma) - \\Gamma(1+1/\\sigma)^{2}\\right]",
  ld = function(y, th) {
    b <- th$mean / gamma(1 + 1 / th$sigma)
    log(th$sigma) - log(b) + (th$sigma - 1) * (log(y) - log(b)) - (y / b)^th$sigma
  },
  grid = function(th) seq(0.2, 9, length.out = 40),
  note = "The `WEI3` of gamlss, `weibull1` reparametrized so that a linear
  predictor acts on the mean; `weibull2` is deliberately absent, following
  the gamlss numbering."
)

DISTRIBS$gumbel <- list(
  title = "Gumbel",
  ctor = "gumbel_distrib()",
  obj = function() gumbel_distrib(),
  theta = list(mu = 1, sigma = 2),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{\\sigma}\\exp\\!\\left\\{-z - e^{-z}\\right\\}, \\qquad z = \\frac{y-\\mu}{\\sigma}",
  moments = "\\mathbb{E}[Y] = \\mu + \\gamma\\sigma, \\qquad \\operatorname{Var}(Y) = \\pi^{2}\\sigma^{2}/6",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    -log(th$sigma) - z - exp(-z)
  },
  grid = function(th) seq(-4, 12, length.out = 40),
  note = "The extreme-value family for maxima, with fixed skewness
  $12\\sqrt{6}\\,\\zeta(3)/\\pi^{3}$ and excess kurtosis $12/5$. The density
  is skewed, so $\\mathbb{E}[\\ell^{(\\mu\\sigma)}]$ does not vanish: location
  and scale are not orthogonal here, unlike in a symmetric location-scale
  family."
)

DISTRIBS$skewnormal1 <- list(
  title = "Skew normal, direct",
  ctor = "skewnormal1_distrib()",
  obj = function() skewnormal1_distrib(),
  theta = list(mu = 1, sigma = 2, alpha = 3),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{2}{\\sigma}\\,\\phi(z)\\,\\Phi(\\alpha z), \\qquad z = \\frac{y-\\mu}{\\sigma}",
  moments = "\\mathbb{E}[Y] = \\mu + \\sigma\\delta\\sqrt{2/\\pi}, \\quad \\operatorname{Var}(Y) = \\sigma^{2}\\left(1 - 2\\delta^{2}/\\pi\\right), \\quad \\delta = \\alpha/\\sqrt{1+\\alpha^{2}}",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    log(2) - log(th$sigma) + dnorm(z, log = TRUE) + pnorm(th$alpha * z, log.p = TRUE)
  },
  grid = function(th) seq(-2, 8, length.out = 40),
  note = "Azzalini's direct parametrization. Every derivative is written in
  the inverse Mills ratio, formed on the log scale by
  `numericals7::mills_ratio()`; the expected information is singular at
  $\\alpha = 0$, a property of the parametrization, and the skewness the
  family reaches is bounded by $0.9953$, which is the reason the skew t
  exists."
)

DISTRIBS$skewnormal2 <- list(
  title = "Skew normal, centered",
  ctor = "skewnormal2_distrib()",
  obj = function() skewnormal2_distrib(),
  theta = list(mu = 1, sigma = 2, gamma1 = 0.5),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{2}{\\omega}\\,\\phi\\!\\left(\\frac{y-\\xi}{\\omega}\\right)\\Phi\\!\\left(\\alpha\\,\\frac{y-\\xi}{\\omega}\\right), \\qquad (\\xi, \\omega, \\alpha) = \\mathrm{DP}(\\mu, \\sigma, \\gamma_1)",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}",
  ld = function(y, th) {
    b <- sqrt(2 / pi)
    r <- sign(th$gamma1) * (2 * abs(th$gamma1) / (4 - pi))^(1 / 3)
    xi <- th$mu - th$sigma * r
    omega <- th$sigma * sqrt(1 + r^2)
    alpha <- r / sqrt(b^2 + (b^2 - 1) * r^2)
    z <- (y - xi) / omega
    log(2) - log(omega) + dnorm(z, log = TRUE) + pnorm(alpha * z, log.p = TRUE)
  },
  grid = function(th) seq(-4, 7, length.out = 40),
  note = "The centered parametrization: the parameters are the mean, the
  standard deviation and the skewness, carried onto the direct ones through
  $r = \\{2\\gamma_1/(4-\\pi)\\}^{1/3}$, with $\\xi = \\mu - \\sigma r$,
  $\\omega = \\sigma\\sqrt{1+r^{2}}$ and
  $\\alpha = r/\\sqrt{b^{2} + (b^{2}-1)r^{2}}$, $b = \\sqrt{2/\\pi}$. The map
  is hand-written because a generic chain rule computes it as a difference of
  large numbers. The skewness travels on `bounded(lwr, upr)` with
  $\\lvert\\gamma_1\\rvert < 0.9953$, the supremum of the skewness a skew
  normal can attain: a wider chart (a rhobit link, say, reaching all of
  $(-1,1)$) would let the optimizer propose skewness values for which no
  member of the family exists and the map to $(\\xi, \\omega, \\alpha)$ has
  no solution."
)

DISTRIBS$skewt <- list(
  title = "Skew t",
  ctor = "skewt_distrib()",
  obj = function() skewt_distrib(),
  theta = list(mu = 1, sigma = 2, alpha = 2, nu = 6),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{2}{\\sigma}\\,t_{\\nu}(z)\\;T_{\\nu+1}\\!\\left(\\alpha z\\sqrt{\\frac{\\nu+1}{\\nu+z^{2}}}\\right), \\qquad z = \\frac{y-\\mu}{\\sigma}",
  moments = "",
  ld = function(y, th) {
    z <- (y - th$mu) / th$sigma
    log(2) - log(th$sigma) + dt(z, th$nu, log = TRUE) +
      pt(th$alpha * z * sqrt((th$nu + 1) / (th$nu + z^2)), th$nu + 1, log.p = TRUE)
  },
  grid = function(th) seq(-4, 9, length.out = 40),
  note = "The mean exists for $\\nu > 1$ and the variance for $\\nu > 2$.
  The four-parameter family a location-scale-shape framework wants:
  unbounded skewness where the skew normal saturates. The $(\\mu, \\sigma,
  \\alpha)$ block is closed to fourth order; every component involving $\\nu$
  carries $T_{\\nu+1}$, whose derivative in the degrees of freedom has no
  elementary form, and comes from one verified five-point stencil."
)

DISTRIBS$exponential <- list(
  title = "Exponential",
  ctor = "exponential_distrib()",
  obj = function() exponential_distrib(),
  theta = list(mu = 3),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{\\mu}\\,e^{-y/\\mu}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu^{2}",
  ld = function(y, th) -log(th$mu) - y / th$mu,
  grid = function(th) seq(0.1, 15, length.out = 40),
  note = "In its mean. Not reachable from the Gamma by `fixed()`: unit shape
  is the relation $\\sigma^{2} = \\mu^{2}$ between two parameters, and
  `fixed()` holds a parameter at a value. It agrees with
  `fixed(weibull1_distrib(), sigma = 1)` to machine precision, which is a
  test of both."
)

DISTRIBS$geometric <- list(
  title = "Geometric",
  ctor = "geometric_distrib()",
  obj = function() geometric_distrib(),
  theta = list(mu = 3),
  support = "y \\in \\{0, 1, \\dots\\}",
  pdf_latex = "P(Y=y) = \\frac{1}{1+\\mu}\\left(\\frac{\\mu}{1+\\mu}\\right)^{y}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu(1+\\mu)",
  ld = function(y, th) y * (log(th$mu) - log1p(th$mu)) - log1p(th$mu),
  grid = function(th) 0:15,
  note = "In its mean, with success probability $1/(1+\\mu)$. It agrees with
  `fixed(negbin2_distrib(), theta = 1)` to machine precision."
)

DISTRIBS$chisq <- list(
  title = "Chi-squared",
  ctor = "chisq_distrib()",
  obj = function() chisq_distrib(),
  theta = list(mu = 5),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{y^{\\mu/2-1}e^{-y/2}}{2^{\\mu/2}\\,\\Gamma(\\mu/2)}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = 2\\mu",
  ld = function(y, th) (th$mu / 2 - 1) * log(y) - y / 2 -
    (th$mu / 2) * log(2) - lgamma(th$mu / 2),
  grid = function(th) seq(0.3, 18, length.out = 40),
  note = "The degrees of freedom are the mean and need not be an integer.
  Like the exponential, it is a relation ($\\sigma^{2} = 2\\mu$) and not a
  fixed value, so it is a family of its own."
)

DISTRIBS$betabinom1 <- list(
  title = "Beta-binomial, mean and dispersion",
  ctor = "betabinom1_distrib(size = 12)",
  obj = function() betabinom1_distrib(size = 12),
  theta = list(mu = 0.4, sigma = 0.5),
  support = "y \\in \\{0, \\dots, n\\}",
  pdf_latex = "P(Y=y) = \\binom{n}{y}\\frac{B\\!\\left(y + \\mu/\\sigma,\\; n - y + (1-\\mu)/\\sigma\\right)}{B\\!\\left(\\mu/\\sigma,\\; (1-\\mu)/\\sigma\\right)}",
  moments = "\\mathbb{E}[Y] = n\\mu, \\qquad \\operatorname{Var}(Y) = n\\mu(1-\\mu)\\,\\frac{1+n\\sigma}{1+\\sigma}",
  ld = function(y, th) {
    a <- th$mu / th$sigma; b <- (1 - th$mu) / th$sigma
    lchoose(12, y) + lbeta(y + a, 12 - y + b) - lbeta(a, b)
  },
  grid = function(th) 0:12,
  note = "The `BB` of gamlss, with `size` a constant of the family like a
  binomial's. It is `betabinom2` at shapes $(\\mu/\\sigma, (1-\\mu)/\\sigma)$
  and inherits its third and fourth orders through the chain rule; the
  expected information is an exact finite sum over the support."
)

DISTRIBS$betabinom2 <- list(
  title = "Beta-binomial, shapes",
  ctor = "betabinom2_distrib(size = 12)",
  obj = function() betabinom2_distrib(size = 12),
  theta = list(alpha = 2, beta = 3),
  support = "y \\in \\{0, \\dots, n\\}",
  pdf_latex = "P(Y=y) = \\binom{n}{y}\\frac{B(y+\\alpha,\\; n-y+\\beta)}{B(\\alpha, \\beta)}",
  moments = "\\mathbb{E}[Y] = \\frac{n\\alpha}{\\alpha+\\beta}, \\qquad \\operatorname{Var}(Y) = \\frac{n\\alpha\\beta\\,(\\alpha+\\beta+n)}{(\\alpha+\\beta)^{2}(\\alpha+\\beta+1)}",
  ld = function(y, th) lchoose(12, y) + lbeta(y + th$alpha, 12 - y + th$beta) -
    lbeta(th$alpha, th$beta),
  grid = function(th) 0:12,
  note = "The classical shape parametrization: a binomial whose success
  probability is Beta$(\\alpha, \\beta)$, integrated out."
)

DISTRIBS$gpd <- list(
  title = "Generalized Pareto",
  ctor = "gpd_distrib()",
  obj = function() gpd_distrib(),
  theta = list(sigma = 2, xi = 0.3),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{\\sigma}\\left(1 + \\frac{\\xi y}{\\sigma}\\right)^{-1/\\xi - 1}",
  moments = "\\mathbb{E}[Y] = \\frac{\\sigma}{1-\\xi} \\;(\\xi < 1), \\qquad \\operatorname{Var}(Y) = \\frac{\\sigma^{2}}{(1-\\xi)^{2}(1-2\\xi)} \\;(\\xi < 1/2)",
  ld = function(y, th) -log(th$sigma) - (1 / th$xi + 1) * log1p(th$xi * y / th$sigma),
  grid = function(th) seq(0.1, 15, length.out = 40),
  note = "The excess-over-threshold family. The information exists exactly
  for $\\xi > -1/2$; the compiled kernel evaluates
  $W = \\log(1+u)\\,(z/u)$ with a series below $|\\xi z| = 0.2$, where the
  direct expression cancels catastrophically."
)

DISTRIBS$gengamma1 <- list(
  title = "Generalized gamma, Stacy",
  ctor = "gengamma1_distrib()",
  obj = function() gengamma1_distrib(),
  theta = list(a = 2, d = 3, p = 1.5),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{p\\,y^{d-1}}{a^{d}\\,\\Gamma(d/p)}\\exp\\!\\left\\{-\\left(\\frac{y}{a}\\right)^{p}\\right\\}",
  moments = "\\mathbb{E}[Y] = a\\,\\frac{\\Gamma((d+1)/p)}{\\Gamma(d/p)}, \\qquad \\operatorname{Var}(Y) = a^{2}\\left[\\frac{\\Gamma((d+2)/p)}{\\Gamma(d/p)} - \\frac{\\Gamma((d+1)/p)^{2}}{\\Gamma(d/p)^{2}}\\right]",
  ld = function(y, th) log(th$p) + (th$d - 1) * log(y) - (y / th$a)^th$p -
    th$d * log(th$a) - lgamma(th$d / th$p),
  grid = function(th) seq(0.2, 8, length.out = 40),
  note = "Stacy's three-parameter family, containing the gamma ($p = 1$) and
  the Weibull ($d = p$). Its distribution function saturates well inside the
  support in double precision, so round trips are tested at points obtained
  from the quantile function."
)

DISTRIBS$gengamma2 <- list(
  title = "Generalized gamma, mean form",
  ctor = "gengamma2_distrib()",
  obj = function() gengamma2_distrib(),
  theta = list(mean = 3, d = 3, p = 1.5),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{p\\,y^{d-1}}{a^{d}\\,\\Gamma(d/p)}\\,e^{-(y/a)^{p}}, \\qquad a = \\mu\\,\\frac{\\Gamma(d/p)}{\\Gamma((d+1)/p)}",
  moments = "\\mathbb{E}[Y] = \\mu",
  ld = function(y, th) {
    a <- th$mean * gamma(th$d / th$p) / gamma((th$d + 1) / th$p)
    log(th$p) + (th$d - 1) * log(y) - (y / a)^th$p - th$d * log(a) - lgamma(th$d / th$p)
  },
  grid = function(th) seq(0.2, 8, length.out = 40),
  note = "`gengamma1` with the scale replaced by the mean, through
  `reparametrize()` and a hand-written table of map partials."
)

DISTRIBS$vonmises1 <- list(
  title = "von Mises, concentration",
  ctor = "vonmises1_distrib()",
  obj = function() vonmises1_distrib(),
  theta = list(mu = 0.5, kappa = 2),
  support = "y \\in [-\\pi, \\pi)",
  pdf_latex = "f(y) = \\frac{e^{\\kappa\\cos(y-\\mu)}}{2\\pi I_{0}(\\kappa)}",
  moments = "\\mathbb{E}[\\cos(Y-\\mu)] = A(\\kappa) = I_{1}(\\kappa)/I_{0}(\\kappa)",
  ld = function(y, th) th$kappa * cos(y - th$mu) - log(2 * pi) -
    (log(besselI(th$kappa, 0, expon.scaled = TRUE)) + th$kappa),
  grid = function(th) seq(-3, 3, length.out = 40),
  note = "The circular family; $\\mu$ is the directional mean.
  $\\log I_{0}$ is evaluated by
  `numericals7::log_bessel_i()`, finite past the point where the scaled
  Bessel underflows to an exact zero; the mean direction is carried on a
  bounded link, since an unbounded chart would make the likelihood periodic
  with infinitely many maxima."
)

DISTRIBS$vonmises2 <- list(
  title = "von Mises, mean resultant length",
  ctor = "vonmises2_distrib()",
  obj = function() vonmises2_distrib(),
  theta = list(mu = 0.5, rho = 0.6),
  support = "y \\in [-\\pi, \\pi)",
  pdf_latex = "f(y) = \\frac{e^{\\kappa\\cos(y-\\mu)}}{2\\pi I_{0}(\\kappa)}, \\qquad \\kappa = A^{-1}(\\rho)",
  moments = "\\mathbb{E}[\\cos(Y-\\mu)] = \\rho",
  ld = function(y, th) {
    k <- uniroot(function(k)
      besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE) - th$rho,
      c(1e-8, 500), tol = 1e-14)$root
    k * cos(y - th$mu) - log(2 * pi) - (log(besselI(k, 0, expon.scaled = TRUE)) + k)
  },
  grid = function(th) seq(-3, 3, length.out = 40),
  note = "The directional mean is $\\mu$; the second parameter is the
  mean resultant length
  $\\rho \\in (0, 1)$, which a linear predictor reaches through a logit
  link. The map $\\kappa = A^{-1}(\\rho)$ has no elementary inverse and is
  differentiated by the inverse function rule of
  `numericals7::bessel_i_ratio_inverse()`."
)

DISTRIBS$invgauss2 <- list(
  title = "Inverse Gaussian, mean and shape",
  ctor = "invgauss2_distrib()",
  obj = function() invgauss2_distrib(),
  theta = list(mu = 2, lambda = 3),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\sqrt{\\frac{\\lambda}{2\\pi y^{3}}}\\exp\\!\\left\\{-\\frac{\\lambda(y-\\mu)^{2}}{2\\mu^{2}y}\\right\\}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu^{3}/\\lambda",
  ld = function(y, th) 0.5 * (log(th$lambda) - log(2 * pi) - 3 * log(y)) -
    th$lambda * (y - th$mu)^2 / (2 * th$mu^2 * y),
  grid = function(th) seq(0.2, 8, length.out = 40),
  note = "The textbook $(\\mu, \\lambda)$ form; `invgauss1` carries the
  dispersion $\\varphi = 1/\\lambda$ instead."
)

DISTRIBS$lognormal2 <- list(
  title = "Lognormal, mean and variance",
  ctor = "lognormal2_distrib()",
  obj = function() lognormal2_distrib(),
  theta = list(mean = 3, var = 4),
  support = "y \\in (0, \\infty)",
  pdf_latex = "f(y) = \\frac{1}{y\\sqrt{2\\pi s^{2}}}\\exp\\!\\left\\{-\\frac{(\\log y - m)^{2}}{2s^{2}}\\right\\}, \\quad s^{2} = \\log\\!\\left(1+\\frac{v}{\\mu^{2}}\\right)\\!, \\; m = \\log\\mu - \\frac{s^{2}}{2}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = v",
  ld = function(y, th) {
    s2 <- log1p(th$var / th$mean^2); m <- log(th$mean) - s2 / 2
    -log(y) - 0.5 * log(2 * pi * s2) - (log(y) - m)^2 / (2 * s2)
  },
  grid = function(th) seq(0.2, 12, length.out = 40),
  note = "The parameters are the moments of $Y$ itself rather than of
  $\\log Y$, which is what a model for a positive response usually wants a
  linear predictor to act on."
)

DISTRIBS$beta2 <- list(
  title = "Beta, shapes",
  ctor = "beta2_distrib()",
  obj = function() beta2_distrib(),
  theta = list(alpha = 2, beta = 3),
  support = "y \\in (0, 1)",
  pdf_latex = "f(y) = \\frac{y^{\\alpha-1}(1-y)^{\\beta-1}}{B(\\alpha, \\beta)}",
  moments = "\\mathbb{E}[Y] = \\frac{\\alpha}{\\alpha+\\beta}, \\qquad \\operatorname{Var}(Y) = \\frac{\\alpha\\beta}{(\\alpha+\\beta)^{2}(\\alpha+\\beta+1)}",
  ld = function(y, th) (th$alpha - 1) * log(y) + (th$beta - 1) * log1p(-y) -
    lbeta(th$alpha, th$beta),
  grid = function(th) seq(0.03, 0.97, length.out = 40),
  note = "The classical shape form of the mean-precision `beta1`, matching
  `dbeta(y, alpha, beta)`."
)

DISTRIBS$student_t2 <- list(
  title = "Student t, standard deviation",
  ctor = "student_t2_distrib()",
  obj = function() student_t2_distrib(),
  theta = list(mu = 1, sigma = 2, nu = 6),
  support = "y \\in \\mathbb{R}",
  pdf_latex = "f(y) = \\frac{1}{s_0}\\,t_{\\nu}\\!\\left(\\frac{y-\\mu}{s_0}\\right), \\qquad s_0 = \\sigma\\sqrt{\\frac{\\nu-2}{\\nu}}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\sigma^{2}",
  ld = function(y, th) {
    s0 <- th$sigma * sqrt((th$nu - 2) / th$nu)
    dt((y - th$mu) / s0, th$nu, log = TRUE) - log(s0)
  },
  grid = function(th) seq(-6, 8, length.out = 40),
  note = "`student_t1` with the scale replaced by the standard deviation of
  $Y$, which requires $\\nu > 2$: the family exists only where its second
  moment does, and this is the price of the readable parameter."
)

DISTRIBS$pig1 <- list(
  title = "Poisson-inverse Gaussian, dispersion",
  ctor = "pig1_distrib()",
  obj = function() pig1_distrib(),
  theta = list(mu = 3, sigma = 1.2),
  support = "y \\in \\{0, 1, \\dots\\}",
  pdf_latex = "P(Y=y) = \\sqrt{\\frac{2\\alpha}{\\pi}}\\,\\frac{\\mu^{y}e^{1/\\sigma}}{(\\alpha\\sigma)^{y}\\,y!}\\,K_{y-1/2}(\\alpha), \\qquad \\alpha = \\sqrt{\\frac{1}{\\sigma^{2}} + \\frac{2\\mu}{\\sigma}}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu + \\sigma\\mu^{2}",
  ld = function(y, th) {
    al <- sqrt(1 / th$sigma^2 + 2 * th$mu / th$sigma)
    0.5 * (log(2) + log(al) - log(pi)) + y * log(th$mu) + 1 / th$sigma -
      y * (log(al) + log(th$sigma)) - lfactorial(y) +
      log(besselK(al, y - 0.5, expon.scaled = TRUE)) - al
  },
  grid = function(th) 0:15,
  note = "A Poisson whose rate is inverse Gaussian, integrated out: heavier
  in the tail than the negative binomial at the same variance. The
  parameters are $\\mu$ and $\\sigma$ alone; $\\alpha$ is not a third
  parameter but the abbreviation defined in the display, written out
  because every derivative is organized around it. The Bessel order is
  half-integer, so $K_{y-1/2}$ is a finite sum and the derivatives are
  hand-written closed forms."
)

DISTRIBS$pig2 <- list(
  title = "Poisson-inverse Gaussian, Bessel argument",
  ctor = "pig2_distrib()",
  obj = function() pig2_distrib(),
  theta = list(mu = 3, alpha = 2),
  support = "y \\in \\{0, 1, \\dots\\}",
  pdf_latex = "P(Y=y) = \\sqrt{\\frac{2\\alpha}{\\pi}}\\,\\frac{\\mu^{y}e^{1/\\sigma}}{(\\alpha\\sigma)^{y}\\,y!}\\,K_{y-1/2}(\\alpha), \\qquad \\sigma = \\frac{1}{\\sqrt{\\mu^{2}+\\alpha^{2}} - \\mu}",
  moments = "\\mathbb{E}[Y] = \\mu, \\qquad \\operatorname{Var}(Y) = \\mu + \\sigma\\mu^{2}",
  ld = function(y, th) {
    s <- 1 / (sqrt(th$mu^2 + th$alpha^2) - th$mu)
    0.5 * (log(2) + log(th$alpha) - log(pi)) + y * log(th$mu) + 1 / s -
      y * (log(th$alpha) + log(s)) - lfactorial(y) +
      log(besselK(th$alpha, y - 0.5, expon.scaled = TRUE)) - th$alpha
  },
  grid = function(th) 0:15,
  note = "`pig1` reparametrized by the Bessel argument $\\alpha$, in which
  the derivatives of the normalizing constant are rational. The parameters
  are $\\mu$ and $\\alpha$ alone; here $\\sigma$ is the abbreviation,
  defined in the display, that carries the density over from `pig1`."
)


# The entries render in list order; sorting by constructor keeps the numbered
# parametrizations of one family adjacent, whatever order they were written in.
DISTRIBS <- DISTRIBS[order(vapply(DISTRIBS, function(r)
  sub("_distrib$", "", sub("\\(.*$", "", r$ctor)), character(1)))]

# ---------------------------------------------------------------------------
# Catalog rendering
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

catalog_table <- function() {
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
  d <- rec$obj()
  cat("| parameter | interpretation | domain | default link |\n")
  cat("|:---|:---|:---|:---|\n")
  for (par in d@params) {
    cat(sprintf("| `%s` | %s | %s | `%s` |\n",
                par,
                d@params_interpretation[[par]],
                fmt_bounds(d@params_bounds[[par]], FALSE),
                d@link_params[[par]]@link_name))
  }
  cat("\n")
  if (nzchar(rec$moments)) cat(sprintf("$$%s$$\n\n", rec$moments))
  cat(gsub("\n\\s+", "\n", rec$note), "\n\n")
}

# Silent consistency gate for Chapter 3: every printed density against the
# implementation, the Bartlett identities, the link-scale machinery, and the
# package's own validator. Stops the render on any disagreement; never
# rendered.
assert_distributions_ok <- function() {
  # every univariate constructor has a catalog record, so a new family cannot
  # ship without an entry
  rec_ctors <- vapply(DISTRIBS, function(r) sub("\\(.*$", "", r$ctor), character(1))
  missing <- setdiff(ALL_FAMILY_CTORS, rec_ctors)
  if (length(missing)) {
    stop("Families without a catalog record: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
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

  # The cdf-derivative section: each family checked against the route the
  # implementation does not use, so neither comparison is tautological.
  for (id in c("gaussian", "logistic", "gamma", "beta")) {
    rec <- DISTRIBS[[id]]; d <- rec$obj(); th <- rec$theta
    q <- rec$grid(th)[c(10, 20, 30)]
    a <- distrib_grad_cdf(d, q, th, log = FALSE)
    for (p in d@params) {
      b <- vapply(q, function(qq) stats::integrate(
        function(y) distrib_pdf(d, y, th) * distrib_gradient(d, y, th)[[p]],
        d@bounds[1], qq, rel.tol = 1e-10)$value, numeric(1))
      if (rel_err(a[[p]], b) > 1e-6) {
        stop(sprintf("cdf gradient of '%s' in %s disagrees with the partial score.", id, p),
             call. = FALSE)
      }
    }
  }
  for (id in c("poisson", "negbin")) {
    rec <- DISTRIBS[[id]]; d <- rec$obj(); th <- rec$theta
    q <- c(1, 4, 9)
    a <- distrib_grad_cdf(d, q, th, log = FALSE)
    for (p in d@params) {
      j <- match(p, d@params); hh <- 1e-5 * max(1, abs(th[[j]]))
      tp <- tm <- th; tp[[j]] <- th[[j]] + hh; tm[[j]] <- th[[j]] - hh
      b <- (distrib_cdf(d, q, tp) - distrib_cdf(d, q, tm)) / (2 * hh)
      if (rel_err(a[[p]], b) > 1e-5) {
        stop(sprintf("cdf gradient of '%s' in %s disagrees with finite differences.", id, p),
             call. = FALSE)
      }
    }
  }
  # dF/dmu = -f for a location-scale family
  gg <- gaussian1_distrib()
  qq <- c(0, 1.2, 3)
  if (rel_err(distrib_grad_cdf(gg, qq, list(mu = 1.2, sigma = 1.7), log = FALSE)$mu,
              -dnorm(qq, 1.2, 1.7)) > 1e-12) {
    stop("The location-scale closed form for dF/dmu no longer holds.", call. = FALSE)
  }

  # The non-regular section makes three specific numerical claims; pin them.
  # (a) the shipped Laplace has a closed-form E[H], so approx is ignored and all
  #     strategies return -1/sigma^2;
  # (b) on a bare Laplace without closed forms, bartlett recovers -1/sigma^2 while
  #     integrate and mc, which both average the a.e.-zero observed l_mumu,
  #     return zero.
  lap <- laplace_distrib()
  th <- list(mu = 0, sigma = 2)
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
    list(mu = sign(r) / b, sigma = (abs(r) / b - 1) / b)
  }
  S7::method(distrib_hessian, BareLap) <- function(distrib, y, theta,
                                                   scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]; b <- theta[[2]]; n <- length(y)
    list(mu_mu = rep(0, n), sigma_sigma = (b - 2 * abs(r)) / b^3, mu_sigma = -sign(r) / b^2)
  }
  bare <- BareLap(
    distrib_name = "bare laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "location", sigma = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = identity_link(), sigma = log_link()),
    params_smooth = c(mu = FALSE, sigma = TRUE)
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
# Certificate: check_distrib() on the whole catalog
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


# ---------------------------------------------------------------------------
# The full family list, the numbering convention and its gates (section 3.1).
#
# Everything here is generated from the package's own constructors, so the
# tables cannot drift from the code: a family added to the package appears,
# a renamed one changes, and the gate below fails if a numbered name promises
# a parametrization that does not exist.
# ---------------------------------------------------------------------------

# every univariate constructor the package exports, by its exported name
ALL_FAMILY_CTORS <- c(
  "gaussian1_distrib", "gaussian2_distrib", "gaussian3_distrib",
  "cauchy_distrib", "logistic_distrib", "student_t1_distrib",
  "student_t2_distrib", "laplace_distrib", "laplace2_distrib",
  "enet_distrib", "pseudohuber_distrib",
  "skewnormal1_distrib", "skewnormal2_distrib", "skewt_distrib",
  "gumbel_distrib", "gamma1_distrib", "gamma2_distrib",
  "invgauss1_distrib", "invgauss2_distrib", "lognormal1_distrib",
  "lognormal2_distrib", "weibull1_distrib", "weibull3_distrib",
  "exponential_distrib", "chisq_distrib", "beta1_distrib", "beta2_distrib",
  "bernoulli_distrib", "binomial_distrib", "poisson_distrib",
  "negbin1_distrib", "negbin2_distrib", "geometric_distrib",
  "betabinom1_distrib", "betabinom2_distrib", "gpd_distrib",
  "pig1_distrib", "pig2_distrib",
  "vonmises1_distrib", "vonmises2_distrib", "gengamma1_distrib",
  "gengamma2_distrib"
)

.family_object <- function(ctor) {
  f <- get(ctor, envir = asNamespace("distributions7"))
  if ("size" %in% names(formals(f))) f(size = 10) else f()
}

full_catalog_table <- function() {
  rows <- lapply(ALL_FAMILY_CTORS, function(ctor) {
    d <- .family_object(ctor)
    is_disc <- S7::S7_inherits(d, discrete_distrib)
    data.frame(
      constructor = sprintf("`%s()`", ctor),
      kind = if (is_disc) "discrete" else "continuous",
      parameters = paste(sprintf("`%s`", d@params), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

numbered_families_table <- function() {
  nm <- sub("_distrib$", "", ALL_FAMILY_CTORS)
  stem <- sub("[0-9]+$", "", nm)
  digit <- sub("^[a-z_]+", "", nm)
  keep <- stem %in% stem[digit != ""]
  df <- data.frame(stem = stem[keep], name = nm[keep],
                   stringsAsFactors = FALSE)
  agg <- aggregate(name ~ stem, df, function(x) paste(sprintf("`%s`", x),
                                                      collapse = ", "))
  names(agg) <- c("family", "parametrizations")
  agg[order(agg$family), ]
}

# The claims of the numbering and moments subsections.
.certify_numbering <- function() {
  out <- character()

  # every listed constructor exists and returns a distribution
  for (ctor in ALL_FAMILY_CTORS) {
    ok <- tryCatch({
      d <- .family_object(ctor)
      S7::S7_inherits(d, distrib)
    }, error = function(e) FALSE)
    if (!ok) out <- c(out, sprintf("%s does not construct a distribution", ctor))
  }

  # the list is the package's export list, not a stale copy: every exported
  # *_distrib constructor that is univariate appears exactly once
  exported <- grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE)
  exported <- setdiff(exported, c(
    # the multivariate families, which have their own section
    "mvgaussian_distrib", "mvstudent_t_distrib",
    "dirichlet_distrib", "multinomial_distrib",
    # exports that end in _distrib without being family constructors
    "check_distrib", "fit_distrib",
    "continuous_distrib", "discrete_distrib", "multivariate_distrib"
  ))
  missing <- setdiff(exported, ALL_FAMILY_CTORS)
  extra <- setdiff(ALL_FAMILY_CTORS, exported)
  if (length(missing)) {
    out <- c(out, sprintf("families missing from the catalog: %s",
                          paste(missing, collapse = ", ")))
  }
  if (length(extra)) {
    out <- c(out, sprintf("catalog names not exported: %s",
                          paste(extra, collapse = ", ")))
  }

  # two parametrizations of one family are the same law: equal densities on
  # a grid, through the map each pair documents
  y <- c(0.3, 1.1, 2.7, 5)
  a <- distrib_pdf(distributions7::gaussian1_distrib(), y,
                   list(mu = 1, sigma = 2), log = TRUE)
  b <- distrib_pdf(distributions7::gaussian2_distrib(), y,
                   list(mu = 1, sigma2 = 4), log = TRUE)
  cc <- distrib_pdf(distributions7::gaussian3_distrib(), y,
                    list(mu = 1, tau = 0.25), log = TRUE)
  if (max(abs(a - b), abs(a - cc)) > 1e-12) {
    out <- c(out, "the three gaussian parametrizations are not the same law")
  }

  out
}

.certify_moments <- function() {
  out <- character()
  # closed moments against the numerical route, which shares no formulas
  cases <- list(
    list(d = distributions7::gamma2_distrib(), th = list(mu = 3, sigma2 = 2)),
    list(d = distributions7::weibull1_distrib(), th = list(mu = 2, sigma = 1.5)),
    list(d = distributions7::negbin2_distrib(), th = list(mu = 4, theta = 1.3))
  )
  for (cs in cases) {
    m <- mean(cs$d, cs$th)
    v <- distributions7::variance(cs$d, cs$th)
    m_num <- distributions7::expectation(cs$d, function(y, theta) y, cs$th)
    v_num <- distributions7::expectation(cs$d, function(y, theta) (y - m)^2, cs$th)
    if (abs(m - m_num) > 1e-6 * max(1, abs(m)) ||
        abs(v - v_num) > 1e-6 * max(1, abs(v))) {
      out <- c(out, sprintf("closed moments of %s disagree with quadrature",
                            cs$d@distrib_name))
    }
  }
  out
}

.certify_reparametrize_twin <- function() {
  out <- character()
  build <- function(md) {
    distributions7::reparametrize(
      distributions7::gaussian1_distrib(),
      map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
      params = c("mu", "sigma2"),
      bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
      links = list(mu = linkfunctions7::identity_link(),
                   sigma2 = linkfunctions7::log_link()),
      map_derivs = md
    )
  }
  md <- function(psi) {
    v <- psi[[2]]
    list(
      list("1" = rep_len(1, length(v))),
      list("2" = 0.5 / sqrt(v), "2,2" = -0.25 / v^1.5,
           "2,2,2" = 0.375 / v^2.5, "2,2,2,2" = -0.9375 / v^3.5)
    )
  }
  h <- distributions7::gaussian2_distrib()
  y <- c(-0.5, 0.8, 2.2)
  th <- list(mu = 1.2, sigma2 = 3.5)
  gap_for <- function(r, fn) {
    ga <- fn(r, y, th)
    gb <- fn(h, y, th)
    max(vapply(names(gb), function(nm) max(abs(ga[[nm]] - gb[[nm]])),
               numeric(1)))
  }
  # the explicit tables the chunk shows are exact
  r <- build(md)
  for (fn in list(distributions7::distrib_gradient,
                  distributions7::distrib_hessian)) {
    if (gap_for(r, fn) > 1e-10) {
      out <- c(out, "reparametrize() and the hand-written gaussian2 disagree")
    }
  }
  # the stencil fallback serves at its own accuracy
  r0 <- build(NULL)
  if (gap_for(r0, distributions7::distrib_gradient) > 1e-5) {
    out <- c(out, "the stencil-fallback gradient is worse than its floor")
  }
  out
}

assert_numbering_ok <- function() {
  problems <- c(.certify_numbering(), .certify_moments())
  if (length(problems)) {
    stop("Section 3.1's family catalog disagrees with the package:
  ",
         paste(problems, collapse = "
  "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_reparametrize_ok <- function() {
  problems <- .certify_reparametrize_twin()
  if (length(problems)) {
    stop("Section 3.2's reparametrization claims disagree with the package:
  ",
         paste(problems, collapse = "
  "), call. = FALSE)
  }
  invisible(TRUE)
}
