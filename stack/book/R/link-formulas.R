# The link catalogue: one record per link class.
#
# Each record carries, side by side:
#
#   *_latex   the formula as it is printed in Chapter 2;
#   g, dg, h, dh   an R transcription of exactly that printed formula, written
#                  here and nowhere else -- in particular NOT copied from the
#                  package, and deliberately naive where the package is clever
#                  about overflow.
#
# Chapter 2 renders the LaTeX from this file and, in the same pass, evaluates the
# transcription against the method linkfunctions7 actually dispatches. Because
# both come from this one record, the text cannot drift away from what is tested:
# if the printed formula is wrong, the certification table says DISAGREE.
#
# A third, genuinely independent check compares the package against numDeriv.
# That one shares no code with either side.

# --- small helpers used by the transcriptions ------------------------------
.zero <- function(x) rep(0, length(x))
.one <- function(x) rep(1, length(x))
.sigma <- function(x) 1 / (1 + exp(-x))
.falling <- function(a, k) prod(a - seq_len(k) + 1) # a(a-1)...(a-k+1)

LINKS <- list()

# ---------------------------------------------------------------- identity --
LINKS$identity <- list(
  title = "Identity",
  ctor = "identity_link()",
  obj = function() identity_link(),
  theta_dom = "\\theta \\in \\mathbb{R}",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\theta",
  h_latex = "h(\\eta) = \\eta",
  dg_latex = c("1", "0", "0", "0"),
  dh_latex = c("1", "0", "0", "0"),
  g = function(th) th,
  dg = list(.one, .zero, .zero, .zero),
  h = function(e) e,
  dh = list(.one, .zero, .zero, .zero),
  grid_theta = seq(-3, 3, length.out = 25),
  grid_eta = seq(-3, 3, length.out = 25),
  note = "The trivial case, kept as an object so that a parameter that needs no
  transformation is still described by the same interface as one that does."
)

# --------------------------------------------------------------------- log --
LINKS$log <- list(
  title = "Log",
  ctor = "log_link()",
  obj = function() log_link(),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log\\theta",
  h_latex = "h(\\eta) = e^{\\eta}",
  dg_latex = c("\\dfrac{1}{\\theta}", "-\\dfrac{1}{\\theta^{2}}",
               "\\dfrac{2}{\\theta^{3}}", "-\\dfrac{6}{\\theta^{4}}"),
  dh_latex = c("e^{\\eta}", "e^{\\eta}", "e^{\\eta}", "e^{\\eta}"),
  g = function(th) log(th),
  dg = list(function(th) 1 / th, function(th) -1 / th^2,
            function(th) 2 / th^3, function(th) -6 / th^4),
  h = function(e) exp(e),
  dh = list(exp, exp, exp, exp),
  grid_theta = exp(seq(log(0.05), log(20), length.out = 25)),
  grid_eta = seq(-3, 3, length.out = 25),
  note = "In general $g^{(k)}(\\theta) = (-1)^{k-1}(k-1)!\\,\\theta^{-k}$, and the
  inverse is its own derivative to every order. The implementation floors
  $h$ and its derivatives at `.Machine$double.eps`: without it, a large negative
  $\\eta$ underflows to exactly zero and the subsequent $1/\\theta$ returns `Inf`."
)

# -------------------------------------------------------------------- sqrt --
LINKS$sqrt <- list(
  title = "Square root",
  ctor = "sqrt_link()",
  obj = function() sqrt_link(),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in (0, \\infty)",
  g_latex = "g(\\theta) = \\sqrt{\\theta}",
  h_latex = "h(\\eta) = \\eta^{2}",
  dg_latex = c("\\tfrac{1}{2}\\theta^{-1/2}", "-\\tfrac{1}{4}\\theta^{-3/2}",
               "\\tfrac{3}{8}\\theta^{-5/2}", "-\\tfrac{15}{16}\\theta^{-7/2}"),
  dh_latex = c("2\\eta", "2", "0", "0"),
  g = function(th) sqrt(th),
  dg = list(function(th) 0.5 * th^(-0.5), function(th) -0.25 * th^(-1.5),
            function(th) (3 / 8) * th^(-2.5), function(th) -(15 / 16) * th^(-3.5)),
  h = function(e) e^2,
  dh = list(function(e) 2 * e, function(e) rep(2, length(e)), .zero, .zero),
  grid_theta = seq(0.2, 9, length.out = 25),
  grid_eta = seq(0.2, 3, length.out = 25),
  note = "The inverse is a quadratic, so its third and fourth derivatives vanish
  identically -- the only link in the catalogue where the fourth-order link-scale
  correction disappears without $\\theta$ being an affine function of $\\eta$.
  $\\eta$ is restricted to the positive half-line to keep the map injective."
)

# ----------------------------------------------------------------- inverse --
LINKS$inverse <- list(
  title = "Inverse (reciprocal)",
  ctor = "inverse_link()",
  obj = function() inverse_link(),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in (0, \\infty)",
  g_latex = "g(\\theta) = \\theta^{-1}",
  h_latex = "h(\\eta) = \\eta^{-1}",
  dg_latex = c("-\\dfrac{1}{\\theta^{2}}", "\\dfrac{2}{\\theta^{3}}",
               "-\\dfrac{6}{\\theta^{4}}", "\\dfrac{24}{\\theta^{5}}"),
  dh_latex = c("-\\dfrac{1}{\\eta^{2}}", "\\dfrac{2}{\\eta^{3}}",
               "-\\dfrac{6}{\\eta^{4}}", "\\dfrac{24}{\\eta^{5}}"),
  g = function(th) 1 / th,
  dg = list(function(th) -1 / th^2, function(th) 2 / th^3,
            function(th) -6 / th^4, function(th) 24 / th^5),
  h = function(e) 1 / e,
  dh = list(function(e) -1 / e^2, function(e) 2 / e^3,
            function(e) -6 / e^4, function(e) 24 / e^5),
  grid_theta = seq(0.2, 6, length.out = 25),
  grid_eta = seq(0.2, 5, length.out = 25),
  note = "An involution, $g = h$, so the two directions carry the same formulas:
  $g^{(k)}(\\theta) = (-1)^{k}k!\\,\\theta^{-(k+1)}$. The linear predictor must be
  kept strictly on one side of zero."
)

# -------------------------------------------------------------- inverse_sq --
LINKS$inverse_sq <- list(
  title = "Inverse square",
  ctor = "inverse_sq_link()",
  obj = function() inverse_sq_link(),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in (0, \\infty)",
  g_latex = "g(\\theta) = \\theta^{-2}",
  h_latex = "h(\\eta) = \\eta^{-1/2}",
  dg_latex = c("-\\dfrac{2}{\\theta^{3}}", "\\dfrac{6}{\\theta^{4}}",
               "-\\dfrac{24}{\\theta^{5}}", "\\dfrac{120}{\\theta^{6}}"),
  dh_latex = c("-\\tfrac{1}{2}\\eta^{-3/2}", "\\tfrac{3}{4}\\eta^{-5/2}",
               "-\\tfrac{15}{8}\\eta^{-7/2}", "\\tfrac{105}{16}\\eta^{-9/2}"),
  g = function(th) th^(-2),
  dg = list(function(th) -2 / th^3, function(th) 6 / th^4,
            function(th) -24 / th^5, function(th) 120 / th^6),
  h = function(e) e^(-0.5),
  dh = list(function(e) -0.5 * e^(-1.5), function(e) 0.75 * e^(-2.5),
            function(e) -(15 / 8) * e^(-3.5), function(e) (105 / 16) * e^(-4.5)),
  grid_theta = seq(0.3, 4, length.out = 25),
  grid_eta = seq(0.1, 5, length.out = 25),
  note = "The canonical link of the inverse Gaussian. Both $\\theta$ and $\\eta$
  must stay strictly positive: fractional powers of a negative argument return
  `NaN` rather than raising."
)

# ------------------------------------------------------------------- logit --
LINKS$logit <- list(
  title = "Logit",
  ctor = "logit_link()",
  obj = function() logit_link(),
  theta_dom = "\\theta \\in (0, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log\\dfrac{\\theta}{1-\\theta}",
  h_latex = "h(\\eta) = \\dfrac{1}{1+e^{-\\eta}} \\equiv p",
  dg_latex = c("\\dfrac{1}{\\theta(1-\\theta)}",
               "\\dfrac{2\\theta-1}{\\theta^{2}(1-\\theta)^{2}}",
               "\\dfrac{2}{\\theta^{3}} + \\dfrac{2}{(1-\\theta)^{3}}",
               "-\\dfrac{6}{\\theta^{4}} + \\dfrac{6}{(1-\\theta)^{4}}"),
  dh_latex = c("p(1-p)", "p(1-p)(1-2p)", "p(1-p)\\left(1-6p+6p^{2}\\right)",
               "p(1-p)\\left(1-14p+36p^{2}-24p^{3}\\right)"),
  g = function(th) log(th / (1 - th)),
  dg = list(function(th) 1 / (th * (1 - th)),
            function(th) (2 * th - 1) / (th^2 * (1 - th)^2),
            function(th) 2 / th^3 + 2 / (1 - th)^3,
            function(th) -6 / th^4 + 6 / (1 - th)^4),
  h = function(e) .sigma(e),
  dh = list(function(e) { p <- .sigma(e); p * (1 - p) },
            function(e) { p <- .sigma(e); p * (1 - p) * (1 - 2 * p) },
            function(e) { p <- .sigma(e); p * (1 - p) * (1 - 6 * p + 6 * p^2) },
            function(e) { p <- .sigma(e); p * (1 - p) * (1 - 14 * p + 36 * p^2 - 24 * p^3) }),
  grid_theta = seq(0.05, 0.95, length.out = 25),
  grid_eta = seq(-4, 4, length.out = 25),
  note = "Splitting $g' = \\theta^{-1} + (1-\\theta)^{-1}$ makes the higher orders
  immediate: $g^{(k)} = (-1)^{k-1}(k-1)!\\,\\theta^{-k} + (k-1)!\\,(1-\\theta)^{-k}$.
  On the inverse side every derivative is a polynomial in $p$ times $p(1-p)$;
  the coefficients $1, -2, (1,-6,6), (1,-14,36,-24)$ are the Eulerian-type numbers
  generated by $\\partial_\\eta$ acting on $p(1-p)$."
)

# ------------------------------------------------------------------ probit --
LINKS$probit <- list(
  title = "Probit",
  ctor = "probit_link()",
  obj = function() probit_link(),
  theta_dom = "\\theta \\in (0, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\Phi^{-1}(\\theta)",
  h_latex = "h(\\eta) = \\Phi(\\eta)",
  dg_latex = c("\\dfrac{1}{\\varphi(z)}", "\\dfrac{z}{\\varphi(z)^{2}}",
               "\\dfrac{1+2z^{2}}{\\varphi(z)^{3}}",
               "\\dfrac{7z+6z^{3}}{\\varphi(z)^{4}}"),
  dh_latex = c("\\varphi(\\eta)", "-\\eta\\,\\varphi(\\eta)",
               "(\\eta^{2}-1)\\,\\varphi(\\eta)", "(3\\eta-\\eta^{3})\\,\\varphi(\\eta)"),
  g = function(th) qnorm(th),
  dg = list(function(th) { z <- qnorm(th); 1 / dnorm(z) },
            function(th) { z <- qnorm(th); z / dnorm(z)^2 },
            function(th) { z <- qnorm(th); (1 + 2 * z^2) / dnorm(z)^3 },
            function(th) { z <- qnorm(th); (7 * z + 6 * z^3) / dnorm(z)^4 }),
  h = function(e) pnorm(e),
  dh = list(function(e) dnorm(e), function(e) -e * dnorm(e),
            function(e) (e^2 - 1) * dnorm(e), function(e) (3 * e - e^3) * dnorm(e)),
  grid_theta = seq(0.05, 0.95, length.out = 25),
  grid_eta = seq(-3, 3, length.out = 25),
  note = "The inverse derivatives are the probabilists' Hermite polynomials:
  $h^{(k)}(\\eta) = (-1)^{k-1}He_{k-1}(\\eta)\\varphi(\\eta)$ with
  $He_0 = 1$, $He_1 = \\eta$, $He_2 = \\eta^{2}-1$, $He_3 = \\eta^{3}-3\\eta$.
  In the forward direction there is no closed form free of $\\Phi^{-1}$, so the
  implementation evaluates $z = \\Phi^{-1}(\\theta)$ once and writes everything in
  terms of $z$ and $\\varphi(z)$."
)

# ----------------------------------------------------------------- cloglog --
LINKS$cloglog <- list(
  title = "Complementary log-log",
  ctor = "cloglog_link()",
  obj = function() cloglog_link(),
  theta_dom = "\\theta \\in (0, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log\\!\\left(-\\log(1-\\theta)\\right)",
  h_latex = "h(\\eta) = 1 - \\exp\\!\\left(-e^{\\eta}\\right)",
  dg_latex = c("-\\dfrac{1}{vL}", "-\\dfrac{L+1}{v^{2}L^{2}}",
               "-\\dfrac{2L^{2}+3L+2}{v^{3}L^{3}}",
               "-\\dfrac{6L^{3}+11L^{2}+12L+6}{v^{4}L^{4}}"),
  dh_latex = c("z\\,e^{-z}", "(z - z^{2})\\,e^{-z}",
               "(z - 3z^{2} + z^{3})\\,e^{-z}",
               "(z - 7z^{2} + 6z^{3} - z^{4})\\,e^{-z}"),
  g = function(th) log(-log(1 - th)),
  dg = list(function(th) { v <- 1 - th; L <- log(v); -1 / (v * L) },
            function(th) { v <- 1 - th; L <- log(v); -(L + 1) / (v^2 * L^2) },
            function(th) { v <- 1 - th; L <- log(v); -(2 * L^2 + 3 * L + 2) / (v^3 * L^3) },
            function(th) { v <- 1 - th; L <- log(v)
                           -(6 * L^3 + 11 * L^2 + 12 * L + 6) / (v^4 * L^4) }),
  h = function(e) 1 - exp(-exp(e)),
  dh = list(function(e) { z <- exp(e); z * exp(-z) },
            function(e) { z <- exp(e); (z - z^2) * exp(-z) },
            function(e) { z <- exp(e); (z - 3 * z^2 + z^3) * exp(-z) },
            function(e) { z <- exp(e); (z - 7 * z^2 + 6 * z^3 - z^4) * exp(-z) }),
  grid_theta = seq(0.05, 0.95, length.out = 25),
  grid_eta = seq(-3, 1.2, length.out = 25),
  note = "Here $v = 1-\\theta$, $L = \\log v < 0$ and $z = e^{\\eta}$. The
  asymmetric link of the extreme-value family. The implementation evaluates the
  inverse derivatives as $\\sum_k c_k \\exp(k\\eta - z)$ rather than in the
  factored form $(\\cdot)\\,e^{-z}$ shown above: for large $\\eta$ the factored
  form computes $\\infty \\times 0$ and returns `NaN`, whereas the combined
  exponent underflows cleanly to zero. Wherever both are finite the two forms
  are algebraically identical."
)

# ------------------------------------------------------------------ loglog --
LINKS$loglog <- list(
  title = "Log-log",
  ctor = "loglog_link()",
  obj = function() loglog_link(),
  theta_dom = "\\theta \\in (0, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = -\\log(-\\log\\theta)",
  h_latex = "h(\\eta) = \\exp\\!\\left(-e^{-\\eta}\\right)",
  dg_latex = c("-\\dfrac{1}{\\theta\\,l}", "\\dfrac{1+l}{\\theta^{2}l^{2}}",
               "-\\dfrac{2+3l+2l^{2}}{\\theta^{3}l^{3}}",
               "\\dfrac{6+12l+11l^{2}+6l^{3}}{\\theta^{4}l^{4}}"),
  dh_latex = c("z\\,e^{-z}", "(z^{2}-z)\\,e^{-z}",
               "(z^{3}-3z^{2}+z)\\,e^{-z}", "(z^{4}-6z^{3}+7z^{2}-z)\\,e^{-z}"),
  g = function(th) -log(-log(th)),
  dg = list(function(th) { l <- log(th); -1 / (th * l) },
            function(th) { l <- log(th); (1 + l) / (th^2 * l^2) },
            function(th) { l <- log(th); -(2 + 3 * l + 2 * l^2) / (th^3 * l^3) },
            function(th) { l <- log(th); (6 + 12 * l + 11 * l^2 + 6 * l^3) / (th^4 * l^4) }),
  h = function(e) exp(-exp(-e)),
  dh = list(function(e) { z <- exp(-e); z * exp(-z) },
            function(e) { z <- exp(-e); (z^2 - z) * exp(-z) },
            function(e) { z <- exp(-e); (z^3 - 3 * z^2 + z) * exp(-z) },
            function(e) { z <- exp(-e); (z^4 - 6 * z^3 + 7 * z^2 - z) * exp(-z) }),
  grid_theta = seq(0.05, 0.95, length.out = 25),
  grid_eta = seq(-1.2, 3, length.out = 25),
  note = "With $l = \\log\\theta < 0$ and $z = e^{-\\eta}$. The mirror image of
  the complementary log-log: $g_{\\mathrm{loglog}}(\\theta) =
  -g_{\\mathrm{cloglog}}(1-\\theta)$, which is why the two sets of numerators are
  the same polynomials with alternating signs."
)

# ----------------------------------------------------------------- cauchit --
LINKS$cauchit <- list(
  title = "Cauchit",
  ctor = "cauchit_link()",
  obj = function() cauchit_link(),
  theta_dom = "\\theta \\in (0, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\tan\\!\\left(\\pi(\\theta - \\tfrac{1}{2})\\right)",
  h_latex = "h(\\eta) = \\tfrac{1}{2} + \\tfrac{1}{\\pi}\\arctan\\eta",
  dg_latex = c("\\pi(1+z^{2})", "2\\pi^{2}z(1+z^{2})",
               "2\\pi^{3}(1+z^{2})(1+3z^{2})", "8\\pi^{4}z(1+z^{2})(2+3z^{2})"),
  dh_latex = c("\\dfrac{1}{\\pi(1+\\eta^{2})}", "\\dfrac{-2\\eta}{\\pi(1+\\eta^{2})^{2}}",
               "\\dfrac{2(3\\eta^{2}-1)}{\\pi(1+\\eta^{2})^{3}}",
               "\\dfrac{24\\eta(1-\\eta^{2})}{\\pi(1+\\eta^{2})^{4}}"),
  g = function(th) tan(pi * (th - 0.5)),
  dg = list(function(th) { z <- tan(pi * (th - 0.5)); pi * (1 + z^2) },
            function(th) { z <- tan(pi * (th - 0.5)); 2 * pi^2 * z * (1 + z^2) },
            function(th) { z <- tan(pi * (th - 0.5)); 2 * pi^3 * (1 + z^2) * (1 + 3 * z^2) },
            function(th) { z <- tan(pi * (th - 0.5)); 8 * pi^4 * z * (1 + z^2) * (2 + 3 * z^2) }),
  h = function(e) 0.5 + atan(e) / pi,
  dh = list(function(e) 1 / (pi * (1 + e^2)),
            function(e) -2 * e / (pi * (1 + e^2)^2),
            function(e) 2 * (3 * e^2 - 1) / (pi * (1 + e^2)^3),
            function(e) 24 * e * (1 - e^2) / (pi * (1 + e^2)^4)),
  grid_theta = seq(0.1, 0.9, length.out = 25),
  grid_eta = seq(-4, 4, length.out = 25),
  note = "With $z = g(\\theta)$. The heavy-tailed alternative to probit: the
  inverse is the Cauchy cdf, so $h'$ is the Cauchy density and the higher
  derivatives are its successive derivatives. Because $h'$ decays only like
  $\\eta^{-2}$, an extreme linear predictor still moves the probability, which is
  the point of choosing it."
)

# ------------------------------------------------------------------ rhobit --
LINKS$rhobit <- list(
  title = "Rhobit (Fisher's z)",
  ctor = "rhobit_link()",
  obj = function() rhobit_link(),
  theta_dom = "\\theta \\in (-1, 1)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\operatorname{artanh}\\theta = \\tfrac{1}{2}\\log\\dfrac{1+\\theta}{1-\\theta}",
  h_latex = "h(\\eta) = \\tanh\\eta \\equiv t",
  dg_latex = c("\\dfrac{1}{1-\\theta^{2}}", "\\dfrac{2\\theta}{(1-\\theta^{2})^{2}}",
               "\\dfrac{2+6\\theta^{2}}{(1-\\theta^{2})^{3}}",
               "\\dfrac{24\\theta(1+\\theta^{2})}{(1-\\theta^{2})^{4}}"),
  dh_latex = c("1-t^{2}", "-2t(1-t^{2})", "-2+8t^{2}-6t^{4}", "16t-40t^{3}+24t^{5}"),
  g = function(th) atanh(th),
  dg = list(function(th) 1 / (1 - th^2),
            function(th) 2 * th / (1 - th^2)^2,
            function(th) (2 + 6 * th^2) / (1 - th^2)^3,
            function(th) 24 * th * (1 + th^2) / (1 - th^2)^4),
  h = function(e) tanh(e),
  dh = list(function(e) { t <- tanh(e); 1 - t^2 },
            function(e) { t <- tanh(e); -2 * t * (1 - t^2) },
            function(e) { t <- tanh(e); -2 + 8 * t^2 - 6 * t^4 },
            function(e) { t <- tanh(e); 16 * t - 40 * t^3 + 24 * t^5 }),
  grid_theta = seq(-0.9, 0.9, length.out = 25),
  grid_eta = seq(-2.5, 2.5, length.out = 25),
  note = "The link for a correlation. It is the logit of $(\\theta+1)/2$ up to a
  factor of two, which is why the inverse derivatives are the logistic ones
  rewritten in $t = 2p-1$."
)

# ------------------------------------------------------------------- power --
LINKS$power <- list(
  title = "Power",
  ctor = "power_link(lambda = 1.5)",
  obj = function() power_link(1.5),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in (0, \\infty)",
  g_latex = "g(\\theta) = \\theta^{\\lambda}, \\qquad \\lambda \\neq 0",
  h_latex = "h(\\eta) = \\eta^{1/\\lambda}",
  dg_latex = c("\\lambda\\,\\theta^{\\lambda-1}",
               "\\lambda(\\lambda-1)\\,\\theta^{\\lambda-2}",
               "\\lambda(\\lambda-1)(\\lambda-2)\\,\\theta^{\\lambda-3}",
               "\\lambda(\\lambda-1)(\\lambda-2)(\\lambda-3)\\,\\theta^{\\lambda-4}"),
  dh_latex = c("k\\,\\eta^{k-1}", "k(k-1)\\,\\eta^{k-2}",
               "k(k-1)(k-2)\\,\\eta^{k-3}", "k(k-1)(k-2)(k-3)\\,\\eta^{k-4}"),
  g = function(th) th^1.5,
  dg = local({
    lam <- 1.5
    lapply(1:4, function(k) function(th) .falling(lam, k) * th^(lam - k))
  }),
  h = function(e) e^(1 / 1.5),
  dh = local({
    kk <- 1 / 1.5
    lapply(1:4, function(j) function(e) .falling(kk, j) * e^(kk - j))
  }),
  grid_theta = seq(0.3, 4, length.out = 25),
  grid_eta = seq(0.3, 6, length.out = 25),
  note = "Shown at $\\lambda = 1.5$; $k = 1/\\lambda$ throughout. The falling
  factorial makes every order a one-liner. Note the trap this family sets in R:
  when $\\lambda = 2$ the second derivative is the constant
  $\\lambda(\\lambda-1)\\theta^{0}$, and `NA^0` evaluates to `1`, so a missing
  parameter would silently become a number. The package propagates missingness
  explicitly through `na_from()`."
)

# ---------------------------------------------------------------- softplus --
LINKS$softplus <- list(
  title = "Softplus",
  ctor = "softplus_link(a = 1)",
  obj = function() softplus_link(1),
  theta_dom = "\\theta \\in (0, \\infty)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\tfrac{1}{a}\\log\\!\\left(e^{a\\theta}-1\\right)",
  h_latex = "h(\\eta) = \\tfrac{1}{a}\\log\\!\\left(1+e^{a\\eta}\\right)",
  dg_latex = c("\\dfrac{v+1}{v}", "-\\dfrac{a(v+1)}{v^{2}}",
               "\\dfrac{a^{2}(v+1)(v+2)}{v^{3}}",
               "-\\dfrac{a^{3}(v+1)(v^{2}+6v+6)}{v^{4}}"),
  dh_latex = c("p", "a\\,p(1-p)", "a^{2}p(1-p)(1-2p)",
               "a^{3}p(1-p)(1-6p+6p^{2})"),
  g = function(th) log(expm1(1 * th)) / 1,
  dg = list(function(th) { v <- expm1(th); (v + 1) / v },
            function(th) { v <- expm1(th); -1 * (v + 1) / v^2 },
            function(th) { v <- expm1(th); 1 * (v + 1) * (v + 2) / v^3 },
            function(th) { v <- expm1(th); -1 * (v + 1) * (v^2 + 6 * v + 6) / v^4 }),
  h = function(e) log1p(exp(e)),
  dh = list(function(e) .sigma(e),
            function(e) { p <- .sigma(e); p * (1 - p) },
            function(e) { p <- .sigma(e); p * (1 - p) * (1 - 2 * p) },
            function(e) { p <- .sigma(e); p * (1 - p) * (1 - 6 * p + 6 * p^2) }),
  grid_theta = seq(0.1, 4, length.out = 25),
  grid_eta = seq(-3, 4, length.out = 25),
  note = "Shown at $a = 1$, with $v = e^{a\\theta}-1$ and $p = \\sigma(a\\eta)$.
  A smooth positivity constraint that, unlike the log link, is asymptotically the
  identity: $h(\\eta) \\to \\eta$ as $\\eta \\to \\infty$, so a large parameter is
  not compressed. The inverse derivatives are the logistic ones scaled by
  $a^{k-1}$."
)

# --------------------------------------------------- bounded: doubly bound --
LINKS$bounded_both <- list(
  title = "Bounded, two-sided",
  ctor = "bounded_link(lwr = -2, upr = 5)",
  obj = function() bounded_link(lwr = -2, upr = 5),
  theta_dom = "\\theta \\in (\\ell, u)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log\\dfrac{q}{1-q}, \\qquad q = \\dfrac{\\theta-\\ell}{W},\\; W = u-\\ell",
  h_latex = "h(\\eta) = \\ell + W\\,p, \\qquad p = \\sigma(\\eta)",
  dg_latex = c("\\dfrac{1}{W\\,q(1-q)}", "\\dfrac{2q-1}{W^{2}q^{2}(1-q)^{2}}",
               "\\dfrac{1}{W^{3}}\\left(\\dfrac{2}{q^{3}}+\\dfrac{2}{(1-q)^{3}}\\right)",
               "\\dfrac{1}{W^{4}}\\left(-\\dfrac{6}{q^{4}}+\\dfrac{6}{(1-q)^{4}}\\right)"),
  dh_latex = c("W\\,p(1-p)", "W\\,p(1-p)(1-2p)",
               "W\\,p(1-p)(1-6p+6p^{2})",
               "W\\,p(1-p)(1-14p+36p^{2}-24p^{3})"),
  g = function(th) { W <- 7; q <- (th + 2) / W; log(q / (1 - q)) },
  dg = list(function(th) { W <- 7; q <- (th + 2) / W; 1 / (W * q * (1 - q)) },
            function(th) { W <- 7; q <- (th + 2) / W; (2 * q - 1) / (W^2 * q^2 * (1 - q)^2) },
            function(th) { W <- 7; q <- (th + 2) / W; (2 / q^3 + 2 / (1 - q)^3) / W^3 },
            function(th) { W <- 7; q <- (th + 2) / W; (-6 / q^4 + 6 / (1 - q)^4) / W^4 }),
  h = function(e) -2 + 7 * .sigma(e),
  dh = list(function(e) { p <- .sigma(e); 7 * p * (1 - p) },
            function(e) { p <- .sigma(e); 7 * p * (1 - p) * (1 - 2 * p) },
            function(e) { p <- .sigma(e); 7 * p * (1 - p) * (1 - 6 * p + 6 * p^2) },
            function(e) { p <- .sigma(e); 7 * p * (1 - p) * (1 - 14 * p + 36 * p^2 - 24 * p^3) }),
  grid_theta = seq(-1.6, 4.6, length.out = 25),
  grid_eta = seq(-4, 4, length.out = 25),
  note = "Shown at $\\ell = -2$, $u = 5$, so $W = 7$. This is the logit link
  conjugated by the affine map that sends $(\\ell, u)$ to $(0,1)$; every formula
  is the logit one with $\\theta \\mapsto q$ and a factor $W^{-k}$ in the forward
  direction, $W$ in the inverse."
)

# ---------------------------------------------------- bounded: lower bound --
LINKS$bounded_lower <- list(
  title = "Bounded below",
  ctor = "bounded_link(lwr = 2)",
  obj = function() bounded_link(lwr = 2),
  theta_dom = "\\theta \\in (\\ell, \\infty)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log(\\theta - \\ell)",
  h_latex = "h(\\eta) = \\ell + e^{\\eta}",
  dg_latex = c("\\dfrac{1}{\\theta-\\ell}", "-\\dfrac{1}{(\\theta-\\ell)^{2}}",
               "\\dfrac{2}{(\\theta-\\ell)^{3}}", "-\\dfrac{6}{(\\theta-\\ell)^{4}}"),
  dh_latex = c("e^{\\eta}", "e^{\\eta}", "e^{\\eta}", "e^{\\eta}"),
  g = function(th) log(th - 2),
  dg = list(function(th) 1 / (th - 2), function(th) -1 / (th - 2)^2,
            function(th) 2 / (th - 2)^3, function(th) -6 / (th - 2)^4),
  h = function(e) 2 + exp(e),
  dh = list(exp, exp, exp, exp),
  grid_theta = 2 + exp(seq(log(0.05), log(20), length.out = 25)),
  grid_eta = seq(-3, 3, length.out = 25),
  note = "Shown at $\\ell = 2$: the log link translated. The translation is
  invisible to the inverse derivatives, since $\\partial_\\eta(\\ell + e^{\\eta})$
  does not see the constant."
)

# ---------------------------------------------------- bounded: upper bound --
LINKS$bounded_upper <- list(
  title = "Bounded above",
  ctor = "bounded_link(upr = 3)",
  obj = function() bounded_link(upr = 3),
  theta_dom = "\\theta \\in (-\\infty, u)",
  eta_dom = "\\eta \\in \\mathbb{R}",
  g_latex = "g(\\theta) = \\log(u - \\theta)",
  h_latex = "h(\\eta) = u - e^{\\eta}",
  dg_latex = c("-\\dfrac{1}{u-\\theta}", "-\\dfrac{1}{(u-\\theta)^{2}}",
               "-\\dfrac{2}{(u-\\theta)^{3}}", "-\\dfrac{6}{(u-\\theta)^{4}}"),
  dh_latex = c("-e^{\\eta}", "-e^{\\eta}", "-e^{\\eta}", "-e^{\\eta}"),
  g = function(th) log(3 - th),
  dg = list(function(th) -1 / (3 - th), function(th) -1 / (3 - th)^2,
            function(th) -2 / (3 - th)^3, function(th) -6 / (3 - th)^4),
  h = function(e) 3 - exp(e),
  dh = list(function(e) -exp(e), function(e) -exp(e),
            function(e) -exp(e), function(e) -exp(e)),
  grid_theta = 3 - exp(seq(log(0.05), log(20), length.out = 25)),
  grid_eta = seq(-3, 3, length.out = 25),
  note = "Shown at $u = 3$. Every derivative is negative: the link is decreasing,
  which is admissible -- what a link must be is a $C^4$ *bijection* onto
  $\\mathbb{R}$ with a non-vanishing derivative, not an increasing one. Note the
  sign pattern differs from the log link: all four forward derivatives share the
  same sign here, because each differentiation contributes both a factor $-1$
  from the chain rule and the alternating sign of the power rule."
)


# ---------------------------------------------------------------------------
# Certification machinery for the catalogue.
#
# Three certificates, answering three different questions, assembled into one
# table per link:
#
#   book       the transcription above, against the method the package
#              dispatches. Answers: is the printed formula the one that runs?
#   chain      order k of the package, against ONE numerical differentiation of
#              order k-1 of the package. Independent of the book entirely.
#   ift        the inverse-link derivatives predicted from the forward ones by
#              the inverse function theorem, against the inverse-link
#              derivatives the package implements separately. Pure algebra, so
#              exact to machine precision, and it fails if either family is
#              wrong.
# ---------------------------------------------------------------------------

QUANTITIES <- c("$g$", "$g'$", "$g''$", "$g'''$", "$g''''$",
                "$h$", "$h'$", "$h''$", "$h'''$", "$h''''$")

certify_link <- function(rec) {
  lk <- rec$obj()
  th <- rec$grid_theta
  et <- rec$grid_eta
  # numDeriv's Richardson stencil reaches about 8e-4 * |x|; a stencil point
  # falling outside the domain comes back NaN and would report a failure for
  # code that is exact, so the numerical column stays away from the ends.
  th_in <- th[3:(length(th) - 2)]
  et_in <- et[3:(length(et) - 2)]

  book <- c(
    compare_on_grid(rec$g, function(t) linkfun(lk, t), th),
    vapply(1:4, function(k) compare_on_grid(rec$dg[[k]], function(t) linkderiv(lk, t, k), th),
           numeric(1)),
    compare_on_grid(rec$h, function(e) linkinv(lk, e), et),
    vapply(1:4, function(k) compare_on_grid(rec$dh[[k]], function(e) linkinvderiv(lk, e, k), et),
           numeric(1))
  )

  chain <- c(
    NA_real_,
    vapply(1:4, function(k) {
      tryCatch(rel_err(linkderiv(lk, th_in, k),
                       num_grad1(function(t) linkderiv(lk, t, k - 1L), th_in)),
               error = function(e) NA_real_)
    }, numeric(1)),
    NA_real_,
    vapply(1:4, function(k) {
      tryCatch(rel_err(linkinvderiv(lk, et_in, k),
                       num_grad1(function(t) linkinvderiv(lk, t, k - 1L), et_in)),
               error = function(e) NA_real_)
    }, numeric(1))
  )

  # theta = h(eta), so the forward derivatives are evaluated where the theorem
  # says they should be.
  ph <- linkinv(lk, et)
  g1 <- linkderiv(lk, ph, 1); g2 <- linkderiv(lk, ph, 2)
  g3 <- linkderiv(lk, ph, 3); g4 <- linkderiv(lk, ph, 4)
  predicted <- list(
    1 / g1,
    -g2 / g1^3,
    (3 * g2^2 - g1 * g3) / g1^5,
    (-15 * g2^3 + 10 * g1 * g2 * g3 - g1^2 * g4) / g1^7
  )
  ift <- c(rep(NA_real_, 6), vapply(1:4, function(k) {
    tryCatch(rel_err(predicted[[k]], linkinvderiv(lk, et, k)), error = function(e) NA_real_)
  }, numeric(1)))

  data.frame(quantity = QUANTITIES, book = book, chain = chain, ift = ift,
             stringsAsFactors = FALSE)
}

# Format one certification table for printing. `tol` is applied per column:
# the algebraic columns must agree to machine precision, the numerical one only
# to what one Richardson extrapolation can deliver.
format_link_certificate <- function(cert, caption,
                                    tol = c(book = 1e-10, chain = 1e-5, ift = 1e-10)) {
  worst <- pmax(
    ifelse(is.na(cert$book), 0, cert$book / tol[["book"]]),
    ifelse(is.na(cert$chain), 0, cert$chain / tol[["chain"]]),
    ifelse(is.na(cert$ift), 0, cert$ift / tol[["ift"]])
  )
  out <- data.frame(
    quantity = cert$quantity,
    `book vs package` = fmt_err(cert$book),
    `vs numerical` = fmt_err(cert$chain),
    `via inverse fn thm` = fmt_err(cert$ift),
    verdict = ifelse(worst < 1, "agree", "DISAGREE"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  knitr::kable(out, caption = caption, align = "lrrrl", row.names = FALSE)
}

fmt_err <- function(x) ifelse(is.na(x), "", formatC(x, format = "e", digits = 1))

# Render one catalogue entry as markdown. Called from a results='asis' chunk, so
# everything it cats becomes part of the document.
render_link_entry <- function(id, rec) {
  cat(sprintf("\n### %s {#sec-link-%s}\n\n", rec$title, gsub("_", "-", id)))
  cat(sprintf("``` r\n%s\n```\n\n", rec$ctor))
  cat(sprintf("Domain: $%s$, so that $%s$.\n\n",
              sub("^\\\\", "\\\\", rec$theta_dom), rec$eta_dom))
  cat(sprintf("$$%s, \\qquad %s$$\n\n", rec$g_latex, rec$h_latex))

  tab <- data.frame(
    k = 1:4,
    forward = paste0("$", rec$dg_latex, "$"),
    inverse = paste0("$", rec$dh_latex, "$"),
    stringsAsFactors = FALSE
  )
  names(tab) <- c("$k$", "$g^{(k)}(\\theta)$", "$h^{(k)}(\\eta)$")
  print(knitr::kable(tab, align = "crr", row.names = FALSE,
                     caption = sprintf("%s: derivatives to fourth order", rec$title)))
  cat("\n\n")
  cat(gsub("\n\\s+", "\n", rec$note), "\n\n")
}

# Silent consistency gate. Run from a hidden chunk; if any formula printed in
# the chapter disagrees with the implementation, or the implementation with
# itself, this stops the render. The reader never sees it -- the check is ours.
assert_links_ok <- function(tol_algebraic = 1e-10, tol_numeric = 1e-5) {
  for (id in names(LINKS)) {
    cert <- certify_link(LINKS[[id]])
    bad <- (!is.na(cert$book) & cert$book > tol_algebraic) |
           (!is.na(cert$chain) & cert$chain > tol_numeric) |
           (!is.na(cert$ift) & cert$ift > tol_algebraic) |
           is.na(cert$book)
    if (any(bad)) {
      stop(sprintf("Link '%s': printed formulas disagree with the package (%s).",
                   id, paste(cert$quantity[bad], collapse = ", ")), call. = FALSE)
    }
  }
  invisible(TRUE)
}
