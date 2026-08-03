# statmodels7

Read this first. It is the orientation file for the whole project, and it is the only
place where the accumulated knowledge lives: Claude Code's memory is keyed by directory
path, and these packages were moved here on 2026-07-22, so nothing from the earlier
sessions carries over automatically.

`statmodels7/` is a plain directory, not a repository — the package repositories sit
inside it. This file, `logo/make-logos.R` and the package plans (`piano_*.txt`, the
decision records drafted before a package's code exists) are therefore versioned in
the portal repository (`site/`, which is `statmodels7.github.io`) and copied back up
here, so a fresh clone can recover them. If you edit any of them, copy the change
into `site/` and commit it there; `site/sync-stack-files.sh` does exactly that.

---

## 0. How to talk about it

Get this right before writing any prose for a README, a site or a paper. The
positioning is **not** "a framework with exact derivatives and optimisation on the
unconstrained scale" — that is a description of the machinery, and it buries the point.

The point is that **these things get rewritten in every package that needs them**.
Almost every R modelling package carries its own distributions and link functions,
written as internal helper functions — a `switch` on a character string, closures
private to the package that owns them. They are not objects, so nothing outside can
reuse them, extend them, or ask them for anything their author did not happen to need.
And since each package writes only what it needs, that is usually the density and the
score, sometimes the Hessian, and nothing beyond.

A Gamma is a fixed mathematical object. It should be written once, correctly, with
everything a modelling routine could want already computed — third and fourth
derivatives, derivatives with respect to the response, derivatives with respect to the
unconstrained parameters — and then everyone builds on top. That is what this toolkit is:
infrastructure to be reused, not another modelling package.

Extension matters just as much in the telling: a toolkit that only worked for the fourteen
distributions it ships with would have solved nothing.

The same argument is what `optimizers7` is for, and it is worth having the sentence
ready. Almost every R package that fits a model carries its own optimiser, written
inside the function that needs it — a loop, a `while (!converged)`, a tolerance compared
against whatever quantity the author had to hand. The stopping rule, which decides what
the whole thing means by *finished*, is a number buried three levels down. Lead with
that, not with the list of algorithms.

Giovanni corrected an earlier draft that led with the link scale and the optimisation
story. Lead with the reuse argument.

## 1. What this is

`statmodels7` is an umbrella for a statistical modelling toolkit built entirely on the
**S7** object system. The naming convention is simply that **every package name ends
in 7**, the 7 being S7 — `linkfunctions7`, `distributions7`, `statmodels7`. Until
2026-08-03 the convention also required a plural noun; Giovanni dropped that when
naming `basis7` ("basis7 è più bello di bases7"), so the plural is history, not a
rule — do not "fix" singular names, and do not cite the plural reading in new prose.

**Built so far**

| package | what it provides |
|---|---|
| `linkfunctions7` | 16 link classes (14 constructors) with exact analytical derivatives to 4th order, both directions, plus numerical fallbacks for user-defined links |
| `distributions7` | 14 distributions with exact score, information and 3rd/4th derivatives, plus wrappers, transformations, MLE |
| `optimizers7` | 11 algorithms as objects — newton, bfgs, lbfgs, cg, bb, gd, adam, nelder_mead, compass, bundle, multistart — with composable stopping rules, self-reporting safeguards, box bounds removed by reparametrisation, starting values that need not be written out, and multistart parallel by default |
| `basis7` | bases as objects: evaluation, derivatives of any order, the integral anchored at the lower endpoint, and exact Gram matrices against a choice of measure. B-splines, Fourier and Legendre; one `TransformedBasis` wrapper for orthonormalisation, constraints and the Demmler-Reinsch construction; `tensor_basis()` for several variables, with `basis_contract()` computing what a fit needs without forming the product; numerical fallbacks make an evaluation-only basis complete |

**Planned** — `modelterms7`, `penalties7`, and eventually the `statmodels7`
package itself, which assembles everything into models. That last one is the destination:
a GAMLSS-like framework but far more organised, where

- every distribution parameter can be modelled, not just the mean;
- optimizers are pluggable, and the user chooses frequentist or Bayesian;
- users can define new model terms and drop them straight into formulas —
  `bs()`, `ps()`, `seg()`, `jump()`, `jseg()`, `ridge()`, `lasso()`, `glasso()`;
- richer specifications compose naturally, e.g. `gas(p = 1, q = 2, by = ~(1|id))` for a
  generalized autoregressive score model on panel data with population-average
  coefficients and random effects.

Everything below `statmodels7` exists to make that possible. The design pressure
throughout has been: *exact derivatives, on the unconstrained scale, for arbitrary
user-defined components.*

---

## 2. Where things are

```
C:\Users\giova\Desktop\labstatr\statmodels7\
    linkfunctions7\
    distributions7\
    optimizers7\
    basis7\
    book\            the Quarto book (see §9); `quarto render` from inside it
    logo\            hex logos: make-logos.R draws them, run from this directory
    site\            the portal, its own repository (statmodels7.github.io)
    articoli\        papers Giovanni drops in as design input (LaTeX sources).
                     LOCAL ONLY — copyrighted material, never committed to any
                     public repository. Entries: Rügamer, "Scalable Higher-Order
                     Tensor Product Spline Models" (AISTATS 2024), the factorized
                     tensor-basis reference for piano_basis7 §3b; Bach & Klein,
                     "Bayesian Effect Selection in Additive Models with an
                     Application to Time-to-Event Data", whose §2 (Demmler-Reinsch
                     basis for P-splines) drives piano_basis7 §3c
    piano_*.txt      package plans, drafted before code; carried in site/stack/
```

`book/` is versioned through `site/` (source under `stack/book/`, rendered HTML
under `book/`, both carried by `sync-stack-files.sh`), and the rendered book is
**published at `statmodels7.github.io/book/`** (since 2026-07-26). After editing
the book: `quarto render` in `book/` **from PowerShell** (see 3), then
`sh sync-stack-files.sh` in `site/` and commit there — the render is a manual step
because it executes R against the working tree, and it takes about twenty minutes
now that there are four package chapters.

GitHub: `github.com/statmodels7/{linkfunctions7,distributions7,optimizers7,basis7}`, all on `master`. The repositories were
transferred from `giovannitinervia9/*` on 2026-07-22; GitHub keeps redirects, so
`install_github("giovannitinervia9/distributions7")` still resolves.

Websites, all live and rebuilt by a `pkgdown.yaml` workflow on every push:

| | |
|---|---|
| portal | `statmodels7.github.io` — plain HTML, own repository, deploys from `main` |
| linkfunctions7 | `statmodels7.github.io/linkfunctions7` — pkgdown, from `gh-pages` |
| distributions7 | `statmodels7.github.io/distributions7` — pkgdown, from `gh-pages` |
| optimizers7 | `statmodels7.github.io/optimizers7` — pkgdown, from `gh-pages` |
| basis7 | `statmodels7.github.io/basis7` — pkgdown, from `gh-pages` |

Predecessors, kept for reference only: `labstatr\distrib` (S3, ~8500 lines) and
`labstatr\linkfunctions`. `distributions7` is a rewrite of `distrib`, not a port —
several formulas in the original are wrong (see §7).

---

## 3. Toolchain

**R is not on PATH.** Use the absolute path:

```
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
```

(The 4.5.1 directory exists but is empty.) Rtools/gcc 14.2 is available and Rcpp
compiles fine.

```bash
# run the test suite (loads from source via pkgload, no install needed)
Rscript.exe -e "setwd('<pkg>'); testthat::test_local()"

# regenerate documentation
Rscript.exe -e "setwd('<pkg>'); roxygen2::roxygenise()"

# recompile the Rcpp interface after adding a new .cpp
Rscript.exe -e "setwd('<pkg>'); Rcpp::compileAttributes()"

# install
Rscript.exe -e "setwd('<pkg>'); install.packages('.', repos=NULL, type='source')"

# full check (installs into its own directory, unaffected by library locks)
R.exe CMD build <pkg> --no-manual
# --as-cran, because that is what the CI action runs and warnings fail it
R.exe CMD check <pkg>_0.0.0.9000.tar.gz --no-manual --as-cran
```

**Check with `--as-cran`, and read WHICH warning.** The CI action runs
`--as-cran` and treats any warning as a failure, so a local check without it
cannot see a whole class of defect: adding a test that called
`pkgdown::check_pkgdown()` without declaring pkgdown in `Suggests` gave
*"'::' import not declared"* on all five platforms after a local
`Status: OK`.

The converse also holds and is worth knowing before chasing a phantom.
`--as-cran` locally reports a **CRAN incoming feasibility** warning that CI
never will, because that check needs the network and the action disables it.
`distributions7` shows one here for its `Remotes` field while its CI is green
and its workflow is identical to the others'. A local warning is a red push only
if it is not that one.

Vignettes need pandoc: `$env:RSTUDIO_PANDOC="C:\Users\giova\AppData\Local\Programs\Quarto\bin\tools"`.

### Windows traps that cost time

- **`install.packages` fails with "cannot remove earlier installation, is it in use?"**
  when an RStudio session has the package loaded — it holds the DLL. Ask Giovanni to run
  `detach("package:distributions7", unload = TRUE)` or restart his R session.
  `testthat::test_local()` is unaffected: pkgload compiles and loads from source.
- **PowerShell's `Set-Location` does not change the process's current directory**, only
  the provider location. To release a directory before moving it you also need
  `[System.IO.Directory]::SetCurrentDirectory(...)` — and even then, the agent session's
  own working directory stays locked. The workaround that always works is to move a
  directory's *contents* file by file rather than the directory itself.
- PowerShell here-strings and `-e` scripts with `$`, `^`, backticks get mangled. Write
  R scripts to a file and run the file. `^` is PowerShell's escape character, so
  `git rev-parse 'HEAD^{tree}'` must be run from bash.
- **`C:\Users\giova\.Rprofile` is a 5-byte BOM-only file**, and every R launched from
  Git Bash dies on it with `Error: 1:1: unexpected input`. Run R and anything that
  shells out to R — `quarto render` above all — **from PowerShell**. Quarto reports it
  as *"Problem with running R ... Please check your installation of R"*, which points
  nowhere near the cause. Deleting the file would fix it properly; it has not been
  done because it is in Giovanni's home directory, not the project.
- **Escapes collapse one level when a script is passed through the shell.** A Python
  heredoc written with `\\n` arrives as a newline, `\\code` as `\code` (with a
  SyntaxWarning that is the only clue), and `\\rVert` as a carriage return inside the
  file. A `sed 's/\\pi/.../g'` loses the backslash and matches bare `pi`, which
  silently rewrote *piecewise*, *stopping*, *epigraph* and five other words in a book
  chapter. **Anything containing a backslash goes into a file with the Write tool and
  is run from there** — the same rule as the PowerShell one above, for the same reason.
- **Run `Rcpp::compileAttributes()` BEFORE `roxygen2::roxygenise()`** after adding a
  `.cpp`. Running roxygenise alone dropped `RcppExports.R` from the `Collate` field,
  and the package then failed to install with *"files ... missing from 'Collate'"*.

---

## 4. Architecture

### linkfunctions7

An S7 class `link` (properties: `link_name`, `link_bounds`, `link_params`) with one
subclass per link. Every link implements ten generics:

```
linkfun / linkinv                        g(theta), g^-1(eta)
dlinkfun  … d4linkfun                    derivatives of g wrt theta
dlinkinv  … d4linkinv                    derivatives of g^-1 wrt eta
```

`linkderiv(x, theta, order)` and `linkinvderiv(x, eta, order)` are convenience routers
over those. **They dispatch twice** — once on themselves, then on the order-specific
generic — so performance-critical code calls `dlinkinv()` and friends directly;
`distributions7::inverse_link_derivs` does exactly this.

`check_link(x)` is the built-in validator: invertibility both ways, monotonicity, the
inverse function theorem, and all derivative orders against `numDeriv`.

**Numerical fallbacks** (added 2026-07-30). Only `linkfun` and `linkinv` are
compulsory; the base class supplies the eight derivative generics a link does not
implement, so a user-defined link runs immediately — the same bargain distributions7
offers a density-only distribution. The rule they obey is the one `check_link()`
already enforced: each fallback finds the highest order `m` implemented analytically
and applies **one** central stencil of order `k - m` to it, never a chain of first
differences. It is worth several digits: on a log link defined with nothing, the
fourth derivative is good to ~2e-5, and supplying just the analytic first and second
takes it to ~2e-8.

Consequence for the validator, and the reason it changed: an order that comes from a
fallback would be compared against a numerical differentiation of the order below —
the same arithmetic twice, agreeing however wrong the link is. `check_link()`
therefore leaves such orders `NA`, meaning **not checked**, and prints
`[numerical]` or `[PASSED to order k, n numerical]` instead of a bare `[PASSED]`.
`link_fallback_orders(x)` answers the same question programmatically, and the result
of `check_link()` carries it as the attribute `"analytic_orders"` — an attribute
rather than a seventh list element, so that every element of the result stays a
logical vector and callers can still reduce over it.

14 constructors: identity, log, logit, probit, cloglog, loglog, cauchit, rhobit, sqrt,
inverse, inverse_sq, power(lambda), softplus(a), bounded(lwr, upr). The last one returns
a lower-, upper- or doubly-bounded link depending on which arguments are given, and the
identity link when neither is — so 16 distinct classes in total.

### distributions7

Base class `distrib` → `continuous_distrib` / `discrete_distrib` → one class per
distribution. Each distribution file defines the S7 class, registers methods, and exports
a constructor wrapper (`gaussian_distrib()`, `negbin_distrib()`, …) taking link functions
from linkfunctions7.

**Parameters** travel as a named list `theta`. Every generic normalises it through the
internal `align_theta()`, which reorders by name, strips stray names off the values, and
validates against `params_bounds` treated as **open** intervals.

**Generics**

```
distrib_pdf / cdf / quantile / rng
distrib_gradient / hessian / expected_hessian      scale = "parameter" | "link"
distrib_deriv3 / deriv4                            + expected =, approx =, nsim =
distrib_grad_y / hess_y                            derivatives wrt the response
distrib_cross_y                                    mixed d2l/dy dtheta_i, scale = as above
expectation, mean, variance, std_dev, skewness, kurtosis, moment
```

**A distribution needs only `distrib_pdf`.** Everything else has a numerical fallback:
cdf by quadrature, quantile by root-finding, RNG by Generalized Ratio-of-Uniforms,
derivatives by finite differences. See `vignettes/defining-a-distribution.Rmd`.

**14 distributions**: gaussian, cauchy, logistic, student_t, laplace, pseudohuber, gamma,
invgauss, lognormal, beta, bernoulli, binomial, poisson, negbin. All have closed-form
observed derivatives to 4th order (Rcpp kernels in `src/*_hd.cpp`).

**Wrappers**: `zero_inflated()`, `zero_adjusted()`, `truncated()`, `transformation()`
(12 transformers: log, exp, sqrt, inverse, power, box-cox, yeo-johnson, softplus, asinh,
logit, expit, affine), and `fixed()` (2026-08-03), which holds parameters at known
values -- the only wrapper that REMOVES parameters. It derives nothing: every method
splices the fixed values into theta and delegates, and the derivatives are the parent's
components among the free parameters, subset BY NAME from the same enumeration
(`deriv_names`/`hess_names` on the free set produce the parent's strings because the
free set preserves the parent's order -- never re-parse a component name). `n_params = 0`
is legal (a fully known prior); fixing a wrapper's own parameter is legal and useful
(`fixed(zero_inflated(d), zi = 0.3)`); `fixed()` of `fixed()` collapses. Built as the
prerequisite of `penalties7`: ridge is `fixed(gaussian_distrib(), mu = 0)` with the
scale free. One S7 trap found writing its `print`: a method registered on a BASE
generic is an S3 method, so the parent's is reached with `NextMethod()` --
`S7::super()` only works inside S7 generics and `S7::method()` only retrieves from
them. All propagate `params_smooth`; all validate their parent in the
constructor (see §7). `zero_adjusted()` of a **continuous** parent is a **mixed**
distribution — a density plus an atom — and declares that atom through
`distrib_atoms()`, which `check_distrib()` consults.

`truncated(distrib, lower, upper)` restricts to an interval, either side optional, both
endpoints **included**. It is the odd one out: it adds **no parameter** (the endpoints
are constants, like a binomial's `size`), and instead adds a θ-dependent normalising
constant `Z(θ) = F(u) − F(ℓ⁻)`. With `m_i = E_T[s_i]` and `M_ij = E_T[H_ij + s_i s_j]`,
both expectations under the *truncated* law,

```
d_i   l_T = s_i(y) - m_i                      (the parent's score, recentred)
d_ij  l_T = H_ij(y) - M_ij + m_i m_j
E[d_ij l_T] = -Cov_T(s_i, s_j)                (2nd Bartlett identity)
```

`m` and `M` have no closed form in general and go through `expectation()`, so truncated
derivatives are much dearer than the parent's.

**Orders 3 and 4 are closed form for every wrapper** (`R/wrapper_derivatives.R`,
2026-07-27). Every wrapper's log-likelihood is the parent's log-density plus, or instead
of, `log L` for some θ-dependent `L`, so two partition sums cover all of them:

- `d^I f / f = Σ_π Π_B l^(B)` — the complete Bell polynomial, i.e. the Bartlett lemma
  read backwards;
- `d^I log L = Σ_π (-1)^{|π|-1}(|π|-1)! Π_B (d^B L / L)` — the moment-to-cumulant
  relation. Only the **ratios** are needed, never `L`'s own derivatives.

Orders 1–2 fall out as special cases and reproduce the hand-written closed forms exactly
(`w0(1-w0)`, `C_ij`, `M_ij - m_i m_j`), so the two derivations check each other. `L0` is
**affine in zi**, which kills every block containing two or more zi's. For truncation each
distinct block costs one quadrature, memoised across the partition sum.

Two traps found while writing it, both worth keeping:
- **`$` on a list does partial matching.** `d4[["mu_mu_mu"]]` is NULL but `d4$mu_mu_mu`
  silently resolves to `mu_mu_mu_mu` when that prefix is unique. Always `[[ ]]` in tests.
- **`deriv_index_list()`'s order-2 case is ordered for `hess_names()`** (diagonal first)
  while `deriv_names()` is lexicographic; pairing them would label `"mu_sigma"` with the
  index `(sigma, sigma)`. Orders 3–4 agree, so the bug is invisible until someone reuses
  the helper. `order_indices()` in `wrapper_derivatives.R` generates its own.

**Link scale.** `scale = "link"` gives derivatives with respect to the unconstrained
parameters, to 4th order, via Faà di Bruno with a diagonal Jacobian and partial Bell
polynomials. The transformation is written **once in the generic body**, intercepting
`S7_dispatch()`; methods always return the parameter scale.
⚠️ **Breaking for user-written derivative methods**: they must declare
`scale = c("parameter", "link"), ...` in their signature or dispatch fails — S7 requires
a method's formals to include the generic's named arguments.

**Expected derivatives.** When there is no closed form, `approx` selects the strategy:
`"bartlett"` (default at order 2), `"integrate"`, `"mc"`. The unifying idea is that
**OPG is the order-2 Bartlett identity**; the general identity is implemented by
enumerating set partitions, so any order follows from products of lower orders.

**CDF derivatives** (`R/cdf_derivatives.R`, 2026-07-27): `distrib_grad_cdf()` /
`distrib_hess_cdf()`, with `lower.tail` and `log` arguments. This is what makes
**censored likelihoods** possible — right, left and interval — and what a quantile
residual's delta-method SE needs. The governing identity is one exchange of derivative
and integral, the region not depending on θ:

```
d^I F(q) / F(q) = E[ d^I f / f | Y <= q ]
```

a *partial* expectation of the same complete-Bell quantity the wrappers use. Two
implementations follow from it: for a **discrete** parent the expectation is a finite
sum, so the identity is exact and is used directly; for a **continuous** one it is a
semi-infinite integral, and differencing the (analytic) cdf is both cheaper and more
accurate, so that is the fallback. **Closed-form coverage (2026-07-27).** 12 of 14 distributions have a closed-form
gradient; the exceptions are **gamma** and **beta**, and that is a mathematical fact,
not an omission — see below.

| family | route | note |
|---|---|---|
| gaussian, cauchy, logistic, laplace | closed, both orders | location-scale: `dF/dmu = -f`, `dF/dsigma = -z f`; 2nd order via `f` and `l_y` |
| lognormal | closed (grad) | location-scale on the log scale; `dF/dmu = -y f(y)`, `dF/ds2 = -y f z/(2 sigma)` |
| invgauss | closed (grad) | its cdf is elementary in `Phi`; the `exp(2/(phi mu))` factor must be combined with `Phi(b)` **on the log scale** or it overflows |
| student_t, pseudohuber | closed in `(mu, sigma)`, FD in `nu` | for pseudohuber this is a real **accuracy** gain, not just speed: its cdf is itself a quadrature, so differencing it is poor (1e-6) while `-f` is exact |
| poisson | closed | the sum telescopes: `dF(k)/dmu = -f(k)` |
| binomial, bernoulli | closed | `dF(k)/dp = -n dbinom(k, n-1, p)` |
| negbin | closed in `mu`, exact sum in `theta` | `dF(k)/dmu = -f(k)(k+theta)/(theta+mu)`, → Poisson as `theta → inf` |
| **gamma, beta** | FD of the analytic cdf | **no elementary closed form exists** |

On gamma and beta: the shape direction is a derivative of the incomplete gamma/beta
*in its parameter*, which is hypergeometric — there is no elementary expression. For
the gamma the **rate** direction is elementary (`dF/dbeta = y f(y)/beta`), but the
package parametrises by `(mu, sigma2)` and both of those involve the shape, so neither
is closed form. What does survive is an exact identity, worth keeping as a test:

```
dF/dmu + (2 sigma2/mu) dF/dsigma2 = -y f(y) / mu
```

the shape direction cancelling from the combination. The FD fallback delivers ~2e-10
against the exact partial-expectation integral, so this is not a practical loss.

⚠️ **Test independence.** Because the implementations differ by kind, a test must use
the route the implementation does *not*: continuous-with-fallback against the
partial-expectation integral, discrete against FD of the cdf, closed forms against the
partial expectation. Checking discrete against the partial sum gives 0.00e+00 and
proves nothing — it is the same sum twice.

**Mixed response-parameter derivatives** (`R/cross_derivatives.R`, 2026-08-03):
`distrib_cross_y()` returns `d2 l / dy dtheta_i`, one component per parameter --
the off-diagonal block of the joint Hessian in `(y, theta)` whose diagonals are
`hess_y` and `hessian`. Built for `penalties7` (a penalty is a negative
log-density at the coefficients; joint estimation of coefficients and
hyperparameters needs this block, and so does a profiled objective's gradient
via the implicit function theorem -- `glmm_prova.R` wrote it by hand). Link
scale is the gradient's own first-order diagonal chain rule, reused verbatim:
the y-derivative does not interact with a reparametrisation of theta. The
fallback differentiates `distrib_grad_y` (one FD layer on an analytic response
gradient; when that is itself the fallback, the two differences act on
DIFFERENT variables and commute into the four-point mixed stencil -- cross-
variable composition is not the forbidden same-variable nesting). Closed forms:
gaussian and Student t. Truncation delegates exactly (log Z has no y); fixed()
subsets. Continuous only, like grad_y. Not yet in `check_distrib()` -- for a
user distribution both sides would be the same FD arithmetic; add the check
when a user-facing closed-form route exists to compare against.

**Truncation uses the cdf derivatives where they are exact.** `d^B Z = d^B F(U) -
d^B F(L^-)` replaces one quadrature per component with two calls on the parent, taking
the truncated Gaussian Hessian from ~8 ms to 0.85 ms. But the route is **gated**: it is
taken only when the parent has a genuine closed form, or is a lattice family (whose cdf
derivatives are an exact sum). The reason is a regression it caused when applied
unconditionally — with a FD-based parent the Hessian picks up ~1e-8 of noise, and
`numerical_deriv4()` differentiates that Hessian, so the *reference* of the fourth-order
check degrades and the check fails on correct code. Detecting "does the parent have a
closed form" uses the documented S7 trick: `attr(m, "signature")[[1]]` is the class the
method was registered on. `E[H_T] = -Cov_T(score)` still needs one quadrature per
component; the cdf derivatives cannot separate `E_T[l^(ij)]` from `E_T[l^(i)l^(j)]`.

**User-facing tools**: `check_distrib()` (thirteen numerical checks on a continuous
distribution, twelve on a discrete one; catches deliberately broken components),
`fit_distrib()` (MLE on the link scale, Fisher scoring / Newton / BFGS, delta-method SEs,
CIs built on the link scale and mapped back so they respect domains), plus `simulate()`
and `plot()` on the fit. Extractors: `coef()`, `vcov()`, `confint()` — all three
taking `scale = c("parameter", "link")` — and `logLik()`.

**`confint()` was missing until 2026-08-03**, and so was the interval in the
`print()` method's link-scale block, although `ci_eta` had been computed, stored
and documented from the start. The link-scale interval is the one actually
computed and the parameter-scale table is its image under `g^{-1}`, so printing
both is what makes the mapping legible; `confint()` recomputes at any `level`
from the stored estimate and standard error, without refitting.

**`fit_distrib()` delegates its optimisation to `optimizers7`** (2026-08-03), and
distributions7 now Imports it. Until then the function did the exact thing the
toolkit's positioning condemns in other packages: a hand-written scoring loop with
its own step halving, the convergence test spelled out inline, and `stats::optim`
for the BFGS branch — an *external* optimiser while the sibling package sat unused.
Three things made the translation nearly literal:

- **Fisher scoring is `newton()` with the expected information passed as `he`.**
  It needs no implementation of its own; `"fisher"` and `"newton"` differ only in
  that one argument.
- the inline test `max|U| < tol || |Δl| < tol(|l|+tol)` *is*
  `crit_any(crit_grad(tol), crit_rel_obj(tol))`;
- `method` now also accepts **any optimiser object**, used as given.

Two things improved rather than merely moved: the line search requires sufficient
decrease instead of mere non-decrease, and a non-PD Hessian is repaired by flooring
its eigenvalues instead of abandoning the start, which `solve()` used to force.
The BFGS fallback is kept for `"fisher"` and `"newton"` only — silently replacing
an optimiser the caller *chose* would report a fit obtained by a different method.

⚠️ **A catch-all for numerical failures must not swallow a configuration error.**
The restart loop's `tryCatch(..., error = function(e) NULL)` exists to absorb a
divergent quadrature in a numerically approximated expected information. It was
also absorbing optimizers7's *refusal* of a stopping rule the method cannot
evaluate, and reporting it as `"Optimisation failed from every starting value;
supply 'start'"` — naming the wrong cause entirely. `check_criterion()` is exported
by optimizers7 for exactly this and is now called before the loop. The general
shape: when a `tryCatch` is there for one class of failure, check what else reaches
it. This was found by writing the test for the new feature, not by using it.

---

## 5. Working preferences

- **No `Co-Authored-By: Claude` trailer** on commits. Giovanni is the sole author of what
  he publishes.
- **No code generators that require Mathematica.** Wolfram output is hand-transcribed
  into Rcpp and validated numerically. He rejected a generator explicitly: *"non ha senso"*
  that the package only builds because Mathematica happens to be installed.
- **Commit messages explain why, not what.** The reasoning in the history is deliberate —
  it is what stops a future reader from "simplifying" a decision that was measured.
- **Verify, do not assume.** Several times a suspected bug turned out to be a bad test,
  and several times a "cosmetic" finding turned out to be a real defect. Measure first.
- Report failures faithfully, with the numbers.

---

## 6. Current state

| | |
|---|---|
| `linkfunctions7` | 886 tests, `R CMD check` OK, CI green |
| `distributions7` | 1656 tests, `R CMD check` OK (2026-08-03, local), CI green |
| `optimizers7` | 660 tests, `R CMD check` OK with vignettes, CI green on all five platforms (2026-07-31). Published the same day; the Rcpp kernels had until then only ever been compiled by one compiler on one machine. |
| `basis7` | 682 tests, `R CMD check --as-cran` clean apart from the two environment notes, CI green (2026-08-03). Version `0.3.1`, a `NEWS.md` from the first commit and a vignette. Phases 1 to 4 of `piano_basis7.txt` are done; phase 5 is the handoff to `penalties7` and `modelterms7`. |

All three repositories run `R-CMD-check` on macOS, Windows and three Linux/R
combinations (devel, release, oldrel-1) plus a coverage workflow, all green. That matrix
matters for the two that ship Rcpp kernels, `distributions7` and `optimizers7`: until the
CI existed those had only ever been compiled by one compiler on one machine.

Vignettes: `defining-a-distribution`, `fitting-a-model`,
`derivatives-and-the-link-scale` (distributions7), `link-functions`
(linkfunctions7) and `extending-optimizers7` (optimizers7). All three packages have a
README with badges.

`R CMD check --as-cran` on distributions7 is clean apart from submission metadata:
new submission, version `0.0.0.9000`, the `Remotes:` field, and its two upstream
packages not being on CRAN. **The CRAN blocker is that BOTH linkfunctions7 and
optimizers7 must be published first** — since 2026-08-03 distributions7 imports
optimizers7 as well (see §4). The two are independent of each other, so they can
be submitted in parallel rather than in sequence.
Both names, and `modelterms7` / `basis7` / `penalties7`, are free on CRAN.
linkfunctions7 is now in the same state: two notes, both the dev version string and
pandoc missing from the check environment.

**Everything is documented, exported or not** (2026-07-30). Both packages had full
coverage of their exported surface and none at all of their internals — 19 objects in
linkfunctions7, 81 in distributions7. All of them now carry roxygen with
`@keywords internal` and *without* `@noRd`, so a page is generated and reachable
through `?name` and the site while staying out of the index. That is the convention
distributions7 had already chosen for its distribution classes; it is now the whole
toolkit's. `.onLoad` is the one exception and keeps `@noRd`, a man page for an R load
hook being noise CRAN would query.

Two things that only bite at submission and that `R CMD check` does not raise locally:
**`\value` on every exported topic** and **an executable example on every exported
function**. Both were missing throughout linkfunctions7 and are the two most common
reasons a first submission comes back. Closed there; worth checking before any future
package is submitted.

---

## 7. Hard-won knowledge

This section is the reason the file exists. Each item cost real time to discover.

### Mathematics

- **The old `distrib` package's NB2 excess kurtosis is wrong.** It has
  `6/θ + (θ+μ)/(μθ)`; the correct expression is `6/θ + θ/(μ(θ+μ))`. Verified numerically.
  Treat every formula inherited from `distrib` as unverified.
- **Wolfram's expected derivatives for the logistic are unusable**: they contain `PolyLog`,
  the imaginary unit, `ConditionalExpression`, and depend on `mu` — impossible for a
  location-scale family. The observed ones were rederived in the sigmoid instead:
  with `z = (y-mu)/sigma`, `t = plogis(z)`, `u = 1-t`, `l = -log(sigma) + g(z)` and
  `g1 = 1-2t`, `g2 = -2tu`, `g3 = g2*g1`, `g4 = g2*(1-6tu)`; then
  `d/dmu[A/sigma^k] = -A'/sigma^(k+1)` and `d/dsigma[A/sigma^k] = -(zA' + kA)/sigma^(k+1)`
  generate every component as a polynomial in z. Overflow-free to |z| = 4000.
  Seven of the nine expected components are known exactly (`E[l_μμσ] = 1/(2σ³)`,
  `E[l_σσσ] = (π²+2)/(2σ³)`, `E[l_μ⁴] = 1/(15σ⁴)`, rest zero by symmetry); the two open
  ones need `∫ wᵏ sech⁴w tanh²w dw` for k = 2, 4.
- **Non-regular models.** The Laplace location has a kink: its observed Hessian is
  degenerate while the Fisher information is `-1/b²`. Precision matters here, because a
  sloppy version of this note put a false claim in the book (caught by Giovanni,
  2026-07-26). The measured facts:
  - the **shipped** `laplace_distrib()` implements `distrib_expected_hessian` in closed
    form, so `approx` is **ignored** and bartlett/integrate/mc all print `-1/b²`;
  - on a **bare** Laplace (density + analytic score + analytic Hessian, no expected
    method), `"bartlett"` recovers `-1/b²`; `"integrate"` **and** `"mc"` both return
    exactly 0, because both average the observed `l_μμ`, which is 0 a.e. — mc is only
    "noisy" when the observed Hessian is itself the FD fallback smearing the kink;
  - the deeper point: averaging the second derivative gives `E[∂²l]`, which for this
    family genuinely *is* zero; what fails is its identification with `-I(θ)` (the
    second Bartlett identity). Only the score-based route survives.
  The base class carries `params_smooth` to record which parameters are differentiable.
  `fit_distrib` works on the Laplace *because* Fisher scoring uses the expected
  information. Terminology trap fixed in docs/vignette the same day: the score-variance
  identity `I = E[ss']` is the **second** Bartlett identity, not the first (the first is
  `E[s] = 0`); the information for a non-regular family is *defined* as `Var(score)`.
  When prose describes what code returns, **run the code first** — the book's hidden
  gates now pin these three claims.
- **Pseudo-Huber Bessel terms** are degree-homogeneous, so the exponentially scaled
  `bessel_k(x, nu, 2)` is exact and avoids overflow (verified to nu = 2000).

### The two zero wrappers

The derivative formulas in `zero_inflated.R` and `zero_adjusted.R` were re-derived and
confirmed correct (2026-07-26), including the neat collapse in the ZI mixed block, where
`-f'(0)/L0 - (1-ζ)f'(0)(1-f0)/L0²` simplifies to `-f'(0)/L0²` because the bracket is
exactly 1. What was wrong was everything around them.

- **They are different models, not two spellings of one.** Zero-inflation *adds* to the
  mass the parent already puts at zero, so `P(Y=0) = ζ + (1-ζ)f(0) > f(0)`: it can only
  ever produce more zeros, and no single zero can be attributed to a mechanism.
  Zero-adjustment *replaces* it — the parent is truncated away from zero — so `π` is
  free to be smaller than `f(0)` too, and the likelihood factorises into a binary part
  and a positive part. Hence: inflation needs a discrete parent **with mass at zero**
  (a continuous one has `P(Y=0)=0` and nothing to inflate — that request is
  `zero_adjusted()`); adjustment takes either, and for a continuous parent needs no
  truncation at all.
- **Neither wrapper can be stacked.** `zero_adjusted(zero_inflated(pois))` was accepted
  and is not a model: truncating at zero cancels `(1-ζ)` between the numerator and the
  truncation constant, so `ζ` leaves the likelihood entirely and its score is
  *identically* zero. `zero_inflated(zero_adjusted(pois))` keeps only `ζ + (1-ζ)π`.
  Both verified numerically; the constructors now refuse.
- **The identifiability rule is a counting rule.** A lattice distribution on `k` points
  has `k-1` free probabilities and the wrapper spends `n_params + 1` of them, so
  `k >= n_params + 2` is necessary — the same bound for inflation and for adjustment.
  It rules out exactly the Bernoulli and `binomial_distrib(size = 1)`. Zero-adjusting a
  Bernoulli leaves the truncated part on `{1}`, and `mu` disappears: the pmf is
  literally the same for `mu = 0.2` and `mu = 0.9`. None of this is visible at run time
  — the pmf sums to one, `check_distrib()` passes, `fit_distrib()` converges somewhere
  on the ridge — so the constructor is the only place it can be caught.
- **`params_smooth` was silently dropped** by all three wrappers, so
  `zero_adjusted(laplace_distrib())` claimed `mu` was smooth and the finite-difference
  guard in `check_distrib()` switched itself off. One line each.
- Parsing a Hessian component name with `strsplit(nm, "_")` and taking the first and
  last piece is wrong for a parameter whose own name contains an underscore. Use the
  internal `hess_pairs(params)`, which inverts `hess_names()` from the parameter vector.
- **That fix did not propagate, and the three sites it missed were the worst ones**
  (found 2026-07-30). `numerical_deriv3()`, `numerical_deriv4()` and
  `expected_by_bartlett()` still recovered their multi-index by splitting the name, so
  a distribution with a parameter called `log_scale` could not obtain a third
  derivative, a fourth, or an expected Hessian by the *default* `approx` — all three
  stopped with an error. Those are precisely the **fallbacks**, the machinery a
  distribution gets when it does *not* implement its own derivatives, so the failure
  was reachable only from a user-defined distribution, which is the whole point of the
  package. None of the fourteen shipped has such a name, which is why nothing noticed.
  `deriv_indices(params, order)` now generates the indices from the same enumeration
  that generates the names, and `deriv_names()` is defined in terms of it, so the two
  cannot disagree. `order_indices()` in `wrapper_derivatives.R` was a third copy of
  that enumeration and delegates too. **`deriv_index_list()` in `link_scale.R` is
  deliberately not merged in** — its order-2 case is ordered for `hess_names()`,
  diagonal first, not lexicographically, and merging them attaches `"mu_sigma"` to the
  index `(sigma, sigma)`. All four now say which is which in their documentation.
  General lesson: when a defect is a *shape of mistake* rather than one line, grep for
  the shape across the package before closing it.

### Truncation

- **"Subtract the mass at ℓ" is not the discrete case, it is the *atom* case.** The
  lower endpoint is included in the truncated support, so `F(ℓ⁻) = F(ℓ) − P(Y = ℓ)`.
  Branching that correction on `discrete_distrib` looks right and is wrong: the cdf of
  `zero_adjusted(gamma)` already contains the point mass at zero, so `F(0) ≠ F(0⁻)`
  even though the object is a `continuous_distrib`. Truncating it from above then drops
  exactly that mass out of `Z`, and the density integrates to something other than one
  while every formula still reads correctly. `parent_mass_at()` asks `distrib_atoms()`
  instead of asking the class. Caught only because `check_distrib()` is atom-aware —
  6 of 7 truncation cases passed and this one failed 5 checks.
- **Truncating zero away from a zero wrapper is the stacking defect again.** For
  `y ∈ [ℓ,u]` with 0 excluded, the `(1−ζ)` factor cancels between numerator and `Z`, so
  ζ leaves the likelihood exactly as in `zero_adjusted(zero_inflated(D))`. Refused.
  Truncating a zero wrapper *elsewhere* is fine and the atom is carried through.
- The support of the truncated law must **not depend on θ** — that is what licenses
  differentiating `Z` under the integral sign and makes truncation at *fixed* points a
  regular problem. Truncating at a θ-dependent point would not be one.
- Truncation composes with itself harmlessly, so the constructor **collapses** nested
  truncation to the intersection rather than nesting: nesting is correct but pays the
  quadrature cost twice.

### Numerics

- **More quadrature knots can be worse.** `expectation()` splits the integral at
  quantiles. Nine knots integrate densities that three cannot, but they also turn loud
  failures into *silent wrong answers*: at Gamma shape 0.02, nine knots return −54.9 for a
  quantity that is exactly zero. Three knots (0.1, 0.5, 0.9) refuse instead, and are more
  accurate in the ordinary range because each extra panel adds its own error.
  **An error is worth more than a plausible number.**
- **Near a non-zero boundary, the spacing of doubles is absolute.** This is why
  Beta(0.9, 0.1) puts 2.5% of its mass within one ulp of 1 and no double-precision sampler
  can resolve it — `rbeta` has the same limitation. Near a boundary at *zero* there is no
  problem at all, denormals giving relative precision to 1e-308.
- **A link maps onto an OPEN interval, and in double precision it did not**
  (fixed 2026-07-31). `linkinv(logit_link(), 37)` was exactly 1, `linkinv(bounded_link(lwr
  = 2), -40)` exactly 2, and the round trip through `linkfun()` then `Inf`;
  `distributions7` validates against open intervals, so a link could hand it a value its
  own validator rejects. Reported as a `bounded_link()` defect and it was not: probing
  every link over |eta| ≤ 800, **nine of eighteen** reached a bound or went non-finite.
  Fixed once in the **generic body** of `linkinv()` — the same interception
  `distributions7` uses for its link scale — so every link inherits it, user-written ones
  included, and every method stays the plain formula.
  Three things in it worth keeping. The bump is `b + |b| * eps` because R has no
  `nextafter` and the note above says a single additive constant cannot serve both a
  bound at 2 and one at 1e-300. **The clamp fires only on EXACT equality with a bound**:
  saturation lands precisely on the endpoint, whereas `inverse_link()` at eta = -40
  returns -0.025 because 1/eta is a bijection from (0, ∞) and -40 is not an admissible
  predictor for it — clamping that would turn a complaint into a small positive number.
  My first version used `>=` and did exactly that.
  And the test for it found a second defect at the other end: `cloglog`'s `linkfun`
  computed `log(-log(1 - theta))`, which for small theta rounds to `log(0) = -Inf` while
  the true value (-176.66 at theta = 1.9e-77) is perfectly representable. `log1p(-theta)`
  returns it, and the four forward derivatives that divide by that logarithm are finite
  again. **`log1p` and `expm1` are not micro-optimisations here; they are the difference
  between a number and an infinity.**
- **`numDeriv`'s Richardson stencil reaches ~8e-4·|x|.** Any grid that comes closer than
  that to a domain boundary is differentiated using points outside the domain, which come
  back NaN. This made `check_link` report failures for links that were exact to 1e-11.
- **Golden-section on a compactified scale is not accurate enough** for locating a mode.
  Its tolerance is expressed in the compactified variable, whose derivative can be
  enormous: with the tangent map `dy/dt ~ y²`, the default tolerance put the "mode" of a
  density centred at 1000 off by 125 standard deviations. `find_pdf_anchor` now refines a
  bracket on a grid and stops on the width measured in y.
- **In R, `NA^0` is 1.** So `theta^(lambda - 2)` silently turns a missing parameter into a
  number as soon as lambda is 2. Any derivative that reduces to a constant must propagate
  missingness explicitly (`const_like()` / `na_from()` in linkfunctions7).
- **A guard clamp is a number, and it has to be derived rather than reached for**
  (2026-07-30). The exponential links floored `exp(eta)` at `.Machine$double.eps`, which
  is not a lower bound on anything relevant — it is the resolution of a double *near 1*,
  not the smallest one that exists. The constraint that actually binds is that the log
  link's fourth forward derivative `-6/theta^4` must not overflow, giving
  `(24/double.xmax)^0.25 ≈ 1.9e-77`. The old value was 61 orders of magnitude higher
  than necessary and corrupted theta for **every** `eta < -36`: `linkinv(log_link(), -40)`
  returned `2.2e-16` instead of `4.2e-18`, and the round trip came back `-36.04`. It also
  did not protect what it claimed to, since the forward derivatives take theta from the
  caller and were never floored at all. Same clamp, same mistake, in `cloglog` and the
  lower-bounded link. Note the near-miss: `double.xmin` looks like the obvious
  replacement and is wrong too — it keeps theta exact further down but lets `1/theta^2`
  overflow. Derive the bound from the expression that binds; do not pick a famous
  constant.
- **`expm1(z)` is `Inf` past z = 709, and a derivative dividing by `val^k` fails long
  before that.** The softplus link was written in `val = expm1(a*theta)`, so
  `d4linkfun` was `NaN` from `theta = 177` with `a = 1`, and from `theta = 24` with
  `a = 30` — an entirely ordinary value for a steep softplus. Rewriting in
  `u = -expm1(-a*theta)` via `log(expm1(z)) = z + log(-expm1(-z))` removes it: the
  identity is exact at both ends, and every derivative follows by dividing through by a
  power of `exp(z)`, which turns each growing quantity into a decaying one. Verify such
  a rewrite by pinning it to the expression it replaces over the range where that one
  still works — agreement to 2e-16 there, plus finiteness beyond, is a much stronger
  statement than agreement with `numDeriv` alone.
- **`chol()` is not a rank test** (2026-08-03). `basis7` asked whether a Gram
  matrix or a penalty was usable by calling `chol()` and treating the error as
  the answer. On a matrix with an *exactly* zero eigenvalue the pivot that
  should be zero comes out positive or negative according to rounding, so the
  same matrix was accepted here and refused on the CI: `orthonorm_basis()` and
  `dr_basis()` disagreed with themselves across machines. The verdict now comes
  from `min(ev) <= tol * max(ev)` on the eigenvalues, which is a statement about
  the matrix rather than about the arithmetic, and the eigendecomposition is
  free next to what both callers already do (`chol_pd()` in
  `R/transformed_basis.R`). Anywhere a `tryCatch(chol(...))` decides a branch,
  the same doubt applies.
  How it surfaced is the second half of the lesson, and it is the **second**
  time in this toolkit that **only the coverage job** saw a defect while
  `R CMD check` stayed green on all five platforms — the first was the S7
  identity comparison above. The two have nothing in common mechanically, so
  the rule is not about `covr`: a green R-CMD-check matrix does not clear a
  change whose outcome depends on floating-point luck, because five platforms
  agreeing is five draws from the same coin. Read a red coverage job.

### Random number generation

- **GRoU (`rng_grou`) is the continuous fallback**, not inverse transform, which costs one
  `uniroot` per draw when the quantile is itself numerical. Two devices make it safe:
  the kernel is **recentred at the mode** (without which a density at mu = 1000 gives a
  degenerate bounding box) and **normalised to a maximum of 1**.
- It refuses far less than expected: bimodal densities, Student t with half a degree of
  freedom, and Pareto with infinite mean are all fine. The only refusal is a density that
  **diverges at an edge**.
- That case is transformed away. If `f(y) ~ |y-a|^(α-1)` near edge `a`, then
  `X = |Y-a|^(1/λ)` has density `~ x^(λα-1)`, bounded once `λα > 1`. **α does not have to
  be searched for**: walking towards the edge in decades lifts the log-density by
  `(1-α)·log(10)` per step, so the probe that *detects* the divergence also *measures* it,
  to four decimals. Gamma shape 0.4 went from 27 ms/draw to 0.8 µs, still exact.
- Divergence at **both** edges needs a map behaving like a different power at each end:
  `T(u) = u^p / (u^p + (1-u)^q)`, `Y = a + (b-a)T(U)`.
- **Discrete distributions needed no new algorithm.** The cdf of a lattice variable is a
  step function, so inverting it is exact and solves nothing. The cost was one R-level
  call per draw; vectorising the table lookup took it from 12.6 µs to 0.1.

### S7

- `S7::method(gen, cls) <- fn` **mutates the generic in place**, so assigning to a local
  copy of the generic still registers globally.
- `identical()` on a method object does **not** work for asking "is this the inherited
  fallback?" — S7 wraps it. Use `attr(m, "signature")[[1]]`, which is the class the method
  was registered on.
- **...and then do not compare that class with `identical()` either** (2026-07-30).
  `identical()` on an S7 class is *object identity*, so it is `FALSE` for a class
  re-created from the same definition — which is what happens whenever a package's code
  is re-evaluated rather than loaded, i.e. under `covr`. Compare `attr(cls, "name")` and
  `attr(cls, "package")`, keeping identity as a fast path; both packages now have a
  helper for it (`is_base_link_class()`, `is_class()`).
  The failure this caused is worth remembering for its shape rather than its size. In
  linkfunctions7 the base numerical fallback was mistaken for an analytic method, so
  every fallback differentiated the order below it and the chain degenerated into the
  nested first differences the design exists to forbid: the log link's fourth derivative
  came back wrong **by a factor of 900**. The local suite passed, `R CMD check --as-cran`
  passed, and the five-platform matrix passed. **Only the coverage job failed**, because
  `covr` is the one context that re-evaluates the code. Two lessons: a green
  `R-CMD-check` does not clear a change that depends on how code was loaded, and the
  reason the failure was legible at all is that the new tests assert on *accuracy*
  rather than on finiteness.
- **S7 objects do not survive a pkgload/installed namespace mismatch** (2026-08-03).
  With the package loaded from source via pkgload, `multistart`'s PSOCK workers
  load the INSTALLED copy; an S7 optimizer built in the dev namespace, shipped to
  a worker, dispatches against the installed namespace's methods and the inner
  runs come back wrong *quietly* -- one iteration, converged FALSE, no error.
  The guard is `exists(".__DEVTOOLS__", asNamespace(pkg))`: under pkgload the
  starts run sequentially, and the parallel path belongs to installed sessions
  (CI exercises it). General shape: any code that serialises S7 objects to
  another R process must ensure both sides run the same namespace.
- Reading an S7 property costs ~2.2 µs against 0.87 for a plain attribute. In hot paths,
  read properties once into locals.
- S7 dispatch itself is ~5 µs against ~0.9 for a plain call. Irreducible.
- **Roxygen block placement matters.** Inserting a helper function between a class's
  roxygen block and its `new_class()` call silently reattaches the documentation to the
  helper. If the helper is `@noRd`, the class ends up undocumented and every `\link{}` to
  it breaks. Internal helpers go *after* the class definition.
  This had already happened to `check_distrib()`: `fd_is_reliable()` sat between its
  roxygen block and the function, so the next `roxygenise()` deleted `check_distrib.Rd`
  and wrote `fd_is_reliable.Rd` instead. The stale Rd had been surviving in `man/`
  because nobody had regenerated since. If a documentation trap is not visible, it is
  because roxygen has not run — run it after touching any file with roxygen in it.
- **An `@include` block not terminated by `NULL` is the same hazard by another route**
  (2026-07-30). It attaches to the next object in the file, which is harmless for as
  long as that object is undocumented — and silently swallows its documentation the
  moment it is not. `cdf_derivatives.R` and `numerical_fallbacks.R` were both like
  this, so `cdf_tail_scale` and `find_pdf_anchor` were written, roxygenised, and
  produced no page at all. It surfaced only as a `checking Rd cross-references`
  warning about links to topics that did not exist, never as anything at the point of
  the mistake. Most files here already end the block with `NULL`; the ones that do not
  are fine only by accident. Also: a `\link{}` to a function carrying `@noRd` is a
  check warning, so documenting a helper that others link to means dropping its
  `@noRd`, not just writing the block.

### Documentation build

- **Regenerate a README with `devtools::build_readme()`**, or
  `rmarkdown::render(output_format = "github_document")`. Not `knitr::knit()`: it runs
  the code correctly but leaves the YAML header at the top of the .md as literal text,
  which GitHub renders as a rule, a line reading `output: github_document`, and another
  rule. pkgdown then carries that malformed markdown into the home page, so the site
  quietly keeps showing stale content even though the workflow reports success — which
  is exactly how this went unnoticed for a while. `linkfunctions7` has a pre-commit hook
  that wants README.Rmd and README.md staged together; it is right to, so satisfy it
  rather than passing `--no-verify`.
- **Two man pages whose topics differ only in case are ONE FILE on Windows.** A page's
  name comes from its topic, so roxygen wrote `Adam.Rd` and then `adam.Rd` over it, and
  seven class pages in `optimizers7` did not exist here while existing on Linux: `?Adam`
  gave the constructor's page and the built package differed by platform. Nothing
  complained — `R CMD check` passed throughout. It bit a *second* time during the
  repair, when roxygen wrote `Adam-class.Rd` and then deleted `Adam.Rd` as obsolete,
  which on this filesystem is `adam.Rd`.
  The classes now live at `X-class` with the plain name as an `@aliases`. But the
  naming is not the fix — a test is: `tests/testthat/test-docs.R` asks *which object in
  the namespace has no help topic*, which is the question that found it, and asks it on
  the machine where the evidence is destroyed. A second check catches the collision at
  its source, where both files can exist. Worth copying into any package whose classes
  are the CamelCase of their constructors.
- **`\value` on every exported topic and an executable example on every exported
  function** are the two commonest reasons a first CRAN submission comes back, and
  `R CMD check` raises neither locally. In `optimizers7` the six topics missing an
  example were the abstract classes and the generics — exactly where an example is most
  useful and least likely to be written.
- **Exporting anything means adding it to `_pkgdown.yml`, and pkgdown enforces
  that in CI rather than locally.** A topic missing from the reference index
  fails the site build — a red push, minutes later, for one line of YAML.
  `pkgdown::check_pkgdown()` asks the question in a second, and both
  `linkfunctions7` and `optimizers7` now ask it from `tests/testthat/test-docs.R`
  rather than relying on anyone remembering: I ran it for one package the same
  afternoon and forgot the other, which is what a habit does. That test found a
  missing `\value` the moment it existed, and then caught its own repair —
  writing the word `\value` inside `\code{}` makes roxygen emit an Rd that will
  not parse, which pkgdown reports as *the topic is missing from the index*,
  because it cannot read the keywords of a file it cannot read. One malformed
  page, two error messages, neither naming the cause.

### Keeping the three packages in step

- **Compare the packages against each other, file by file, from time to time.**
  `git ls-files | grep -v '^R/\|^man/\|^src/\|^tests/'` in each and diff the
  lists: it takes a minute and on 2026-07-31 it found that `optimizers7` had no
  `.Rproj`, no logo, no favicons, and — the one that mattered — none of the
  `check_pkgdown()` test that `linkfunctions7` carries. **That test had been
  failing in CI for two pushes**, for one line of missing YAML, while
  `R-CMD-check` and coverage stayed green on all five platforms; a red job that
  is always the same red job stops being read. A habit that lives in one package
  out of three is not a habit, so when a guard is added anywhere, add it
  everywhere the same afternoon.
- **`git mv` does not happen by itself on Windows.** Seven `man/*.Rd` were
  tracked as `Adam.Rd`, `Bfgs.Rd`, `MultiStart.Rd`… while existing on disk in
  lower case: the `X-class` repair renamed them here, the filesystem could not
  tell the difference, and git kept the old names, so the Linux checkout and
  this one differed by exactly those seven files. Nothing complained, because R
  reads the topic from `\name{}` inside the file and not from the file name.
  The check is `comm -23` of `git ls-files man` against `ls man`, and the fix is
  `git mv -f`.
- **Favicons do not need `pkgdown::build_favicons()`**, which posts the logo to
  an external service. Four PNGs, an `.ico`, the SVG and a copied
  `site.webmanifest` are seven lines of `magick`, and the result is what pkgdown
  looks for.

### GitHub Actions

- **A workflow added in a push does not run on that push** -- in an EXISTING
  repository. GitHub registers the file and triggers nothing, so the first run
  has to be started by hand. Keep `workflow_dispatch:` in every workflow's `on:`
  block so that is possible without an empty commit.
  It does NOT apply to the first push of a NEW repository: `optimizers7` was
  created and pushed on 2026-07-31 and all three workflows ran immediately on
  that push. Worth knowing before waiting for something that has already
  happened, or firing a second run on top of it.
- **A new repository does not serve GitHub Pages until Pages is ENABLED**, even
  though the pkgdown workflow pushes to `gh-pages` and reports success. basis7
  was green on all three workflows, `gh-pages` held the built site, and
  `statmodels7.github.io/basis7/` returned 404 for hours; `gh api
  repos/statmodels7/basis7/pages` answered 404 too, which is the diagnostic.
  The fix is one call:
  `gh api -X POST repos/<org>/<pkg>/pages -f "source[branch]=gh-pages" -f "source[path]=/"`.
  Worth doing immediately after `gh repo create`, because the only thing that
  reports the problem is `R CMD check --as-cran` flagging the package's own URL
  as a possibly invalid one -- and that check is disabled in CI, so nothing at
  all would have said so.
- The Codecov action fetches its own uploader and verifies its GPG signature; when the
  keyserver does not answer that fails the whole job for a reason unrelated to the
  package. `use_pypi: true` avoids that path, and `fail_ci_if_error: false` keeps a
  coverage outage from reddening a badge on a healthy package.

### Testing and measurement

- **A continuous `ks.test` is invalid on lattice data** (ties). It reported p = 0.000 for a
  correct discrete generator. Use a chi-square on the counts instead.
- **Benchmark with the minimum of several batches, not a mean.** A timer that
  self-calibrates its repetition count from one noisy measurement gave numbers that swung
  by 2× between runs and briefly suggested a regression that did not exist.
- **Inject defects to test a validator.** `check_distrib` caught broken gradients,
  Hessians, third derivatives, response derivatives and RNGs — but *not* a wrong cdf,
  because when the quantile is the numerical fallback it is derived from the cdf, so the
  round-trip stays self-consistent however wrong the cdf is. A cdf shifted by a constant
  passed all twelve checks. Always also confirm that an *un*broken reference still passes,
  so the checks cannot be trivially failing.
- **Inject the *missing* case too, not only the wrong one** (2026-07-30). `check_link()`
  caught a derivative that was 5% out and reported `[PASSED]` on all four orders for a
  link that implemented only its first: an unimplemented derivative raised, which left
  `NA` and broke the loop, and the summary reduced with `na.rm = TRUE`. Absent and
  wrong are different failure modes and a validator needs a probe for each — the more
  so here, since `check_link()` exists for whoever is writing a link of their own, and
  an omission is exactly what they are likeliest to have. Pair the two probes in the
  tests: one asserting the missing case now fails, one asserting the 5%-wrong case
  still does, so the first cannot be satisfied by weakening the check.
- **Test grids sit inside the domain; the defects sit at the edges.** Every derivative
  test in linkfunctions7 evaluated well inside the domain, deliberately, so that
  numerical differentiation stayed there — and that is precisely why the softplus
  overflow and the exponential clamp survived 711 passing tests. A separate file that
  goes to the tails costs little and is where the bodies are buried
  (`tests/testthat/test-extremes.R`).
- **A validator that does not know about mixed distributions rejects correct code.**
  `check_distrib()` treated `zero_adjusted(gamma_distrib())` as purely continuous and
  reported four failures on a wrapper that is exactly right: the density integrates to
  `1 - π`, the quantile round-trip cannot hold inside the jump of the cdf, and a central
  difference in `y` across the atom has no derivative to converge to. That is worse than
  a missing check — a user validating their own mixed distribution had no way to tell a
  real defect from this one. The fix is the `distrib_atoms()` generic: the default
  returns nothing (right for both ordinary cases), and the four affected checks ask for
  the atoms and adjust. All four wrappers now pass all thirteen checks, and a gradient
  made 5% wrong is still caught — the guard did not blunt anything.
- **Sweep parameters randomly**, not one fixed value per distribution. That is what
  surfaced the divergence of `expectation()` for Gamma shape below 0.5 and a spurious
  negative eigenvalue in the pseudo-Huber information.
- `purrr::transpose()` **keys rows by the names found inside the elements.** Parameters
  that have been through a link function carry their own name, so a three-parameter theta
  collapsed onto its first column and `fit_distrib`'s default method was broken for the
  pseudo-Huber while Newton and BFGS worked. purrr has been removed entirely.
- **A test that fails on one CI platform is not necessarily about that platform.**
  R-devel went red while the other four passed; the cause was that `check_distrib()`
  called `set.seed(NULL)`, which reseeds from the clock and discards the caller's seed.
  The check was simply random, so it failed roughly one run in several hundred, anywhere.
  Two lessons: never `set.seed(NULL)` in a function a test will call, and when one
  platform disagrees, suspect non-determinism before suspecting the platform.
- **Comparing an analytical derivative against finite differences is invalid at a kink.**
  The Laplace has no derivative at `y = mu`, so a central difference straddling it returns
  a wrong number — the *reference* is wrong, not the value under test. The guard, in
  `fd_is_reliable()`, recomputes the reference with the step halved and drops observations
  where the two disagree; it runs only when `params_smooth` declares a non-smooth
  parameter, and never drops every observation. Note the subtlety: the two estimates must
  be compared **relative to their own magnitude**, not with a denominator floored at 1 —
  near a kink both are tiny yet differ by a factor of two, which a floor of 1 flattens
  into apparent agreement. Always confirm such a guard has not blunted the check by
  re-injecting a deliberate error (a gradient 5% wrong must still be caught).
- **Separate what a component PROMISES from how WELL it does it.** `check_optimizer()`
  runs twelve checks on what an optimiser *reports* — that `value` is the objective at
  `par`, that `converged` follows the rule and is never inferred from the run ending,
  that bounds hold strictly — and then a battery of standard problems whose gaps are
  printed as information, not as pass or fail. Conflating them would have made the
  function useless: gradient descent does not solve Rosenbrock in 500 iterations and is
  not broken, it is slow. The separation earned itself immediately. `bundle()` passed
  every numbered check while diverging on three problems — it did not lie, it merely did
  not work — and the battery is where that showed. A single verdict would have said
  either "fine" or "broken" and neither is true.
- **The validator found two defects the same afternoon it was written**, which is the
  whole argument for writing one. First, `t0` in the bundle method was read as a bare
  multiplier, so the first step was as long as the gradient was big — 233 units on
  Rosenbrock — and the run spent 3000 iterations taking rejected steps. Second, and
  worse, its stationarity measure was the predicted decrease `-v = t||p||^2 + alpha`,
  which carries a factor of `t`, and `t` halves at every null step: it could be driven
  under any tolerance **by the trust parameter shrinking** rather than by the point
  becoming stationary, reporting `converged = TRUE` after zero serious and zero null
  steps. General shape: a convergence measure must not contain a quantity the algorithm
  is free to shrink for reasons of its own.
- **The way to decide whether Rcpp is worth it is to measure, and the answer was
  mostly no** (2026-07-31). Giovanni asked, before publishing `optimizers7`,
  whether `statmodels7` would really be able to hand the optimisers
  log-likelihoods with gradients and Hessians assembled in C++ from model terms.
  Measured on a Gamma regression, an R objective against the identical thing
  compiled: **0.88x at n = 200, 0.55x at 2000, 1.39x at 20 000, 2.09x at
  200 000** — slower below about twenty thousand observations and never better
  than twice as fast. The R callback into the loop costs about 0.5 microseconds
  and is unmeasurable against a real objective; a two-rule criterion costs 24
  microseconds an iteration, which is 0.5 to 5 per cent of a fit.
  The decisive number is structural rather than incidental. **Assembling a
  gradient or a Hessian from term blocks is `crossprod(Z, s)` and
  `crossprod(Z, W * Z)`**, and `src/Makevars` links against R's own BLAS, so the
  two languages call literally the same routine: at n = 50 000, p = 200 the
  Hessian assembly costs 1.33 s in both, and the gradient assembly saves 5 ms on
  it — 0.4 per cent. Whatever `modelterms7` composes, the composition is BLAS.
  Compiled kernels earn their keep on **irregular elementwise arithmetic**,
  which is where `distributions7` already puts them. The conclusion for the
  model layer is *vectorised R over compiled kernels*, and the real order of
  magnitude is algorithmic — Fisher scoring exploiting block structure and sparse
  penalties — not linguistic.
- **A feature with one caller is that caller's feature.** `adam(resample=)` and
  `finite_sum()` let the optimiser draw its own minibatches, and between them
  cost a second objective class, a rule for which criteria it permitted, a token
  in the criterion machinery, a branch in the compiled loop and an argument to
  select the path. An objective that resamples is a closure; the caller who
  knows what an observation is writes it in one line and can also say what a
  stratum is. Removing it took 378 lines out. Same judgement removed
  `cpp_objective()`, which additionally *could not carry data* — the pointer type
  `double(*)(const arma::vec&)` has no closure, so any real objective needed C++
  globals, which is not reentrant.
  Watch for the second-order consequence: `need_value` had made the objective
  evaluation conditional on the criterion reading it, and the token answered
  that question. Removing the token made the condition uncomputable, and a stale
  value would have let a rule comparing `f_new` with `f_old` compare a number
  with itself and fire at once.
- **Writing the documentation is a test of the design.** The `optimizers7` vignette
  works through a user-written optimiser, and writing it showed that a user *could not*
  write a conforming one: honouring `bounds` and refusing an unevaluable stopping rule
  both needed internal functions. The package promised extensibility and did not supply
  the parts. Three exports later it does. The same vignette also disproved two claims
  in its own prose — measure the example before describing it.

---

## 8. Open items

- The two logistic integrals above (cosmetic — `approx = "integrate"` already delivers
  them within Monte Carlo noise).
- Publishing `linkfunctions7` **and `optimizers7`** to CRAN, which together unblock
  `distributions7`. They do not depend on each other, so they go in parallel.
  Before either: version `0.1.0` in place of `0.0.0.9000`, a `NEWS.md` (neither has
  one), `cran-comments.md`, and a run on win-builder.
- **A URL of the form `https://statmodels7.github.io/<pkg>` redirects 301** to the
  same address with a trailing slash, and `--as-cran` reports it as possibly
  invalid. All six occurrences across the three READMEs were fixed on 2026-08-03.
  CI never sees this: the URL check lives inside *checking CRAN incoming
  feasibility*, which the workflow disables. CRAN would return it at submission.
- ~~A nonmonotone line search for `bb()`~~ **done 2026-07-31.**
  `nonmonotone(c1, shrink, memory, max_step)` is a third `line_search` subclass;
  the loop keeps a `std::deque` of the last `memory` values and Armijo tests
  against their maximum. `memory = 0` is `armijo()` value for value, which is
  what makes the comparison a comparison of that one number: Rosenbrock 68
  iterations and 77 evaluations against 82 and 186, eleven uphill steps accepted
  out of sixty-seven where Armijo accepts none.
  What the measurement exposed was a **separate defect in `bb()`**, and its
  shape is the one worth keeping. When the secant pair reports no usable
  curvature, alpha has to come from outside, and **both constants one reaches
  for are absorbing states**: freezing the previous alpha traps a short step in
  its own shortness (bb2 spent 873 of 945 iterations there), and restarting at
  `alpha0` does exactly the same wherever `alpha0` is itself too short — on a
  boxed quadratic seen through its reparametrisation, 1395 refusals in 1521
  iterations. `alpha_max`, which is what the SPG literature does, hands the
  search a direction of length 1e10 that thirty halvings do not rescale, and
  that run stopped a unit from the solution *reporting success*. What works is
  `alpha = 1/max|g|`, a step of order one in the parameters: it cannot freeze,
  since it does not depend on the alpha it replaces, and cannot explode, since
  it scales with the gradient. Measured over Rosenbrock, Beale, Powell and the
  boxed quadratic it is the only one of the four never worst. The curvature test
  is also now **relative**, `s'y > curv_tol*||s||*||y||`, as `bfgs()` already had
  it. General lesson: a policy chosen on one problem is a policy chosen on one
  problem — the boxed case disagreed with Rosenbrock about all three constants.
- Next packages: `modelterms7`, `basis7`, `penalties7`.
- **`infer_npar()` cannot decide for a vectorised objective, and that is
  arithmetic rather than a gap.** `start_zeros()` with no `npar` and no bounds
  probes the objective, and both plausible guesses about R are wrong: recycling
  warns only when the shorter length is not a **divisor** of the longer, so
  `sum((p - c(1,2,3))^2)` accepts a length-one vector *in silence* and returns a
  perfectly finite 14; and its gradient `2*(p - c(1,2,3))` returns six
  components for a length-six argument, so requiring `length(gr(x)) ==
  length(x)` passes it too. What decides is an objective with a fixed width
  built in — `X %*% beta` of the wrong length is an error — which is every real
  model and no toy. So `statmodels7`, which always knows `npar`, will never meet
  this; a user at the console with a toy objective will, and is told to say
  `npar` rather than being guessed at.
- **A censored-likelihood front end.** `distrib_grad_cdf()` supplies the pieces, but
  nothing yet assembles them: a `fit_distrib(..., censored = )` taking a status vector,
  or a `Surv()`-like object. That is the step that turns the capability into a feature.
- **Second-order closed forms** for the families that only have the gradient
  (lognormal, invgauss, student t, pseudohuber). Mechanical but fiddly; the FD fallback
  is ~1e-8 and adequate.
- **Gamma and beta shape derivatives** would need the series representation of
  `d/da P(a,x)` (a convergent alternating sum, badly cancelling for large x) or a
  Meijer-G evaluation. Probably *not* worth it: the FD of pgamma/pbeta already gives
  ~2e-10, and a hand-rolled series would likely be worse in some regime. Revisit only
  if profiling shows it matters.
- **Third and fourth cdf derivatives** would let truncation drop quadrature at orders
  3-4 as well.
- No `NEWS.md` on either package.
- **`expected_by_bartlett()` recomputes too much.** Inside the integrand it calls
  `observed_deriv()` for a whole order and keeps one component, once per block per
  partition per component — and the integrand runs at every quadrature node. Memoising
  per integrand call looks worth it at orders 3 and 4, but `"bartlett"` is only the
  default at order 2, where the score is all it needs, so measure before touching it.
- **Logos** exist but are competent rather than designed — drawn by a script, not by
  someone with visual judgement. Worth redoing with a designer if the identity matters.
  `logo/make-logos.R` regenerates them; the curves are the real `plogis()` and `dgamma()`.
  The palette is deliberately taken from Giovanni's earlier `mvreg` sticker — chalkboard
  green `#3D6B4C`, rust border `#9C3E11`, chalk `#F7F4D4`, monospace wordmark — so the
  toolkit reads as one family with his existing work. Keep it if you add a package.
  A decoration that carries no information does not belong: an earlier version marked
  the mode of the density with a dot, and he was right to have it removed.
- The portal is hand-written HTML. If the toolkit grows past a handful of packages it may
  deserve a generator, but not yet.
- The book now covers links, distributions, the transformation wrappers and fitting
  (chapter 5, added 2026-07-30: the link scale as the place to optimise, Fisher
  scoring versus Newton and why the congruence corollary makes the expected
  information the right matrix to invert, the line search, the delta method as the same
  congruence applied to the inverse, and intervals built on the link scale and mapped
  back; revised 2026-08-03, when `fit_distrib()` started delegating to
  `optimizers7` — §5.2 now says the optimisation is not written there, the step
  rule is the sufficient-decrease condition of §4.1 rather than plain halving,
  the non-PD Hessian is repaired by §4.2 rather than abandoned, and a chunk shows
  `method = lbfgs(...)`), and §3.4 on the numerical fallbacks — the ratio-of-uniforms theorem, the
  mode recentring, the divergence transform at one edge and at two, the
  discrete cumulative table, and the two warnings that belong with them. Chapter 4
  (added 2026-07-31) covers `optimizers7`: descent directions and what a line search
  must guarantee, with Zoutendijk's theorem; Newton's Hessian repairs and why
  the eigenvalue floor is what keeps that theorem's angle bound alive; the secant
  equation and the BFGS update, then conjugate gradients
  and Barzilai-Borwein as the bottom of a ladder ordered by how much curvature a
  method stores, from p(p+1)/2 down to one number; the subdifferential, Fermat's
  rule, and the bundle subproblem with its dual; and the box
  reparametrisation and Adam; and 4.5 (added 2026-07-31) on where a run begins --
  why a constant is only sensible on the unconstrained scale, how the number of
  parameters is settled, the Latin hypercube, the count of distinct optima as
  the thing several starts actually measure, and why the parallelism is at the
  level of processes. Its gate is `assert_starting_ok()`, injection-checked
  against a printed h(0) wrong by 0.5 and against independent draws replacing
  the hypercube. Chapter 5 (added 2026-08-03) covers `basis7`: the expansion
  and its design matrix, the anchored integral and why the constant is fixed,
  the Vandermonde stencil; the Cox-de Boor recurrence with local support and
  the partition of unity, the Fourier phase shift, the Legendre recurrences;
  the Gram matrix as the matrix of a roughness penalty, its exact integration
  knot interval by knot interval, and the three measures the inner product can
  be taken against; orthonormalisation, constraints and Demmler-Reinsch as one
  congruence; and tensor products, separability, and the two coefficient
  shapes a contraction accepts. Its five gates are in
  `book/R/basis-certificates.R`, injection-checked seventeen times. What the
  book still lacks is anything on **censored likelihoods**, which waits on the
  front end that does not exist yet, and it will need a chapter per package as
  the others arrive.

---

## 9. The book

`book/` is a **Quarto book** (`quarto render` from inside it; Quarto CLI 1.8.24 is on
PATH, and no new R packages are needed — knitr and numDeriv suffice). It is the
mathematical companion Giovanni asked for on 2026-07-26: every formula the toolkit
implements, explained and attributed, with the argument in full left to the
citation. English, on purpose — it is the public document of a toolkit aiming at
CRAN.

**Structure: one chapter per package** (Giovanni, 2026-07-30, explicit). The book
is about `statmodels7`, not about distributions, and the earlier arrangement — a
chapter each for links, distributions, transformations, fitting — read as a book
on distributions with a chapter on links attached. Now:

```
Preface / 1 Introduction / 2 The linkfunctions7 package /
3 The distributions7 package  (3.1 Distributions, 3.2 Transformations,
                               3.3 Fitting, 3.4 Fallbacks) /
4 The optimizers7 package     (4.1 Descent, 4.2 Curvature, 4.3 Non-smooth,
                               4.4 Constraints, 4.5 Starting points) /
5 The basis7 package          (5.1 Expansions, 5.2 Families, 5.3 Inner
                               products, 5.4 Transformations, 5.5 Several
                               variables) /
A Notation / B References
```

When `modelterms7` and the rest arrive they arrive as chapters, and nothing
already written moves. Quarto numbers chapters per *file*, so a package chapter is
one `.qmd` that pulls its sections in with `{{< include >}}` from `_`-prefixed
files that are **not** listed in `_quarto.yml`; that keeps 3.1, 3.2 numbering
without a two-thousand-line file. `toc-depth` is 4, `number-depth` 3.

⚠️ **Demoting headings when reorganising:** shift `#` levels only outside fenced
code (where `#` is an R comment) and outside `:::` divs — a `##` inside
`::: {#thm-x}` is the *theorem's title*, which Quarto reads as the block's name,
not a section heading.

**Terminology, fixed 2026-08-03 (Giovanni, explicit): "lattice" is BANNED** as a
name for discrete distributions, everywhere in documentation, book and code
comments -- the word is **"discrete"** (or "support points" for the points
themselves). The one legitimate survivor is the *partition lattice* in the book,
which really is an order-theoretic lattice. Swept from distributions7 and the
book on 2026-08-03; grep for it before writing new prose.

**Documentation style, fixed 2026-08-03 (Giovanni, explicit).** Package
documentation is **aseptic and synthetic**: it states what a function does and
briefly which method it uses, for EVERY function and S7 method. No anecdotes,
no measured war stories (the three-knot Gamma tale, iteration counts from named
runs), no punchy subsection titles, no stock phrases, and no strings of short
sentences -- flowing prose. Methods whose behaviour a generic's or constructor's
page fully specifies may share that page; methods with behaviour of their own
(defaults, closed-form overrides, numerical fallbacks, wrapper re-derivations)
get their own internal page. He also said plainly that my writing style did not
please him: no eloquence, in docs OR book titles -- every heading names its
subject ("The BFGS update", not "Below a matrix").

**The word is "toolkit", never "stack"** (Giovanni, 2026-07-30) — in the book, the
portal and both package READMEs. Careful with the verb: "wrappers do not stack"
is the other sense and must survive, as must the `sec-no-stacking` label.

⚠️ **`quarto render` needs R on PATH.** Quarto resolves R from the registry and
finds `R-4.5.1`, the empty directory of §3, then stops with *"Unable to locate an
installed version of R"*. Prefix the call:
`PATH="/c/Program Files/R/R-4.6.0/bin:$PATH" quarto render`.

⚠️ **`freeze` is `false`, deliberately, and must stay that way** (2026-07-30).
Quarto's freeze cache is keyed on the `.qmd` and does not see the R files a
chapter `source()`s — which is exactly where the printed formulas and their
transcriptions live, and where the packages themselves sit one directory up. With
`freeze: auto`, editing a formula in `R/link-formulas.R` and re-rendering left the
old text in the output, and, far worse, **the consistency gate for that chapter did
not re-run either**: a full render finished in seconds and reported success without
having checked anything. The verification story of this section is only true when
every chapter actually executes. The cost is a render of twenty-odd minutes; it is
already a manual step, so when it is run it should be real.

**Editorial line (Giovanni, 2026-07-26, explicit):** the book takes it for granted
that the formulas shown are the ones the packages run. The book-vs-package
consistency checking is *ours* and must not be discussed in the text — no
certification tables, no "verdict" columns, no meta-narrative about verification.
And the prose must be flowing, discursive, book-like — not strings of short
sentences — with the mathematics explained step by step.

**Editorial line, revised (Giovanni, 2026-07-31, explicit) — this supersedes the
"nothing is cited away" rule the preface used to state.** **No proofs.** Show the
relevant formulas and enough reasoning to make them intelligible — where a
quantity comes from, what it depends on, which cancellation the implementation
relies on — and then **cite** the paper or textbook that gives the argument in
full. `references.bib` and appendix `A2-references.qmd` (a `::: {#refs} :::`
block) were added for this; `bibliography:` and `link-citations: true` sit in
`_quarto.yml`, and pandoc's default author-date style is used, no CSL file.
`sync-stack-files.sh` carries `references.bib`. Keep the bibliography free of
uncited entries — pandoc silently ignores them, so they rot unnoticed; the check
is a shell loop over `grep -o "@[a-z]*[0-9]\{4\}"` against `^@[a-z]*{`.

**Prose register, fixed 2026-08-03 (Giovanni, explicit, after two milder passes
did not satisfy him): the book reads as a plain textbook.** Two constructions
are banned everywhere: the announce-then-reveal opening ("Notice first what the
definition does NOT say: ...") -- state the fact directly instead -- and the
compressed two-beat contrast ("The point of X is not Y. It buys Z.") -- say
"X matters because it yields Z, which Y alone does not". No second person, no
suspense colons, no one-line punch sentences; every paragraph opens by saying
its subject. The full restyling was done 2026-08-03, chapters 2 and 4 rewritten
wholesale. It caught a live defect: the 4.4 box example still passed the removed
bounds= argument, silently swallowed by ..., so its table had been computed
UNCONSTRAINED. When rewriting prose around an executable chunk, check the chunk
against the current API.

The book is also **didactic**, and that rules out two things it used to do:
**no punchy titles** ("The catch, stated plainly", "Below a matter", "Where the
derivative fails" → say what the section is about), and **no over-segmented
prose** — sentences joined into paragraphs that read continuously, not strings of
three-word declaratives. All seventeen proof blocks were removed on 2026-07-31
and the preface rewritten; the headings were swept at the same time.

⚠️ **The book's own source had escape damage, committed and live.**
`_04b-curvature.qmd` §"Below a matrix" had `\beta`, `\alpha`, `\frac`, `\top`,
`\times`, `\ne` collapsed into BELL, BACKSPACE, FORMFEED and TAB characters —
the §3 trap, landed in a chapter written through a shell heredoc, and it had been
rendering as garbage at `statmodels7.github.io/book/` for a day. `cat -A` on
every `.qmd` looking for `^[GHIL]` is the check; `^M` alone is only CRLF and is
harmless. **Some book files are CRLF and some are LF**, so a Python edit script
must normalise on read and restore on write, or every marker containing `\n`
silently fails to match.

**Notation, fixed 2026-07-26 (Giovanni).** Derivatives of the log-likelihood use
**parenthesized superscripts** — `l^(i)`, `l^(ij)`, `l^(ijk)` — never subscripts. A
subscript on `l` means what it conventionally means in likelihood work: the
contribution of one observation, or the model the log-likelihood belongs to (`l_Y`,
`l_T`). Chapter 4 previously used `s_i`/`H_ij` for the same quantities; unified —
but that unification was **incomplete until 2026-07-30**, four `s_i s_j` surviving in
@thm-truncated and its proof, where they contradicted the notation appendix on the
same page. They survive a careless grep because `s_is_j` has no word boundary between
the two symbols: search for `s_is_j` as well as for `s_i`. Note
this is *not* Einstein/McCullagh convention — there, partial derivatives are covariant
(lower) indices and a bare superscript `l^rs` denotes the *inverse* information. The
parentheses are what disambiguate, as in `f^(k)`.

**Truncation points are `L` and `U`**, not `l` and `u`: the lower point collided with
the log-likelihood symbol, and both appeared in the same comparison table. The uniform
draw in the inverse-transform formula is `V` for the same reason.

The same collision decided chapter 5. `basis7` writes a basis interval
`[l, u]` in its own documentation; the **book writes `[a, b]`**, with `r` and
`s` for basis-function indices, `q` for the intermediate degree of the
Cox-de Boor recurrence, `tau` for the shifted Legendre variable and
`delta_j` for a finite-difference offset. The notation appendix carries a
Bases table and says which symbols are reused across chapters -- `B` is a
curvature model in chapter 4 and a design matrix in chapter 5, `m` a spline
degree here and a truncated expectation there.

**Pipes break tables.** A `|` inside math in a markdown table splits the cell; kable
escapes it to `&#124;` and MathJax prints that verbatim (visible in the transformer
table until 2026-07-26). Always `\lvert ... \rvert`. After every render:
`grep -c '&#124;' _book/chapters/*.html` must be all zeros — this is not covered by the
R gates, which run before rendering. As of 2026-07-30 all zeros, and no `|` in math
sits in a table row; the eight that remain are all in prose, where they are harmless
but still off-convention.

**`date: today` makes the output churn.** Every render stamps a new date into the
HTML, so `sync-stack-files.sh` produces a diff in `site/` even when nothing changed.
Harmless, but do not read a diff after a re-render as evidence that something did.

**How consistency is enforced without being mentioned:** each chapter ends with a
hidden chunk (`include: false`) calling `assert_links_ok()`,
`assert_distributions_ok()`, `assert_transformations_ok()` or `assert_fit_ok()`
(defined in `book/R/*.R`). These evaluate every printed formula against the running
package, the structural identities (inverse function theorem, Bartlett,
expected-vs-observed Hessian), and `check_distrib()` on everything — and **stop the
render** on any disagreement. Verified by injection: a 5% corruption of a printed
derivative or density is caught and kills the build. The reader never sees any of
it — every setup and gate chunk is `include: false`, which matters more than it
looks, since one of the files sourced is named `transformation-certificates.R` and
a visible `source()` line would leak the whole apparatus.

`assert_fit_ok()` (chapter 5, `book/R/fit-certificates.R`) checks the four things
that chapter claims, each against a route the chapter does not itself use: that the
link-scale score really vanishes at the reported optimum, that
`V_theta = diag(h') V_eta diag(h')` with `h'` taken fresh from linkfunctions7, that
no interval escapes its parameter's domain or comes back inverted, and that Fisher
scoring, Newton and BFGS agree on the maximised log-likelihood. Injection-checked
too: at the optimum the score is 4e-13 per observation against a 1e-4 threshold and
2e-1 when displaced by 0.05, and the delta-method check is exact against 1e-1 for a
Jacobian 5% wrong.

Extended 2026-08-03, when `confint()` was added: the interval checks now read
through **`confint()` rather than the `@ci` slot**, because that is what the
chapter displays and a check on a quantity the reader never sees would not cover
the one shown. Two conditions were added with it — that the link-scale interval
is symmetric about the estimate, and that the parameter-scale one is its image
under `linkinv()`, mapped independently of the package's own mapping and sorted,
since a link may decrease. Both are injection-checked: shifting one end of the
link-scale interval by 0.05 and inflating the parameter-scale one by 5% are each
caught, and the gate passes again once they are removed.

`assert_optimizers_ok()` (chapter 4, `book/R/optimizer-certificates.R`) does the
same for six claims: that an accepted step really satisfies the condition claimed
for it, that the BFGS update satisfies the secant equation and stays positive
definite, that the two-loop recursion equals the explicit inverse built from the
same pairs, that the box transform's chain rule matches a numerical derivative of
the composed objective, that the aggregate subgradient vanishes at the median
where an ordinary one has norm 1, and that every method passes `check_optimizer()`.
Injection-checked: a wrong sign in the printed upper-bound transform, a 5% error
in the lower one, and a mis-scaled two-loop recursion are all caught.

⚠️ **A gate that compares the package with itself proves nothing.** The bound
checks originally compared `h'` against a numerical derivative of `h` — both from
the package, so a wrong *printed* formula would have passed. The chapter's three
cases are now transcribed by hand into the gate and compared, and that is what
caught the upper-bounded transform being written `U - exp(-eta)` where the code
computes `U - exp(eta)`. The same discipline as `link-formulas.R`: the printed
LaTeX and an independent transcription live in the same record.

The mechanism that makes this cheap: `book/R/link-formulas.R` and
`distribution-catalogue.R` keep the printed LaTeX and an independent R
transcription of it **in the same record**, and the chapters generate the displayed
formulas from those records. Editing a formula edits both the equation the reader
sees and the thing under test. Do not split them.

Numerical-differentiation rule (also in `check_link()`): **never nest numDeriv** to
reach high orders — by order 3 the reference is noise (the identity link's exactly
zero third derivative comes back as O(1)). Compare order k against ONE numerical
differentiation of the analytical order k−1.

Loading is from **source** via `pkgload` on `../linkfunctions7` and
`../distributions7`, so the book documents the working tree. `execute-dir: project`
keeps the working directory at `book/`. The palette in `assets/theme.scss` is the
toolkit's: chalkboard green `#3D6B4C`, rust `#9C3E11`, chalk `#F7F4D4`.
