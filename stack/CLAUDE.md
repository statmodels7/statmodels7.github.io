# statmodels7

Read this first. It is the orientation file for the whole project, and it is the only
place where the accumulated knowledge lives: Claude Code's memory is keyed by directory
path, and these packages were moved here on 2026-07-22, so nothing from the earlier
sessions carries over automatically.

`statmodels7/` is a plain directory, not a repository — the three repositories sit
inside it. This file and `logo/make-logos.R` are therefore versioned in the portal
repository (`site/`, which is `statmodels7.github.io`) and copied back up here, so a
fresh clone can recover them. If you edit either, copy the change into `site/` and
commit it there; `site/sync-stack-files.sh` does exactly that.

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
unconstrained parameters — and then everyone builds on top. That is what this stack is:
infrastructure to be reused, not another modelling package.

Extension matters just as much in the telling: a stack that only worked for the fourteen
distributions it ships with would have solved nothing.

Giovanni corrected an earlier draft that led with the link scale and the optimisation
story. Lead with the reuse argument.

## 1. What this is

`statmodels7` is an umbrella for a statistical modelling stack built entirely on the
**S7** object system. The name carries the convention: every package is a *plural noun
followed by 7*, which happens to spell `...s7` — `linkfunction**s7**`,
`distribution**s7**`, `statmodel**s7**`. The 7 is S7; the s is the plural, not decoration.

**Built so far**

| package | what it provides |
|---|---|
| `linkfunctions7` | 16 link classes (14 constructors) with exact analytical derivatives to 4th order, both directions |
| `distributions7` | 14 distributions with exact score, information and 3rd/4th derivatives, plus wrappers, transformations, MLE |

**Planned** — `modelterms7`, `basis7`, `penalties7`, and eventually the `statmodels7`
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
    book\            the Quarto book (see §9); `quarto render` from inside it
    logo\            hex logos: make-logos.R draws them, run from this directory
    site\            the portal, its own repository (statmodels7.github.io)
```

⚠️ `book/` has **no repository yet**. Like this file and `logo/`, it is not
covered by any of the three repos, so it needs a home — the natural one is `site/`,
which would also let the rendered book be served from the portal.

GitHub: `github.com/statmodels7/{linkfunctions7,distributions7}`. The repositories were
transferred from `giovannitinervia9/*` on 2026-07-22; GitHub keeps redirects, so
`install_github("giovannitinervia9/distributions7")` still resolves.

Websites, all live and rebuilt by a `pkgdown.yaml` workflow on every push:

| | |
|---|---|
| portal | `statmodels7.github.io` — plain HTML, own repository, deploys from `main` |
| linkfunctions7 | `statmodels7.github.io/linkfunctions7` — pkgdown, from `gh-pages` |
| distributions7 | `statmodels7.github.io/distributions7` — pkgdown, from `gh-pages` |

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
R.exe CMD check <pkg>_0.0.0.9000.tar.gz --no-manual
```

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
logit, expit, affine). All propagate `params_smooth`; all validate their parent in the
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
derivatives are much dearer than the parent's and orders 3–4 use the FD fallback.

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

**User-facing tools**: `check_distrib()` (thirteen numerical checks on a continuous
distribution, twelve on a discrete one; catches deliberately broken components),
`fit_distrib()` (MLE on the link scale, Fisher scoring / Newton / BFGS, delta-method SEs,
CIs built on the link scale and mapped back so they respect domains), plus `simulate()`
and `plot()` on the fit.

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
| `linkfunctions7` | 711 tests, `R CMD check` OK, CI green |
| `distributions7` | 1168 tests, `R CMD check` OK (2026-07-26, local), CI green |

Both repositories run `R-CMD-check` on macOS, Windows and three Linux/R combinations
(devel, release, oldrel-1) plus a coverage workflow, all green. That matrix matters for
`distributions7` in particular: it ships Rcpp kernels, and until the CI existed they had
only ever been compiled by one compiler on one machine.

Vignettes: `defining-a-distribution`, `fitting-a-model`,
`derivatives-and-the-link-scale` (distributions7) and `link-functions`
(linkfunctions7). Both packages have a README with badges.

`R CMD check --as-cran` on distributions7 is clean apart from submission metadata:
new submission, version `0.0.0.9000`, the `Remotes:` field, and linkfunctions7 not being
on CRAN. **The real CRAN blocker is that linkfunctions7 must be published first.**
Both names, and `modelterms7` / `basis7` / `penalties7`, are free on CRAN.

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

### GitHub Actions

- **A workflow added in a push does not run on that push.** GitHub registers the file but
  triggers nothing, so the first run has to be started by hand. Keep `workflow_dispatch:`
  in every workflow's `on:` block so that is possible without an empty commit.
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

---

## 8. Open items

- The two logistic integrals above (cosmetic — `approx = "integrate"` already delivers
  them within Monte Carlo noise).
- Publishing `linkfunctions7` to CRAN, which unblocks `distributions7`.
- Next packages: `modelterms7`, `basis7`, `penalties7`.
- Truncated distributions have no closed-form 3rd/4th derivatives (the FD fallback of
  the analytical Hessian is used). Deriving `d^3 log Z` and `d^4 log Z` from the
  moment-to-cumulant relations would close that, and would also give the zero wrappers
  their missing orders.
- No `NEWS.md` on either package.
- **Logos** exist but are competent rather than designed — drawn by a script, not by
  someone with visual judgement. Worth redoing with a designer if the identity matters.
  `logo/make-logos.R` regenerates them; the curves are the real `plogis()` and `dgamma()`.
  The palette is deliberately taken from Giovanni's earlier `mvreg` sticker — chalkboard
  green `#3D6B4C`, rust border `#9C3E11`, chalk `#F7F4D4`, monospace wordmark — so the
  stack reads as one family with his existing work. Keep it if you add a package.
  A decoration that carries no information does not belong: an earlier version marked
  the mode of the density with a dot, and he was right to have it removed.
- The portal is hand-written HTML. If the stack grows past a handful of packages it may
  deserve a generator, but not yet.
- **`book/` is not versioned anywhere.** It needs a repository — `site/` is the obvious
  one, since it already hosts the stack-wide files and could serve the rendered book
  from the portal. Decide before the book accumulates history worth keeping.
- The book covers links, distributions and the transformation wrappers. `fit_distrib`
  (Fisher scoring, Newton, BFGS, delta-method standard errors, intervals built on the
  link scale and mapped back) has no chapter yet, and neither do the numerical
  fallbacks beyond the summary in §3.7 of the book.

---

## 9. The book

`book/` is a **Quarto book** (`quarto render` from inside it; Quarto CLI 1.8.24 is on
PATH, and no new R packages are needed — knitr and numDeriv suffice). It is the
mathematical companion Giovanni asked for on 2026-07-26: every formula the stack
implements, derived from the definitions with the steps written out. Chapters:
preface, introduction, link functions, distributions, transformations, plus a
notation appendix. English, on purpose — it is the public document of a stack aiming
at CRAN.

**Editorial line (Giovanni, 2026-07-26, explicit):** the book takes it for granted
that the formulas shown are the ones the packages run. The book-vs-package
consistency checking is *ours* and must not be discussed in the text — no
certification tables, no "verdict" columns, no meta-narrative about verification.
And the prose must be flowing, discursive, book-like — not strings of short
sentences — with the mathematics explained step by step.

**How consistency is enforced without being mentioned:** each chapter ends with a
hidden chunk (`include: false`) calling `assert_links_ok()`,
`assert_distributions_ok()` or `assert_transformations_ok()` (defined in
`book/R/*.R`). These evaluate every printed formula against the running package,
the structural identities (inverse function theorem, Bartlett, expected-vs-observed
Hessian), and `check_distrib()` on everything — and **stop the render** on any
disagreement. Verified by injection: a 5% corruption of a printed derivative or
density is caught and kills the build. The reader never sees any of it.

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
stack's: chalkboard green `#3D6B4C`, rust `#9C3E11`, chalk `#F7F4D4`.
