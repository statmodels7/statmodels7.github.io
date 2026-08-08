# statmodels7

Read this first. It is the orientation file for the whole project, and it is the only
place where the accumulated knowledge lives: Claude Code's memory is keyed by directory
path, and these packages were moved here on 2026-07-22, so nothing from the earlier
sessions carries over automatically.

`statmodels7/` is a plain directory, not a repository — the package repositories sit
inside it. This file, the logo sources (`logo/svg/*.svg`, `logo/author-svg.py`,
`logo/make-logos.R`) and the package plans (`piano_*.txt`, the
decision records drafted before a package's code exists) are therefore versioned in
the portal repository (`site/`, which is `statmodels7.github.io`) and copied back up
here, so a fresh clone can recover them. If you edit any of them, copy the change
into `site/` and commit it there; `site/sync-stack-files.sh` does exactly that.

---

## 0. How to talk about it

Get this right before writing any prose for a README, a site or a paper. The
positioning is **not** "a framework with exact derivatives and optimization on the
unconstrained scale" — that is a description of the machinery, and it buries the point.

The point is that **these things get rewritten in every package that needs them**.
Almost every R modeling package carries its own distributions and link functions,
written as internal helper functions — a `switch` on a character string, closures
private to the package that owns them. They are not objects, so nothing outside can
reuse them, extend them, or ask them for anything their author did not happen to need.
And since each package writes only what it needs, that is usually the density and the
score, sometimes the Hessian, and nothing beyond.

A Gamma is a fixed mathematical object. It should be written once, correctly, with
everything a modeling routine could want already computed — third and fourth
derivatives, derivatives with respect to the response, derivatives with respect to the
unconstrained parameters — and then everyone builds on top. That is what this toolkit is:
infrastructure to be reused, not another modeling package.

Extension matters just as much in the telling: a toolkit that only worked for the fourteen
distributions it ships with would have solved nothing.

The same argument is what `optimizers7` is for, and it is worth having the sentence
ready. Almost every R package that fits a model carries its own optimizer, written
inside the function that needs it — a loop, a `while (!converged)`, a tolerance compared
against whatever quantity the author had to hand. The stopping rule, which decides what
the whole thing means by *finished*, is a number buried three levels down. Lead with
that, not with the list of algorithms.

Giovanni corrected an earlier draft that led with the link scale and the optimization
story. Lead with the reuse argument.

## 1. What this is

`statmodels7` is an umbrella for a statistical modeling toolkit built entirely on the
**S7** object system. The naming convention is simply that **every package name ends
in 7**, the 7 being S7 — `linkfunctions7`, `distributions7`, `statmodels7`. Until
2026-08-03 the convention also required a plural noun; Giovanni dropped that when
naming `basis7` ("basis7 è più bello di bases7"), so the plural is history, not a
rule — do not "fix" singular names, and do not cite the plural reading in new prose.

**Built so far**

| package | what it provides |
|---|---|
| `numericals7` | the numerical layer at the ROOT of the toolkit (created 2026-08-05; jets REMOVED 2026-08-06, see the de-jettization item in section 7): the enumerations (`set_partitions`, `tuple_indices`, `compositions` -- ONE copy each), the stencil library (`fd_weights`/`fd_offsets`/`fd_step`/`fd_derivative`), `quad_vec`/`series_vec` (quadrature and series vectorized over the parameters, convergence on the SUM of a row's panel errors), and the special functions (`mills_ratio`, `owen_t`, `bessel_i_ratio` with derivatives and inverse; `log_bessel_i`/`log_bessel_k` with four argument-derivatives, after Plesner-Sorensen-Hauberg ICS 2024 (arXiv:2409.08729), 0.5.0 -- finite wherever the log itself is representable, switching guards tightened on the Wronskian, pure R after an Rcpp transcription measured at 0.9-2.7x). No S7 classes on purpose: these are functions |
| `linkfunctions7` | 16 link classes (14 constructors) with exact analytical derivatives to 4th order, both directions, plus numerical fallbacks for user-defined links; stencils delegate to numericals7 |
| `distributions7` | 41 univariate families (census 2026-08-08: 41 univariate + 4 multivariate = 45 constructors; the "38" predated pig1/pig2 and had rotted on the site, the book and the README at once) -- ONE NAME PER PARAMETRIZATION (2026-08-05, Giovanni): 13 numbered groups (gaussian1/2/3, gamma1/2, negbin1/2 after Cameron-Trivedi, weibull1/3 after gamlss WEI/WEI3 with weibull2 deliberately empty, student_t1/2, skewnormal1/2, vonmises1/2, invgauss1/2, lognormal1/2, beta1/2, betabinom1/2, gengamma1/2, laplace/laplace2 -- the Laplace scale was RENAMED b -> sigma and laplace2 carries the rate lambda = 1/sigma, the lasso-friendly form, both 2026-08-08) -- with exact score and information, closed-form moments where they exist, `reparametrize()` (Faa di Bruno over partitions, map partials as hand-written keyed tables via `map_derivs`, one stencil per partial as the fallback), `folded()`, wrappers, transformations, MLE; **4 multivariate families** -- gaussian and Student t, whose matrix parameter comes from `parameters7`, plus Dirichlet and multinomial, whose simplex parameter does. expectation() and the cdf fallback run on numericals7's batched engines |
| `optimizers7` | 12 algorithms as objects — newton, bfgs, lbfgs, cg, bb, gd, adam, nelder_mead, compass, bundle, prox_grad (2026-08-08: FISTA with backtracking and adaptive restart, the only method written in R rather than C++ since every iteration calls caller-supplied R functions), multistart — with composable stopping rules, self-reporting safeguards, box bounds removed by reparametrization, starting values that need not be written out, and multistart parallel by default |
| `basis7` | bases as objects: evaluation, derivatives of any order, the integral anchored at the lower endpoint, and exact Gram matrices against a choice of measure. B-splines, Fourier and Legendre; one `TransformedBasis` wrapper for orthonormalization, constraints and the Demmler-Reinsch construction; `tensor_basis()` for several variables, with `basis_contract()` computing what a fit needs without forming the product; numerical fallbacks make an evaluation-only basis complete |
| `parameters7` | constrained parameters as maps from an unconstrained vector, exact to 4th order (RENAMED from covstructs7, 2026-08-04): base class `parameter` + SPD branch `matrix_parameter` (rank, log-(pseudo-)determinant, solves); `log_cholesky()`, `matrix_log()`, `correlation_matrix()` (spherical chart), `compound_symmetry()`/`ar1()` (two free values at any p, closed separable logdet, closed inverse), `autoregressive(p, order)` (PACF chart, derivative arrays through Levinson-Durbin in Rcpp, banded precision) (logdet = tr(S) linear, inverse = expm(-S) exact, Frechet derivatives by Daleckii-Krein/Opitz), `diagonal_matrix()`/`scalar_matrix()` (linkfunctions7 links), `scaled_matrix()` (rank-deficient ADMITTED), `simplex()` (ALR, cumulant recursion), `transition_matrix()` (row-wise simplexes); `kron_identity(structure, m)` (0.8.0, 2026-08-08) -- m identical diagonal blocks sharing one free vector, every contract quantity a linear lift (I_m kron A_k, m*logdet, blockwise solve), the FIRST composition wrapper, built for grouped random effects. `piano_parameters7.txt` supersedes piano_covstructs7.txt |

| `penalties7` | penalties as objects (created 2026-08-06 from `piano_penalties7.txt`): rho(D beta; theta) with value (NORMALIZING CONSTANT KEPT -- Giovanni 2026-08-03), exact derivatives in beta and theta, the mixed block (consumes `distrib_cross_y`, closed for all continuous families the same day), kink set, links on the hyperparameters. Three branches: `quadratic_penalty()` (rank/null basis/log-pdet fixed at ONE eigendecomposition, the REML pieces), `distrib_penalty()` (a univariate distributions7 log-density coordinatewise; `ridge_penalty()`/`lasso_penalty()`/`heavy_penalty()` named instances; ridge pinned against its quadratic twin at machine precision), `scad_penalty()`/`mcp_penalty()` (defined by rho', improper by construction). `check_penalty()` with injections. 0.4.0 (2026-08-08) adds `additive_penalty()`, the FOURTH branch: a sum of quadratics with a smoothing parameter on each, which is what an anisotropic tensor smooth needs; the log pseudo-determinant and its two derivatives come from one eigendecomposition of the sum (d/dlk = tr(S+ Pk), d2/dlkdll = -tr(S+ Pk S+ Pl)), and THE RANK IS FIXED AT CONSTRUCTION from the components stacked and normalized -- the null space of a sum of PSD matrices is the intersection of theirs and does not move with the parameters, while a count taken from the assembled S(lambda) falls as they spread apart (the measured trap of section 7, now pinned by a test). 0.3.0 adds `penalty_prox()`/`has_prox()`: one linear solve for the quadratic and structured branches (ANY map), closed forms for the Gaussian and Laplace instances, the closed piecewise operators of SCAD and MCP over their convex regions (a step past t = a-1 / t = gamma is REJECTED, the operator being set-valued there), and a coordinatewise root of (b-v)/t = l^(y)(b) for any other separable parent; a separable penalty under a general map is rejected, that being the generalized-lasso problem rather than a different formula |
| `modelterms7` | model terms as S7 objects (phase 1 shipped 2026-08-08 from `piano_modelterms7.txt`): `model_term` -> `additive_term`/`structural_term` (the second RESERVED for gas() and routed but refused), `linpar()` with a blueprint (terms/xlev/contrasts) so `term_predict()` reapplies rather than rebuilds, `interpret_formula()` with RECOGNITION BY EVALUATION (a call whose value inherits model_term is a term; log(x) stays a covariate; bare covariates collapse into one linpar; the intercept convention is the formula's), `cens()`/`censored_response` (statuses observed/left/right/interval derived from the values), `check_term()` whose subset check DROPLEVELS the subset (a plain row subset keeps unused factor levels and cannot expose a rebuild-from-newdata predict). Phase 2 (same day): `ridge()/lasso()/scad()/mcp()` over formula (intercept removed, linpar blueprint) or matrix input -- a matrix is predicted by re-evaluating its expression in the NEW DATA ONLY (a build-env fallback would silently reuse the build-time rows when the counts coincide), so the intended pattern is a matrix column of the data frame; the penalty is attached at build and `term_smooth` reads `penalty_kinks()` at a probe theta (midpoint-of-bounds, the reparametrize() probe rule); `by =` reserved in all four signatures. Phases 3-4 (same day): `edf()` counting per penalty (exact for linpar, tr[(H+S)^-1 H] on the block for a smooth penalty with S from `penalty_hessian()`, nonzero count for lasso/scad/mcp after Zou-Hastie-Tibshirani 2007), `print` showing a built penalized term's penalty, `plot` at supplied coefficients; `random(~ 1 | g)` = the indicator block with the effect distribution as penalty (ridge by default -- the random intercept IS the ridge -- a parameters7 PRECISION via structured_penalty, or a distributions7 object via distrib_penalty with the joint-mode reading documented). Random SLOPES landed the same day (Giovanni: intercepts AND slopes, correlated or not): `random(~ x | g)` interacts the within-group design with the indicators GROUP-MAJOR (the order I_m kron S assumes), default gaussian unstructured (log_cholesky(d)) or diagonal by `correlated =`, per-group `precision =` replicated by the NEW parameters7::kron_identity(), `distrib =` coordinatewise. 0.6.0 (2026-08-08) adds the STRUCTURAL contract and the SMOOTHS: `term_params()`/`term_links()`/`term_filter()` say what a term that rewrites the likelihood must provide, and `gas()` implements score-driven dynamics -- the filter returns the predictor AND its exact derivative in the term's parameters, propagated alongside the state, because the recursion is the only place it can be computed; the persistence rides a PARTIAL AUTOCORRELATION chart (the stationary region is not a box, the same argument as parameters7's autoregressive()) with Levinson-Durbin and its Jacobian written out. `s()` is the Demmler-Reinsch smooth: a linear column plus the reparametrized basis, penalized by diag(0,1,...,1) -- rank deficient by exactly one because dr_basis makes T'PT the identity -- so edf runs from k down to ONE, a straight line and not a constant. `te()` is the tensor product with the summed marginal penalties (⚠️ tensor_basis varies the FIRST margin fastest, so the Kronecker product is taken over the REVERSED blocks), ANISOTROPIC by default since 0.7.0, one lambda per margin through penalties7's `additive_penalty()`, with `anisotropic = FALSE` summing them under one parameter instead. `nl()` (0.8.0) is the nonlinear parametric term, and Giovanni's observation is what shapes it: while f is smooth the design block IS the Jacobian, so it stays an additive_term and only needs `term_refresh(built, coef)` to recompute the block as the parameters move, with `term_value()` reporting f itself -- the pair a Gauss-Newton step needs, and a test runs that iteration to the truth. TWO INPUT ROUTES WITH DIFFERENT REACH, which the documentation states rather than blurs: a FORMULA is read symbolically (`stats::deriv()`, falling back to a central difference where deriv cannot read the expression -- besselJ is the test case), and an opaque FUNCTION f(x, theta) is always differenced. Links per parameter, and `subformulas` developing a parameter as g^-1(Z gamma), are on the FORMULA route only: substituting inside f requires knowing where the parameter enters, which a function does not say. The gas filter was COMPILED in 0.7.0: measured, the two R callbacks (the score and the curvature, which belong to the model's distribution) are 17-27% of the loop's time and the arithmetic around them the rest, so the port buys 2.3x at p=q=1 and 3.2x at p=q=2, flat in n to 50000; `gas_filter_r()` survives as the twin the kernel is held to at machine precision. Both smooths take `by`: a factor gives one copy per level with the penalty repeated blockwise, a numeric gives a varying coefficient. Phase 5 (2026-08-08): book chapter 8 with `assert_terms_ok()` (injection-checked twice), CI matrix + coverage + pkgdown workflows, Pages enabled, README, logo deployed (make-logos.R now includes it in the loop), portal card activated, eighth member of the meta-package |
| `statmodels7` | the meta-package (2026-08-05, Giovanni asked for a tidyverse-style grouping). Installing it installs the five members and `library(statmodels7)` attaches them, reporting versions. `statmodels7_packages()`, `statmodels7_versions()`, `statmodels7_conflicts()`, `statmodels7_update()`. It is ALSO the destination package below: the modeling code lands here later, so nothing gets renamed |

**Planned** — `modelterms7`, `penalties7`, and the modeling layer of
`statmodels7` itself, which assembles everything into models. That layer is the
destination: a GAMLSS-like framework but far more organized, where

- every distribution parameter can be modeled, not just the mean;
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
    numericals7\
    linkfunctions7\
    distributions7\
    optimizers7\
    basis7\
    parameters7\
    modelterms7\
    statmodels7\     the meta-package, its own repository; note the directory
                     has the same name as the umbrella it sits inside
    book\            the Quarto book (see §9); `quarto render` from inside it
    logo\            hex logos: hand-authored SVG in svg/ (written by
                     author-svg.py); make-logos.R only rasterizes, run from here
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

GitHub: `github.com/statmodels7/{linkfunctions7,distributions7,optimizers7,basis7,parameters7,modelterms7,statmodels7}`, all on `master`. The repositories were
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
| parameters7 | `statmodels7.github.io/parameters7` — pkgdown, from `gh-pages` |
| modelterms7 | `statmodels7.github.io/modelterms7` — pkgdown, from `gh-pages` |
| statmodels7 | `statmodels7.github.io/statmodels7` — pkgdown, from `gh-pages`. Note this is the PACKAGE's site, distinct from the portal above |

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

The asymmetry runs in BOTH directions (2026-08-06): a local
`checking Rd \usage sections ... OK` did not cover what CI flags either --
the `map_derivs`/`reparam_derivs` arguments missing from their Rd passed a
local `--as-cran` on this machine and failed all five CI platforms. After
adding an argument to an exported function, grep its man page for the
argument before pushing; do not rely on the local usage check.

The converse also holds and is worth knowing before chasing a phantom.
`--as-cran` locally reports a **CRAN incoming feasibility** warning that CI
never will, because that check needs the network and the action disables it.
`distributions7` shows one here for its `Remotes` field while its CI is green
and its workflow is identical to the others'. A local warning is a red push only
if it is not that one.

Vignettes need pandoc: `$env:RSTUDIO_PANDOC="C:\Users\giova\AppData\Local\Programs\Quarto\bin\tools"`.

⚠️ **Set that variable when running the test suite too, or the documentation
guard silently skips.** `test-docs.R`'s first check calls
`pkgdown::check_pkgdown()`, which is what catches an exported object missing
from `_pkgdown.yml` -- and it is guarded by `skip_if_not(rmarkdown::pandoc_available())`,
so on this machine it reports `S` and tests nothing. That is precisely how
`optimizer_bounded` reached CI unindexed on 2026-08-08 and reddened pkgdown
minutes after the push. With the variable set, all six packages' guards run
and pass. The rule: after exporting anything, run the docs guard WITH pandoc
before pushing, not the ordinary suite.

### Windows traps that cost time

- **`install.packages` fails with "cannot remove earlier installation, is it in use?"**
  when an RStudio session has the package loaded — it holds the DLL. Ask Giovanni to run
  `detach("package:distributions7", unload = TRUE)` or restart his R session.
  `testthat::test_local()` is unaffected: pkgload compiles and loads from source.
  ⚠️ **A removal that fails PART WAY leaves a corrupt installation**, and that
  is worse than a failed install: `optimizers7` was left as `libs/` alone,
  with no DESCRIPTION and no `R/`, so `packageVersion()` reported it ABSENT
  while the directory existed. The symptom surfaced three steps away — a
  `quarto render` died on *"The package optimizers7 is required"* raised by
  pkgload while loading **distributions7's** imports, naming the consumer
  rather than the corruption. When a package that was working reports as
  absent, look at its installed directory before believing anything else;
  re-installing it is the whole fix.
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
- **Two processes that pkgload the same working tree compile the same `src/`
  and collide** (2026-08-08). A `roxygenise()` batch and a `quarto render`
  launched in parallel both compiled `optimizers7/src/`; the linker lost the
  race on the shared `.o`/`.dll` files and BOTH runs died (`ld returned 1
  exit status` on one side, a dead render on the other). Anything that loads
  the packages from source — roxygenise, `test_local()`, the book render —
  runs SEQUENTIALLY against this tree, never in parallel with another such
  process. After a collision, delete the leftover `src/*.o`/`*.dll` before
  retrying.

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

**Parameters** travel as a named list `theta`. Every generic normalizes it through the
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

**26 distributions**: gaussian, cauchy, logistic, student_t, laplace, pseudohuber, gamma,
invgauss, lognormal, beta, bernoulli, binomial, poisson, negbin -- the original
fourteen, all with closed-form observed derivatives to 4th order (Rcpp kernels
in `src/*_hd.cpp`) -- plus **weibull, gumbel, skewnormal, skewt** (2026-08-04),
and **exponential, geometric, chisq, betabinom, negbin1, gpd, gengamma,
vonmises** (2026-08-05). Of the eight added last, the first three are analytic
to 4th order observed AND expected; the other five carry closed-form score and
information (expected in closed form for the gpd, the gengamma and the von
Mises, an exact finite sum for the beta-binomial, a series against the exact
mass for NB1) and take orders 3-4 from the verified fallback. See section 7
for what each one taught.
Since later that day the new families and the laplace carry closed-form
observed 3rd/4th orders too (`weibull_hd.cpp`, `gumbel_hd.cpp`,
`skewnormal_hd.cpp`, `laplace_hd.cpp`), so **every univariate family is
analytic to 4th order** except the skew t's nu components, which cannot be.
Weibull and gumbel have closed EXPECTED 3rd/4th orders as well: every
expectation is a derivative of Gamma at 2, i.e. `E[u (log u)^k] =
Gamma^(k)(2)`, assembled from polygammas at 2 by the moment-cumulant
relations; the skew normal's expected values share the obstruction of its
expected information and stay numerical. The four new families keep orders
1-2 in vectorized R -- measured, a port would buy nothing: the gaussian
kernel costs 2.2 ms per gradient at n = 1e5 against weibull's 7.2, but
weibull's own deriv4, already Rcpp, costs 12.2 -- the gap is the
transcendentals per element, which C++ pays identically.

The four new ones, and what is worth knowing about each:

- **weibull** `(mu = scale, sigma = shape)`, gamlss's `WEI`. `mu` is NOT the
  mean -- that is `mu*Gamma(1+1/sigma)`, and a mean parametrization would make
  every derivative a derivative of the gamma function and of its inverse.
- **gumbel** `(mu, sigma)`, for maxima. Location-scale with a FIXED shape: its
  skewness `12*sqrt(6)*zeta(3)/pi^3 = 1.1395` and excess kurtosis `12/5` are
  constants. Note that `E[l_mu_sigma]` does NOT vanish -- the density is
  skewed, so location and scale are not orthogonal, unlike in a symmetric
  location-scale family.
- The two share an expectation. In the Weibull `u = (Y/mu)^sigma` is standard
  exponential whatever the parameters, and in the Gumbel so is `w = exp(-Z)`,
  so EVERY expectation either needs is a derivative of `Gamma` at 2:
  `E[u] = 1`, `E[u log u] = 1-gamma`, `E[u (log u)^2] = (1-gamma)^2 + pi^2/6 - 1`.
  Both expected informations are therefore closed form in one line. The two
  families are also one another: `exp(-Gumbel)` is Weibull.
- **skewnormal** `(mu, sigma, alpha)`, Azzalini. Everything is written in the
  inverse Mills ratio `R = phi/Phi` and `R' = -R(t+R)`; `mills_ratio()` forms it
  ON THE LOG SCALE, because below about `t = -38` both phi and Phi underflow
  while the ratio is finite and close to `-t`. The cdf uses **Owen's T**
  (`owen_t()`), one bounded 1-D quadrature per point, which beats the base
  class's semi-infinite integral of the density. Two things to know: the
  expected information is **singular at alpha = 0**, which is a property of the
  parametrization and not a defect, and the skewness the family can reach is
  bounded by 0.9953 -- which is the reason the skew t exists.
- **skewt** `(mu, sigma, alpha, nu)`, the four-parameter family a
  location-scale-shape framework actually wants. Score and Hessian are closed
  form in `(mu, sigma, alpha)`; everything involving `nu` is NOT, because the
  density carries `T_{nu+1}` and the derivative of a Student t cdf in its
  degrees of freedom has no elementary form -- the same obstruction as the
  gamma and beta cdf in their shape. Those components come from ONE five-point
  stencil (`fd5_first`, `fd5_second`) applied to an analytic quantity: to the
  log-density for the `nu` score and `nu_nu`, and to the CLOSED-FORM score for
  the mixed `(par, nu)` entries. `method = "newton"` is much the cheaper way to
  fit it, since the expected information would be quadrature per component.

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
are constants, like a binomial's `size`), and instead adds a θ-dependent normalizing
constant `Z(θ) = F(u) − F(ℓ⁻)`. With `m_i = E_T[s_i]` and `M_ij = E_T[H_ij + s_i s_j]`,
both expectations under the *truncated* law,

```
d_i   l_T = s_i(y) - m_i                      (the parent's score, recentered)
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
distinct block costs one quadrature, memoized across the partition sum.

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
package parametrizes by `(mu, sigma2)` and both of those involve the shape, so neither
is closed form. What does survive is an exact identity, worth keeping as a test:

```
dF/dmu + (2 sigma2/mu) dF/dsigma2 = -y f(y) / mu
```

the shape direction canceling from the combination. The FD fallback delivers ~2e-10
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
the y-derivative does not interact with a reparametrization of theta. The
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

**`fit_distrib()` delegates its optimization to `optimizers7`** (2026-08-03), and
distributions7 now Imports it. Until then the function did the exact thing the
toolkit's positioning condemns in other packages: a hand-written scoring loop with
its own step halving, the convergence test spelled out inline, and `stats::optim`
for the BFGS branch — an *external* optimizer while the sibling package sat unused.
Three things made the translation nearly literal:

- **Fisher scoring is `newton()` with the expected information passed as `he`.**
  It needs no implementation of its own; `"fisher"` and `"newton"` differ only in
  that one argument.
- the inline test `max|U| < tol || |Δl| < tol(|l|+tol)` *is*
  `crit_any(crit_grad(tol), crit_rel_obj(tol))` — **replaced on 2026-08-04 by
  `crit_grad(tol)` alone** (Giovanni, explicit). The OR was the weaker of the
  two rules in practice: the objective is flat near the maximum, so
  `crit_rel_obj` fired first and the run stopped at whatever gradient it
  happened to have. Dropping it costs iterations (the 2-d gaussian of the book
  goes from a handful to 78) and buys a stationary point: the summed score
  per observation lands at 1e-16 where before it was platform-dependent
  between 1e-15 and 1e-8. A test that asserted digits of the *point* rather
  than of the *objective* was what exposed the difference, and it had already
  gone red on macOS once;
- `method` now also accepts **any optimizer object**, used as given.

**The objective handed to the optimizer is `-l(eta)/n`** (2026-08-04, Giovanni's
proposal), with the gradient and the Hessian divided by `n` too. The maximizer is
unchanged and so is every Newton step, since the factor cancels in `H^-1 g`; what
changes is what a *threshold* means. `crit_grad(tol)` on a summed score asks a
sample of ten million for an accuracy per observation ten million times finer than
it asks of a sample of ten, which is why an absolute `1e-10` was unreachable on
macOS and green everywhere else. My first fix was `crit_grad(tol * n)` and it was
worse than it looks: **a criterion the caller supplies gets no such scaling**, so
`method = lbfgs(criterion = crit_grad())` — which Giovanni actually ran — kept the
whole problem. Scaling the objective fixes both at once, and it is the only place
that can. `ll_hat` and `I_eta` are recomputed **unscaled** at the optimum, so
`logLik`, AIC/BIC and every standard error are untouched (verified: gap to the
closed-form maximum 0, `se(mu)` against `sigma/sqrt(n)` to 1.2e-16). The book's
§3.3 states the averaged rule and `.certify_objective_scale()` pins it, injection
-checked against standard errors left on the averaged scale (a factor of 70.7 at
n = 5000) and against a 5% error.

Two things improved rather than merely moved: the line search requires sufficient
decrease instead of mere non-decrease, and a non-PD Hessian is repaired by flooring
its eigenvalues instead of abandoning the start, which `solve()` used to force.
The BFGS fallback is kept for `"fisher"` and `"newton"` only — silently replacing
an optimizer the caller *chose* would report a fit obtained by a different method.

⚠️ **A catch-all for numerical failures must not swallow a configuration error.**
The restart loop's `tryCatch(..., error = function(e) NULL)` exists to absorb a
divergent quadrature in a numerically approximated expected information. It was
also absorbing optimizers7's *refusal* of a stopping rule the method cannot
evaluate, and reporting it as `"Optimization failed from every starting value;
supply 'start'"` — naming the wrong cause entirely. `check_criterion()` is exported
by optimizers7 for exactly this and is now called before the loop. The general
shape: when a `tryCatch` is there for one class of failure, check what else reaches
it. This was found by writing the test for the new feature, not by using it.

**`fisher_scoring()`** (2026-08-04, Giovanni, after correcting an earlier
attempt). `fit_distrib()` takes ONE argument saying how to optimize, and it
takes either an `optimizers7` optimizer or `fisher_scoring(approx =, nsim =,
criterion =, maxit =)`. The loose `approx`/`nsim` arguments are gone: how the
expected information is approximated is a property of Fisher scoring, and had
no business sitting next to optimizers that never look at it. Fisher scoring
is not an algorithm of its own -- it is `newton()` with one matrix replaced --
which is exactly why it is an object and not an optimizer:

| `method =` | what it does |
|---|---|
| `fisher_scoring()` (default) | Newton with the **expected** information |
| `optimizers7::newton()` | Newton with the **observed** Hessian |
| `optimizers7::lbfgs()` etc. | whatever that optimizer does |
| `"fisher"` / `"newton"` / `"bfgs"` | kept as short names |

A strategy chosen where it would be ignored is refused, through
`has_exact_expected_hessian()` -- see section 7 for the shadowing bug that made
that predicate answer backwards.

**`maxit` and `tol` left `fit_distrib()`'s signature on 2026-08-04**
(Giovanni: they *"dovrebbero essere tutti governati dagli optimizers"*), and
he was right for a sharper reason than tidiness -- with `method = <an
optimizer>` they were **silently ignored**, so a call setting both the
optimizer's `maxit` and the fit's got no complaint and no effect from the
second. The budget and the stopping rule now live on the method: on an
optimizer object, or on `fisher_scoring(criterion =, maxit =)`, and
otherwise at `crit_grad()`'s and the optimizer's own defaults. The internal
BFGS fallback inherits whatever the chosen method set. See section 7.

**Multivariate distributions** (2026-08-03/04). Base class `multivariate_distrib`
sits beside `continuous_distrib`/`discrete_distrib` and **not under** it: the
one-dimensional defaults registered there — a cdf by quadrature, a quantile by
root finding — have no counterpart, and a cdf on R^p is an orthant probability
while a quantile needs an ordering of R^p. Both are refused rather than
approximated.

The design decision that makes everything else free: **the matrix parameter is
flattened into scalars**. A p-variate family's parameters are `mu1..mup`
followed by the `free_names` of a `parameters7` structure, every one of them a
scalar with an identity link, so `align_theta`, `deriv_names`, `hess_names`, the
link scale and `fit_distrib` need no special case at all. The constraint lives
in the structure, where it belongs — it is a constraint on the matrix as a
whole and a scalar link cannot express it.

- `mvgaussian_distrib(n_dim, sigma=, omega=)` — **one** constructor,
  the two arguments mutually exclusive, not two constructors. Score
  `-½∂_k log|Σ| + ½ w'A_k w` with `w = Σ⁻¹r`; the expected information is
  `-½tr(Σ⁻¹A_kΣ⁻¹A_l)` and needs **no `A_kl`**, the mixed mean-matrix block being
  exactly 0 since `E[w] = 0`. For the precision form `∂Σ/∂η = -ΣA_kΣ`.
- `mvstudent_t_distrib(n_dim, sigma, link_nu)` — the gaussian's score
  with every data term multiplied by `c = (ν+p)/(ν+q)`, plus a `ν` component;
  the only multivariate family whose link scale is **not** its parameter scale.
  `mv_sigma()` is the **scale** matrix and `variance()` the covariance
  (`νΣ/(ν−2)`, infinite for `ν ≤ 2`) — keeping the two apart is what lets it be
  fitted where the second moment does not exist. Marginals keep the **same** ν.
- `dirichlet_distrib(n_dim, mean = simplex(n_dim), link_phi = log_link())` —
  the first family here that is **not elliptical**: no location and scale to
  separate, the support a set of dimension `p-1`, the covariance singular by
  construction. Parameters `mean_alr1..alr(p-1)` and `phi`; shapes are
  `alpha = phi*mu`. Score `phi*A'g` and `psi(phi) + mu'g` with
  `g_j = log y_j - psi(alpha_j)`. Its expected information is **closed form**
  for one reason used twice: differentiating `sum(mu) = 1` once makes every
  column of `A` sum to zero and again makes every `param_d2` vector do so, and
  `E[g] = 0` kills the rest. Marginals are **beta with the same `phi`**, the
  concentration being shared — which is why `mv_marginal()` is a real method
  here rather than a refusal.
- `multinomial_distrib(n_dim, size, probs = simplex(n_dim))` — the first
  multivariate **discrete** family. Expected information
  `-n*sum(A_jk A_jl / p_j)`, the second-derivative term vanishing by the same
  zero sum. Marginals are binomial.
- **`mv_support()`** (new generic; base class refuses) returns the finite
  support as a matrix, built by `compositions(n, k)`, the weak compositions.
  What it buys is that **every expectation is an exact sum**: the mass over
  the support comes back 1 to 7.8e-16 and the closed-form information agrees
  with the summed observed Hessian to 3e-16, where an importance-sampling
  check would only ever compare against Monte Carlo error. A normalization
  wrong by a thousandth is caught by the sum and would not be by the sample.
- **`mv_reference_draw()`** (new generic) supplies the proposal
  `check_distrib()` integrates the density against. The base method is the
  inflated gaussian it always used; the Dirichlet overrides with the uniform
  on the simplex (constant density `Gamma(p)`).
  ⚠️ The base proposal does **not** fail loudly on the Dirichlet:
  `chol()` accepts the singular covariance — the §7 note that *`chol()` is
  not a rank test*, met a third time — and the estimate of an integral that
  is 1 comes back **2.0e-08**. A silently wrong number, not an error.
- `check_distrib()` runs a **nine-check multivariate battery** (both Bartlett
  identities, the moments, the normalization by the route above). Two checks
  are **emitted only when they apply**, following the univariate convention of
  twelve checks on a discrete family against thirteen on a continuous one:
  the response derivatives are skipped for a family with an enumerable support
  (discrete, so there is none) and for one that has not registered
  `distrib_grad_y` (the base class refuses it *by design*, so not registering
  it is a choice, not a gap). `has_mv_support()` and `has_mv_grad_y()` are the
  predicates, and they ask the **method's owning class**, since the
  multivariate branch sits beside `continuous_distrib`/`discrete_distrib`
  rather than under either and there is no class to test.
- `mv_marginal()` is a generic; the base class **refuses**, because a numerical
  marginal would be a different object with the same name. It exists because a
  panel of a pairs plot *is* a marginal.
- `plot()` on a multivariate distribution or fit draws a panel matrix: marginal
  density on the diagonal (plus a kernel estimate when there are data), contours
  below, correlation above. Refused above three coordinates, with `which=` /
  `mv_which=` offered instead.
- `n_obs()` is the row count, not the length. A parameter may **not** vary by
  observation here: the structure describes one matrix for the whole sample.
  ⚠️ **And `length(y)` was used for it once, in the link-scale branch of
  `distrib_expected_hessian` (`R/generics.R`)**, where a zero vector is built
  to stand for a parameter that does not vary. For a matrix response
  `length(y)` counts entries, so that vector came out `n*p` long and recycled
  against the `p`-long components; the first-order term of the order-2 chain
  rule appears only on the diagonal, so **every diagonal entry of the
  information was inflated by `p` and every standard error of a multivariate
  fit came back a factor of `sqrt(p)` too small** — measured ratios 0.70, 0.56
  and 0.48 at p = 2, 3, 4 against the exact `sd_j/sqrt(2n)`, now 1.000. The
  arithmetic was confirmed by simulation (0.02589 sampled against 0.02506
  theoretical) before touching anything. A quiet factor of `sqrt(p)` is what a
  wrong length looks like in a language that recycles: **anywhere a matrix
  response meets code written for a vector, `length()` is a defect and
  `n_obs()` is the question.** It also mattered for the fitting: the corrupted
  expected information made Fisher scoring take poor steps from a distant
  start, and the iris fit that used to run 500 iterations without arriving now
  converges from the origin in 42.
- **The structure's free names are PREFIXED by the matrix they describe**
  (2026-08-04, Giovanni): `sigma_log_L1` for a covariance or a scale matrix,
  `omega_log_L1` for a precision. A free name says how the matrix is BUILT and
  not WHICH matrix it is, so the same structure on the two sides of a model
  gave identical parameter names for two genuinely different models. The
  distribution applies the prefix, not the structure -- the structure does not
  know which side it has been handed to. `mv_prefixed_names()` in
  `multivariate.R`. **The names themselves changed on 2026-08-04**: a free
  name now says which TRANSFORM produced the coordinate, so an AR(2)
  covariance reports `sigma_log_scale`, `sigma_z_pacf1`, `sigma_z_pacf2`.
  See section 7 -- the old names promised bounded quantities and reported
  free values.
- **A structured matrix reports its own quantities as a block**
  (`mv_param_block()`, 2026-08-04). `parameters7::param_readable()` declares
  what a family is about with the Jacobian from its free vector, and
  `mv_derived()` appends it to the standard deviations and correlations, its
  Jacobian widened by placing the structure's columns in the stretch of
  `distrib@params` they occupy. Two things to know: the label says
  `(precision)` when the distribution is inverted, since the structure does
  not know which side it is on; and a multivariate family written OUTSIDE the
  package has no `@param` property at all, so the property is asked for with
  `S7::prop_names()` rather than assumed -- which the existing user-defined
  test caught immediately.
- **A fit reports the quantities a reader reads, not the coordinates**
  (`mv_summary()`, `mv_derived()`, `R/mv_summary.R`). Nobody reads the
  logarithm of a diagonal entry of a Cholesky factor; what a fitted
  multivariate gaussian is about is `Sigma = D R D`, so the standard deviations
  and correlations get the estimate, the standard error and the interval, and
  `print()` shows them in blocks. Three parts worth keeping:
  - the standard errors are the **delta method**, `J V J'`, with `J` closed
    form: `d s_j/d eta = A_k[j,j]/(2 s_j)` and
    `d rho_jk/d eta = A[j,k]/(s_j s_k) - (rho/2)(A[j,j]/S_jj + A[k,k]/S_kk)`,
    where `A_k` is the structure's own `param_d1` -- nothing new is
    computed;
  - each **interval is built on the scale that keeps the quantity in its own
    set** and mapped back, exactly as `fit_distrib()` does for a univariate
    parameter: log for a standard deviation, Fisher's z for a correlation. On
    the raw scale a correlation's interval routinely exceeds 1;
  - a **precision** parametrization reports the same standard deviations and
    correlations (they are properties of the law) and ADDS the readings that
    are its own: the conditional standard deviations `1/sqrt(Omega_jj)` and the
    partial correlations `-Omega_jk/sqrt(Omega_jj Omega_kk)`. At `p = 2` the
    partial correlation IS the correlation, so it is not printed twice.
  The base-class method returns the distinct entries of `mv_sigma()` named
  `sigma_v1_v2`, with a numerical Jacobian, so a family that decomposes into
  nothing still reports its matrix on its own scale rather than a coordinate.
  The **t** names its diagonal quantities `scale_sd_` -- the matrix is a scale
  matrix and those are not standard deviations of the response -- while its
  correlations need no qualification, a positive multiple of a matrix having
  the same correlations.

---

## 5. Working preferences

- **No `Co-Authored-By: Claude` trailer** on commits. Giovanni is the sole author of what
  he publishes.
- **No code generators that require Mathematica.** Wolfram output is hand-transcribed
  into Rcpp and validated numerically. He rejected a generator explicitly: *"non ha senso"*
  that the package only builds because Mathematica happens to be installed.
- **Commit messages are short and factual** (Giovanni, 2026-08-04, explicit, and
  he had the whole history rewritten to match). A subject line saying what was
  done, in the imperative, under about 65 characters; a body of a few plain
  lines listing the changes, and the measured numbers where they matter. No
  rhetorical subjects, no two-beat contrasts, no capitals for emphasis, no
  narrative of how the work went. "Add multistart with Latin hypercube starts",
  not "Multi-start, and the answer to Q6 turns out to be the opposite of the
  premise". This supersedes the earlier rule that messages should explain why
  rather than what: the reasoning belongs in this file and in the code
  comments, where it is read, not in the history.
- **Do not wait for CI, or for `R CMD build`/`check --as-cran`, between
  queued changes** (Giovanni, 2026-08-05, twice). When several changes are
  lined up, run `testthat::test_local()` during the work and leave the build,
  the `--as-cran` check and the CI to ONE pass at the end of the batch:
  waiting on a five-platform matrix, or on a build that recompiles every
  kernel, between steps is dead time. The local suite catches almost
  everything; what `--as-cran` adds is the documentation checks, and those do
  not change between one family and the next. The rule that a red coverage job
  must be read still stands -- it is about reading the result, not about
  blocking on it.
- **No jets, explicit derivatives, Rcpp when it measures faster** (Giovanni,
  2026-08-06, emphatic: "abbasso il jet, viva implementazioni esplicite in
  rcpp"). Jet-based differentiation is BANNED from production paths: the PIG
  kernels measured its fixed composition overhead at 2x-36x the hand-written
  closed forms. Derivatives are written out (Faa di Bruno per component),
  validated against an independent route -- one numDeriv pass per order, or a
  mechanical twin kept in tests only, like the pig*_hd_jet_cpp kernels. When
  porting, benchmark R against Rcpp on the real workload and prefer Rcpp
  "anche se basso" -- any measured advantage decides for it. External
  packages such as gamlss.dist are not dependencies even in Suggests;
  references are implemented in-test (recursions, series, closed identities).
- **Verify, do not assume.** Several times a suspected bug turned out to be a bad test,
  and several times a "cosmetic" finding turned out to be a real defect. Measure first.
- Report failures faithfully, with the numbers.

---

## 6. Current state

| | |
|---|---|
| `numericals7` | 0.7.0 (2026-08-06): enumerations (0.1.0), stencils (0.2.0), quad_vec/series_vec (0.3.0), special functions (0.4.0), log-Bessel (0.5.0), jets REMOVED (0.6.0), log-Bessel kernels compiled with R twins in tests (0.7.0). `R CMD check --as-cran` clean apart from the environment notes, CI green. Logo tracked, favicons in place |
| `linkfunctions7` | 891 tests, `R CMD check` OK, CI green; stencils delegate to numericals7; the transcendental links' derivatives compiled 2026-08-06 (src/link_kernels.cpp, first compiled code) |
| `distributions7` | 2993 tests, `R CMD check --as-cran` clean apart from the submission notes (2026-08-05, local), CI green |
| `optimizers7` | 710 tests, `R CMD check --as-cran` OK with vignettes, two notes (2026-08-04). Published 2026-07-31; the Rcpp kernels had until then only ever been compiled by one compiler on one machine. |
| `parameters7` | renamed from covstructs7 on 2026-08-04 (clean cut; GitHub redirects the repo, the PAGES URL DOES NOT redirect). Version `0.8.0` (0.8.0: kron_identity, the first composition wrapper -- identical blocks sharing a free vector; 0.7.0: log-Cholesky Leibniz assembly compiled, 30x at p = 8 order 4; 0.6.0: AR derivative arrays in Rcpp, first compiled code). Earlier 0.3.0: base/matrix split, orders 3-4 everywhere, simplex, transition_matrix, matrix_log, phase 2's correlation_matrix/compound_symmetry/ar1/autoregressive, free names tagged by their transform, and `param_readable()`. Composition wrappers: kron_identity done (identical blocks); D R D', general block diagonals and sums still open. |
| `basis7` | 683 tests, `R CMD check --as-cran` clean apart from the two environment notes, CI green (2026-08-03). Version `0.3.1`, a `NEWS.md` from the first commit and a vignette. Phases 1 to 4 of `piano_basis7.txt` are done; phase 5 is the handoff to `penalties7` and `modelterms7`. |

| `penalties7` | 105 tests, version `0.4.0` (2026-08-08: the proximal operator, then the additive branch), created 2026-08-06; ALL FIVE PHASES of `piano_penalties7.txt` done by 2026-08-07: 0.2.0 adds `structured_penalty()` (a parameters7 matrix_parameter as the PRECISION, free vector = theta, identity links, every derivative from param_d1/param_d2 and the logdet contract; at a zero log-Cholesky free vector it IS the plain ridge, pinned at machine precision) and the book gained chapter 8 with three injection-checked gates (`book/R/penalty-certificates.R`). Repo pushed, Pages enabled |
| `modelterms7` | 272 tests, `R CMD check --as-cran` clean apart from the environment notes and the local-only `Remotes` warning, all five phases of `piano_modelterms7.txt` done 2026-08-08 and EVERY reserved term implemented -- `gas()`, `s()`, `te()`, `nl()`; first compiled code (src/gas_filter.cpp). Version `0.8.0`, `NEWS.md` from the first commit. CI, coverage, pkgdown and Pages live. Still reserved: `nl()`, whose block is the Jacobian refreshed per iteration |
| `statmodels7` | 32 tests, `R CMD check --as-cran` with one deliberate NOTE (see below) and the submission warning, created 2026-08-05. Version `0.1.0` with a `NEWS.md` from the first commit. |

All six repositories run `R-CMD-check` on macOS, Windows and three Linux/R
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
reasons a first submission comes back.

**And they were still missing throughout distributions7 until 2026-08-04**, which
is the point of the item rather than an aside: the gap was closed in the package
where it was found and nowhere else, exactly the failure mode below calls *a habit
that lives in one package out of three is not a habit*. The audit found 13 exported
topics with no `\value`, 61 with no executable example, and 101 internal pages with
no `\value` -- 174 in all. The guard is now `tests/testthat/test-docs.R`, **copied
into all five packages**, which asks four questions: does every exported object
appear in `_pkgdown.yml`, does every object in the namespace have a help topic, do
any two topics collide when case is ignored, and does every exported topic carry a
`\value` and an example.

⚠️ **"Copied into all five" was true of the file name and not of the file**
(corrected 2026-08-05). Three of the five carried weaker versions that could
not ask the second question at all, and seven undocumented objects were sitting
in exactly those three. The file is now generated from one source and is
byte-identical apart from the package name; section 7 has the census. It also
excludes a package landing page from the `\value` requirement, which every
package now has: `?linkfunctions7`, `?distributions7`, `?optimizers7`,
`?basis7`, `?parameters7` all reach a page carrying the title, the description,
the logo and the URLs from `DESCRIPTION`. Two packages had one and three did
not; the block is the standard `#' @keywords internal` before `"_PACKAGE"`, and
it is deliberately **not** listed in any `_pkgdown.yml` — pkgdown drops a
keyword-internal topic from the reference index without reporting it missing,
so listing it in one package and not the others was the asymmetry rather than
the fix.

Two things learned closing it. `S7::method(coef, cls) <- fn` at the top level of a
package **creates a binding `coef` in the namespace**, because the replacement form
expands to an assignment; those are S7's shims over the base generics and belong in
the test's exclusion list next to `print`/`plot`/`summary`/`mean`/`simulate`. And
adding a test that calls `pkgdown::check_pkgdown()` **needs pkgdown in `Suggests`**
-- the same *"'::' import not declared"* warning §3 already records, walked into a
second time, and then a THIRD time the same day for `numDeriv`, used in the
tests of the skew families. The rule is mechanical and worth applying without
thinking about it: **any package a test names with `::` goes in `Suggests`,
added in the same edit as the test.**

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
  free to be smaller than `f(0)` too, and the likelihood factorizes into a binary part
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
  again. **`log1p` and `expm1` are not micro-optimizations here; they are the difference
  between a number and an infinity.**
- **`numDeriv`'s Richardson stencil reaches ~8e-4·|x|.** Any grid that comes closer than
  that to a domain boundary is differentiated using points outside the domain, which come
  back NaN. This made `check_link` report failures for links that were exact to 1e-11.
- **Golden-section on a compactified scale is not accurate enough** for locating a mode.
  Its tolerance is expressed in the compactified variable, whose derivative can be
  enormous: with the tangent map `dy/dt ~ y²`, the default tolerance put the "mode" of a
  density centered at 1000 off by 125 standard deviations. `find_pdf_anchor` now refines a
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
  A **second instance of the same shape**, found while planning `parameters7`
  (2026-08-03): counting eigenvalues above a relative tolerance is not a rank
  test either. On a tensor-product penalty `l1*P1 + l2*P2` of true rank 28 out
  of 32, the count reads 28 while the two are comparable and **24** once the
  ratio reaches 1e10 -- the small contributions sink below the tolerance and
  are read as zeros -- while the null-space residual `max|M N| / max|M|` never
  moves off 3e-16 across the whole range. (Re-measured 2026-08-04 at a
  tolerance of 1e-10, which is what the book prints and its gate checks; an
  earlier draft of this note quoted 26 at 1e8 and 16 at 1e12, from a different
  tolerance.) Smoothing parameters ten orders of
  magnitude apart are an ordinary fitted model, not a pathology. So a rank
  comes from the **stacked, individually normalized components** (the null
  space of a sum of PSD matrices is the intersection of the null spaces) and
  membership is tested through the stored null basis, never by counting
  eigenvalues of the assembled matrix. General form of both instances: a
  membership question is about the matrix, not about whichever arithmetic was
  performed on it.
  How it surfaced is the second half of the lesson, and it is the **second**
  time in this toolkit that **only the coverage job** saw a defect while
  `R CMD check` stayed green on all five platforms — the first was the S7
  identity comparison above. The two have nothing in common mechanically, so
  the rule is not about `covr`: a green R-CMD-check matrix does not clear a
  change whose outcome depends on floating-point luck, because five platforms
  agreeing is five draws from the same coin. Read a red coverage job.

### Phase 2 of parameters7, and two claims the checks refused

Done 2026-08-04. `correlation_matrix()` uses the SPHERICAL chart
(Rapisarda-Brigo-Mercurio): rows of the Cholesky factor are unit-sphere
points in angular coordinates, so the unit diagonal and the definiteness are
structural and the angles are free through `bounded_link(0, pi)`. It was
chosen over the canonical partial-correlation recursion for one reason that
decides everything here: the recursion carries square roots and conditional
terms whose fourth derivatives are unwritable, while the spherical factor is
a product of sines and cosines each depending on ONE free value, so the
Leibniz machinery already in the package applies verbatim. Its
log-determinant is `2*sum(log(sin(theta)))` -- one term per free value, hence
separable, hence every mixed derivative exactly zero.

`compound_symmetry()` and `ar1()` are two free values at any p. Both have a
separable closed log-determinant and a closed inverse (Sherman-Morrison;
tridiagonal). The bound that matters: compound symmetry is definite only for
`rho > -1/(p-1)`, so its correlation link is `bounded_link(-1/(p-1), 1)` and
NOT `rhobit_link()` -- rhobit would hand a consumer an indefinite matrix at
an ordinary free value. AR(1) has no such dimension-dependent bound and does
use rhobit.

Two things the checks refused, both worth keeping:

- **`off[, c(2, 1)]` on a ONE-ROW index matrix collapses to a vector**, which
  R then reads as LINEAR positions rather than as a row-column pair. The
  AR(1) tridiagonal inverse therefore wrote `-rho` into `[1, 1]` at p = 2 and
  nowhere else -- invisible at p >= 3, where the index matrix has more than
  one row. `drop = FALSE` is the fix, and **sweeping the dimension down to
  its minimum is what found it**: the test now runs `check_parameter()` over
  p = 2:6 for all three families. The smallest legal case is where a matrix
  degenerates, and it is the case least likely to be tried by hand.
- **"a cross-row component of a correlation derivative is zero" was my
  claim, and it is false.** It holds for the FACTOR, whose rows are
  independent; it does not hold for `R = LL'`, whose entry (i, j) is the
  inner product of rows i and j. What is true is that such a component is
  SUPPORTED on entries (i, j) and (j, i) -- and is zero even there when the
  two angles sit beyond the columns the rows share, which is why the first
  correction ("non-zero at exactly those two entries") was refused too. The
  test asserts the support and then asserts, separately, that some cross
  component really is non-zero, so the claim cannot be satisfied vacuously.
  General shape: a structural claim about a product is not a claim about its
  factors, and writing the test before believing the sentence is what
  separates them.

### The non-smooth block: what the measurements decided (2026-08-08)

GAP 1 of `piano_modelterms7.txt`, closed. `penalty_prox()` in penalties7
and `prox_grad()` in optimizers7, then three comparisons on one lasso
objective, every route agreeing on the optimum to machine precision:

- **LLA is not needed.** The direct SCAD proximal operator beats the
  iterated weighted lasso on total inner iterations at every conditioning
  tried -- 26 against 45, 123 against 300, 370 against 390 -- at the same
  objective to 1e-16. A route that solves a sequence of surrogates does
  not beat one that solves the thing.
- **Coordinate descent is faster, and is NOT an optimizer.** 1.1x to 5.3x
  faster than prox_grad across condition numbers 2.8 to 7900, at a tighter
  stationarity (1e-10 against 1e-8). But it needs the columns of X and the
  running residual -- the MODEL, not fn and gr -- so it cannot live behind
  optimizers7's black-box interface. It belongs to statmodels7's fitting
  layer as a specialized path for a linear predictor with a separable
  penalty; prox_grad is the general route, for any smooth loss and any
  operator.
- **Acceleration is for ill-conditioned problems.** At a condition number
  of 3 the plain iteration wins narrowly (39 against 24); at 55 it is 4153
  against 126; at 480 the plain method does not converge in 50000
  iterations while the accelerated one takes 334.

⚠️ **An iteration COUNT on an ill-conditioned problem is platform
arithmetic, and a ratio calibrated here fails elsewhere.** The plain method
takes 4153 iterations on this machine and **652 on macOS** from the same
start on the same problem: the products are summed differently, the
trajectory diverges, and the count moves by a factor of six. The
accelerated one is far steadier (126 against 133), which is itself part of
the story. A test asserting the measured ratio went red on macOS alone;
what it asserts now is the structural claim with room (acceleration at
least halves the count) and it PRINTS both counts, per the `fit_report()`
lesson -- an assertion that prints only TRUE/FALSE cannot be diagnosed on a
platform you do not have.

Two defects the measurements exposed, both in the first version of the
loop, and both of a shape section 7 already records in other words:

- **a restart consumed an iteration of the budget.** The accelerated
  method then looked SLOWER than the plain one for no reason but the
  bookkeeping, and the trace had gaps. A retry that makes no progress must
  not spend the budget that measures progress.
- **the stationarity measure was read at the EXTRAPOLATED point, not at
  the iterate.** With momentum the two differ, so the mapping never
  vanished: the run circled the answer at 3.7e-9 for 20000 iterations
  reporting failure, while the plain variant reached exactly 0 in 49. Read
  at the iterate it agrees with the KKT violation to the digit (4.45e-07
  against 4.45e-07, 7.12e-09 against 7.12e-09), which is what says the
  reported quantity is the right one. It costs one extra gradient per
  iteration under acceleration, since the point a step is taken from is
  not the point being reported.

⚠️ **An S7 generic declared without a body gets a `...` in its usage, and
only CI calls it undocumented.** `S7::new_generic("f", "x")` generates
`function(x, ...)`, so `--as-cran` reports *"Undocumented arguments in Rd
file"* -- on all five CI platforms, while the local check said `checking Rd
\usage sections ... OK`. The sibling `optimizer_provides` passes an explicit
`function(optimizer) S7::S7_dispatch()` and has no dots; matching that shape
is the fix, and is better than documenting a `...` no method reads. Third
instance of the §3 asymmetry, and the second one this year on a `\usage`
section.

⚠️ **A test that names a package UP the dependency graph closes a cycle,
and only CI says so.** The first version put the "penalty_prox drives a
proximal run" test in optimizers7, which meant `penalties7` in its
`Suggests` -- and penalties7 depends on distributions7, which depends on
optimizers7. `R CMD check --as-cran` was clean here, because penalties7
is installed on this machine and a local check never resolves anything;
on CI **all three workflows failed at dependency setup**, before a single
test ran, with *"Can't find package called penalties7"*. The test moved
to penalties7, where the direction is right, and it is the better home
anyway: the question it asks is whether the OPERATOR works, not whether
the optimizer does. General shape: a package may only name packages
BELOW it, in Suggests as much as in Imports, and the local check cannot
see the violation.

⚠️ The attainable mapping obeys the rounding floor of §7's tolerance
item: with acceleration it stalls at ~4e-9 on an objective of order 2, so
`crit_grad(1e-10)` is unreachable there and asking for it reports failure
at the answer. The tests ask 1e-8.

**`optimizer_bounded()`** was added with it: a new generic, TRUE by
default, FALSE for `prox_grad`, which `check_optimizer()` consults before
testing box bounds. A proximal method takes its constraint INSIDE the
operator, where it composes with the term already there, so bounds beside
the objective would be a second and conflicting route to the same thing --
and the check would otherwise fail a correct optimizer, the shape §7
records as *separate what a component promises from how well it does it*.

### reparametrize() against a family written in full (measured 2026-08-08)

Giovanni asked whether reparametrize() should build new families or whether
they should be written out, and the measurement decides for writing out.
The same law three ways -- the shipped hand-written gaussian2 (Rcpp
kernels), a reparametrize() twin with hand map_derivs tables, and the same
twin on the stencil fallback -- at n = 1e5, min-of-5: the tables route costs
1.5x (gradient), 2.8x (Hessian), 4.0x (order 3), 6.6x (order 4); the
stencil route 4.9x to 8.1x. The shipped weibull3 against its hand-written
parent weibull1 shows the same shape, 1.3x to 4.7x. All routes agree to
1e-16 (tables) and 1e-11 (stencils). The absolute costs stay in
milliseconds, so the five shipped reparametrized families (weibull3,
student_t2, lognormal2, gengamma2, invgauss2) are fine as they are; NEW
families are written in full (laplace2 was, same day), and reparametrize()
remains the user-facing escape hatch rather than the package's own
construction route.

### The Rcpp review of 2026-08-06 (#83), the numbers that decided each package

Giovanni asked for every package to be benchmarked R against Rcpp, porting
wherever Rcpp measures ahead "anche se basso". The evidence, package by
package (min-of-5 timings, R 4.6.0/gcc 14.2 on this machine):

- **numericals7, log-Bessel: PORTED (0.7.0).** Mixed workload of 1e6 points
  across every branch: log K R 1.00 s vs Rcpp 0.34 s (2.9x), log I 3.35 vs
  3.09 (1.1x -- the series branch is lgamma-bound). The R implementations
  stay as internal twins (`.log_bessel_i_r/.log_bessel_k_r`) compared
  against the compiled route in a test; the u_k table is injected at load.
  The rest of the package: mills_ratio/bessel_i_ratio are one-liners over
  C-backed dnorm/pnorm/besselI, owen_t delegates to quad_vec, and the
  batched engines are matrix-op-bound -- nothing to port.
- **linkfunctions7: PORTED (2026-08-06, sera).** Thirteen scalar kernels
  (src/link_kernels.cpp) serve the derivative methods of logit, probit,
  cloglog, loglog, cauchit, rhobit, softplus; the doubly bounded link and
  the softplus inverse REUSE the logit kernel (scaled by the width,
  shifted one order). Measured through S7 dispatch: d4linkinv(logit)
  40 -> 25 ms at n = 1e6. identity/log/sqrt/inverse/power stay in R,
  their bodies being single C-level primitives. The old "transcendentals
  cost the same in C++" argument is half-true: the transcendental is
  shared but R's vector TEMPORARIES are not, and an order-4 polynomial in
  m*(1-m) allocates a dozen of them. Validator: check_link's numDeriv
  battery plus the extremes suite, no twins needed.
- **parameters7: AR ported in #82; log_cholesky assembly PORTED (0.7.0,
  2026-08-06 sera).** Every derivative of the Cholesky factor is a
  SINGLE-ENTRY matrix, so each Leibniz term is one row, one column or one
  cell -- chol_leibniz_cpp assembles all four orders without a matrix
  product, param_d1/d2 included. Measured at p = 8: order 4 from 9.39 s
  to 0.28 s (30x; the remainder is allocating the 82,251 result matrices
  the contract returns), order 2 from 6 ms to 1.2 ms. The dense R
  assembly stays as the twin `.chol_leibniz_r`, compared at machine
  precision in the tests. leibniz_gram survives for correlation_matrix,
  whose factor derivatives are NOT single-entry.
- **basis7: NOTHING TO PORT.** bspline_design delegates to splines2's
  compiled C++; Fourier/Legendre are vectorized recurrences over
  trig/polynomial ops.
- **distributions7 / optimizers7: already compiled where it counts** (the
  *_hd.cpp kernels; descent.cpp), and the 2026-07-31 measurement stands
  for the assembly layer: crossprod IS BLAS in both languages.

Method note: at n = 1e4 a single call is under the Windows timer's
resolution, so every number above is min-of-5 over repetition loops sized
to ~2e6 elements; the earlier NaN-ridden table from raw system.time was
discarded. Same lesson as the benchmark item in Testing and measurement.

### The de-jettization (2026-08-06), and what replaced each jet

The jets were built for `autoregressive(p, order = q)` (2026-08-04): its
chart HAS to be the partial autocorrelations -- the stationary region in
the coefficients is not a box, so no collection of scalar links covers it,
and Levinson-Durbin carries the PACF onto the coefficients bijectively
(Barndorff-Nielsen & Schou). The recursion's fourth derivative expanded by
hand is pages of algebra, so the first implementation carried jets through
it. Giovanni then measured what generic jet composition costs at the R
level -- the PIG kernels put it at 2x-36x the hand-written closed forms --
and banned the machinery outright (see section 5). What replaced it:

- **reparametrize() takes `map_derivs`**: keyed tables of the map's
  non-zero partials ("1", "1,2", "2,2,3,3"; a missing key is an exact
  zero), hand-derived for all six shipped second parametrizations in
  `R/reparam_maps.R` (fdb1/fdb2 are the written-out univariate and
  bivariate Faa di Bruno templates). Without tables, ONE stencil per
  partial on the analytic map (`reparam_stencil_derivs`) -- usable, and the
  documented reason to write tables for a family fitted in earnest.
- **skewnormal2's CP-to-DP map reduces to r = cbrt(2 gamma1/(4-pi))**: xi =
  mu - sigma*r, omega = sigma*sqrt(1+r^2), alpha = r/sqrt(b^2+(b^2-1)r^2),
  so `md_skewnormal2` is three univariate chains (a' = b^2 D^{-3/2} and
  friends). The power rule written through r/gamma keeps every derivative
  real on both signs of the skewness.
- **autoregressive propagates plain derivative arrays through
  Levinson-Durbin in Rcpp** (`src/ar_taylor.cpp`, parameters7's first
  compiled code): each tracked quantity holds its value and full symmetric
  tensors to order 4 in the q+1 free values, combined by the product rule
  written out per order (16 terms at order 4). Seeds are the link inverses
  with four derivatives, on one variable's diagonal each. The q = 1 case is
  pinned against products of link derivatives, formulas sharing no code
  with the kernel; `ar_prediction` and `param_logdet` need no derivatives
  and now call `linkinv` directly.
- **The jet twins survive in tests only** (`pig*_hd_jet_cpp`, self-contained
  C++): a mechanical transcription is an independent reference for the
  hand-written kernels, which is the one role the ban leaves it.

The general lesson survives the machinery: when a family's map is a
recursion of sums and products, propagate derivative arrays through it with
explicit per-order product rules; when the map is a formula, write the
derivative out. What is banned is the generic composition OBJECT with its
per-element overhead, not the mathematics.

### Renaming a PACKAGE reaches further still (covstructs7 -> parameters7)

Done 2026-08-04 (Giovanni: the old name was too narrow once simplex-valued
and transition-matrix parameters arrived, and "structure"/"struct_" named
nothing). Decisions and mechanics worth keeping:

- **`gh repo rename` keeps git and web redirects but the PAGES SITE DOES NOT
  redirect**: statmodels7.github.io/covstructs7/ died the moment the repo
  became parameters7. Every link in the portal, READMEs and book had to be
  swept by hand; `install_github("statmodels7/covstructs7")` still resolves.
- **A clean cut beats aliases** for a package not on CRAN with one consumer:
  struct_* -> param_*, covstruct -> parameter, and distributions7 swept in the
  same afternoon -- pushing distributions7 BEFORE its sweep would have broken
  its CI, because the renamed repo installs a package named parameters7.
- **The base/branch split is what made the rename worth it**: `parameter`
  keeps only the map; `matrix_parameter` owns dimension/rank/null
  basis/logdet/solve/factor, so a simplex has no logdet BY CONSTRUCTION
  rather than by run-time refusal. Registrations that read `@rank` moved down
  with it; ad-hoc subclasses in tests had to re-parent.
- **The Daleckii-Krein/Opitz machinery earns a note**: Frechet derivatives of
  expm contract direction chains against divided differences of exp; computed
  by the quotient recursion those cancel catastrophically under near-repeated
  eigenvalues, while the Opitz route (expm of a small upper bidiagonal, read
  off the corner) is exact there -- validated at eigenvalue gaps of 1e-9.
  The multilinear form sums over ALL orderings of the directions: summing the
  DISTINCT orderings requires the prod(mult!) factor, and the first draft of
  the multivariate gaussian's P-tensor recursion double-counted mixed terms
  the same way (error ~1e2, caught only by the stencil validation). The
  correct expansion of the derivative of an inverse is the sum over ORDERED
  BLOCK PARTITIONS, `P_t = sum (-1)^q P A_B1 P ... A_Bq P`, implemented
  recursively rather than transcribed.
- **The multivariate gaussian's d3/d4 are now closed** on top of
  param_d3/param_d4 (0.14 s against the fallback's 3.9 s at p = 4, order 3),
  both parametrizations validated at 4e-10/5e-7. The multivariate t stays on
  the disciplined fallback: nu blocks closure exactly as in the univariate
  skew t.

### Seven families in one sitting, and what each one taught

2026-08-05, Giovanni's list. Beyond the formulas, these are the items worth
keeping.

- **A simple family is a test of the machinery, not a duplicate of it.** That
  was Giovanni's argument for the exponential and the geometric and it was
  right: agreeing to 2e-15 with `fixed()` over a bigger family validates
  `fixed()` against something derived rather than delegated. Reach for that
  pairing whenever a new family happens to be a special case of an old one.
- **`ifelse()` returns a result the length of its TEST.** A scalar shape in
  the generalized Pareto's quantile collapsed a vector to one number, and the
  symptom was a fit converging to (0.5, -0.5) instead of (1.5, 0.3) --- three
  steps away, the random draws having been wrong. Anywhere a distribution
  branches on a PARAMETER rather than on the data, `ifelse` is a defect;
  recycle the test to the answer's length or branch explicitly.
- **A Monte Carlo check can be the weaker reference.** The generalized
  Pareto's closed-form information disagreed with a 4e5-draw Monte Carlo by 9%
  at `xi = -0.3`, and the formula was right: the second derivative blows up at
  the upper endpoint, so the mean converges slowly. Integrating on the
  PROBABILITY scale instead --- `E[h] = int_0^1 h(Q(u)) du` --- turns the
  endpoint into an ordinary point and gave 1e-11. When a closed form and a
  Monte Carlo disagree, ask which one has the heavier tail before believing
  the Monte Carlo.
- **The condition for an expectation to exist is the condition for an integral
  to converge, and it is worth deriving rather than quoting.** The generalized
  Pareto's information exists exactly for `xi > -1/2` because on the
  probability scale the integrand grows like `(1-u)^(-2|xi|)`. `NA` is
  returned below that, which is the honest answer; the quadrature really does
  diverge there.
- **A bounded chart for a circular parameter is a choice with a cost.** The
  von Mises carries its mean direction on `bounded_link(-pi, pi)`, which keeps
  it identified but stops a fit walking across the boundary. Leaving it
  unbounded makes the likelihood periodic and every maximum one of infinitely
  many --- the same shape as the folded normal's sign.
- **`expon.scaled` is what makes a Bessel constant usable.** `besselI(900, 0)`
  is `Inf`; the scaled form with the exponent added back is finite at 5000,
  and the ratio `I_1/I_0` needs no exponent at all. Same lesson as `log1p` and
  `expm1` above: the difference between a number and an infinity.
- **A distribution function saturates well inside its support.** The
  generalized gamma at a small shape has `F(y) = 1` exactly in double
  precision at `y = 4`, so no quantile comes back from there. A round-trip
  test must be asked at points obtained FROM the quantile function, not at
  arbitrary ones.

### An argument's POSITION is part of its interface

Giovanni reported *"un problema con il caricamento dei metodi S7"* from

```r
fit <- fit_distrib(d, y, lbfgs())
#> Error in as.vector(x, "list") :
#>   cannot coerce type 'object' to vector of type 'list'
```

There was nothing wrong with S7 or with the meta-package. `fit_distrib()` is
`(distrib, y, start, method, level, n_start)`, so an optimizer written
positionally lands in **`start`**, and `align_theta()` refuses to coerce it
several frames down. `method = lbfgs()` fits in 22 iterations.

Two things to keep. The message named neither the argument nor the mistake,
which is the failure mode section 7 already records twice in other words: an
error that surfaces far from its cause is worse than no error. `start` is now
checked in the function itself, and an `optimizers7::optimizer` or a
`FisherScoring` is told apart from any other wrong type so that the message
can say which argument was meant.

And the general form: **the order of arguments is part of an interface, and
the commonly-passed one being fourth is a trap the documentation cannot
close.** Reordering was rejected -- it would silently change the meaning of
every existing positional call -- so the guard is the answer. Anywhere two
adjacent arguments take unrelated types and the later is the one people
actually set, check the earlier one by type and say so.

### The meta-package, and the note it keeps

`statmodels7` was created on 2026-08-05 (Giovanni: a tidyverse-style
grouping, so that one command installs and attaches everything). Three
decisions worth keeping.

**Members in `Imports`, not `Depends`.** `Depends` attaches them through R's
own mechanism, in the order the field happens to list them, with no message;
`Imports` pulls them at installation — which is what makes one install
command install the toolkit — and leaves the attaching to `.onAttach`, where
it can be reported. The price is a NOTE, *"All declared Imports should be
used"*, correct as far as the heuristic can see, since nothing in a
meta-package calls a member by name. **It is left standing deliberately.**
The two ways to silence it are both worse: importing an arbitrary symbol
from each member puts five unused bindings in the namespace and says
something false about what the package uses, and moving to `Depends` gives
up the control the package exists to exercise. CI fails only on warnings.

**The member list is read, not written.** `statmodels7_packages()` parses the
package's own `Imports` and keeps the names ending in `7`, so a member added
there is a member everywhere. The test re-parses the field a second way, so
a mistake in the parsing cannot confirm itself.

**Measured:** the five members export **no name in common**, so
`statmodels7_conflicts()` is empty and the attach order does not matter.
(Outside the toolkit there is one: `testthat` masks
`distributions7::expectation`. The report's scope is the members, correctly.)

Two things the check caught, both worth generalizing:

- **A shared test file carries a hidden dependency.** `test-docs.R` calls
  `rmarkdown::pandoc_available()`, and the five older packages satisfy that
  only **incidentally**, having `rmarkdown` in `Suggests` for their
  vignettes. A package without vignettes gets *"'::' import not declared"*,
  a WARNING, which CI treats as a failure. This is the fourth time the
  `::`-goes-in-`Suggests` rule has been walked into; the new thing is that a
  file copied between packages brings its dependencies with it, and they may
  be satisfied by accident at the destination.
- **`Title` must be title case, and a lowercase package name in it is not.**
  *"The statmodels7 Toolkit"* is flagged; the siblings' style — a
  description with no package name — avoids it.

### The portal is hand-written HTML, and nothing executes it

The rename of 2026-08-04 swept the packages, the book and the READMEs. It did
not reach `site/index.html`, which still showed `struct_matrix`,
`struct_dmatrix`, `struct_logdet`, `struct_dlogdet` and `scaled_struct` — five
names that had not existed for a day — plus a `statmodels7` card marked
*planned* with `pointer-events: none`. Nothing failed, because nothing runs
the portal.

This is the same shape as the non-executed README block below and as the
`\value`-less pages above: **a surface no process exercises rots silently**,
and the only guard is to include it in the sweep by hand. After any rename,
grep `site/index.html` as well as the packages.

### A name names the coordinate, not the quantity

Giovanni, 2026-08-04, reading a fit: *"dal fatto che sigma_scale sigma_pacf1 e
sigma_pacf2 abbiano link identità non si capisce qualcosa: i nomi
sembrerebbero quelli di parametri vincolati ma in realtà il legame è
identità!"* He was right, and the printed evidence is the whole argument:

```
sigma_pacf1   0.9679   [0.9125, 1.0233]
```

A partial autocorrelation with an interval **past one**. The number is
`atanh(0.7478)`, and 0.7478 is the `cor_v1_v2` printed four lines above.

The cause is structural rather than cosmetic. A `parameters7` family's free
vector is flattened by `distributions7` into scalars carrying **identity**
links, because the free vector is already unconstrained — that design is
right and stays. What was wrong is that four families named the coordinate
after the **constrained** quantity it produces, so the name promised a
bounded thing and the value reported was on the free scale. The other six
families were already honest, in two different ways: `log_L1` names the
transform, `z2.1`/`S2.1`/`alr1` name the chart. So the package held one
defect and two good conventions at once.

**The rule now**: where a link carries a quantity onto the free scale the
name records that link (`log_scale`, `z_rho`, `logit_rho`, `z_pacf1`,
`log_d1`, `sqrt_d1` under a square-root link); where the coordinate is
already unrestricted the name is the plain quantity (`L2.1`, `S2.1`).
`link_tag()` in `R/naming.R` maps a link CLASS — not its `link_name`, since
a parametric link names itself `"bounded(lwr=-0.25, upr=1)"`, which cannot
appear in an identifier. The property the spelling reflects is worth
asserting separately from the spelling: a labeled coordinate ranges over
the whole line, tested by sweeping each free value to ±20 and checking the
matrix stays finite and non-negative-definite. At ±20 a correlation reaches
its boundary in double precision (min eigenvalue 1e-33 of the max), so
strict definiteness is asserted only at ±6 — an open set approached far
enough is a closed one.

**And what a coordinate is not is what a reader reads.** `param_readable()`
lets a family declare its quantities with the Jacobian from the free vector
and the scale each interval belongs on; `mv_summary()` prints them as a
block. For an AR(2) that is the marginal variance, the partial
autocorrelations and the **coefficients**, which appear nowhere in Σ at all
and were previously reachable only through the recursion's internals. Their
Jacobian is free: the Levinson-Durbin recursion already propagates its
derivative arrays, so the first-order block is exactly what the delta
method needs. The
coefficients are intervalled on the identity scale deliberately — the
stationary region is not a box, so no scalar transform expresses it.

Two checks in the output are checks of the code: `pacf1` must equal
`cor_v1_v2` to the digit (the lag-one correlation of an AR **is** its first
partial autocorrelation) and `phi_q` must equal `pacf_q`. Both hold, and
both are asserted.

### An argument that is accepted and ignored is worse than one that errors

`fit_distrib()` carried `maxit` and `tol` while every optimizer carries its
own, and with `method = <an optimizer>` the two were **silently discarded**.
Giovanni's call was `fit_distrib(..., method = bfgs(maxit = 1000), maxit = 100)`:
1000 governed, the 100 was read by nobody, and nothing said so. Both are
gone from the signature — Fisher scoring is Newton with one matrix replaced,
so its budget and its stopping rule live on `fisher_scoring()` like any
other optimizer's, and the internal BFGS fallback inherits them rather than
inventing its own.

The consequence worth knowing for next time: **three places restated the
default**, and my first grep found only two of them, because it searched for
`formals(distributions7::fit_distrib)` while the third said
`formals(fit_distrib)`. Grep for the SHAPE of the reference, not for one
spelling of it. The rule from the tolerance episode still stands and is now
uniform: a threshold that mirrors a default reads that default
(`eval(formals(optimizers7::crit_grad)$tol)`), never a copy of the number.

### The same guard in three strengths is not one guard

`test-docs.R` was in all five packages, and in **three different versions**
(2026-08-04). distributions7 and optimizers7 asked about every object in the
namespace, reading `man/` from disk. linkfunctions7 asked only about
`getNamespaceExports()`, so an undocumented internal was invisible to it.
basis7 and parameters7 asked about every object but through
`tools::Rd_db()`, which under pkgload skips with *"package not installed"* —
that is, always, since `test_local()` is how the suite runs.

CLAUDE.md §6 has claimed since 2026-07-30 that everything is documented,
exported or not. The census found **seven objects with no page**:
`eta_bounds` and `is_base_link_class` in linkfunctions7,
`diag_dlog`/`diag_higher`/`diag_logdet_higher`/`mlog_higher`/`scaled_dlog` in
parameters7. All seven sat in the packages whose guard could not see them.
The file is now byte-identical across the five apart from the package name,
generated from one source, and its header says so.

One correction went the other way. The strong version demanded `\value` on
**every** page, including a package landing page, which is not a function
and by convention has none; basis7 and parameters7 generate one through
`"_PACKAGE"` and the other three do not, so the guard was over-strict rather
than those two being wrong. The exclusion is by what the page is, not by
whether it happens to exist.

### A fenced block that is not executed is documentation that cannot fail

parameters7's README had a ` ``` r ` block — plain, not `{r}` — showing
`mvgaussian_distrib(2, struct_sigma = ...)` and a hardcoded output of
`log_L1`, `log_L2`, `L2.1`. The argument was renamed to `sigma` and the
names gained a `sigma_` prefix, and neither reached it: the block is never
run, so nothing could go red. It had been wrong through two renames. Where
a package genuinely cannot execute an example — parameters7 must not depend
on its own consumer — the block still has to be swept by hand at every
rename, and the same applies to a summary table: the README's "what is in
the box" still listed four families out of eleven and derivatives to second
order.

The README also had `ar1 <- Ar1(...)` as a user-defined example, written
before `ar1()` existed as a real constructor. `devtools::build_readme()`
caught it — the next chunk passed the FUNCTION to `check_parameter()` — but
only because a later chunk happened to use the variable. **Grep a package's
own exports against the variable names in its documentation** after adding
any constructor.

### Renaming a parameter reaches further than the package

Prefixing the multivariate free names with `sigma_`/`omega_` (2026-08-04)
touched 97 literals across the package, its tests, the README and the book. Two
of the places a mechanical sweep does NOT reach are worth knowing in advance:

- **Hessian component names are CONCATENATIONS of parameter names**, so
  `eh[["mu1_log_L1"]]` had to become `eh[["mu1_sigma_log_L1"]]`. A regex with a
  word boundary correctly refuses to touch `log_L1` inside `mu1_log_L1` --
  which is what makes it safe, and also what makes it miss these.
- **The book's gate built its lookup from `param@free_names`**, the
  STRUCTURE's names, not the distribution's. After the rename every component
  it looked up was NULL, and the failure surfaced as
  *"non-numeric argument to mathematical function"* several frames away from
  the cause. Anywhere a consumer reconstructs a parameter name from the
  structure rather than reading `distrib@params`, a rename will break it
  silently or confusingly.

### An argument named after a class shadows it

`has_exact_expected_hessian()` took an argument called `distrib`, and the base
class of the package is also called `distrib`. Inside the function the
argument won, so the line meant to ask *is this method registered on the base
class?* compared the owning class with the DISTRIBUTION OBJECT instead. It was
therefore false for every family whose expected information comes from the
base class -- which is precisely the set the predicate exists to identify --
and the function answered "closed form" for all of them (2026-08-04).

The visible consequence was the refusal firing backwards: `fisher_scoring()`
with an `approx` was refused on the skew normal, which needs it, and accepted
on families that would ignore it. The three sibling checks in the same
expression, against `continuous_distrib`, `discrete_distrib` and
`multivariate_distrib`, were fine, because nothing shadows those -- so the
multivariate t behaved correctly while the skew families did not, which is
what makes this kind of bug slow to see.

**The argument is now `x`.** General shape: in a package whose classes and
whose arguments are drawn from the same vocabulary, an argument named after a
class is a silent rebinding. Grep for a formal whose name matches an exported
class before trusting a comparison against that class.

### A starting value is not a detail

Found on `mvgaussian_distrib(4)` fitted to `iris[, 1:4]` (2026-08-04, Giovanni
reported it as "the fitting takes really long -- maybe because we have not put
in Rcpp what should be in Rcpp?"). It was not Rcpp. Measured, from the same
data, the same derivatives, the same optimizer:

| start | iterations | -logLik | time |
|---|---|---|---|
| zeros on the link scale | 500, DID NOT CONVERGE | 836.17 | 2.31 s |
| sample mean + covariance | **1**, converged | **379.914630** | 0.00 s |

The gradient at the closed-form maximum is 2.5e-12, so nothing was wrong with
the mathematics. `fit_distrib()` was starting at `eta = 0`, which for a
multivariate gaussian is a zero mean and a unit covariance, and on data whose
scale is nothing like that the run never arrives.

⚠️ **Half of that measurement was a second defect, and the table above is the
"before" of two fixes at once.** With the `sqrt(p)` information bug of §4
repaired, the same fit started at the origin **does** arrive -- 42 iterations,
0.28 s, the same maximum to 6e-14. The corrupted expected information was
making Fisher scoring take poor steps from far away, so the origin looked
fatal when it was merely expensive. What the starting value buys is therefore
**work, not an answer**, and the book's §3.3 and its gate were both weakened to
say that. General shape: when one change fixes a symptom, the number that
motivated it has to be **re-measured** afterwards, because a second cause can
be hiding inside the first one's evidence.

- **`distrib_start(distrib, y, n_start)`** is the fix: a generic returning a
  starting value computed from the DATA. The default is the old random draw
  from the parameter domains, so a family that says nothing loses nothing; the
  multivariate gaussian returns the sample mean and covariance -- its own MLE
  for an unstructured matrix -- and the t returns those with `nu = 8`. For a
  precision structure the covariance is inverted first, and for a structure
  that cannot represent the sample covariance exactly, `param_free()` refuses
  and a short least-squares fit stands in. That last fallback is allowed to be
  approximate BECAUSE it is a starting value; `param_free()` itself must stay
  exact-or-refuse, which is the contract parameters7 documents.
- **The restart loop kept the LAST start's result, not the best.** With five
  random starts and none converging, the fit reported whichever start came
  last -- a mean of 2.03 where the data's mean is 5.84. The rule is now:
  a converged run beats a non-converged one, and among runs of equal status
  the lower objective wins.
- **AIC and BIC differed between the covariance and the precision
  parametrization** of the same model (1594 against 812.9). That was the same
  bug seen from the other end: neither run had converged, so the two were
  comparing accidents. They now agree to the digit, which is the invariance a
  reparametrization has to satisfy and a useful thing to assert.
- **optimizers7's gradient check fires at a stationary start.** The difference
  of `fn` along a direction of no slope is its own truncation error, so
  comparing it with a gradient of the same order compares two kinds of
  nothing -- and it warned at precisely the caller who had done the best
  possible thing, handing in the exact MLE. It now skips when the gradient
  norm is negligible against the size of the objective.

Read together with the item above on finite differences inside a score: both
are cases where the arithmetic was right and the run still reported failure.

### A gradient tolerance is bounded by the value of the objective

The last red macOS job (2026-08-04) was three fits reporting
`converged = FALSE`, and the number that explained it came from CI itself:
the assertions were changed to print what the run did, and macOS answered
*"BFGS, 14 iterations, 'the line search found no acceptable step',
score/n 2.086e-10"* against a tolerance of 1e-10. **An assertion that prints
only TRUE/FALSE cannot be diagnosed on a platform you do not have**;
`fit_report()` in `tests/testthat/helper-distributions.R` is the fix and is
worth copying anywhere a convergence flag is asserted.

The mechanism is arithmetic, not platform. A line search accepts a step only
when the objective decreases by a definite amount, and near the optimum that
decrease is about `|g|^2/(2*lambda)`. Once it drops below the rounding of the
objective itself, `eps*|f|`, no step in any direction can be verified and the
search refuses all of them, so the smallest gradient a run can reach is
around

```
|g|_floor ~ sqrt(2 * lambda * eps * |f*|)
```

and it **grows with the value at the solution**. The decisive measurement is
one line: adding a constant to the objective moves neither the minimizer nor
the gradient, and on conjugate gradients applied to Rosenbrock it takes the
attainable gradient from 1.9e-9 at `f* = 0` to 4.4e-8 at `f* = 1`, 2.8e-6 at
1e3 and 6.5e-5 at 1e6 -- a square-root law, visible with no reference to any
platform. BFGS is unaffected on those problems because it arrives before it
hits the wall, which is exactly why the defect looked like a BFGS/Fisher
disagreement.

⚠️ **Every problem in `optimizers7::test_problems()` has `f* = 0`**, where the
rounding floor is itself zero, so the package's own battery could never have
shown this. A likelihood is the opposite case: `-l/n` is of order one at the
maximum, and the floor is then ~1e-8. Measured across families, methods and
seeds (72 runs) it is usually 1e-15 or 1e-16 but reaches **1.06e-8**, and
macOS reported 1.345e-8 on the multivariate t -- the same number.

Both defaults were therefore raised to **1e-6**, two orders above the worst
observed floor: `crit_grad()` in optimizers7 and `tol` in `fit_distrib()`.
Nothing statistical is lost, a score of 1e-6 per observation placing the
estimate a small fraction of a standard error from the maximum. Two
consequences worth knowing:

- **the optimizers no longer restate the constants.** `gd`, `cg`, `bb`,
  `newton`, `bfgs` and `lbfgs` all had `crit_any(crit_grad(1e-8),
  crit_rel_obj(1e-12))` written out, so changing `crit_grad()`'s default would
  not have reached any of them. They now say `crit_any(crit_grad(),
  crit_rel_obj())`;
- **a test's tolerance on an estimate cannot be tighter than the rule that
  stopped the run.** Several asserted a parameter to 1e-6 or 1e-7 while the
  fit promised only what its stopping rule promised, and passed by the same
  luck that made the CI green four times out of five. They now follow the
  tolerance, and say so. The same applies to the book's gates: one of them
  restated `1e-10` and went red the moment the default moved, so it now reads
  `formals(distributions7::fit_distrib)$tol` — **a threshold that mirrors a
  default belongs to that default, not to a copy of it.**

**What `fit_distrib()` kept from the optimizer, and what it threw away**
(2026-08-04, Giovanni asked). It has delegated to `optimizers7::minimize()`
since 2026-08-03 and already carried the iteration count, the evaluation
counts, the rule that fired and the note. It discarded `elapsed` and
`gradient`, both of which are worth having: the time is now **accumulated over
every start and every fallback**, since what a caller wants is what the fit
cost rather than what the surviving run cost, and the gradient of `-l/n` is
*by construction* the score per observation the rule tested, so `@score` says
how close to stationary a run ended. `print()` shows the time always and the
score only when the run did not converge, which is when it answers something.
Consequence to know: the book prints a fit in two chunks, so its HTML now
carries a timing that changes at every render — the same harmless churn as
`date: today`, and not evidence that anything changed.

⚠️ **The first version printed the time only when it was positive, and
Windows went red for it.** On a coarse clock a fast fit measures *exactly*
zero seconds, so what the object reported depended on the platform's timer
resolution — the same class of defect as the tolerance this whole episode was
about, reintroduced in the fix for it. Zero is the reading, not the absence of
one: the line is printed whenever the time is finite. **A guard of the form
`if (x > 0)` on a measured quantity is a claim that zero cannot be measured**,
and for a duration that claim is false.

### Completing the higher derivatives (2026-08-04, sera)

What closing the d3/d4 gap on the five families taught, beyond the formulas:

- **The generic d3/d4 fallback differentiates `distrib_hessian`, so on the
  skew t it nested differences in nu** -- the components whose Hessian entry
  is itself a stencil (`(i,nu,nu)`, `(nu,nu,nu)` and their order-4 analogues)
  were a difference of a difference, the shape the toolkit forbids. The
  skew t now registers its own d3/d4: the generic construction is kept for
  every component whose Hessian entry is closed (there one stencil lands on
  an analytic quantity and is legal), and the nested ones are replaced by
  single higher-order stencils on the closed score or the log-density
  (`fd5_third`, `fd5_fourth` beside `fd5_first/fd5_second`).
- **The fourth difference gets its own step.** Rounding grows as `h^-4`: at
  the family's base step (1e-3 relative) the per-observation noise on
  `nu_nu_nu_nu` is ~1e-2 relative, at ten times that it is negligible and the
  `h^2` truncation (~6e-4) is what remains. Swept over nu in {3, 6, 15};
  same lesson as the score's step -- measured, not chosen.
- **`approx = "integrate"` fails on the Gumbel's expected d3 with
  "non-finite function value"**: in the left tail the observed derivative
  carries `w = e^{-z}` which overflows while the density underflows, and the
  quadrature meets Inf*0. The right independent reference there is the Monte
  Carlo mean of the observed kernel (all 18 components validated at |z| < 1.4
  on 1e6 draws). A quadrature failing is not evidence against a formula.
- **Tiny components sitting beside huge ones cannot be checked relatively.**
  At t = -40 the skew normal's `mu_mu_mu` is -0.016 next to components of
  1200, and kernel, fallback and Richardson disagree at the percent level --
  all three references fight over crumbs. The asymptotic form (here
  `psi''' ~ -2/t^3`) is what settles it, and it sided with the kernel.
- **Multivariate d3/d4 stay on the fallback, and analytic ones are blocked on
  parameters7's API**: the components with three or four structure indices
  need `struct_d3matrix`/`struct_d4matrix`, which the parameters7 contract
  (exact to 2nd order) does not offer. Cost of the fallback, measured at
  n = 50: order 3 takes 0.14 s at p = 2, 0.55 at p = 3, 3.9 at p = 4 (560
  components, two Hessians each); order 4 takes 2.7 s already at p = 3.
  Extending parameters7 is a design decision for Giovanni, not a patch.
- The multivariate base class now **refuses `grad_y`/`hess_y`/`cross_y`**
  (the univariate fallbacks difference along a line and would return numbers
  of the wrong shape); the gaussian and the t override with closed forms. The
  t's `hess_y` is the reweighted gaussian expression
  `-c Sigma^-1 + (2c/(nu+q)) w w'`, one `p x p` matrix per row (a
  `p x p x n` array), where the gaussian returns a single `p x p` because its
  is constant -- the two return shapes are documented on each page.
- **A method file registering on a generic defined in another file needs that
  file in its `@include`**: `multivariate.R` gained methods on
  `distrib_cross_y`, whose generic lives in `cross_derivatives.R`, and
  pkgload died with "object not found" until the Collate order said so.

### Finite differences inside a score

Three things learned adding the skew t (2026-08-04), all of which apply to any
family whose score is not entirely closed form.

- **A stopping rule on the gradient cannot be met below the error of the
  finite difference inside the gradient.** The skew t's `nu` score is one
  stencil on the log-density; measured, its error on the SUMMED score is
  1e-11 to 1e-9 depending on n and nu. `fit_distrib()` defaults to
  `crit_grad(1e-10)`, which is at or under that floor, so the run reaches the
  maximum -- the three closed-form components go to 1e-13 -- and then spends
  its whole iteration budget reporting `converged = FALSE` at the maximum.
  `tol = 1e-8` is the honest ask for such a family, and the constructor's page
  says so. Note this is a NEW consequence of dropping `crit_rel_obj` from the
  criterion the same day: the old OR would have stopped on the objective and
  reported success.
- **The step is measured, not chosen.** Swept over nu from 2 to 30 and n from
  500 to 4000: a three-point stencil bottoms out at ~1e-8 (truncation, which
  is a BIAS and so does not cancel over the sum), while a five-point one at a
  relative step of 1e-3 reaches 1e-11 to 1e-9. Going smaller makes it worse --
  the five-point stencil amplifies rounding by 18/(12h). Two extra density
  evaluations bought an order of magnitude, which was the difference between a
  fit that converges and one that does not.
- **Never let the numerical Hessian differentiate a numerical gradient.** The
  base fallback differentiates `distrib_gradient` once; if that gradient is
  itself a difference in the same variable, the result is the nested
  differencing the toolkit forbids everywhere. So the skew t registers its own
  Hessian: closed form where it can, one stencil on the log-density for
  `nu_nu`, and one stencil on the CLOSED-FORM score for the mixed `(par, nu)`
  entries. A mixed stencil differences two DIFFERENT variables and is therefore
  a single stencil, not a difference of a difference.
- Corollary for testing: a gradient that contains a difference must be checked
  against **Richardson extrapolation** (`numDeriv::grad`), never against a plain
  central difference, which would be the same arithmetic twice.

### A run that reached a point is not a run that failed

`fit_distrib()`'s restart loop discarded a BFGS fallback that ran but did not
converge, so `res` stayed NULL and the caller got
*"Optimization failed from every starting value; supply 'start'"* for a fit
that existed and was correct (found on the skew t at the default tolerance,
2026-08-04). The fallback is now kept as a last resort when the chosen method
raised, and the loop keeps the last result it obtained. General shape: an
error message that names a cause is worse than useless when the cause is not
the real one -- the same lesson as the `check_criterion()` item in section 4.

**The same shape, one layer down, in `optimizers7`** (2026-08-04). A descent
run asked its stopping rule only *after* an accepted step, so a line search
that could find no acceptable step broke the loop with `converged = FALSE` and
*"the line search found no acceptable step"* -- **including when the start was
already the answer**, which is now the ordinary case, since `distrib_start()`
hands a multivariate gaussian its closed-form estimate. Fisher and Newton
reported one iteration and convergence while BFGS reported failure on the same
fit, and that disagreement was the last red job on macOS. `src/descent.cpp`
now asks the criterion before giving up, with `have_old = false`, so `f_old`
and `x_old` arrive NULL: a rule reading a change in the objective returns
FALSE by construction and only the state at the point -- its gradient, its
stationarity -- can end the run. The counterexample is in the tests beside it:
a mis-stated gradient makes every direction ascend far from the optimum, and
that run must still fail. **A method that cannot move is converged only if the
point says so, never because the move failed.**

### The two non-elliptical families, and what they cost the validator

Adding the Dirichlet and the multinomial (2026-08-05) broke nothing in the
derivative machinery -- the flattening of §4 carried them unchanged -- and
broke three of the nine checks in `check_distrib()`'s multivariate battery,
on code that was correct. Each failure was the battery assuming something
only an elliptical family satisfies. What that cost, in order:

- **The normalization check assumed a proposal exists in `R^p`.** For the
  Dirichlet the support has measure zero there, so a gaussian proposal
  estimates 2.0e-08 for an integral that is 1 -- and does not raise, because
  `chol()` accepts a covariance whose smallest eigenvalue is 8.7e-19. That is
  the *`chol()` is not a rank test* note met for the **third** time, in a
  third package-independent context. The fix is a generic,
  `mv_reference_draw()`, whose base method is the old gaussian and whose
  Dirichlet method is the uniform on the simplex; a discrete family never
  reaches it, its normalization being an exact sum over `mv_support()`.
- **The response-derivative check assumed a refusal is a defect.** The
  multivariate base class refuses `distrib_grad_y` *by design*, so a family
  that has not registered one has made a choice; and a discrete family has no
  such derivative at all. Both are now not emitted rather than failed --
  the convention the univariate battery already follows with twelve checks on
  a discrete family against thirteen on a continuous one. A third status
  (`"n/a"`) was considered and rejected: several consumers, including three
  book gates, read `status != "OK"` as failure, and a new value would have
  made every one of them wrong at once.
- **`mean()` and `variance()` were simply missing**, which surfaced three
  frames away as *"Can't find method for `expectation(<DirichletDistrib>)`"*
  -- the base class falling through to the univariate quadrature.

Two smaller things worth keeping:

- **`@include` naming a file pulls that file EARLIER in Collate**, which can
  drag it ahead of a file it needs. Naming `mvgaussian_distrib.R` (to reach
  `mv_location`) pulled it ahead of `moments.R` and the package failed to
  load with *"object 'variance' not found"*. The right fix was not to
  hand-edit Collate but to move the three generics -- `mv_location`,
  `mv_sigma`, `mv_marginal` -- out of `mvgaussian_distrib.R` and `mv_plot.R`
  into `multivariate.R`, **where they belonged all along**: they are generics
  of the multivariate layer, not of the gaussian. A Collate cycle is usually
  telling you a definition is in the wrong file.
- **The roxygen-placement trap fired again**, and this time it had been live
  and committed. `mv_param_block` was inserted between the
  `mv_derived.multivariate_distrib` block and the method it documents, with
  no blank line, so the two comment runs fused: the man page kept the
  method's `@name` and got the *helper's* `\usage`, and `R CMD check`
  reported it only as *"Documented arguments not in \usage"*. It is invisible
  to the tests, since the topic exists and is indexed. **`--as-cran` is what
  found it; the docs guard cannot.**

### Multivariate distributions

- **A rank read off an assembled matrix is not a property of the family.** The
  measured case, on a tensor-product penalty `l1*P1 + l2*P2` of true rank 28 out
  of 32: counting eigenvalues above a relative tolerance of 1e-10 gives 28 while
  the two are comparable and **24** once the ratio reaches 1e10, while the
  residual `max|M N| / max|M|` against the stored null basis never leaves 3e-16
  across the whole range. (An earlier note here said 26 at 1e8 and 16 at 1e12;
  those came from a different tolerance. The numbers above are the ones the book
  now prints and its gate now checks, so they are the ones to quote.) Note also
  that the 1-D pair -- second and first differences on 8 coefficients -- does
  **not** reproduce the effect at any ratio up to 1e14. A demonstration of this
  needs the tensor case; a smaller example silently demonstrates nothing.
- **A `parameter` refuses `param_solve()` when it is rank deficient**, and that
  is deliberate rather than an omission: what a consumer of an improper prior
  needs is the quadratic form and the log pseudo-determinant, since the matrix it
  actually inverts is `X'X + lambda P`, which is non-singular. I wrote a book
  section claiming a Moore-Penrose solve before checking, and the code was right.
- **The matrix parameter is flattened into scalars.** Everything else in
  distributions7 then applies unchanged -- see section 4. The corollary worth
  keeping is negative: a multivariate distribution must NOT inherit from
  `continuous_distrib`, because the defaults registered there (cdf by quadrature,
  quantile by root finding) are one-dimensional and would silently produce
  numbers.
- **The gaussian's expected information needs no second derivative of the
  structure.** `E[l^(kl)] = -0.5*tr(Sigma^-1 A_k Sigma^-1 A_l)`: both
  second-derivative terms cancel under expectation, and the mixed mean-matrix
  block is exactly 0 because `E[w] = 0`. That is why Fisher scoring is cheaper
  than Newton here and not merely differently accurate.
- **A multivariate t is the gaussian score reweighted.** Every data term carries
  `c = (nu+p)/(nu+q)`; only the `nu` component is new. Marginals keep the SAME
  `nu` -- it follows from the scale-mixture representation, since conditioning on
  the mixing variable leaves a gaussian.
- **The scale matrix and the covariance are different objects and must stay
  apart.** `mv_sigma()` returns what the parametrization carries, `variance()`
  the moment. Conflating them makes the family undescribable at `nu <= 2`, which
  is exactly the regime it exists for.
- **A test harness that evaluates chunks in the global environment must not
  use short names of its own.** A scratch script looping `for (s in starts)`
  clobbered the structure each chapter chunk had just built, and reported the
  chapter as broken. The chapter was fine.

### Random number generation

- **GRoU (`rng_grou`) is the continuous fallback**, not inverse transform, which costs one
  `uniroot` per draw when the quantile is itself numerical. Two devices make it safe:
  the kernel is **recentered at the mode** (without which a density at mu = 1000 gives a
  degenerate bounding box) and **normalized to a maximum of 1**.
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
  call per draw; vectorizing the table lookup took it from 12.6 µs to 0.1.

### The batched engines (numericals7, 2026-08-05)

- **A quadrature budget allocated by panel length cannot meet an integrable
  endpoint singularity at any depth.** Near x^(s-1) the error is concentrated
  in the innermost panel however deep the bisection goes, so a per-length
  share demands of it an accuracy no depth can reach: the first quad_vec
  refused a beta(0.8, 1.2) normalization that stats::integrate handles
  routinely. Convergence is judged on the SUM of a row's panel errors, as
  QUADPACK judges it, and the worst panels are bisected.
- **A row whose max-depth panels alone exceed the budget must fail AT ONCE.**
  Without that clause the loop keeps splitting the smooth panels' rounding
  noise -- errors all of one tiny magnitude, so half the list clears the
  0.5*max threshold at every pass and the panel count doubles until memory
  runs out. The observed symptom was a hang, not an error.
- **A series tail guard needs three conditions, not two.** Block sum small
  and last term small are both true BEFORE the mode of a hump-shaped term:
  dpois(0:63, 300) sums to nearly nothing and ends on ~1e-52. The third
  condition -- the terms are not growing across the block -- is what tells a
  premature block from a finished tail.
- **The old scalar engine was the weaker twin.** On 300 gamma means the
  looped stats::integrate reference sat 2.5e-6 from the closed form (its
  default rel.tol) while quad_vec sits at 2.8e-12; when a twin comparison
  fails, ask WHICH side is wrong before touching the new one. Same lesson as
  the GPD Monte Carlo.
- **R's scaled besselI underflows to an EXACT ZERO between kappa = 1e5 and
  1e6**, so the ratio I1/I0 came back NaN there -- the old vonmises2 note
  saying "past about 1e15" was wrong by ten orders. bessel_i_ratio switches
  to the asymptotic expansion 1 - 1/(2k) - 1/(8k^2) - 1/(8k^3) at kappa =
  1e4, where the next term is below the resolution of a double, so the whole
  range is evaluable. Same shape as expon.scaled and log1p: the difference
  between a number and a NaN.
- **expectation()'s integrand contract is now elementwise in y AND theta
  jointly**: every parameter combination shares one batched quad_vec /
  series_vec call, so f receives vector theta components recycled against y.
  Every internal f (distrib_* calls, moment powers) already satisfied this;
  a user f that branches on a scalar theta does not. A combination whose
  quadrature fails raises an error naming it, preserving the old
  integrate-error semantics that expected_by_integrate's tryCatch relies on.
- **A text cut keyed on "the first brace after the definition" eats the next
  function when the body has no braces.** Removing the one-line vm_dA took
  the section comment and the roxygen header of the pdf method behind it,
  fused as in the roxygen-placement trap; the deleted Rd in git status was
  the only visible symptom, and the suite stayed green. After any scripted
  deletion, read the git diff around every cut, and treat a deleted man page
  as a mangled roxygen block until proved otherwise.

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
  (CI exercises it). General shape: any code that serializes S7 objects to
  another R process must ensure both sides run the same namespace.
- **A method registered on a BASE generic needs `S7::methods_register()` in
  `.onLoad`, and pkgload hides its absence** (2026-08-08). `S7::method(print,
  cls) <- fn` registers with the S3 dispatch table under pkgload, so
  `test_local()` is green; from an INSTALLED package the method is simply
  absent and the object prints as the raw S7 property dump. penalties7 shipped
  that way and nothing noticed until modelterms7's `R CMD check` (which
  installs) failed on its own print tests. Both packages now carry the
  `.onLoad`; any future package with methods on base generics needs it, and a
  print test is what makes the omission visible.
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
  The comparison is now over **six** packages, and the meta-package is the one
  most likely to drift: it has no `src/`, no vignettes and almost no `R/`, so
  a habit that lives in the others' boilerplate has nothing to attach to here.
  Its missing `rmarkdown` in `Suggests` was exactly that.
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
- **`gh run list --commit <sha>` needs the FULL forty-character sha** and
  matches nothing at all against an abbreviated one -- silently, with exit
  status 0 and an empty list (2026-08-04). A script that looped over six
  repositories and printed a per-repository verdict therefore reported
  *"TUTTO VERDE"* from **zero runs found**, which is the worst possible
  failure for a check whose entire job is to say whether something is green.
  Any script that summarizes a query must **assert that the query returned
  something** before summarizing it; an empty result is a third answer, not a
  pass. `git rev-parse HEAD` gives the full sha, and `--json` output is easier
  to count than to read.
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
  runs twelve checks on what an optimizer *reports* — that `value` is the objective at
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
  whether `statmodels7` would really be able to hand the optimizers
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
  model layer is *vectorized R over compiled kernels*, and the real order of
  magnitude is algorithmic — Fisher scoring exploiting block structure and sparse
  penalties — not linguistic.
- **A feature with one caller is that caller's feature.** `adam(resample=)` and
  `finite_sum()` let the optimizer draw its own minibatches, and between them
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
  works through a user-written optimizer, and writing it showed that a user *could not*
  write a conforming one: honoring `bounds` and refusing an unevaluable stopping rule
  both needed internal functions. The package promised extensibility and did not supply
  the parts. Three exports later it does. The same vignette also disproved two claims
  in its own prose — measure the example before describing it.

---

## 8. Open items

- The two logistic integrals above (cosmetic — `approx = "integrate"` already delivers
  them within Monte Carlo noise).
- Publishing to CRAN now starts one step earlier: **`numericals7` first** (dependency-free), then `linkfunctions7` **and `optimizers7`**, which together unblock
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
  boxed quadratic seen through its reparametrization, 1395 refusals in 1521
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
- ~~`modelterms7` owes `nl()`~~ **done 2026-08-08**. modelterms7 is feature complete against `piano_modelterms7.txt`; what it owes now is whatever the model layer asks of it. Two contract questions are already visible, both the model layer's (Giovanni, 2026-08-08). A structural term's parameters are not penalizable, so `gas(by =, penalty =)` waits on how a term's own parameters reach a penalty. And `gas()` drives ONE distribution parameter, where the general Creal-Koopman-Lucas form has a VECTOR level with matrix loadings so the scale can respond to the location's score: the arithmetic generalizes mechanically and `term_filter()`'s signature survives (only the jacobian's shape changes, n x npar to n x k x npar), but a term written inside ONE parameter's formula has no place to declare that it spans several -- that construct is statmod()'s. The chart would have to change too: the PACF construction is scalar, and a matrix B's stationary region is a bound on the companion matrix's spectral radius, not a box.
- ~~GAP 1 of `piano_modelterms7.txt`~~ **done 2026-08-08**; the next batch is statmodels7's own fitting layer, which owes a coordinate-descent path for the separable case (measured 1.1x-5.3x faster than prox_grad, and not expressible as an optimizer). kron_identity() closed the identical-blocks case of the composition item; D R D', general block diagonals of distinct structures and sums remain open.
- Next packages: **`parameters7` phase 2** (`piano_covstructs7.txt`:
  correlations via canonical partial correlations, D R D', compound
  symmetry/AR(1), block diagonal, composition -- phase 1 shipped 2026-08-03),
  then `penalties7`, then `modelterms7`.
- **parameters7 orders 3-4?** Analytic multivariate d3/d4 need
  `struct_d3matrix`/`struct_d4matrix`, which the parameters7 contract (exact
  to 2nd order) does not offer; for log-Cholesky they are cheap (the diagonal
  is the only nonlinearity). The numerical fallback costs 3.9 s at p = 4 for
  order 3 alone (560 components, two Hessians each), so the extension becomes
  worth it the day something consumes multivariate third derivatives --
  a decision for Giovanni, recorded 2026-08-04.
- **More distributions for `distributions7`** (assessment 2026-08-04, Giovanni
  asked whether it was worth it). What the fourteen do NOT cover, in the order
  the model layer will want them:
  - ~~**Weibull** and **Gumbel/extreme value**~~ **done 2026-08-04.**
  - ~~**Skew normal / skew t**~~ **done 2026-08-04.**
  - ~~**exponential, geometric, chi-squared**~~ **done 2026-08-05**, all in
    their mean. Giovanni asked for them as a TEST BENCH for the wrappers, and
    that is what they turned out to be worth: `exponential_distrib()` agrees
    with `fixed(weibull_distrib(), sigma = 1)` and `geometric_distrib()` with
    `fixed(negbin_distrib(), theta = 1)` to at most 2e-15 on every quantity,
    both scales, observed and expected -- two independent implementations of
    one object, so no tolerance has to be chosen. None of the three is a Gamma
    with a parameter fixed: unit shape is the RELATION `sigma2 = mu^2` and the
    chi-squared `sigma2 = 2*mu`, and `fixed()` holds a parameter at a value,
    not a relation between two.
  - ~~**generalized Pareto**~~, ~~**beta-binomial**~~, ~~**NB1**~~,
    ~~**generalized gamma**~~, ~~**von Mises**~~ **done 2026-08-05.**
  - ~~**Dirichlet**~~ and ~~**multinomial**~~ **done 2026-08-05**, both on a
    `parameters7::simplex()`. See §4 for what they cost the validator.
  - **Skipped deliberately**: anything obtained by a wrapper already
    (`truncated`, `zero_*`, `transformation` cover lognormal-like families,
    hurdle models and Box-Cox families), and anything whose only interest is a
    reparametrization of a family already present. The book's
    `sec-transform-reach` now shows which standard names those are -- the
    inverse Gamma, the Rayleigh, the log-logistic, the Pareto -- with the
    distinction that decides it: `fixed()` holds a parameter at a VALUE and
    cannot impose a RELATION between parameters, which is why the exponential
    and the chi-squared are families of their own and the Rayleigh is not.
- ~~**`penalties7` design decisions agreed 2026-08-03**~~ **done 2026-08-06**:
  `piano_penalties7.txt` written first, phases 1-3 shipped the same day.
  The decisions below are preserved for the record: a penalty is
  rho(D beta; theta) -- a linear map, a scalar function, parameters. THREE
  branches: (i) quadratic with matrix P (covers ridge, spline Grams, fused
  quadratic; the correlated gaussian prior, later through parameters7);
  (ii) separable, built from a univariate distributions7 object applied
  coordinatewise to D beta -- ridge IS fixed(gaussian_distrib(), mu = 0), the
  heavy-tailed prior IS fixed(student_t_distrib(), mu = 0); (iii) families
  defined by their derivative (SCAD, MCP), which are improper priors and
  cannot come from a distribution. The NORMALIZING CONSTANT IS KEPT (Giovanni,
  explicit): dropping it makes the hyperparameter estimation degenerate
  (lambda -> 0), and estimating nu of a t prior requires it; even so, a
  spline's smoothing parameter needs REML/marginal likelihood, not joint
  maximization, so the class must serve both routes. Proximal operators and a
  prox_grad() in optimizers7: DEFERRED, explicitly. Prerequisites already
  built: fixed() wrapper and distrib_cross_y() in distributions7.
- **`infer_npar()` cannot decide for a vectorized objective, and that is
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
- **Audit analitico 2026-08-06, e la chiusura completa che ne e' seguita.**
  Matrice famiglia x generica su 41 famiglie univariate, owner del metodo via
  `attr(m, "signature")[[1]]`. Alla fine della giornata:
  - **ordini 1-2, 3-4: chiusi OVUNQUE.** Ogni famiglia univariata e' analitica
    al quarto ordine tranne le componenti in nu della skew t, che non possono
    esserlo. Chiuse in giornata vonmises1/2, gengamma1, negbin1, betabinom1,
    gpd.
  - **derivate della cdf: chiuse dove esistono.** Le riparametrizzate portano
    le forme del genitore con la regola della catena sulla mappa, la via
    CONDIZIONATA a che il genitore ne abbia davvero una
    (`has_exact_cdf_deriv`), altrimenti nulla finge. Nuove: gumbel,
    exponential, weibull1, gpd, lognormal1 al secondo ordine, e il blocco
    location-scale di skewnormal1/skewt/student_t1/pseudohuber. Misurato
    contro l'integrale di aspettativa parziale: 1e-14/1e-16 dove le differenze
    finite arrivano a 9e-11/1.5e-9, e 1.03x-1.60x piu' veloce a n = 2e5.
  - **cross_y: da 2 famiglie a TUTTE le 30 continue.** Un test cammina il
    namespace e fallisce se una famiglia continua resta sul fallback.
  - **grad_y/hess_y**: scritte per gpd e gengamma1, delegate per invgauss2 e
    vonmises2.
  - **expected info** a fallback per pig1/pig2 (somma esatta sul supporto,
    accettabile) e skewnormal1/skewt (ostruzione nota).
  Per le discrete il fallback cdf e' la somma parziale esatta, quindi non e'
  una lacuna. Ostruzioni che restano tali: gamma/beta/chisq/gengamma nella
  direzione di forma per la cdf, skew t in nu.

  **Le tre tecniche che hanno fatto il lavoro**, in ordine di quanto hanno
  reso:
  - **far girare un assemblaggio generale agli ordini dove esiste gia' un
    kernel compilato.** gengamma1 riproduce gradiente e Hessiano a 2.6e-16,
    NB1 a 7e-15, betabinom1 esattamente: questo, non un confronto con
    differenze finite, e' cio' che autorizza a fidarsi agli ordini dove non
    c'e' nulla da confrontare.
  - **cercare l'identita' prima della formula.** Se la risposta entra solo
    attraverso z = (y-mu)/sigma allora d2l/dy dmu = -l_yy e
    d2l/dy dsigma = -z l_yy - l_y/sigma: una riga che ha chiuso nove famiglie,
    verificata contro le due forme chiuse che esistevano gia'.
  - **la delega**: una famiglia scritta come mappa di un'altra non deve
    riderivare nulla. betabinom1 e' betabinom2 alle forme, e ne eredita gli
    ordini 3-4 via `chain_derivatives`; cross_y e' una catena del PRIMO
    ordine, perche' la derivata nella risposta non interagisce con una
    riparametrizzazione dei parametri.

  **E due volte il riferimento e' stato il lato debole**, come gia' successo
  col Monte Carlo della GPD: (i) il kernel C++ di `xi_xi` della gpd somma
  -7.998 + 7.996 + 0.0021 per ottenere 9.92e-05 a xi*z piccolo, cioe' un
  pavimento di 1.79e-11, e Richardson su di esso amplifica a 1.65e-4 --
  l'assemblaggio nuovo e' meglio condizionato, e il kernel e' stato
  RISCRITTO lo stesso giorno (#89): W = log1p(u)*(z/u) esatto ovunque (la
  vecchia guardia su |xi| con serie corta in xi troncava (xi z)^3/4,
  quattro decimali a xi = 1e-8, z = 1e7), e le due derivate passano alla
  serie sotto |xi z| = 0.2 -- il corner torna a 2.2e-16 e il test
  kernel-contro-assemblaggio stringe da 1e-9 a 1e-13; (ii) la componente
  in xi di `cross_y` della gpd vale ESATTAMENTE zero a y = sigma, dove nessuna
  misura puntuale relativa esiste, e il confronto va fatto sulla scala
  dell'intero vettore.
- **A censored-likelihood front end.** `distrib_grad_cdf()` supplies the pieces, but
  nothing yet assembles them: a `fit_distrib(..., censored = )` taking a status vector,
  or a `Surv()`-like object. That is the step that turns the capability into a feature.
- **Second-order closed forms** for the families that only have the gradient
  (lognormal, invgauss, student t, pseudohuber). Mechanical but fiddly; the FD fallback
  is ~1e-8 and adequate.
- **Gamma and beta shape derivatives** would need the series representation of
  `d/da P(a,x)` (a convergent alternating sum, badly canceling for large x) or a
  Meijer-G evaluation. Probably *not* worth it: the FD of pgamma/pbeta already gives
  ~2e-10, and a hand-rolled series would likely be worse in some regime. Revisit only
  if profiling shows it matters.
- **Third and fourth cdf derivatives** would let truncation drop quadrature at orders
  3-4 as well.
- No `NEWS.md` on linkfunctions7, distributions7 or optimizers7. `basis7` and
  `parameters7` have one from their first commit, which is the right habit.
- **`expected_by_bartlett()` recomputes too much.** Inside the integrand it calls
  `observed_deriv()` for a whole order and keeps one component, once per block per
  partition per component — and the integrand runs at every quadrature node. Memoizing
  per integrand call looks worth it at orders 3 and 4, but `"bartlett"` is only the
  default at order 2, where the score is all it needs, so measure before touching it.
- **The logos were redesigned on 2026-08-07** (Giovanni: the old set was
  "blando", all thin single-weight chalk lines). The stickers are now
  hand-authored SVG under `logo/svg/`, one per package, written by
  `logo/author-svg.py`; `logo/make-logos.R` only rasterizes them into each
  repo's `man/figures/logo.png`, the `pkgdown/favicon/` set and a contact
  sheet. The system: white chalk `#F7F4D4` with real weight hierarchy and
  translucent fills; the sanguine chalk `#DD7644` survives ONLY at the center
  of the statmodels7 honeycomb (Giovanni stripped it from every other glyph).
  Each glyph is the package's own mathematics, sampled exactly: adaptive
  quadrature panels plus blackboard symbols (numericals7), the bare bold logit
  (linkfunctions7, reduced at his request to the curve alone), one filled
  gamma density (distributions7), contours with an exact-Cauchy-step zigzag in
  dashed chalk (optimizers7), five wide cubic B-spline bells with the central
  one filled (basis7 — his call: few and wide, not many and tall), concentric
  circles of the ALR chart pressing into the open 2-simplex (parameters7),
  ridge/lasso/SCAD told apart by stroke weight (penalties7), a seven-cell
  honeycomb with the members' lowercase names and numericals at the center as
  the root (statmodels7), and a scatter of open circles with a chalk fit line
  and a border-less confidence band (modelterms7, drawn ahead of its package;
  PNG stays in `logo/png/` until the repo exists). Rules that survived the
  redesign: palette from the `mvreg` sticker (chalkboard green `#3D6B4C`, rust
  border `#9C3E11`, chalk, monospace Courier wordmark); a glyph must not
  repeat another package's shape; no decoration without information; logos
  must be git-tracked (`git ls-files man/figures` is the check). One geometry
  trap worth keeping: the honeycomb's cells sit in the vertex directions
  (30°+60k) of a pointy-top hexagon, so they face each other vertex-to-vertex
  and need center distance > 2r — at spacing √3·r their strokes cross.
  The SVG sources and `author-svg.py` are carried by `site/sync-stack-files.sh`
  like CLAUDE.md itself.
- The portal is hand-written HTML. If the toolkit grows past a handful of packages it may
  deserve a generator, but not yet.
- The book now covers links, distributions, the transformation wrappers and fitting
  (§3.3, written 2026-07-30 as a chapter of its own and folded in when the book
  went to one chapter per package: the link scale as the place to optimize, Fisher
  scoring versus Newton and why the congruence corollary makes the expected
  information the right matrix to invert, the line search, the delta method as the same
  congruence applied to the inverse, and intervals built on the link scale and mapped
  back; revised 2026-08-03, when `fit_distrib()` started delegating to
  `optimizers7` — §3.3 now says the optimization is not written there, the step
  rule is the sufficient-decrease condition of §4.1 rather than plain halving,
  the non-PD Hessian is repaired by §4.2 rather than abandoned, and a chunk shows
  `method = lbfgs(...)`; revised again 2026-08-04 for the averaged objective and
  for `distrib_start()`), and §3.4 on the numerical fallbacks — the ratio-of-uniforms theorem, the
  mode recentering, the divergence transform at one edge and at two, the
  discrete cumulative table, and the two warnings that belong with them. Chapter 4
  (added 2026-07-31) covers `optimizers7`: descent directions and what a line search
  must guarantee, with Zoutendijk's theorem; Newton's Hessian repairs and why
  the eigenvalue floor is what keeps that theorem's angle bound alive; the secant
  equation and the BFGS update, then conjugate gradients
  and Barzilai-Borwein as the bottom of a ladder ordered by how much curvature a
  method stores, from p(p+1)/2 down to one number; the subdifferential, Fermat's
  rule, and the bundle subproblem with its dual; and the box
  reparametrization and Adam; and 4.5 (added 2026-07-31) on where a run begins --
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
  be taken against; orthonormalization, constraints and Demmler-Reinsch as one
  congruence; and tensor products, separability, and the two coefficient
  shapes a contraction accepts. Its five gates are in
  `book/R/basis-certificates.R`, injection-checked seventeen times.
  Section 3.5 (added 2026-08-04) covers the multivariate families: the
  flattening of a matrix parameter into scalars, the gaussian's score and
  Hessian in terms of the structure's `A_k` and `A_kl`, the collapse of the
  expected information onto `-0.5*tr(Sigma^-1 A_k Sigma^-1 A_l)`, the Student t
  as a scale mixture with its weight `c = (nu+p)/(nu+q)`, marginals and what a
  pairs plot shows, and what is refused (the cdf, the quantile). Chapter 6
  (added 2026-08-04) covers `parameters7`: the map and its free vector, why the
  construction is not a link -- the Jacobian is dense, and that is the whole
  reason -- the four things a consumer asks and the trace identity behind the
  log-determinant; the log-Cholesky map with its two derivative formulas and
  its log-determinant linear in the free vector, the diagonal families that
  reuse `linkfunctions7` links, and the fixed matrix carried by one scale;
  then rank, the pseudo-determinant, `lambda = r/(b'Pb)`, why an eigenvalue
  count is not a rank test, and why the solve is refused rather than replaced
  by a pseudo-inverse. Their gates are `book/R/multivariate-certificates.R`
  and `book/R/parameter-certificates.R`, injection-checked three times each.
  Section 3.5 gained a part on **what a fit reports** (2026-08-04): the delta
  method carrying the variance matrix onto the standard deviations and
  correlations, their Jacobian written out, and why each interval is built on
  log or on Fisher's z. Its gate transcribes that Jacobian by hand and also
  checks it against `numDeriv`, which shares no code with either --
  and the third parameters7 injection is worth keeping: the free-vector
  ordering could not be tested at `p = 3`, where reading the below-diagonal
  entries by row and by column give the SAME sequence, so the gate carries a
  `p = 4` case as well. What the
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
                               3.3 Fitting, 3.4 Fallbacks,
                               3.5 Several dimensions) /
4 The optimizers7 package     (4.1 Descent, 4.2 Curvature, 4.3 Non-smooth,
                               4.4 Constraints, 4.5 Starting points) /
5 The parameters7 package     (5.1 Matrix parameters, 5.2 The families,
                               5.3 Rank and the null space) /
6 The basis7 package          (6.1 Expansions, 6.2 Families, 6.3 Inner
                               products, 6.4 Transformations, 6.5 Several
                               variables) /
7 The penalties7 package      (7.1 Definition, 7.2 The three constructions,
                               7.3 Marginal pieces, 7.4 In practice) /
8 The modelterms7 package     (8.1 The term contract, 8.2 The terms,
                               8.3 The formula interpreter, 8.4 Degrees of
                               freedom and the fitted views, 8.5 In
                               practice) /
9 The numericals7 package     (9.1 Stencils, 9.2 Quadrature and series,
                               9.3 Special functions) /
A Notation / B References
```

**The chapter order was changed on 2026-08-08** (Giovanni: the order
followed package creation and did not have to). It is now pedagogical by
dependency: parameters7 right after its consumers (the multivariate section
of 3 and the fitting of 4), penalties7 beside its two suppliers (5 and 6),
numericals7 closing as the machinery. The section files were renamed with
their chapters (`_05a-parameters.qmd` etc.); the `{#sec-...}` anchors are
unchanged, so no cross-reference moved. `piano_libro.txt` carries the
export-versus-book coverage census (62 of distributions7's 122 exports,
the whole crit_* catalog and every check_* reader were unmentioned) and
the chapter-by-chapter "In practice" expansion program, to be consumed in
batches; the penalties7 section of 2026-08-08 is the template.

When `modelterms7` and the rest arrive they arrive as chapters, and nothing
already written moves. Quarto numbers chapters per *file*, so a package chapter is
one `.qmd` that pulls its sections in with `{{< include >}}` from `_`-prefixed
files that are **not** listed in `_quarto.yml`; that keeps 3.1, 3.2 numbering
without a two-thousand-line file. `toc-depth` is 4, `number-depth` 3.

⚠️ **Demoting headings when reorganizing:** shift `#` levels only outside fenced
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
sentences -- flowing prose. Methods whose behavior a generic's or constructor's
page fully specifies may share that page; methods with behavior of their own
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

**Prose register, revised 2026-08-07 (Giovanni, explicit): STYLE 1 of
`piano_refactoring.txt` — the dry classical manual — governs every prose
surface of the toolkit** ("mi sembra quello più accademico e chiaro"). Both
editorial passes are DONE and their records are `piano_refactoring.txt`
(book + portal, 2026-08-07) and `piano_refactoring_docs.txt` (package
documentation, 2026-08-08), carried in `site/stack/` like this file.

The rules of the register:
- the grammatical subject is the function or the quantity, never a
  personified package: "knows", "refuses", "says so", "prefers", "honest"
  with a software subject are out. The standing replacements: refuse →
  reject (or "signals an error"), "knows X" → "carries/has/records X",
  "the honest answer" → say what is returned (`NaN`, `NA`, an error);
- the last sentence of a paragraph is a fact or a consequence, never a
  maxim — the gnomic closer ("An explicit failure is preferable to a
  plausible wrong number", "both know it", "a wrong answer wearing the
  shape of a right one") is the single most recognizable LLM tell;
- no compressed two-beat contrasts, no "The first... The second..."
  mechanical enumerations, "rather than" only for real alternatives and at
  most once per paragraph; parentheses over em dashes; "worth
  noting/stating/recording" is filler and goes;
- section titles name their subject, without appended relative clauses;
- numeric examples live in parentheses or executed chunks, not in the
  narrative, and no counts or versions are hand-typed where an executable
  source exists (the site's fake `library()` output and three stale
  "fourteen/four sibling" counts were exactly that rot).

WHICH SURFACES the register covers, learned by finding tells in each: the
book, the portal, READMEs (edit the .Rmd and re-knit), roxygen (`#'` lines),
vignettes, NEWS.md, the `desc:` blocks of `_pkgdown.yml`, and **error
message strings** — a `stop()` text is user-visible prose ("It is refused
rather than fitted" became "rejected"). Code comments are exempt: they
follow the code's idiom and are not a reader-facing surface. Also exempt:
titles of cited papers (reported as published, "modelling" included),
"lattice" in the geometric grid sense, and "stacked matrix". Before
renaming a check label or an error text, grep the tests for the literal —
none matched this time, but a label is load-bearing by construction.

Spelling is American everywhere: parametrization/reparametrization (the
parameteris*/reparameteris* forms were swept 2026-08-07), recenter (not
recentre). Notation follows the book: the scalar step length is $\alpha$
and $s$ is the secant vector — a collision the notation appendix itself
had recorded without noticing.

The check greps for any NEW prose, before committing it: `refus`, `knows`,
`honest`, `wearing`, `worth `, `says so`, `parameteris`, `recentre`, a
paragraph-final aphorism, and hand-typed counts of families/packages/links.
This tightens, and does not replace, the register below.

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
must normalize on read and restore on write, or every marker containing `\n`
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

`assert_fit_ok()` (§3.3, `book/R/fit-certificates.R`) checks the four things
that chapter claims, each against a route the chapter does not itself use: that the
link-scale score really vanishes at the reported optimum, that
`V_theta = diag(h') V_eta diag(h')` with `h'` taken fresh from linkfunctions7, that
no interval escapes its parameter's domain or comes back inverted, and that Fisher
scoring, Newton and BFGS agree on the maximized log-likelihood. Injection-checked
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

Extended again 2026-08-04 with `.certify_starting_values()` and
`.certify_objective_scale()`. The first pins that `distrib_start()` returns the
sample mean and second moment for a multivariate gaussian and that the fit
confirms them in one iteration; its counterexample **had to be rewritten**,
because with the expected information corrected the iris fit started at the
origin now reaches the same maximum in 42 iterations rather than failing at
500. So the section's claim was weakened to what remains true — a start buys
work, not an answer — and the gate measures the ratio of iterations. **A
counterexample that stops holding is a statement in the prose that stopped
being true**, and the gate is what said so. The second checks the averaged
objective: the maximum against a closed form, `se(mu)` and `se(sigma)` against
`sigma/sqrt(n)` and `sigma/sqrt(2n)`, and the averaged score at the reported
optimum. Injection-checked at a factor of `sqrt(n)` (70.7) and at 5%.

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
`distribution-catalog.R` keep the printed LaTeX and an independent R
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
