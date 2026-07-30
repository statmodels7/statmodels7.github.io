# The statmodels7 book

The accountability document for the toolkit: every formula the packages implement,
derived from the definition, and checked against the running code while the book
renders.

## Build

Needs the Quarto CLI and R on the path.

``` bash
cd book
quarto render
```

Output goes to `_book/`. Open `_book/index.html`.

After rendering, check that no pipe escaped into a formula:

``` bash
grep -c '&#124;' _book/chapters/*.html
```

Every count must be zero. A `|` written inside a math expression that sits in a
markdown table breaks the cell, and kable escapes it to `&#124;`, which MathJax
prints verbatim. Use `\lvert ... \rvert` for absolute values everywhere.

The two packages are loaded **from source** with `pkgload`, from `../linkfunctions7`
and `../distributions7`, so the book always documents the working tree rather than
an installed snapshot. Nothing needs to be installed first.

## Layout

**One chapter per package.** The sections of a package chapter live in their own
files, prefixed with `_` so Quarto does not promote them to chapters, and are
pulled in with `{{< include >}}` by the chapter that owns them. That is what keeps
the numbering at 3.1, 3.2, ... without putting two thousand lines in one file, and
it is how the next package will be added.

```
book/
  _quarto.yml          project and format configuration
  index.qmd            preface: what the book is and the rule it follows
  chapters/
    01-introduction.qmd        what the toolkit is for
    02-linkfunctions7.qmd      the inverse function theorem, 16 links
    03-distributions7.qmd      the distributions7 chapter; includes the four below
    _03a-distributions.qmd       likelihood theory, Bartlett, link scale, 14 families
    _03b-transformations.qmd     change of variables, zero wrappers, truncation
    _03c-fitting.qmd             maximum likelihood, Fisher scoring, delta method
    _03d-fallbacks.qmd           what a density-only distribution gets for free
    A1-notation.qmd            notation, session info
  R/
    _setup.R                       shared setup and the certification helpers
    link-formulas.R                the 16 link records + their certificates
    distribution-catalogue.R       the 14 distribution records + certificates
    transformation-certificates.R  transformer and zero-wrapper certificates
    fit-certificates.R             the fitting chapter's certificates
  assets/theme.scss      the toolkit's palette
```

## Notation

Derivatives of the log-likelihood carry parenthesized superscripts —
`\ell^{(i)}`, `\ell^{(ij)}`, `\ell^{(ijk)}` — never subscripts. A subscript on
`\ell` is reserved for what subscripts conventionally mean in likelihood work:
the contribution of one observation, or the model a log-likelihood belongs to
(`\ell_Y`, `\ell_T`). Truncation points are `L` and `U`, not `\ell` and `u`.

## The one rule

Every formula printed in the book is checked against the implementation while the
book renders, by the hidden `assert_*_ok()` gates at the end of each chapter. The
`R/` files keep the printed LaTeX and a machine-readable transcription of it **in
the same record**, so the text cannot drift away from what is tested without the
render failing. None of this is mentioned in the book itself: the reader is
entitled to assume the formulas shown are the ones that run, and the checking is
ours to do.

If you edit a formula in `R/link-formulas.R` or `R/distribution-catalogue.R`, you
are editing both the equation the reader sees and the thing that gets tested.
That is deliberate.

## Adding a chapter or a section

A **chapter** is a package. Add the `.qmd` under `chapters/`, list it in
`_quarto.yml`, and start it with

``` r
source("R/_setup.R")
```

in an `include: false` chunk. The working directory is always `book/`
(`execute-dir: project`).

A **section** of an existing package chapter is a `_`-prefixed file starting at
`##`, added to that chapter with `{{< include >}}` and *not* listed in
`_quarto.yml`. It needs its own `source("R/_setup.R")` chunk with a label unique
across the whole book.

## Rendering

`quarto render`, from this directory, with R on `PATH` — Quarto resolves R from
the registry, finds an empty `R-4.5.1`, and stops with "Unable to locate an
installed version of R" otherwise.

`freeze` is **off**, deliberately. Quarto keys its cache on the `.qmd` and does
not see the files a chapter `source()`s, which is exactly where the formulas and
their transcriptions live; with `freeze: auto` a render finished in seconds
carrying stale text and, worse, without re-running a single gate. Renders take
twenty minutes and mean something. Do not turn it back on.

After rendering, `grep -c '&#124;' _book/chapters/*.html` must be all zeros: a
`|` inside math in a table splits the cell, so absolute values are written
`\lvert ... \rvert`.
