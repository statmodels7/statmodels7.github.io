# The statmodels7 book

The accountability document for the stack: every formula the packages implement,
derived from the definition, and checked against the running code while the book
renders.

## Build

Needs the Quarto CLI and R on the path.

``` bash
cd book
quarto render
```

Output goes to `_book/`. Open `_book/index.html`.

The two packages are loaded **from source** with `pkgload`, from `../linkfunctions7`
and `../distributions7`, so the book always documents the working tree rather than
an installed snapshot. Nothing needs to be installed first.

## Layout

```
book/
  _quarto.yml          project and format configuration
  index.qmd            preface: what the book is and the rule it follows
  chapters/
    01-introduction.qmd      what the framework is for
    02-link-functions.qmd    theory, the inverse function theorem, 16 links
    03-distributions.qmd     likelihood theory, Bartlett, link scale, 14 distributions
    04-transformations.qmd   change of variables, zero-inflation, zero-adjustment
    A1-notation.qmd          notation, session info
  R/
    _setup.R                       shared setup and the certification helpers
    link-formulas.R                the 16 link records + their certificates
    distribution-catalogue.R       the 14 distribution records + certificates
    transformation-certificates.R  transformer and zero-wrapper certificates
  assets/theme.scss      the stack's palette
```

## The one rule

Every formula printed in the book is checked against the implementation while the
book renders. The `R/` files keep the printed LaTeX and a machine-readable
transcription of it **in the same record**, so the text cannot drift away from
what is tested without the certification table reporting `DISAGREE`.

If you edit a formula in `R/link-formulas.R` or `R/distribution-catalogue.R`, you
are editing both the equation the reader sees and the thing that gets tested.
That is deliberate.

## Adding a chapter

Add the `.qmd` under `chapters/`, list it in `_quarto.yml`, and start it with

``` r
source("R/_setup.R")
```

in an `include: false` chunk. The working directory is always `book/`
(`execute-dir: project`).
