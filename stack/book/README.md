# The statmodels7 book

*Statistical models built from reusable parts.* The book shows how to fit
models with the toolkit and how to extend it with a distribution, a link, a
model term or an optimizer of one's own.

## Build

The book needs the Quarto CLI and R on the path.

``` bash
quarto render
```

The command runs from the book directory, and the output goes to `_book/`.

The packages are loaded from source with `pkgload`, from the sibling
directories `../numericals7` to `../statmodels7`, so the book always uses the
working tree and not an installed copy. Nothing needs to be installed first.

## Layout

```
_quarto.yml          project and format configuration
index.qmd            preface
chapters/            one file per chapter, in four parts:
                       Using statmod()
                       Special model terms
                       Choosing
                       Extending the toolkit
                     and the bibliography (A2-references.qmd)
R/                   _setup.R, which loads the packages, and one
                     certificate file per chapter
assets/theme.scss    the palette of the toolkit
references.bib       the bibliography
```

## The rule

Every number in the text comes from code that runs during the render. Each
chapter ends with a hidden chunk that calls the `assert_*_ok()` function of its
certificate file in `R/`. That function checks the claims of the chapter by a
route that the chapter does not take (a reference package, a closed form, or
the optimality conditions of a problem), and the render stops when a claim no
longer holds.

## Rendering

Quarto must find R 4.6.0 on `PATH`. The `freeze` option is off on purpose:
Quarto keys its cache on the `.qmd` file and does not see the files that a
chapter sources, which is where the certificates live, so a frozen render
would skip the checks.
