# Shared setup for every chapter of the book.
#
# The book documents the *working tree*, not an installed snapshot: the packages
# are loaded from source with pkgload. If a formula in the text and the code in
# the repository ever disagree, the certification chunks below say so while the
# book renders, which is the whole point of the exercise.

suppressPackageStartupMessages({
  library(S7)
  library(numDeriv)
})

# _quarto.yml sets execute-dir: project, so the working directory is always
# statmodels7/book and the two packages sit one level up.
.toolkit_root <- normalizePath("..", mustWork = TRUE)

# Load both packages from source, quietly.
suppressMessages({
  pkgload::load_all(file.path(.toolkit_root, "numericals7"), quiet = TRUE, export_all = FALSE)
  pkgload::load_all(file.path(.toolkit_root, "linkfunctions7"), quiet = TRUE, export_all = FALSE)
  pkgload::load_all(file.path(.toolkit_root, "distributions7"), quiet = TRUE, export_all = FALSE)
  pkgload::load_all(file.path(.toolkit_root, "optimizers7"), quiet = TRUE, export_all = FALSE)
  pkgload::load_all(file.path(.toolkit_root, "parameters7"), quiet = TRUE, export_all = FALSE)
  pkgload::load_all(file.path(.toolkit_root, "basis7"), quiet = TRUE, export_all = FALSE)
})

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 4.2,
  dpi = 150,
  out.width = "100%",
  warning = FALSE,
  message = FALSE
)

options(digits = 7, width = 82)


# ---------------------------------------------------------------------------
# Certification helpers
#
# Three of them, answering three different questions:
#
#   certify_formula()  does the formula printed in this book agree with the
#                      function the package actually runs?
#   certify_numeric()  does the package's analytical derivative agree with a
#                      derivative computed numerically, with no shared code?
#   certify_table()    collect a set of such comparisons into one table.
#
# All of them report a *relative* error, max |a - b| / max(|b|, 1), and a
# verdict against a stated tolerance. Nothing here is allowed to fail silently:
# an error inside a comparison is caught and shown as the verdict.
# ---------------------------------------------------------------------------

rel_err <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1))
}

# Compare two functions of one argument over a grid.
compare_on_grid <- function(f_book, f_pkg, grid) {
  tryCatch(rel_err(f_book(grid), f_pkg(grid)), error = function(e) NA_real_)
}

verdict <- function(err, tol) {
  ifelse(is.na(err), "ERROR", ifelse(err < tol, "agree", "DISAGREE"))
}

# Render a certification data.frame as a compact table.
certification_table <- function(df, caption, tol = 1e-8) {
  df$verdict <- verdict(df$error, tol)
  df$error <- formatC(df$error, format = "e", digits = 2)
  knitr::kable(df, caption = caption, align = "lrl", row.names = FALSE)
}

# One numerical derivative of a scalar function, evaluated pointwise.
#
# Deliberately ONE layer. Nesting numDeriv::grad to reach a fourth derivative
# does not work: each layer multiplies the extrapolation error, and by order
# three the reference is pure noise -- for the identity link, whose third
# derivative is exactly zero, nested differentiation returns a number of order
# one. A comparison against such a reference reports a failure for code that is
# exact, so it must never appear in this book. Every high-order comparison here
# differentiates the analytical derivative of the order below, exactly once,
# which is also what check_link() does.
num_grad1 <- function(f, x) {
  vapply(x, function(xi) numDeriv::grad(func = f, x = xi), numeric(1))
}

# A grid of n points strictly inside (lo, hi), pulled away from the boundaries.
# numDeriv's Richardson stencil reaches about 8e-4 * |x|; any grid closer than
# that to a boundary is differentiated using points outside the domain, which
# come back NaN and would report a failure for code that is exact.
inner_grid <- function(lo, hi, n = 25, pad = 0.02) {
  if (!is.finite(lo) && !is.finite(hi)) return(seq(-3, 3, length.out = n))
  if (!is.finite(hi)) return(lo + exp(seq(log(0.05), log(6), length.out = n)))
  if (!is.finite(lo)) return(hi - exp(seq(log(0.05), log(6), length.out = n)))
  w <- hi - lo
  seq(lo + pad * w, hi - pad * w, length.out = n)
}
