# The catalog's figures: each family drawn at several values of one of its
# parameters.
#
# The parameter varied is declared per family in the catalog's `curves` field,
# because which one changes the SHAPE is a fact about the family and not
# something to be guessed: for a location-scale family it is the scale, for the
# Student t the degrees of freedom, for the generalized Pareto the tail index.
# assert_curves_ok() refuses a record without the field, so a family cannot be
# added without saying what its picture should show.


# The parameter settings a record asks for: everything at `theta`, one
# parameter taken through the declared values.
distrib_curve_thetas <- function(rec) {
  pname <- names(rec$curves)[1L]
  lapply(rec$curves[[1L]], function(v) {
    th <- rec$theta
    th[[pname]] <- v
    th
  })
}

# The window to draw in. The record's `grid` exists for the density gate and
# is fixed per family, so on a family whose parameter concentrates the mass it
# is far wider than the picture wants; the range is therefore trimmed to where
# some setting still puts a visible amount of density. Only distrib_pdf() is
# called, never the quantile function, which for several families is itself a
# root-finding fallback.
distrib_curve_data <- function(rec) {
  d <- rec$obj()
  thetas <- distrib_curve_thetas(rec)
  discrete <- S7::S7_inherits(d, distributions7::discrete_distrib)

  rng <- range(unlist(lapply(thetas, rec$grid)), finite = TRUE)

  eval_at <- function(y) {
    lapply(thetas, function(th) {
      v <- suppressWarnings(distributions7::distrib_pdf(d, y, th))
      v[!is.finite(v)] <- NA_real_
      v
    })
  }
  make_y <- function(r) {
    if (discrete) seq.int(floor(r[1L]), ceiling(r[2L]))
    else seq(r[1L], r[2L], length.out = 400L)
  }

  # The window is the record's grid, trimmed below. Pushing it outward when an
  # endpoint still carries density was tried and is wrong: the test fires on a
  # heavy tail, which is exactly where the extra width is empty, and it turned
  # the Cauchy, the Student t and the Laplace into a spike in the middle of a
  # blank panel. Where a setting genuinely runs past the grid the remedy is to
  # choose settings the family's own grid covers.
  y <- make_y(rng)
  dens <- eval_at(y)

  peak <- do.call(pmax, c(dens, list(na.rm = TRUE)))
  keep <- which(is.finite(peak) & peak > 1e-4 * max(peak, na.rm = TRUE))
  if (length(keep) > 1L) {
    span <- seq.int(min(keep), max(keep))
    y <- y[span]
    dens <- lapply(dens, function(v) v[span])
  }
  list(y = y, dens = dens, thetas = thetas, discrete = discrete,
       param = names(rec$curves)[1L], values = rec$curves[[1L]])
}

# One panel, drawn by the package's own plot method.
#
# The overlay -- the colors, the line types, the key and its side, the offset
# stems of a discrete family -- belongs to distributions7, where a reader of
# the package gets it too, and the book adds only what it alone knows: the
# window, which comes from the record's grid rather than from quantiles (the
# 99.5% quantile of a Cauchy is two orders of magnitude past anything worth
# drawing), and an empty title, the section heading having named the family
# already.
plot_distrib_entry <- function(rec) {
  cd <- distrib_curve_data(rec)
  theta <- rec$theta
  theta[[cd$param]] <- cd$values

  op <- graphics::par(mar = c(4, 4, 0.6, 0.6), cex = 0.9)
  on.exit(graphics::par(op), add = TRUE)

  plot(rec$obj(), theta,
       xlim = range(cd$y), main = "", xlab = "y",
       ylab = if (cd$discrete) "probability" else "density",
       bty = "n", lwd = 2)
  invisible(cd)
}

# The sentence under the picture: which parameter moves and where the others
# are held, so the figure is readable without counting back to the code.
distrib_curve_caption <- function(rec) {
  what <- if (S7::S7_inherits(rec$obj(), distributions7::discrete_distrib)) {
    "The probability mass"
  } else {
    "The density"
  }
  n <- length(rec$curves[[1L]])
  count <- c("two", "three", "four", "five")[max(1L, min(4L, n - 1L))]
  held <- setdiff(names(rec$theta), names(rec$curves)[1L])
  if (!length(held)) {
    return(sprintf("%s at %s values of `%s`.", what, count,
                   names(rec$curves)[1L]))
  }
  sprintf("%s at %s values of `%s`, with %s held fixed.",
          what, count, names(rec$curves)[1L],
          paste(sprintf("`%s = %s`", held,
                        format(unlist(rec$theta[held]), trim = TRUE)),
                collapse = ", "))
}

# Silent gate. A picture that is empty, flat or non-finite says nothing and
# would not be noticed in forty-two of them, so each is checked for the three
# ways it can be uninformative: a setting outside the parameter's domain, a
# density that is not finite and positive somewhere, and a set of settings
# that draw the same curve.
assert_curves_ok <- function() {
  missing <- names(DISTRIBS)[!vapply(DISTRIBS,
                                     function(r) is.list(r$curves) &&
                                       length(r$curves) == 1L &&
                                       length(r$curves[[1L]]) >= 2L,
                                     logical(1))]
  if (length(missing)) {
    stop("Families without a 'curves' field: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  for (id in names(DISTRIBS)) {
    rec <- DISTRIBS[[id]]
    d <- rec$obj()
    pname <- names(rec$curves)[1L]

    if (!pname %in% d@params) {
      stop(sprintf("Distribution '%s': curves vary '%s', which is not a parameter.",
                   id, pname), call. = FALSE)
    }
    b <- d@params_bounds[[pname]]
    v <- rec$curves[[1L]]
    if (any(!is.finite(v)) || any(v <= b[1L]) || any(v >= b[2L])) {
      stop(sprintf("Distribution '%s': a value of '%s' is outside (%s, %s).",
                   id, pname, format(b[1L]), format(b[2L])), call. = FALSE)
    }

    cd <- distrib_curve_data(rec)
    for (i in seq_along(cd$dens)) {
      f <- cd$dens[[i]]
      if (!any(is.finite(f) & f > 0)) {
        stop(sprintf("Distribution '%s': the density at %s = %s is nowhere finite and positive.",
                     id, pname, format(v[i])), call. = FALSE)
      }
    }
    # the settings must actually differ on the page
    spread <- max(vapply(seq_along(cd$dens)[-1L], function(i) {
      a <- cd$dens[[1L]]; b2 <- cd$dens[[i]]
      ok <- is.finite(a) & is.finite(b2)
      if (!any(ok)) return(0)
      max(abs(a[ok] - b2[ok]))
    }, numeric(1)))
    if (!is.finite(spread) || spread < 1e-6) {
      stop(sprintf("Distribution '%s': the curves for '%s' are indistinguishable.",
                   id, pname), call. = FALSE)
    }
  }
  invisible(TRUE)
}
