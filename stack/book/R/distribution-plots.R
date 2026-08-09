# The catalog's figures: one panel per parameter, each drawing the family at
# three values of that parameter with the others held at the record's theta.
#
# The values are declared per family in the catalog's `curves` field rather
# than derived, because what is worth showing is a fact about the family: a
# location wants positions either side of its default, a scale a small, a
# middle and a large one, a shape the values that change the shape. Every
# parameter gets a panel, so none of them goes unshown, and
# assert_curves_ok() rejects a record whose `curves` are not exactly the
# family's parameters in the family's own order.


# The parameter settings one panel asks for: everything at `theta`, the named
# parameter taken through its declared values.
distrib_curve_thetas <- function(rec, pname) {
  lapply(rec$curves[[pname]], function(v) {
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
distrib_curve_data <- function(rec, pname = names(rec$curves)[1L]) {
  d <- rec$obj()
  thetas <- distrib_curve_thetas(rec, pname)
  discrete <- S7::S7_inherits(d, distributions7::discrete_distrib)

  rng <- range(unlist(lapply(thetas, rec$grid)), finite = TRUE)
  # A parameter free on the whole line moves the density ALONG the axis, and
  # the record's grid is sized for one position of it, so the candidate range
  # is widened by how far the values travel. Every other kind of parameter
  # reshapes the density in place and needs no room.
  b <- d@params_bounds[[pname]]
  if (is.infinite(b[1L]) && is.infinite(b[2L])) {
    travel <- diff(range(rec$curves[[pname]]))
    rng <- rng + c(-1, 1) * travel / 2
  }

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

  keep <- if (discrete) {
    # For a discrete family the criterion is the mass itself rather than a
    # ratio of densities: the window is the smallest one holding 99% of
    # every setting, so a long thin tail does not fill the panel with stems
    # that carry nothing. negbin is the only family this bites, its grid
    # running to 40 while the mass is under 15.
    cuts <- vapply(dens, function(v) {
      v[!is.finite(v)] <- 0
      tot <- sum(v)
      if (!(tot > 0)) return(c(1, length(v)))
      c(which(cumsum(v) >= 0.005 * tot)[1L],
        which(cumsum(v) >= 0.995 * tot)[1L])
    }, numeric(2))
    seq.int(min(cuts[1L, ]), max(cuts[2L, ]))
  } else {
    peak <- do.call(pmax, c(dens, list(na.rm = TRUE)))
    which(is.finite(peak) & peak > 1e-4 * max(peak, na.rm = TRUE))
  }
  if (length(keep) > 1L) {
    span <- seq.int(min(keep), max(keep))
    y <- y[span]
    dens <- lapply(dens, function(v) v[span])
  }
  list(y = y, dens = dens, thetas = thetas, discrete = discrete,
       param = pname, values = rec$curves[[pname]])
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
  d <- rec$obj()
  pars <- names(rec$curves)
  op <- graphics::par(mfrow = c(1, length(pars)),
                      mar = c(4, 4, 1.6, 0.6), cex = 0.75)
  on.exit(graphics::par(op), add = TRUE)

  out <- lapply(pars, function(pname) {
    cd <- distrib_curve_data(rec, pname)
    theta <- rec$theta
    theta[[pname]] <- cd$values
    # bty = "l" rather than "n": with no box R draws each axis only between
    # its extreme ticks, and a discrete family whose support runs past the
    # last tick then has stems standing on nothing
    plot(d, theta,
         xlim = range(cd$y), main = "", xlab = "y",
         ylab = if (cd$discrete) "probability" else "density",
         bty = "l", lwd = 2)
    cd
  })
  names(out) <- pars
  invisible(out)
}

# The sentence under the picture: which parameter moves and where the others
# are held, so the figure is readable without counting back to the code.
distrib_curve_caption <- function(rec) {
  what <- if (S7::S7_inherits(rec$obj(), distributions7::discrete_distrib)) {
    "the probability mass"
  } else {
    "the density"
  }
  pars <- names(rec$curves)
  fmt <- function(v) paste(format(v, trim = TRUE), collapse = ", ")
  if (length(pars) == 1L) {
    return(sprintf("%s at `%s` = %s.",
                   sub("^t", "T", what), pars, fmt(rec$curves[[1L]])))
  }
  # the values themselves are in each panel's legend; the caption says what
  # the panels are and where the parameters not varying are held
  sprintf("One panel per parameter: %s at three values of each, the others held at %s.",
          what,
          paste(sprintf("`%s = %s`", pars,
                        format(unlist(rec$theta[pars]), trim = TRUE)),
                collapse = ", "))
}

# Silent gate. A picture that is empty, flat or non-finite says nothing and
# would not be noticed in forty-two of them, so each is checked for the three
# ways it can be uninformative: a setting outside the parameter's domain, a
# density that is not finite and positive somewhere, and a set of settings
# that draw the same curve.
assert_curves_ok <- function() {
  for (id in names(DISTRIBS)) {
    rec <- DISTRIBS[[id]]
    d <- rec$obj()

    if (!is.list(rec$curves) || !length(rec$curves)) {
      stop(sprintf("Distribution '%s': no 'curves' field.", id), call. = FALSE)
    }
    # every parameter gets a panel, so no parameter's effect goes unshown, and
    # the panels appear in the order the family declares its parameters
    if (!identical(names(rec$curves), d@params)) {
      stop(sprintf(paste0("Distribution '%s': 'curves' names (%s) are not the ",
                          "family's parameters (%s). Every parameter needs a ",
                          "panel, in the family's own order."),
                   id, paste(names(rec$curves), collapse = ", "),
                   paste(d@params, collapse = ", ")), call. = FALSE)
    }

    for (pname in d@params) {
      v <- rec$curves[[pname]]
      if (length(v) < 2L) {
        stop(sprintf("Distribution '%s': '%s' needs at least two values.",
                     id, pname), call. = FALSE)
      }
      b <- d@params_bounds[[pname]]
      if (any(!is.finite(v)) || any(v <= b[1L]) || any(v >= b[2L])) {
        stop(sprintf("Distribution '%s': a value of '%s' is outside (%s, %s).",
                     id, pname, format(b[1L]), format(b[2L])), call. = FALSE)
      }

      cd <- distrib_curve_data(rec, pname)
      for (i in seq_along(cd$dens)) {
        fi <- cd$dens[[i]]
        if (!any(is.finite(fi) & fi > 0)) {
          stop(sprintf("Distribution '%s': at %s = %s nothing in the window is finite and positive.",
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
  }
  invisible(TRUE)
}
