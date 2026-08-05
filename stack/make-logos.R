# Hex logos for the statmodels7 stack.
#
# One shape, one palette, one layout; the glyph is what changes. Each glyph is
# the actual mathematics the package is about -- the sigmoid really is
# plogis(), the density really is dgamma() -- sampled here and written out as an
# SVG path, so the curves are correct rather than merely curve-shaped.
#
# Run with:  Rscript logo/make-logos.R
# Writes logo/<pkg>.svg and <pkg>/man/figures/logo.png

suppressMessages({
  library(svglite)
})

# --- shared look -----------------------------------------------------------

# Palette taken from Giovanni's mvreg sticker, so the stack sits alongside his
# earlier work rather than beside it: chalkboard green, a rust border, chalk.
INK    <- "#3D6B4C"   # chalkboard green: hexagon fill
LINE   <- "#F7F4D4"   # chalk: the glyph and the wordmark
ACCENT <- "#F7F4D4"   # the 7, same chalk at a heavier weight
EDGE   <- "#9C3E11"   # rust: hexagon border
FILL   <- "#9BBBA2"   # a lighter green, for shaded area under a curve

W <- 520              # canvas, in the 1:1.1547 ratio a hex sticker wants
H <- 600

# Flat-top hexagon centered on the canvas, the standard R sticker outline.
hex_path <- function(cx, cy, r) {
  ang <- seq(90, 450, by = 60) * pi / 180
  x <- cx + r * cos(ang)
  y <- cy - r * sin(ang)
  paste0("M ", paste(sprintf("%.2f %.2f", x, y), collapse = " L "), " Z")
}

# A polyline through (x, y), mapped from data units into the canvas box.
curve_path <- function(x, y, xlim, ylim, box) {
  px <- box$x + (x - xlim[1]) / diff(xlim) * box$w
  py <- box$y + box$h - (y - ylim[1]) / diff(ylim) * box$h
  paste0("M ", paste(sprintf("%.2f %.2f", px, py), collapse = " L "))
}

svg_header <- function() {
  sprintf(paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" ',
    'viewBox="0 0 %d %d">\n'), W, H, W, H)
}

# Half-width of the hexagon at a given height, so the wordmark can be sized to
# fit rather than sized by eye. A pointy-top hexagon is full width only between
# the two flat sides; below that it closes to a point, and a word set there
# spills over the edge.
hex_halfwidth <- function(y, cy = H / 2, r = 250) {
  dy <- abs(y - cy)
  flat <- r * sin(30 * pi / 180)          # where the sides stop being vertical
  wmax <- r * cos(30 * pi / 180)
  if (dy <= flat) wmax else wmax * (r - dy) / (r - flat)
}

write_logo <- function(file, glyph_svg, name, y_text = 408, size = 46) {
  hex <- hex_path(W / 2, H / 2, 250)

  # Shrink the wordmark until it clears the edge with a margin. Courier advances
  # exactly 0.6 em per character, monospace being predictable that way.
  avail <- 2 * hex_halfwidth(y_text) - 52
  label <- paste0(name, "7")
  while (size > 18 && nchar(label) * 0.6 * size > avail) size <- size - 1

  txt <- paste0(
    sprintf('<path d="%s" fill="%s" stroke="%s" stroke-width="9"/>\n', hex, INK, EDGE),
    glyph_svg,
    # Monospace, as on the mvreg sticker these sit beside.
    sprintf(paste0('<text x="%d" y="%d" font-family="Courier New,Courier,monospace" ',
                   'font-size="%d" font-weight="400" letter-spacing="0.5" ',
                   'text-anchor="middle" fill="%s">%s7</text>\n'),
            W / 2, y_text, size, LINE, name)
  )
  writeLines(paste0(svg_header(), txt, "</svg>\n"), file)
  invisible(file)
}

box <- list(x = 152, y = 142, w = 216, h = 182)

# Shape chosen so the mode sits about a third of the way across: a gamma with a
# mode hard against the left edge reads as a spike rather than as a density.
DENS_SHAPE <- 4
DENS_RATE  <- 1
DENS_XMAX  <- 11

# --- linkfunctions7: the inverse logit, a link function made visible -------
{
  eta <- seq(-6, 6, length.out = 200)
  th  <- plogis(eta)
  p <- curve_path(eta, th, c(-6, 6), c(0, 1), box)

  # the two asymptotes the link maps onto, drawn faintly
  y0 <- box$y + box$h
  y1 <- box$y
  glyph <- paste0(
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            box$x - 10, y0, box$x + box$w + 10, y0, LINE),
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            box$x - 10, y1, box$x + box$w + 10, y1, LINE),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round"/>\n',
            p, LINE),
    # tangent at the origin: the derivative, which is what the package computes
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="3.5" stroke-linecap="round" opacity="0.85"/>\n',
            box$x + box$w / 2 - 52, box$y + box$h / 2 + 45,
            box$x + box$w / 2 + 52, box$y + box$h / 2 - 45, ACCENT)
  )
  write_logo("logo/linkfunctions7.svg", glyph, "linkfunctions")
}

# --- distributions7: a density, deliberately skewed, with its mass ---------
{
  x  <- seq(0, DENS_XMAX, length.out = 220)
  d  <- dgamma(x, shape = DENS_SHAPE, rate = DENS_RATE)
  d  <- d / max(d)
  p  <- curve_path(x, d, c(0, DENS_XMAX), c(0, 1.08), box)

  # filled area, closed along the baseline
  y0 <- box$y + box$h
  px <- box$x + x / DENS_XMAX * box$w
  py <- box$y + box$h - d / 1.08 * box$h
  area <- paste0("M ", sprintf("%.2f %.2f", px[1], y0), " L ",
                 paste(sprintf("%.2f %.2f", px, py), collapse = " L "),
                 sprintf(" L %.2f %.2f Z", px[length(px)], y0))

  glyph <- paste0(
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            box$x - 10, y0, box$x + box$w + 10, y0, LINE),
    sprintf('<path d="%s" fill="%s" opacity="0.30"/>\n', area, FILL),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round" stroke-linejoin="round"/>\n',
            p, LINE)
  )
  write_logo("logo/distributions7.svg", glyph, "distributions")
}

# --- optimizers7: steepest descent zigzagging down a valley ---------------
#
# The zigzag is not drawn freehand. It is steepest descent with an EXACT line
# search on f(x, y) = (x^2 + g y^2)/2, started at the worst point there is,
# (g t, t), for which the whole path can be written down: each step multiplies
# the iterate by (g-1)/(g+1) and flips the sign of y. That is the classic
# picture of why steepest descent is poor -- consecutive directions come out
# orthogonal, so the iterate crosses the floor of the valley instead of
# traveling along it -- and it is the reason cg() exists, so the glyph is the
# real thing rather than something zigzag-shaped.
#
# g = 4 rather than something more dramatic, for two reasons that pull the same
# way: the amplitude of the zigzag relative to the level set it starts on is
# 1/sqrt(g+1), so a very ill-conditioned bowl draws a path too flat to read,
# and the contraction (g-1)/(g+1) = 0.6 leaves enough visible steps to count.
#
# Two things about the drawing rather than the mathematics. The mapping is
# ISOTROPIC -- one number of pixels per unit in both directions -- so the level
# sets keep their true axis ratio of sqrt(g); the other two glyphs map x and y
# independently because a density has no aspect ratio to preserve, and this one
# does. And the whole picture is then rotated as a block, which is a rigid
# motion and so changes nothing about it: a valley running at an angle fills a
# hexagon better than one lying flat.
{
  gam <- 4
  stp <- function(p) {
    g <- c(p[1], gam * p[2])
    p - g * sum(g^2) / (g[1]^2 + gam * (gam * p[2])^2)
  }
  pts <- matrix(NA_real_, 6, 2)
  pts[1, ] <- c(gam * 0.7, 0.7)
  for (k in 2:nrow(pts)) pts[k, ] <- stp(pts[k - 1, ])

  lev0 <- (pts[1, 1]^2 + gam * pts[1, 2]^2) / 2   # the level the start sits on
  cx <- W / 2; cy <- 238                          # glyph center on the canvas
  sc <- 152 / sqrt(2 * 1.85 * lev0)                      # px per unit, both directions
  ang <- -25                                      # degrees, anticlockwise

  xy <- function(x, y) paste0(
    "M ", paste(sprintf("%.2f %.2f", cx + sc * x, cy - sc * y), collapse = " L "))

  # Level sets of the same quadratic: ellipses of axis ratio sqrt(g), at levels
  # falling geometrically so the spacing reads evenly.
  ell <- function(lev, n = 220) {
    a <- seq(0, 2 * pi, length.out = n)
    xy(sqrt(2 * lev) * cos(a), sqrt(2 * lev / gam) * sin(a))
  }
  # The outermost contour is drawn OUTSIDE the level the path starts on, so
  # the picture reads as a descent into the valley rather than a line leaving it.
  levs <- lev0 * c(1.85, 1.18, 0.66, 0.33, 0.13)
  fades <- c(0.28, 0.36, 0.44, 0.52, 0.62)

  glyph <- paste0(
    sprintf('<g transform="rotate(%.1f %.1f %.1f)">\n', ang, cx, cy),
    paste0(mapply(function(lv, op) sprintf(
      '<path d="%s Z" fill="none" stroke="%s" stroke-width="2.2" opacity="%.2f"/>\n',
      ell(lv), LINE, op), levs, fades), collapse = ""),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="3.2" stroke-linecap="round" stroke-linejoin="round"/>\n',
            xy(pts[, 1], pts[, 2]), LINE),
    # the minimum the path is converging on: the only mark on any of these
    # stickers that is a point rather than a curve, and it carries the one thing
    # the picture is about
    sprintf('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s"/>\n', cx, cy, ACCENT),
    "</g>\n"
  )
  write_logo("logo/optimizers7.svg", glyph, "optimizers")
}

# --- basis7: the B-spline bumps that a basis IS ---------------------------
#
# Not bumps drawn to look like splines: the real cubic B-spline basis on the
# unit interval, from the same recurrence the package calls, so the local
# support and the unequal end shapes are the ones a reader would get.
#
# The coloring says the one thing the picture is about. Every function is
# drawn faintly, and one of them -- an interior function, whose support is a
# full four knot intervals -- is drawn in the accent color, because a basis is
# a collection whose members are individually addressable: that is what makes
# it an object rather than a matrix somebody once built.
{
  knots <- seq(0, 1, length.out = 5)[-c(1, 5)]
  xx <- seq(0, 1, length.out = 320)
  bb <- splines::splineDesign(
    knots = c(rep(0, 4), knots, rep(1, 4)), x = xx, ord = 4, outer.ok = TRUE
  )
  hi <- 3L # the interior function drawn in the accent color

  y0 <- box$y + box$h
  paths <- vapply(seq_len(ncol(bb)), function(j) {
    curve_path(xx, bb[, j], c(0, 1), c(0, 1.08), box)
  }, character(1))

  glyph <- paste0(
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            box$x - 10, y0, box$x + box$w + 10, y0, LINE),
    paste0(vapply(setdiff(seq_along(paths), hi), function(j) sprintf(
      '<path d="%s" fill="none" stroke="%s" stroke-width="4" stroke-linecap="round" opacity="0.45"/>\n',
      paths[j], LINE), character(1)), collapse = ""),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round"/>\n',
            paths[hi], ACCENT)
  )
  write_logo("logo/basis7.svg", glyph, "basis")
}

# --- parameters7: a constrained set, and the chart that reaches it --------
#
# The picture is what the package is: a straight, unbounded grid on the free
# scale carried onto a bounded set. The set drawn is the 2-simplex -- the
# probability vectors on three categories -- because it is the one constrained
# set the package parametrizes that can be drawn honestly in two dimensions,
# and the grid is the image of straight lines under the additive log-ratio map
# simplex() implements, sampled here rather than sketched.
#
# The crowding is the whole statement. Equally spaced values on the free scale
# come out bunched against the boundary: a line at free value a meets the
# opposite edge at plogis(a), so the spacing along that edge is the logistic's,
# and no finite value reaches the edge at all. The grid therefore stops short
# and the outline is drawn whole, which is the difference between a
# constrained set and the unconstrained scale that maps onto it.
#
# Three pencils rather than the two the chart names. The simplex does not care
# which category is the reference, so the third family -- the one holding
# p1/p2 fixed -- belongs to the geometry as much as the two the coordinates
# happen to single out, and drawing all three keeps the picture symmetric
# under relabeling.
{
  # The map simplex() carries, written out rather than called, so this script
  # keeps drawing the mathematics without depending on the package.
  alr_inv <- function(e1, e2) {
    z <- c(e1, e2, 0)
    z <- z - max(z)
    p <- exp(z)
    p / sum(p)
  }

  cx <- box$x + box$w / 2
  cy <- box$y + box$h / 2 + 6
  R  <- 128                              # circumradius, in canvas units
  V  <- rbind(
    c(cx,                   cy - R),
    c(cx - R * sqrt(3) / 2, cy + R / 2),
    c(cx + R * sqrt(3) / 2, cy + R / 2)
  )
  proj <- function(p) c(sum(p * V[, 1]), sum(p * V[, 2]))

  path_of <- function(f, n = 90) {
    xy <- vapply(seq(-2.7, 2.7, length.out = n), function(s) proj(f(s)),
                 numeric(2))
    paste0("M ", paste(sprintf("%.2f %.2f", xy[1, ], xy[2, ]), collapse = " L "))
  }

  levels_free <- c(-2.4, -1.2, 0, 1.2, 2.4)
  grid_paths <- unlist(lapply(levels_free, function(a) {
    c(path_of(function(s) alr_inv(a, s)),        # p1/p3 held fixed
      path_of(function(s) alr_inv(s, a)),        # p2/p3 held fixed
      path_of(function(s) alr_inv(s, s - a)))    # p1/p2 held fixed
  }))

  outline <- paste0(
    "M ", paste(sprintf("%.2f %.2f", V[, 1], V[, 2]), collapse = " L "), " Z"
  )

  glyph <- paste0(
    paste0(vapply(grid_paths, function(p) sprintf(
      '<path d="%s" fill="none" stroke="%s" stroke-width="2.4" stroke-linecap="round" opacity="0.42"/>\n',
      p, LINE), character(1)), collapse = ""),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="5.5" stroke-linejoin="round"/>\n',
            outline, ACCENT)
  )
  write_logo("logo/parameters7.svg", glyph, "parameters")
}

# --- statmodels7: the umbrella, both curves layered -----------------------
{
  eta <- seq(-6, 6, length.out = 200)
  s   <- curve_path(eta, plogis(eta), c(-6, 6), c(0, 1), box)
  x   <- seq(0, DENS_XMAX, length.out = 220)
  d   <- dgamma(x, shape = DENS_SHAPE, rate = DENS_RATE); d <- d / max(d)
  dd  <- curve_path(x, d, c(0, DENS_XMAX), c(0, 1.08), box)
  y0  <- box$y + box$h

  glyph <- paste0(
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            box$x - 10, y0, box$x + box$w + 10, y0, LINE),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="5" stroke-linecap="round" opacity="0.9"/>\n',
            dd, FILL),
    sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round"/>\n',
            s, LINE)
  )
  write_logo("logo/statmodels7.svg", glyph, "statmodels")
}

# --- numericals7: the Gauss-Kronrod 7-15 weights, a quadrature made visible --

# The one discrete glyph in the family, deliberately: six siblings are curves,
# and the package is about numbers computed at points. The stems are the REAL
# Kronrod weights at the real nodes (the same constants gauss_kronrod15()
# carries), and the accent dots sit on the seven nodes the embedded Gauss rule
# shares -- one node set, two weight sets, which is the package's whole trick
# for getting an estimate and its error from a single evaluation.
{
  xh <- c(0.991455371120813, 0.949107912342759, 0.864864423359769,
          0.741531185599394, 0.586087235467691, 0.405845151377397,
          0.207784955007898, 0)
  wkh <- c(0.022935322010529, 0.063092092629979, 0.104790010322250,
           0.140653259715525, 0.169004726639267, 0.190350578064785,
           0.204432940075298, 0.209482141084728)
  nodes <- c(-xh[1:7], xh[8], rev(xh[1:7]))
  wk <- c(wkh[1:7], wkh[8], rev(wkh[1:7]))
  # the Gauss nodes are every second one, and carry the accent
  gauss <- seq_along(nodes) %% 2 == 0

  xlim <- c(-1.08, 1.08)
  ylim <- c(0, 0.24)
  px <- function(x) box$x + (x - xlim[1]) / diff(xlim) * box$w
  py <- function(y) box$y + box$h - (y - ylim[1]) / diff(ylim) * box$h

  base <- py(0)
  stems <- paste0(vapply(seq_along(nodes), function(i) {
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="5" stroke-linecap="round"/>\n',
            px(nodes[i]), base, px(nodes[i]), py(wk[i]), LINE)
  }, character(1)), collapse = "")
  dots <- paste0(vapply(which(gauss), function(i) {
    sprintf('<circle cx="%.1f" cy="%.1f" r="8" fill="%s"/>\n',
            px(nodes[i]), py(wk[i]), ACCENT)
  }, character(1)), collapse = "")

  glyph <- paste0(
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" opacity="0.35"/>\n',
            px(-1.05), base, px(1.05), base, LINE),
    stems, dots
  )
  write_logo("logo/numericals7.svg", glyph, "numericals")
}

cat("wrote:\n")
for (f in list.files("logo", pattern = "[.]svg$", full.names = TRUE)) cat("  ", f, "\n")

# --- rasterise into each package's man/figures ---------------------------
if (requireNamespace("magick", quietly = TRUE)) {
  # statmodels7 is the umbrella and now also a package -- the meta-package
  # that installs and attaches the rest -- so its glyph is rasterised too.
  for (p in c("linkfunctions7", "distributions7", "optimizers7", "basis7",
              "parameters7", "statmodels7", "numericals7")) {
    src <- file.path("logo", paste0(p, ".svg"))
    dst_dir <- file.path(p, "man", "figures")
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    img <- magick::image_read_svg(src, width = 480)
    magick::image_write(img, file.path(dst_dir, "logo.png"), format = "png")
    cat("  ", file.path(dst_dir, "logo.png"), "\n")
  }
}
