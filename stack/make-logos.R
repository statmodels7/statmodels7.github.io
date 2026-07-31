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

# Flat-top hexagon centred on the canvas, the standard R sticker outline.
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
# travelling along it -- and it is the reason cg() exists, so the glyph is the
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
  cx <- W / 2; cy <- 238                          # glyph centre on the canvas
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

cat("wrote:\n")
for (f in list.files("logo", pattern = "[.]svg$", full.names = TRUE)) cat("  ", f, "\n")

# --- rasterise into each package's man/figures ---------------------------
if (requireNamespace("magick", quietly = TRUE)) {
  for (p in c("linkfunctions7", "distributions7", "optimizers7")) {
    src <- file.path("logo", paste0(p, ".svg"))
    dst_dir <- file.path(p, "man", "figures")
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    img <- magick::image_read_svg(src, width = 480)
    magick::image_write(img, file.path(dst_dir, "logo.png"), format = "png")
    cat("  ", file.path(dst_dir, "logo.png"), "\n")
  }
}
