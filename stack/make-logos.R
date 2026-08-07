# Rasterizes the hex logos for the statmodels7 toolkit.
#
# The stickers are hand-authored SVG, one per package, under logo/svg/ (drawn
# by logo/author-svg.py; the SVG files are the source of truth and can also be
# edited directly). R's only job here is what R is needed for: turning each
# SVG into the PNG pkgdown wants in man/figures/, the favicon set under
# pkgdown/favicon/, and a contact sheet for review.
#
# Palette (from Giovanni's mvreg sticker): chalkboard green #3D6B4C, chalk
# #F7F4D4, rust border #9C3E11, with a sanguine accent #DD7644 kept only at
# the center of the statmodels7 honeycomb.
#
# Run with:  Rscript logo/make-logos.R   (from logo/ or the umbrella root)

suppressMessages(library(magick))

if (basename(getwd()) != "logo" && dir.exists("logo")) setwd("logo")
root <- normalizePath("..")

# modelterms7 has no repository yet; its sticker waits here until it does.
pkgs <- c("numericals7", "linkfunctions7", "distributions7", "optimizers7",
          "basis7", "parameters7", "penalties7", "statmodels7")

manifest <- '{
  "name": "",
  "short_name": "",
  "icons": [
    {
      "src": "/web-app-manifest-192x192.png",
      "sizes": "192x192",
      "type": "image/png",
      "purpose": "maskable"
    },
    {
      "src": "/web-app-manifest-512x512.png",
      "sizes": "512x512",
      "type": "image/png",
      "purpose": "maskable"
    }
  ],
  "theme_color": "#ffffff",
  "background_color": "#ffffff",
  "display": "standalone"
}'

png_at <- function(svg, width, out) {
  img <- image_read_svg(svg, width = width)
  image_write(img, out, format = "png")
}

dir.create("png", showWarnings = FALSE)
for (pkg in pkgs) {
  svg <- file.path("svg", paste0(pkg, ".svg"))
  stopifnot(file.exists(svg))

  png_at(svg, 480, file.path("png", paste0(pkg, ".png")))

  fig <- file.path(root, pkg, "man", "figures")
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path("png", paste0(pkg, ".png")),
            file.path(fig, "logo.png"), overwrite = TRUE)

  fav <- file.path(root, pkg, "pkgdown", "favicon")
  dir.create(fav, recursive = TRUE, showWarnings = FALSE)
  file.copy(svg, file.path(fav, "favicon.svg"), overwrite = TRUE)
  png_at(svg, 180, file.path(fav, "apple-touch-icon.png"))
  png_at(svg,  96, file.path(fav, "favicon-96x96.png"))
  png_at(svg, 192, file.path(fav, "web-app-manifest-192x192.png"))
  png_at(svg, 512, file.path(fav, "web-app-manifest-512x512.png"))
  image_write(image_read_svg(svg, width = 32),
              file.path(fav, "favicon.ico"), format = "ico")
  writeLines(manifest, file.path(fav, "site.webmanifest"))

  cat("done", pkg, "\n")
}

# modelterms7: PNG only, kept beside the others for the day it has a home.
png_at(file.path("svg", "modelterms7.svg"), 480,
       file.path("png", "modelterms7.png"))

sheet <- image_montage(
  image_join(lapply(sort(list.files("png", full.names = TRUE)), image_read)),
  tile = "3x3", geometry = "+10+10", bg = "white")
image_write(sheet, "contact_sheet.png", format = "png")
cat("contact sheet done\n")
