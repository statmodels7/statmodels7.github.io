#!/bin/sh
# statmodels7/ is a plain directory, so the logo sources and the book have no
# repository of their own. They live here, and this copies the working
# versions in.
#
# THIS REPOSITORY IS PUBLIC. Only what is meant to be read by anyone belongs
# here: the portal, the book and the logo sources. The working files of the
# umbrella directory stay in the umbrella directory.
#
# Run from statmodels7/site:  sh sync-stack-files.sh
set -e
cp ../logo/make-logos.R    stack/make-logos.R
cp ../logo/author-svg.py   stack/author-svg.py
mkdir -p stack/logo-svg
cp ../logo/svg/*.svg       stack/logo-svg/
cp ../logo/svg/statmodels7.svg logo/statmodels7.svg

# The book: source under stack/book (so it is versioned), rendered HTML under
# book/ (so GitHub Pages serves it at /book/). Since 2026-10-08 the published
# book is the one in ../new_book; the earlier book stays in ../book and in the
# history of this repository. The render itself stays a manual step --
# `quarto render` in ../new_book -- because it executes R against the working
# tree and takes minutes.
#
# The rendered book is only replaced when there is one to replace it with. A
# render interrupted part way leaves ../new_book/_book empty, and wiping book/
# against it would take the published book offline while this script still
# reported success.
if [ ! -f ../new_book/_book/index.html ]; then
  echo "../new_book/_book/index.html is missing: render the book first" >&2
  exit 1
fi

rm -rf stack/book book
mkdir -p stack/book
cp -r ../new_book/chapters ../new_book/R ../new_book/assets stack/book/
cp ../new_book/_quarto.yml ../new_book/index.qmd ../new_book/README.md stack/book/
cp ../new_book/references.bib stack/book/
cp -r ../new_book/_book book

# An interrupted render also leaves its intermediates beside the sources, and
# those are not sources: they would be committed as if they were.
find stack/book/chapters \( -name '*.html' -o -name '*.knit.md' \
  -o -name '*.rmarkdown' -o -name '*_files' \) -exec rm -rf {} +

echo "synced; now commit in this repository"
