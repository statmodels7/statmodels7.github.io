#!/bin/sh
# statmodels7/ is a plain directory, so the files that belong to the stack as a
# whole -- the orientation file, the logo generator, and the book -- have no
# repository of their own. They live here, and this copies the working versions in.
#
# Run from statmodels7/site:  sh sync-stack-files.sh
set -e
cp ../CLAUDE.md            stack/CLAUDE.md
cp ../logo/make-logos.R    stack/make-logos.R
cp ../logo/statmodels7.svg logo/statmodels7.svg

# The package plans are decision records drafted before any code exists, so
# the umbrella directory is the only place they live; carry them like CLAUDE.md.
cp ../piano_*.txt          stack/

# The book: source under stack/book (so it is versioned), rendered HTML under
# book/ (so GitHub Pages serves it at /book/). The render itself stays a manual
# step -- `quarto render` in ../book -- because it executes R against the
# working tree and takes minutes.
rm -rf stack/book book
mkdir -p stack/book
cp -r ../book/chapters ../book/R ../book/assets stack/book/
cp ../book/_quarto.yml ../book/index.qmd ../book/README.md stack/book/
cp ../book/references.bib stack/book/
cp -r ../book/_book book

echo "synced; now commit in this repository"
