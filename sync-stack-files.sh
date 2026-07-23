#!/bin/sh
# statmodels7/ is a plain directory, so the files that belong to the stack as a
# whole -- the orientation file and the logo generator -- have no repository of
# their own. They live here, and this copies the working versions in.
#
# Run from statmodels7/site:  sh sync-stack-files.sh
set -e
cp ../CLAUDE.md            stack/CLAUDE.md
cp ../logo/make-logos.R    stack/make-logos.R
cp ../logo/statmodels7.svg logo/statmodels7.svg
echo "synced; now commit in this repository"
