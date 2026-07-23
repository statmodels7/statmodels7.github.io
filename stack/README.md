# Stack-level files

`statmodels7/` on disk is a plain directory holding the three repositories, so
these two files — which belong to the stack as a whole rather than to any one
package — would otherwise not be versioned anywhere.

| file | what it is |
|---|---|
| `CLAUDE.md` | orientation for the whole project: architecture, toolchain, conventions, and the knowledge that cost time to acquire. Read it before working on any package. |
| `make-logos.R` | draws the hex logos. Run from `statmodels7/`, not from here. |

They are copies. The working versions live at `statmodels7/CLAUDE.md` and
`statmodels7/logo/make-logos.R`; `../sync-stack-files.sh` refreshes these.
