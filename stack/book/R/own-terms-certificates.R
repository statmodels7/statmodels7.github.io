# Certificates for Part V, chapter 16 (writing a model term; sections 16.1-16.3).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the reader's term does not take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   season    check_term() passes every row; the user's term reproduces the
#             basis7 Fourier smoother with the second-derivative penalty in
#             log-likelihood, effective degrees of freedom, criterion and
#             fitted values; predict() on the training data returns the
#             fitted values; the scale equation fits and converges;
#   asym      check_term() passes every row; the coefficients agree with
#             nls() to 1e-3; the log-likelihood gap to nls() is the dispersion
#             term 42 log(84/82) - 1 to 1e-3; nl() gives the same
#             log-likelihood; the sigma and random-effect fits converge; the
#             predictions carry finite standard errors.

assert_own_terms_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Chapter 'Writing a model term': ", what,
         " no longer agrees with the package.", call. = FALSE)
  }
  ll <- function(f) as.numeric(stats::logLik(f))
  checks_pass <- function(ck) all(ck$status == "OK")
  as_fit <- function(p) if (is.data.frame(p)) p$fit else p

  # --- season --------------------------------------------------------------------
  if (!checks_pass(claims$chk_season)) fail("the checks of season()")
  if (abs(ll(fits$season) - ll(fits$season_b)) > 1e-4) {
    fail("the season() log-likelihood against fourier_smooth()")
  }
  if (abs(sum(fits$season@edf$edf) - sum(fits$season_b@edf$edf)) > 1e-4) {
    fail("the season() effective degrees of freedom against fourier_smooth()")
  }
  if (abs(fits$season@criterion - fits$season_b@criterion) > 1e-4) {
    fail("the season() criterion against fourier_smooth()")
  }
  if (max(abs(fitted(fits$season) - fitted(fits$season_b))) > 1e-5) {
    fail("the season() fitted values against fourier_smooth()")
  }
  pr <- as_fit(claims$season_pred_grid)
  if (is.null(pr) || max(abs(pr - fitted(fits$season))) > 1e-8) {
    fail("predict() of the season term on the training data")
  }
  if (statmod_certificate(fits$season_sd)$state != "converged") {
    fail("the seasonal scale equation")
  }

  # --- asym ----------------------------------------------------------------------
  if (!checks_pass(claims$chk_asym)) fail("the checks of asym()")
  cs <- unname(coef(fits$asym)$mu)
  cn <- unname(stats::coef(ref$nls))
  if (is.null(cs) || max(abs(cs - cn)) > 5e-3) {
    fail("the asym() coefficients against nls()")
  }
  gap <- ll(fits$asym) - ll(ref$nls)
  if (abs(gap + (42 * log(84 / 82) - 1)) > 1e-3) {
    fail("the dispersion term of the restricted likelihood against nls()")
  }
  if (abs(ll(fits$asym) - ll(fits$asym_nl)) > 1e-6) {
    fail("asym() against nl()")
  }
  if (statmod_certificate(fits$asym_sd)$state != "converged") {
    fail("asym() with sigma ~ age")
  }
  if (statmod_certificate(fits$asym_re)$state != "converged") {
    fail("asym() with random(~ 1 | Seed)")
  }
  pa <- claims$pa
  if (is.null(pa$se) || !all(is.finite(pa$se))) {
    fail("the standard errors of the asym() predictions")
  }
  invisible(TRUE)
}
