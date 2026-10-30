# Certificates for the three sections of chapter 2 on sandwich variances,
# the profile likelihood and the bootstrap.
#
# Same contract as every gate of the book: what the sections state is
# re-derived here by a route the package does not itself take, and the
# render stops on any disagreement. Nothing in this file is visible to the
# reader.
#
#   sandwich   HC0 and HC3 of the schools model are the matrices
#              (X'X)^-1 X' diag(e^2 / (1 - h)^delta) X (X'X)^-1 built here;
#              HC0, HC2, HC3, HC4 are vcovHC()'s and HC1 is vcovHC()'s times
#              (n - k)/(n - k - 1); the state of largest leverage is Alaska,
#              with the largest income and expenditure; CR2 is vcovCL()'s
#              HC2 and is above HC0 on the slope; the HAC matrix is
#              NeweyWest()'s with three lags, the default bandwidth is 3 and
#              the HAC slope se is above the model one; the HC3/bayesian
#              ratios of the motorcycle mean and sigma cross 1, peak within
#              5 ms of an end of the time range and are below 1 over most of
#              the middle;
#   profile    the profile intervals of the Poisson model are glm()'s to
#              1e-5; the simultaneous ones are the 1-d ones at the chi2_2
#              cutoff to 1e-4; profile_taylor() is within 1e-3 and its
#              converged column is NA; at the ends of the REML interval of
#              sigma_b the REML statistic, computed here from the marginal
#              gaussian likelihood with dense matrices, is the chi2_1
#              quantile; the ML interval is lme4's; the REML interval lies to
#              the right of the ML one and is wider;
#   bootstrap  the first 30 case replicates of the GAG slope are glm() fits
#              with the counts of sample.int() from the same seed; the
#              percentile and basic intervals are the type-7 quantiles of the
#              replicates; the bootstrap se is within 10% of HC0; the
#              parametric replicates of sigma_b agree with bootMer()'s in mean
#              within four simulation errors; the cluster bootstrap variance
#              of the slope is within 25% of vcovCL()'s CR0, the model variance
#              of the slope is below 0.8 times it and the one of the intercept
#              within 15%; the subjects' own slopes vary by more than three
#              model standard errors.

assert_beyond_tools_ok <- function(e) {
  fail <- function(what) stop("chapter 2 tools gate: ", what, call. = FALSE)
  near <- function(a, b, tol) all(abs(a - b) <= tol * pmax(abs(b), 1e-300))
  # the checks run in a child of the chapter's environment, so that what
  # they assign does not reach the chapter
  ev <- new.env(parent = e)
  ev$fail <- fail
  ev$near <- near

  ## sandwich -------------------------------------------------------------
  local(envir = ev, {
    X <- model.matrix(lm_ps); r <- residuals(lm_ps); h <- hatvalues(lm_ps)
    XtXi <- solve(crossprod(X))
    hand <- function(delta) XtXi %*% crossprod(X * (r^2 / (1 - h)^delta), X) %*% XtXi
    v0 <- vcov(fit_ps, vcov_type = frequentist_vcov(sandwich = hc0_sandwich()), param = "mu")
    v3 <- vcov(fit_ps, vcov_type = frequentist_vcov(sandwich = hc3_sandwich()), param = "mu")
    if (!near(unname(v0), unname(hand(0)), 1e-8)) fail("HC0 is not the hand sandwich")
    if (!near(unname(v3), unname(hand(2)), 1e-8)) fail("HC3 is not the hand sandwich")
    rat <- se_hc / se_ref
    if (!near(rat[, c("HC0", "HC2", "HC3", "HC4")], 1, 1e-8)) fail("HC0/2/3/4 differ from vcovHC()")
    n <- nrow(schools); k <- length(coef(lm_ps))
    if (!near(rat[, "HC1"], sqrt((n - k) / (n - k - 1)), 1e-8)) fail("HC1 factor is not (n-k)/(n-k-1)")
    st <- names(which.max(h))
    if (st != "Alaska") fail("the state of largest leverage is not Alaska")
    if (schools[st, "Income"] != max(schools$Income) ||
        schools[st, "Expenditure"] != max(schools$Expenditure)) fail("Alaska is not the largest in income and expenditure")
    ref_cl <- sandwich::vcovCL(lm(y ~ x, data = PetersenCL), cluster = ~ firm, type = "HC2")
    v_cl <- vcov(fit_pc, vcov_type = frequentist_vcov(sandwich = cr2), param = "mu")
    if (!near(unname(v_cl), unname(ref_cl), 1e-8)) fail("CR2 differs from vcovCL()")
    if (!(se_pc_cr2[[2]] > se_pc_hc0[[2]])) fail("CR2 slope se is not above HC0")
    ref_hac <- sandwich::NeweyWest(lm(level ~ year, data = huron), lag = 3,
                                   prewhite = FALSE, adjust = FALSE)
    v_hac <- vcov(fit_lh, vcov_type = frequentist_vcov(sandwich = hac), param = "mu")
    if (!near(unname(v_hac), unname(ref_hac), 1e-8)) fail("HAC differs from NeweyWest()")
    if (floor(4 * (nrow(huron) / 100)^(2 / 9)) != 3) fail("default bandwidth is not 3")
    if (!(se_lh_hac[[2]] > se_lh_model[[2]])) fail("HAC slope se is not above the model one")
    if (!(min(r_mu) < 1 && max(r_mu) > 1)) fail("HC3/bayesian of the mean does not cross 1")
    ends <- range(mcycle$times)
    at_end <- function(t) min(abs(t - ends)) < 5
    if (!at_end(t_mu_max)) fail("HC3/bayesian of the mean does not peak at an end")
    if (!at_end(t_sigma_max)) fail("HC3/bayesian of sigma does not peak at an end")
    if (!(min(r_sigma) < 1 && max(r_sigma) > 1)) fail("HC3/bayesian of sigma does not cross 1")
    mid <- grid$times > ends[1] + 10 & grid$times < ends[2] - 10
    if (!(mean(r_mu[mid] < 1) > 0.5 && mean(r_sigma[mid] < 1) > 0.5)) fail("ratios are not below one over most of the middle")
  })

  ## profile --------------------------------------------------------------
  local(envir = ev, {
    if (!(d_prof < 1e-5)) fail("profile intervals differ from glm()")
    if (!(d_sim < 1e-4)) fail("simultaneous intervals differ from the chi2_2 1-d ones")
    if (!(d_taylor < 1e-3)) fail("profile_taylor() is not within 1e-3")
    if (!all(is.na(prof_t@points$converged))) fail("profile_taylor() points are not NA")
    # the REML criterion of y = X b + Z u + e, u ~ N(0, sb^2), e ~ N(0, s^2),
    # computed with dense matrices; sigma profiled out by optimize()
    y <- sleepstudy$Reaction
    X <- cbind(1, sleepstudy$Days)
    Z <- model.matrix(~ 0 + Subject, data = sleepstudy)
    ZZ <- tcrossprod(Z)
    reml <- function(sb, s) {
      V <- s^2 * diag(length(y)) + sb^2 * ZZ
      Vi <- solve(V); XVX <- crossprod(X, Vi %*% X)
      b <- solve(XVX, crossprod(X, Vi %*% y)); r <- y - X %*% b
      -0.5 * (determinant(V)$modulus + determinant(XVX)$modulus + sum(r * (Vi %*% r)))
    }
    prof <- function(sb) optimize(function(s) reml(sb, s), c(10, 60), maximum = TRUE)$objective
    lmax <- optimize(function(lsb) prof(exp(lsb)), log(c(10, 100)), maximum = TRUE)$objective
    stat <- 2 * (lmax - c(prof(ci_s[[1]]), prof(ci_s[[2]])))
    if (!near(stat, qchisq(0.95, 1), 5e-3)) fail(sprintf("REML statistic at the ends is %s", paste(signif(stat, 5), collapse = ", ")))
    if (!near(unname(ci_s_ml), unname(ci_lmer), 1e-3)) fail("ML profile interval differs from lme4")
    if (!(ci_s[[1]] > ci_s_ml[[1]] && ci_s[[2]] > ci_s_ml[[2]] &&
          diff(ci_s) > diff(ci_s_ml))) fail("REML interval is not to the right of and wider than ML")
  })

  ## bootstrap ------------------------------------------------------------
  local(envir = ev, {
    j <- grep("log_age", colnames(reps_g))
    n <- nrow(GAGurine)
    set.seed(11)
    for (b in 1:30) {
      w <- tabulate(sample.int(n, n, replace = TRUE), n)
      g <- glm(GAG ~ log_age, family = Gamma(link = "log"), data = GAGurine, weights = w)
      if (abs(coef(g)[[2]] - reps_g[b, j]) > 1e-5) fail(sprintf("case replicate %d is not the glm fit", b))
    }
    q <- quantile(reps_g[, j], probs = c(0.025, 0.975), type = 7, names = FALSE)
    if (!near(c(ci_perc$lower, ci_perc$upper), q, 1e-10)) fail("percentile interval is not the type-7 quantiles")
    ci_basic <- confint(fit_g, parm = "mu:log_age", interval = basic_interval())
    est <- coef(fit_g)$mu[["log_age"]]
    if (!near(c(ci_basic$lower, ci_basic$upper), 2 * est - rev(q), 1e-10)) fail("basic interval is not 2 est - quantiles")
    if (!(abs(se_boot_g[[2]] / se_hc0_g[[2]] - 1) < 0.1)) fail("bootstrap se is not within 10% of HC0")
    se_d <- sqrt(var(sb_star) / length(sb_star) + var(sb_lmer) / length(sb_lmer))
    if (!(abs(mean(sb_star) - mean(sb_lmer)) < 4 * se_d)) fail("parametric replicates differ from bootMer()")
    if (!(r_cl > 0.8 && r_cl < 1.25)) fail("cluster bootstrap is not within 25% of CR0")
    if (!(v_cl_model[2, 2] / v_cl_boot[2, 2] < 0.8)) fail("model variance of the slope is not below 0.8 times the cluster bootstrap")
    if (!(abs(v_cl_model[1, 1] / v_cl_boot[1, 1] - 1) < 0.15)) fail("model variance of the intercept is not within 15% of the cluster bootstrap")
    # the cause stated in the text: the subjects' own slopes vary
    sl <- sapply(split(sleepstudy, sleepstudy$Subject), function(d) coef(lm(Reaction ~ Days, data = d))[[2]])
    if (!(sd(sl) > 3 * sqrt(v_cl_model[2, 2]))) fail("the subjects' slopes do not vary beyond the model se")
  })
  invisible(TRUE)
}
