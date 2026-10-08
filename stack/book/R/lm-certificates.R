# Certificates for Part II, chapter 3 (linear models and beyond).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   linear     the coefficients of the mean are lm()'s; sigma is
#              sqrt(RSS/(n - p)), lm()'s, and larger than sqrt(RSS/n); the
#              variance block of the mean is lm()'s and the row of sigma's
#              intercept is zero off its diagonal; the z of every row is lm()'s
#              t value, statmod()'s p-value is 2 pnorm(-|z|), the first three
#              rows have p below 1e-10 under both references and the
#              interaction's differ; the intervals are the estimate +- qnorm
#              se; logLik is the gaussian log-likelihood at the REML sigma and
#              lm()'s at sqrt(RSS/n), both with df 5, and AIC is -2 l + 10;
#              the marginal log-likelihood is logLik(lm, REML = TRUE);
#   factors    the separate-lines fit has the same fitted values and its
#              coefficients are the sums the text writes; the sum-contrast
#              fit, through linpar_options() and through linpar(), is lm()'s
#              with the same contrasts, its intercept the average of the two
#              intercepts and Insul1 half their difference; after the
#              insulation the line is lower and less steep;
#   sparse     the dense and the sparse fits agree, their designs are a base
#              matrix and a dgCMatrix, the sparse fit is faster, and the
#              default is sparse, the product of rows and indicator columns
#              exceeding 1e5;
#   spread     the variance regression of cars is gls(varExp)'s, coefficients
#              and REML criterion; sigma at 4 and 25 mph is exp(g0 + g1 v);
#              the constant band is wider at 4 mph and narrower at 25; the
#              likelihood ratio is anova()'s, and both it and the Wald test
#              reject at 5%; the sigma of each insecticide is the sample sd
#              of its counts and gls(varIdent)'s, with the same REML
#              criterion; the pooled sigma is lm()'s; the standard error of
#              each mean is sigma_g / sqrt(12) under the modeled sigma and
#              constant under the pooled one; C, D and E have narrower
#              intervals under the modeled sigma;
#   params     the log-variance coefficients are twice the log-sd ones and
#              the log-precision ones minus twice, with twice the standard
#              errors; the three log-likelihoods and REML criteria agree;
#              both variance fits give lm()'s variance; the identity interval
#              is symmetric; the log-scale standard error is sqrt(2/(n - p));
#              the exact chi-squared interval is closer to the log-scale one
#              than to the identity one.

assert_lm_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part II chapter 3: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  z975 <- stats::qnorm(0.975)
  w <- data$whiteside
  n <- nrow(w)
  p <- 4L

  # --- linear ---------------------------------------------------------------
  f <- fits$lin
  l <- ref$lm_lin
  if (rel(unname(coef(f)$mu), unname(stats::coef(l))) > 1e-10) {
    fail("the coefficients of the mean")
  }
  rss <- sum(stats::residuals(l)^2)
  s_reml <- sqrt(rss / (n - p)); s_ml <- sqrt(rss / n)
  sg <- exp(coef(f)$sigma[[1]])
  if (abs(sg / summary(l)$sigma - 1) > 1e-6 || abs(sg / s_reml - 1) > 1e-6 ||
      !(s_reml > s_ml)) fail("the REML sigma")
  V <- stats::vcov(f)
  if (rel(unname(V[1:4, 1:4]), unname(stats::vcov(l))) > 1e-8) {
    fail("the variance matrix of the mean")
  }
  if (any(V[5, 1:4] != 0)) fail("the zero row of sigma's intercept")
  ci <- claims$ci_lin
  zz <- ci$estimate / ci$se
  tl <- summary(l)$coefficients
  if (rel(zz, unname(tl[, "t value"])) > 1e-6) fail("z against the t value")
  tb <- summary(f)@tables$mu[[1]]$table
  if (rel(tb$p_value, 2 * stats::pnorm(-abs(zz))) > 1e-10) fail("the p-values")
  pn <- 2 * stats::pnorm(-abs(zz)); pt_ <- 2 * stats::pt(-abs(zz), n - p)
  if (any(pn[1:3] > 1e-10) || any(pt_[1:3] > 1e-10) ||
      !(pt_[4] > 1.5 * pn[4])) fail("the p-values under the two references")
  if (rel(ci$lower, ci$estimate - z975 * ci$se) > 1e-10 ||
      rel(ci$upper, ci$estimate + z975 * ci$se) > 1e-10) fail("confint()")
  mu <- stats::fitted(l)
  ll_r <- sum(stats::dnorm(w$Gas, mu, sg, log = TRUE))
  if (abs(as.numeric(logLik(f)) - ll_r) > 1e-8 ||
      abs(as.numeric(stats::logLik(l)) -
          sum(stats::dnorm(w$Gas, mu, s_ml, log = TRUE))) > 1e-8) {
    fail("the two log-likelihoods")
  }
  if (attr(logLik(f), "df") != 5 || attr(stats::logLik(l), "df") != 5 ||
      abs(AIC(f) - (-2 * as.numeric(logLik(f)) + 10)) > 1e-8) fail("AIC")
  if (abs(as.numeric(logLik(f, type = "marginal")) -
          as.numeric(stats::logLik(l, REML = TRUE))) > 1e-6) {
    fail("the REML criterion against lm()")
  }

  # --- factors --------------------------------------------------------------
  b <- coef(f)$mu
  bs <- coef(fits$sep)$mu
  if (max(abs(stats::fitted(f) - stats::fitted(fits$sep))) > 1e-8 ||
      abs(bs[["InsulAfter"]] - (b[["(Intercept)"]] + b[["InsulAfter"]])) >
      1e-8 ||
      abs(bs[["InsulAfter:Temp"]] - (b[["Temp"]] + b[["Temp:InsulAfter"]])) >
      1e-8) fail("the separate-lines fit")
  if (!(bs[["InsulAfter"]] < bs[["InsulBefore"]]) ||
      !(abs(bs[["InsulAfter:Temp"]]) < abs(bs[["InsulBefore:Temp"]]))) {
    fail("the lower and less steep line after the insulation")
  }
  ref_sum <- unname(stats::coef(ref$lm_sum))
  if (rel(unname(coef(fits$sum)$mu), ref_sum) > 1e-10 ||
      rel(unname(coef(fits$lp)$mu), ref_sum) > 1e-10) {
    fail("the sum contrasts")
  }
  if (abs(ref_sum[1] - (bs[["InsulBefore"]] + bs[["InsulAfter"]]) / 2) >
      1e-8 ||
      abs(ref_sum[3] - (bs[["InsulBefore"]] - bs[["InsulAfter"]]) / 2) > 1e-8) {
    fail("the meaning of the sum-contrast coefficients")
  }

  # --- sparse ---------------------------------------------------------------
  if (rel(unname(coef(fits$dense)$mu),
          unname(coef(fits$sparse)$mu)) > 1e-8) {
    fail("the dense and the sparse fits")
  }
  if (!is.matrix(stats::model.matrix(fits$dense)) ||
      !methods::is(stats::model.matrix(fits$sparse), "dgCMatrix") ||
      !methods::is(stats::model.matrix(fits$shops), "dgCMatrix")) {
    fail("the storage of the three designs")
  }
  sh <- data$shops
  if (!(nrow(sh) * (nlevels(sh$shop) - 1) > 1e5)) fail("the default storage")
  if (!(claims$t_sparse[["elapsed"]] < claims$t_dense[["elapsed"]])) {
    fail("the time of the sparse fit")
  }

  # --- spread: cars -----------------------------------------------------------
  c1 <- fits$c1
  g <- ref$gls_c1
  dl <- stats::coef(g$modelStruct$varStruct, unconstrained = FALSE)
  if (rel(unname(coef(c1)$mu), unname(stats::coef(g))) > 1e-4 ||
      abs(coef(c1)$sigma[[1]] - log(g$sigma)) > 1e-4 ||
      abs(coef(c1)$sigma[[2]] - dl[[1]]) > 1e-4) {
    fail("the variance regression against gls()")
  }
  if (abs(as.numeric(logLik(c1, type = "marginal")) -
          as.numeric(stats::logLik(g))) > 1e-6) fail("the REML criterion of cars")
  v4 <- data.frame(speed = c(4, 25))
  s4 <- predict(c1, what = "sigma", newdata = v4)
  gm <- coef(c1)$sigma
  if (rel(s4, exp(gm[[1]] + gm[[2]] * c(4, 25))) > 1e-10) fail("sigma(v)")
  s0 <- predict(fits$c0, what = "sigma", newdata = v4)
  if (!(s0[1] > s4[1]) || !(s0[2] < s4[2])) fail("the constant band")
  an <- stats::anova(ref$gls_c0, g)
  if (abs(claims$lr - an$L.Ratio[2]) > 1e-4 ||
      abs(stats::pchisq(claims$lr, 1, lower.tail = FALSE) -
          an[["p-value"]][2]) > 1e-5) fail("the likelihood ratio of cars")
  cg <- confint(c1, parm = "sigma:speed")
  if (!(2 * stats::pnorm(-abs(cg$estimate / cg$se)) < 0.05) ||
      !(stats::pchisq(claims$lr, 1, lower.tail = FALSE) < 0.05)) {
    fail("the rejection of a constant sigma")
  }

  # --- spread: sprays ---------------------------------------------------------
  sp <- data$sprays
  is <- datasets::InsectSprays
  sdg <- unname(tapply(is$count, is$spray, stats::sd))
  sg1 <- predict(fits$s1, what = "sigma", newdata = sp)
  if (rel(sg1, sdg) > 1e-5) fail("the sd of each insecticide")
  gs <- ref$gls_s1
  rat <- stats::coef(gs$modelStruct$varStruct, unconstrained = FALSE)
  if (rel(sg1, gs$sigma * c(1, unname(rat))) > 1e-4 ||
      abs(as.numeric(logLik(fits$s1, type = "marginal")) -
          as.numeric(stats::logLik(gs))) > 1e-6) {
    fail("the insecticides against gls(varIdent)")
  }
  pooled <- summary(stats::lm(count ~ spray, data = is))$sigma
  if (abs(exp(coef(fits$s0)$sigma[[1]]) / pooled - 1) > 1e-6) {
    fail("the pooled sigma")
  }
  p1 <- claims$p_s1; p0 <- claims$p_s0
  if (rel(p1$se, sdg / sqrt(12)) > 1e-5 || diff(range(p0$se)) > 1e-8) {
    fail("the standard errors of the means")
  }
  w1 <- p1$upper - p1$lower; w0 <- p0$upper - p0$lower
  if (!all(w1[3:5] < w0[3:5])) fail("the narrower intervals of C, D and E")

  # --- params -----------------------------------------------------------------
  a1 <- coef(fits$g1)$sigma
  if (rel(coef(fits$g2)$sigma2, 2 * a1) > 1e-5 ||
      rel(coef(fits$g3)$tau, -2 * a1) > 1e-5) {
    fail("the coefficients of the three parametrizations")
  }
  ll <- vapply(list(fits$g1, fits$g2, fits$g3),
               function(x) as.numeric(logLik(x)), numeric(1))
  lm_ <- vapply(list(fits$g1, fits$g2, fits$g3),
                function(x) as.numeric(logLik(x, type = "marginal")),
                numeric(1))
  if (diff(range(ll)) > 1e-6 || diff(range(lm_)) > 1e-6) {
    fail("the three log-likelihoods")
  }
  se1 <- confint(fits$g1, parm = "sigma:speed")$se
  if (abs(confint(fits$g2, parm = "sigma2:speed")$se / (2 * se1) - 1) > 1e-4 ||
      abs(confint(fits$g3, parm = "tau:speed")$se / (2 * se1) - 1) > 1e-4) {
    fail("the doubled standard errors")
  }
  s2 <- summary(stats::lm(dist ~ speed, data = datasets::cars))$sigma^2
  cl <- claims$ci_v_log; cid <- claims$ci_v_id
  dfc <- nrow(datasets::cars) - 2
  if (abs(exp(cl$estimate) / s2 - 1) > 1e-6 ||
      abs(cid$estimate / s2 - 1) > 1e-6) fail("the two variance fits")
  if (abs((cid$lower + cid$upper) / 2 - cid$estimate) > 1e-8) {
    fail("the symmetric identity interval")
  }
  if (abs(cl$se / sqrt(2 / dfc) - 1) > 1e-4) fail("the log-scale standard error")
  ex <- c(dfc * s2 / stats::qchisq(0.975, dfc), dfc * s2 / stats::qchisq(0.025, dfc))
  d_log <- sum(abs(exp(c(cl$lower, cl$upper)) - ex))
  d_id <- sum(abs(c(cid$lower, cid$upper) - ex))
  if (!(d_log < d_id)) fail("the exact interval")

  invisible(TRUE)
}
