# Certificates for Part III, chapter 10 (penalized regression).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
#
# The claims, section by section:
#
#   basic    the ridge is mgcv's paraPen fit (lambda = sp / sig2, fitted values,
#            REML criterion); the lasso is glmnet's at
#            lambda_g = lambda sigma^2 / sqrt(n (n - 1)); the lasso and the
#            elastic net satisfy their KKT conditions with the penalty's
#            derivative written out here; the lasso keeps Po1 and drops Po2;
#            the matrix input is the formula input; standardized fits do not
#            move when two covariates change units, raw fits do;
#   scad     SCAD is ncvreg's ncvfit at the same conversion; ncvreg's MCP
#            reaches our estimates along a path and from zero; MCP's score
#            is zero beyond gamma lambda / c and equals the derivative of
#            the penalty inside, the derivative written out here; the a grid
#            chooses 3.7;
#   hyper    the path is decreasing and its first point is empty;
#            between the two points of the jump lambda sigma^2 falls further
#            than lambda; a held value has source fixed, a written grid is
#            visited from the largest down and its minimum is chosen; the
#            cyclic search of the elastic net scores fewer pairs;
#   sigma    the lasso of log sigma satisfies its KKT conditions with the
#            score written out here, keeps x3 and x4 only, and the lasso of
#            the mean keeps x1 and x2 only; the constant-sigma model keeps
#            more; nu's lasso satisfies its KKT conditions, the score taken
#            by numDeriv on a sum of stats::dt;
#   mixed    the Boston and Cars93 lasso blocks are glmnet's with the rest as
#            an offset; beside the random intercept the lasso keeps fewer
#            covariates, the ones the chapter names dropping;
#   id       members of an id carry one value; two ridges with one id are the
#            ridge on the union (fitted values, criterion); the likelihood
#            ratio is recomputed; a lasso id visits the single lasso's path
#            and chooses its lambda; across equations the shared path starts
#            at the larger of the two tops;
#   choose   BIC is -2 l + log(n) tau recomputed; the one-se value is read
#            from the path; the count of values within one standard error.

assert_penalized_ok <- function(fits, ref, data, claims) {

  fail <- function(what) {
    stop("Part III chapter 10: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  # an empty or mismatched comparison is a failure, never a vacuous pass
  rel <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (!length(b) || length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
    max(abs(a - b)) / max(1, max(abs(b)))
  }
  gap <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (!length(b) || length(a) != length(b) || anyNA(a) || anyNA(b)) return(Inf)
    max(abs(a - b))
  }
  for (pk in c("glmnet", "ncvreg", "mgcv", "numDeriv")) {
    if (!requireNamespace(pk, quietly = TRUE)) fail(paste("the package", pk))
  }
  hyp <- function(fit, nm, term = NULL) {
    h <- hyper(fit)
    if (!is.null(term)) h <- h[grepl(term, h$term, fixed = TRUE), ]
    h$estimate[h$name == nm]
  }
  nz <- function(b) b != 0

  cr <- data$UScrime
  X <- data$X
  n <- nrow(cr)
  s <- apply(X, 2, stats::sd)
  sig <- function(fit) exp(fit@coefficients$sigma[[1]])

  # --- basic -----------------------------------------------------------------
  gr <- ref$g_ridge
  if (rel(hyp(fits$ridge, "lambda"), gr$sp / gr$sig2) > 1e-3 ||
      rel(fitted(fits$ridge), stats::fitted(gr)) > 1e-5 ||
      abs(fits$ridge@criterion - (-gr$gcv.ubre)) > 1e-3) {
    fail("the ridge against mgcv")
  }
  if (rel(fits$lasso@coefficients$mu, as.numeric(stats::coef(ref$g_lasso))) > 1e-6) {
    fail("the lasso against glmnet")
  }
  # KKT: the derivative of the log-likelihood in beta_j, at sigma's estimate,
  # equals the derivative of the penalty where beta_j is not zero and is at
  # most lambda alpha s_j in absolute value where it is
  score <- function(fit) {
    r <- cr$lrate - X %*% fit@coefficients$mu[-1] - fit@coefficients$mu[[1]]
    as.numeric(crossprod(X, r)) / sig(fit)^2
  }
  kkt <- function(fit, lam, a) {
    b <- fit@coefficients$mu[-1]; g <- score(fit); on <- nz(b)
    dq <- lam * (1 - a) * s^2 * b
    all(abs(g[on] - dq[on] - lam * a * s[on] * sign(b[on])) <= 1e-4 * lam * s[on]) &&
      all(abs(g[!on]) <= lam * a * s[!on] * (1 + 1e-4))
  }
  if (!kkt(fits$lasso, hyp(fits$lasso, "lambda"), 1)) fail("the KKT conditions of the lasso")
  if (!kkt(fits$enet, hyp(fits$enet, "lambda"), hyp(fits$enet, "alpha"))) {
    fail("the KKT conditions of the elastic net")
  }
  bl <- stats::setNames(fits$lasso@coefficients$mu[-1], colnames(X))
  if (!(bl[["Po1"]] != 0 && bl[["Po2"]] == 0)) fail("the lasso's Po1 and Po2")
  if (!identical(claims$keep_l, names(bl)[nz(bl)])) fail("the covariates of the lasso")
  ne <- sum(nz(fits$enet@coefficients$mu[-1]))
  if (!identical(claims$same_e,
                 identical(nz(fits$enet@coefficients$mu[-1]), nz(bl))) ||
      !(fits$enet@edf$edf[2] < ne) || !(hyp(fits$enet, "alpha") < 0.5)) {
    fail("the elastic net's covariates and degrees of freedom")
  }
  if (abs(fits$lasso@edf$edf[2] - sum(nz(bl))) > 1e-8) {
    fail("the degrees of freedom of the lasso")
  }
  if (rel(fits$mat@coefficients$mu, fits$lasso@coefficients$mu) > 1e-8) {
    fail("the matrix input")
  }
  if (rel(fitted(fits$s1), fitted(fits$s2)) > 1e-6) fail("the standardized units")
  if (!(gap(fitted(fits$r1), fitted(fits$r2)) > 1e-3) ||
      !is.finite(gap(fitted(fits$r1), fitted(fits$r2)))) fail("the raw units")

  # --- SCAD and MCP ----------------------------------------------------------
  s_n <- s * sqrt((n - 1) / n)
  b_s <- fits$scad@coefficients$mu[-1]
  if (gap(ref$ncv_scad$beta / s_n, b_s) > 1e-5) fail("SCAD against ncvreg")
  v_m <- fits$mcp@coefficients$mu[-1] * s_n
  np <- ref$ncv_path$beta
  if (gap(np[-1, ncol(np)], v_m) > 1e-5) fail("MCP against ncvreg's path")
  if (gap(ref$ncv_zero$beta, v_m) > 1e-5 ||
      abs(claims$o_zero - claims$o_ours) > 1e-10) {
    fail("MCP against ncvreg from zero")
  }
  # MCP on the scale u = s beta: c = (n - 1) / sigma^2, flat beyond
  # gamma lambda / c, derivative lambda - c |u| / gamma inside
  bm <- fits$mcp@coefficients$mu
  u <- bm[-1] * s
  lm_ <- hyp(fits$mcp, "lambda"); cc <- (n - 1) / sig(fits$mcp)^2
  r <- cr$lrate - bm[[1]] - drop(X %*% bm[-1])
  sc <- as.numeric(crossprod(X, r)) / sig(fits$mcp)^2 / s
  on_m <- nz(u)
  flat <- on_m & abs(u) > 3 * lm_ / cc
  inside <- on_m & !flat
  if (!any(flat) || !any(inside) ||
      max(abs(sc[flat])) > 1e-4 * lm_ ||
      max(abs(sc[inside] - (lm_ - cc * abs(u[inside]) / 3) * sign(u[inside]))) >
        1e-4 * lm_ ||
      !all(abs(sc[!on_m]) <= lm_ * (1 + 1e-4))) {
    fail("the score of MCP against its penalty")
  }
  if (!identical(claims$keep_s, colnames(X)[nz(b_s)])) fail("the covariates of SCAD")
  # the chapter: four MCP coefficients beyond the point, six kept, edf above
  # the count; SCAD has none where its penalty bends and edf equal to the
  # count; Prob has a larger standard error under MCP
  if (sum(flat) != 4 || sum(on_m) != 6 ||
      abs(claims$edf_m - fits$mcp@edf$edf[2]) > 1e-12 || !(claims$edf_m > 6)) {
    fail("the effective degrees of freedom of MCP")
  }
  us <- b_s * s
  cs <- (n - 1) / sig(fits$scad)^2
  ls_ <- hyp(fits$scad, "lambda")
  bend <- nz(us) & abs(us) > ls_ / cs & abs(us) < 3.7 * ls_ / cs
  if (any(bend) || abs(claims$edf_s - sum(nz(b_s))) > 1e-8) {
    fail("the effective degrees of freedom of SCAD")
  }
  se_prob <- function(f) {
    v <- sqrt(diag(vcov(f)))
    v[grep("Prob", names(v))][1]
  }
  if (!isTRUE(se_prob(fits$mcp) > se_prob(fits$scad))) fail("the standard error of Prob")
  four <- match(c("Po1", "Ed", "Ineq", "M"), colnames(X))
  if (max(abs(bm[-1][four] / b_s[four] - 1)) > 0.15) {
    fail("the four coefficients common to SCAD and MCP")
  }
  if (hyp(fits$scad_a, "a") != 3.7) fail("the a chosen by the grid")

  # --- the hyperparameters ---------------------------------------------------
  po <- fits$lasso@history$outer
  if (is.unsorted(rev(po$value))) fail("the order of the path")
  np5 <- claims$n_path
  if (np5[1] != 0 || sum(nz(fits$lasso@coefficients$mu[-1])) == 0) {
    fail("the empty start of the path")
  }
  ij <- claims$i_jump
  le <- claims$lam_eff
  if (!(np5[ij + 1] - np5[ij] >= 2) ||
      !(le[ij] / le[ij + 1] > po$value[ij] / po$value[ij + 1])) {
    fail("the jump of the path")
  }
  hh <- hyper(fits$held)
  if (!(hh$source[1] == "fixed" && hh$estimate[1] == 40)) fail("the held lambda")
  pg <- fits$grid@history$outer
  if (!identical(as.numeric(pg$value), c(80, 40, 20, 10, 5)) ||
      hyp(fits$grid, "lambda") != pg$value[which.min(pg$criterion)]) {
    fail("the written grid")
  }
  p100 <- fits$lasso_100@history$outer
  if (!identical(as.numeric(claims$lams_100), as.numeric(p100$value)) ||
      hyp(fits$lasso_100, "lambda") != claims$lam_100 ||
      claims$lam_100 != p100$value[which.min(p100$criterion)] ||
      !(BIC(fits$lasso_100) < BIC(fits$lasso)) ||
      sum(nz(fits$lasso_100@coefficients$mu[-1])) != sum(nz(bl))) {
    fail("the path of 100 values")
  }
  if (!(nrow(fits$enet_c@history$outer) < nrow(fits$enet@history$outer))) {
    fail("the cyclic search")
  }

  # --- sigma and nu ----------------------------------------------------------
  sm <- data$sim
  Xs <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, sm)[, -1]
  f1 <- fits$sig1
  bm <- f1@coefficients$mu; bs <- f1@coefficients$sigma
  mu <- bm[[1]] + drop(Xs %*% bm[-1]); ls <- bs[[1]] + drop(Xs %*% bs[-1])
  z2 <- (sm$y - mu)^2 * exp(-2 * ls)
  g_s <- as.numeric(crossprod(Xs, z2 - 1))
  lam_s <- hyp(f1, "lambda", "sigma")
  if (!length(lam_s)) lam_s <- hyper(f1)$estimate[hyper(f1)$parameter == "sigma"]
  on <- nz(bs[-1])
  if (!(all(abs(g_s[on] - lam_s * sign(bs[-1][on])) <= 1e-3 * lam_s) &&
        all(abs(g_s[!on]) <= lam_s * (1 + 1e-3)))) {
    fail("the KKT conditions of the lasso of log sigma")
  }
  if (!identical(colnames(Xs)[on], c("x3", "x4")) ||
      !identical(colnames(Xs)[nz(bm[-1])], c("x1", "x2"))) {
    fail("the covariates of the model with sigma")
  }
  if (!(sum(nz(fits$sig0@coefficients$mu[-1])) > 2)) fail("the constant-sigma model")
  fn <- fits$nu
  cn <- fn@coefficients
  ll_nu <- function(g) {
    m <- cn$mu[[1]] + cn$mu[[2]] * sm$x1
    nu <- exp(g[1] + drop(Xs %*% g[-1]))
    sc <- exp(cn$sigma[[1]])
    sum(stats::dt((sm$yt - m) / sc, df = nu, log = TRUE) - log(sc))
  }
  g_nu <- numDeriv::grad(ll_nu, cn$nu)
  lam_nu <- hyper(fn)$estimate[hyper(fn)$parameter == "nu"]
  bn <- cn$nu[-1]; onn <- nz(bn)
  if (!(abs(g_nu[1]) < 1e-3 &&
        all(abs(g_nu[-1][onn] - lam_nu * sign(bn[onn])) <= 1e-3 * max(1, lam_nu)) &&
        all(abs(g_nu[-1][!onn]) <= lam_nu * (1 + 1e-3)))) {
    fail("the KKT conditions of the lasso of log nu")
  }

  # --- beside a smooth -------------------------------------------------------
  bb <- coef(fits$boston)$mu
  il <- startsWith(names(bb), "lasso.")
  if (!any(il)) fail("the lasso block of the Boston fit")
  if (gap(as.numeric(stats::coef(ref$g_b))[-1], bb[il]) > 1e-5) {
    fail("the Boston lasso against glmnet")
  }

  # --- beside a random effect ----------------------------------------------
  ca <- data$Cars93
  nmc <- paste0("lasso.", all.vars(data$cx))
  b_c <- coef(fits$cars)$mu[nmc]
  b_c0 <- coef(fits$cars0)$mu[nmc]
  if (gap(as.numeric(stats::coef(ref$g_c))[-1], b_c) > 1e-5) {
    fail("the Cars93 lasso against glmnet")
  }
  if (!identical(claims$drop_c, all.vars(data$cx)[b_c == 0 & b_c0 != 0]) ||
      !length(claims$drop_c) || !(sum(b_c != 0) < sum(b_c0 != 0))) {
    fail("the covariates dropped beside the random intercept")
  }
  hc <- hyper(fits$cars)
  if (!(hc$source[2] == "reml" && hc$estimate[2] > 0.05 && hc$source[1] == "bic")) {
    fail("the hyperparameters of the Cars93 model")
  }

  # --- shared hyperparameters ------------------------------------------------
  hi <- hyper(fits$id)$estimate
  if (abs(hi[1] - hi[2]) > 1e-10) fail("the shared ridge lambda")
  if (rel(fitted(fits$id), fitted(fits$ridge)) > 1e-5 ||
      abs(fits$id@criterion - fits$ridge@criterion) > 1e-4) {
    fail("two ridges with one id against one ridge")
  }
  lr <- 2 * (as.numeric(logLik(fits$two, type = "marginal")) -
               as.numeric(logLik(fits$id, type = "marginal")))
  if (abs(lr - claims$lr_id) > 1e-8 || !(lr >= 0)) fail("the likelihood ratio of the id")
  hl <- hyper(fits$lid)$estimate
  if (abs(hl[1] - hl[2]) > 1e-10 ||
      rel(hl[1], hyp(fits$lasso, "lambda")) > 1e-8 ||
      rel(fits$lid@history$outer$value, po$value) > 1e-8) {
    fail("the lasso with one id")
  }
  hs <- hyper(fits$sid)$estimate
  if (abs(hs[1] - hs[2]) > 1e-10 || !(BIC(fits$sid) > BIC(fits$sig1))) {
    fail("the id across equations")
  }

  # --- choosing --------------------------------------------------------------
  b_hand <- -2 * as.numeric(logLik(fits$lasso)) + log(n) * sum(fits$lasso@edf$edf)
  if (abs(b_hand - BIC(fits$lasso)) > 1e-6) fail("the BIC of the lasso")
  pc <- claims$path_c
  ok <- is.finite(pc$criterion)
  im <- which.min(ifelse(ok, pc$criterion, Inf))
  if (sum(ok & pc$criterion <= pc$criterion[im] + pc$se[im]) != claims$n_within) {
    fail("the values within one standard error")
  }
  if (hyp(fits$cv, "lambda") != pc$value[im]) fail("the value chosen by cv")
  within <- ok & pc$criterion <= pc$criterion[im] + pc$se[im]
  # since statmodels7 0.196.0 (the folds reapply the full-data terms) the
  # minimum is the third value of the path and the rule picks the second
  if (im != 3L || max(pc$value[within]) != pc$value[2L] ||
      claims$lam_1se != pc$value[2L]) {
    fail("the one-se rule")
  }
  invisible(TRUE)
}
