# Certificates for Part I, chapter 3 (graphs of a fit).
#
# Same contract as every gate of the book: what the chapter states is
# re-derived here by a route the package does not itself take, and the render
# stops on any disagreement. Nothing in this file is visible to the reader.
# Each graph is read through plot_data(), which returns the numbers that the
# graph draws, and compared with the numbers that predict(), residuals() and
# the closed forms of the gaussian and the negative binomial distributions
# give.
#
# The claims, by graph:
#
#   effect       the curve is predict() of the mean at the grid minus the
#                intercept; the partial residuals are the response residuals
#                plus the effect at the observed times; the band is
#                qnorm(0.975) times the square root of x' V x for the
#                prediction matrix x of the term and the block V of vcov();
#                on the response scale the effect of a factor is predict()
#                at the first level of the other factors, the median of the
#                offset variable or the value given in `at`;
#   param        the standard deviation against time is predict() of sigma;
#   quantiles    the curve of probability p is mu + sigma qnorm(p) with mu
#                and sigma from predict(); for the negative binomial it is
#                qnbinom(p, size = theta, mu = mu);
#   density      the density at a time is dnorm() at the predicted mean and
#                standard deviation, and it integrates to one;
#   moments      the mean is mu and the standard deviation is sigma;
#   qq           the points are the sorted quantile residuals against
#                qnorm(ppoints(n)), and the band is qnorm() of the beta
#                quantiles of the order statistics;
#   worm         the points are the qq points minus the abscissa, and the
#                classes of time have equal counts;
#   residuals    the points are the quantile residuals against the covariate;
#   pit          the bars are the counts of pnorm(r) in ten bins over n times
#                the bin width, and the band is the binomial one;
#   calibration  the coverage at p is the fraction of the responses at or
#                below mu + sigma qnorm(p), and the band is
#                p +- qnorm(0.975) sqrt(p (1 - p) / n);
#   ppc          the observed density is stats::density() of the response
#                on the grid of the graph, the simulated densities number
#                n_sim and each is stats::density() of the responses that
#                simulate() draws with the same seed, the observed standard
#                deviation is sd() and the mean of the simulated ones is the
#                standard deviation that the fitted parameters imply, and
#                the standard deviation of the first 20 responses of the
#                data is under a fifth of every simulated one;
#   default      plot() composes four graphs, diagnostics() returns the four
#                named graphs ppc, worm, calibration and residuals, and a
#                graph that reads a covariate signals an error without
#                `data`;
#   table        every exported graph type is in the table of the last
#                section, and the table names nothing else.

assert_graphs_ok <- function(f0, f1, fi, sp, data, ins, species, chapter) {

  fail <- function(what) {
    stop("Part I chapter 3: ", what, " no longer agrees with the package.",
         call. = FALSE)
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
  z <- stats::qnorm(0.975)
  ctx <- function(param = NULL, name = NULL, vars = NULL, d = data) {
    list(param = param, name = name, vars = vars, cols = NULL, data = d)
  }
  nm <- "s(times, bspline_smooth(k = 20))"
  tm <- f0@spec@terms$mu[[nm]]

  # --- effect of a smooth ---------------------------------------------------
  d <- plot_data(plot_effect(), tm, f0, ctx = ctx("mu", nm, "times"))
  b0 <- coef(f0)$mu[["(Intercept)"]]
  grid <- data.frame(times = d$curve$x)
  eff <- predict(f0, what = "mu", newdata = grid) - b0
  if (rel(d$curve$fit, eff) > 1e-8) fail("the curve of the effect")
  obs <- data.frame(times = d$points$x)
  part <- residuals(f0, residual_type = "response") +
    predict(f0, what = "mu", newdata = obs) - b0
  if (rel(d$points$y, part) > 1e-8) fail("the partial residuals")
  des <- statmod_design(f0@spec)
  idx <- des$mu$blocks[[nm]]
  key <- paste0("mu:", des$mu$coef_names[idx])
  Vb <- vcov(f0)[key, key]
  Xg <- modelterms7::term_predict(tm, newdata = grid)
  half <- z * sqrt(rowSums((Xg %*% Vb) * Xg))
  if (rel(d$curve$upper - d$curve$fit, half) > 1e-8 ||
      rel(d$curve$fit - d$curve$lower, half) > 1e-8) {
    fail("the band of the effect")
  }

  # --- effect on the response scale, with an offset ---------------------------
  first <- function(x) {
    factor(levels(x)[1L], levels = levels(x), ordered = is.ordered(x))
  }
  p_i <- ggplot2::autoplot(
    fi, term = "Age", data = ins,
    plot_type = plot_effect(scale = "response", at = list(Holders = 100)))
  ci <- p_i$layers[[1L]]$data
  lev <- levels(ins$Age)
  new_i <- data.frame(District = rep(first(ins$District), length(lev)),
                      Group = rep(first(ins$Group), length(lev)),
                      Age = factor(lev, levels = lev, ordered = TRUE),
                      Holders = 100)
  pr <- predict(fi, what = "mu", newdata = new_i, se = TRUE)
  if (rel(ci$fit, pr$fit) > 1e-8 || rel(ci$lower, pr$lower) > 1e-8 ||
      rel(ci$upper, pr$upper) > 1e-8) {
    fail("the effect on the response scale")
  }
  med <- ggplot2::autoplot(fi, term = "Age", data = ins,
                           plot_type = plot_effect(scale = "response"))
  new_m <- new_i
  new_m$Holders <- stats::median(ins$Holders)
  if (rel(med$layers[[1L]]$data$fit,
          predict(fi, what = "mu", newdata = new_m)) > 1e-8) {
    fail("the median of the offset variable")
  }

  # --- graphs of the fitted distribution --------------------------------------
  gr <- data.frame(times = seq(min(data$times), max(data$times),
                               length.out = 200))
  th <- predict(f1, what = "parameter", newdata = gr)
  dp <- plot_data(plot_param(param = "sigma", over = "times"),
                  f1@spec@distrib, f1, ctx = ctx())
  sg <- predict(f1, what = "sigma", newdata = gr, se = TRUE)
  if (rel(dp$curve$fit, th$sigma) > 1e-8 ||
      rel(dp$curve$lower, sg$lower) > 1e-8 ||
      rel(dp$curve$upper, sg$upper) > 1e-8) {
    fail("the standard deviation against time")
  }
  dq <- plot_data(plot_quantiles(over = "times"), f1@spec@distrib, f1,
                  ctx = ctx())
  cv <- dq$curves
  near <- function(x, g) vapply(x, function(v) which.min(abs(g - v)), 1L)
  j <- near(cv$x, gr$times)
  ref <- th$mu[j] + th$sigma[j] * stats::qnorm(cv$prob)
  if (rel(cv$q, ref) > 1e-8) fail("the quantile curves")
  if (!identical(sort(unique(cv$prob)), c(0.05, 0.25, 0.5, 0.75, 0.95))) {
    fail("the default probabilities of the quantile curves")
  }
  at <- data.frame(times = c(10, 20, 30, 45))
  dd <- plot_data(plot_density(at = at), f1@spec@distrib, f1, ctx = ctx())
  tha <- predict(f1, what = "parameter", newdata = at)
  for (i in seq_len(nrow(at))) {
    r <- dd$dens[dd$dens$index == i, ]
    if (rel(r$dens, stats::dnorm(r$y, tha$mu[i], tha$sigma[i])) > 1e-8) {
      fail("the fitted density")
    }
    area <- sum(diff(r$y) * (r$dens[-1] + r$dens[-nrow(r)]) / 2)
    if (abs(area - 1) > 5e-3) fail("the area of the fitted density")
  }
  dm <- plot_data(plot_moments(over = "times",
                               moments = c("mean", "std_dev")),
                  f1@spec@distrib, f1, ctx = ctx())
  mn <- dm$curve[dm$curve$moment == "mean", ]
  sd_ <- dm$curve[dm$curve$moment == "std_dev", ]
  if (rel(mn$value, th$mu) > 1e-8 || rel(sd_$value, th$sigma) > 1e-8) {
    fail("the moments")
  }

  # --- a discrete distribution ------------------------------------------------
  sg_ <- seq(min(species$log_lake), max(species$log_lake), length.out = 200)
  ths <- predict(sp, what = "parameter", newdata = data.frame(log_lake = sg_))
  dsq <- plot_data(plot_quantiles(over = "log_lake"), sp@spec@distrib, sp,
                   ctx = ctx(d = species))
  cs <- dsq$curves
  i_s <- near(cs$x, sg_)
  if (max(abs(cs$q - stats::qnbinom(cs$prob, size = ths$theta[i_s],
                                    mu = ths$mu[i_s]))) > 1e-8) {
    fail("the quantiles of the negative binomial distribution")
  }
  ats <- data.frame(log_lake = c(2, 5, 8))
  dds <- plot_data(plot_density(at = ats), sp@spec@distrib, sp,
                   ctx = ctx(d = species))
  tsa <- predict(sp, what = "parameter", newdata = ats)
  for (i in seq_len(nrow(ats))) {
    r <- dds$dens[dds$dens$index == i, ]
    if (max(abs(r$dens - stats::dnbinom(r$y, size = tsa$theta[i],
                                        mu = tsa$mu[i]))) > 1e-8) {
      fail("the fitted probabilities of the negative binomial distribution")
    }
  }

  # --- diagnostics ------------------------------------------------------------
  n <- nrow(data)
  r0 <- residuals(f0)
  dq0 <- plot_data(plot_qq(), f0@spec@distrib, f0, ctx = ctx())
  qx <- stats::qnorm(stats::ppoints(n))
  if (rel(dq0$points$x, qx) > 1e-8 || rel(dq0$points$y, sort(r0)) > 1e-8) {
    fail("the points of the quantile-quantile graph")
  }
  i <- seq_len(n)
  lo <- stats::qnorm(stats::qbeta(0.025, i, n + 1 - i))
  hi <- stats::qnorm(stats::qbeta(0.975, i, n + 1 - i))
  if (rel(dq0$band$lower, lo) > 1e-8 || rel(dq0$band$upper, hi) > 1e-8) {
    fail("the band of the quantile-quantile graph")
  }
  dw <- plot_data(plot_worm(), f0@spec@distrib, f0, ctx = ctx())
  if (rel(dw$points$y, sort(r0) - qx) > 1e-8) fail("the worm plot")
  dwo <- plot_data(plot_worm(over = "times"), f0@spec@distrib, f0,
                   ctx = ctx())
  cnt <- as.integer(table(dwo$points$class))
  if (length(cnt) != 4L || max(cnt) - min(cnt) > 1L || sum(cnt) != n) {
    fail("the classes of the worm plot")
  }
  cl <- split(r0, cut(data$times,
                      stats::quantile(data$times, 0:4 / 4),
                      include.lowest = TRUE))
  ref_w <- unlist(lapply(cl, function(r) {
    sort(r) - stats::qnorm(stats::ppoints(length(r)))
  }), use.names = FALSE)
  if (max(abs(sort(dwo$points$y) - sort(ref_w))) > 1e-8) {
    fail("the worm plot within classes")
  }
  dr <- plot_data(plot_residuals(against = "times"), f0@spec@distrib, f0,
                  ctx = ctx())
  if (rel(dr$points$y, r0) > 1e-8 || rel(dr$points$x, data$times) > 1e-8) {
    fail("the residuals against time")
  }
  r1 <- residuals(f1)
  dpi <- plot_data(plot_pit(bins = 10), f1@spec@distrib, f1, ctx = ctx())
  cnt_p <- tabulate(findInterval(stats::pnorm(r1), seq(0, 1, 0.1),
                                 rightmost.closed = TRUE), nbins = 10)
  if (rel(dpi$bars$density, cnt_p / (n * 0.1)) > 1e-8) {
    fail("the histogram of the probability integral transform")
  }
  dc <- plot_data(plot_calibration(), f1@spec@distrib, f1, ctx = ctx())
  th1 <- predict(f1, what = "parameter", newdata = data)
  cov_p <- vapply(dc$points$p, function(p) {
    mean(data$accel <= th1$mu + th1$sigma * stats::qnorm(p))
  }, numeric(1))
  if (rel(dc$points$observed, cov_p) > 1e-8) fail("the observed coverage")
  band_c <- z * sqrt(dc$points$p * (1 - dc$points$p) / n)
  if (rel(dc$points$upper, pmin(dc$points$p + band_c, 1)) > 1e-8 ||
      rel(dc$points$lower, pmax(dc$points$p - band_c, 0)) > 1e-8) {
    fail("the band of the calibration")
  }

  # --- simulated responses ----------------------------------------------------
  y <- data$accel
  dpp <- plot_data(plot_ppc(n_sim = 20, seed = 1), f0@spec@distrib, f0,
                   ctx = ctx())
  rg <- range(y)
  ref_d <- stats::density(y, from = rg[1L] - 0.1 * diff(rg),
                          to = rg[2L] + 0.1 * diff(rg), n = 256)
  if (rel(dpp$observed$y, ref_d$y) > 1e-8 ||
      rel(dpp$observed$x, ref_d$x) > 1e-8) {
    fail("the observed density of the response")
  }
  if (length(unique(dpp$sim$sim)) != 20L) fail("the number of simulated densities")
  sims <- simulate(f0, nsim = 20, seed = 1)
  for (k in 1:20) {
    dk <- stats::density(sims[[k]], from = rg[1L] - 0.1 * diff(rg),
                         to = rg[2L] + 0.1 * diff(rg), n = 256)
    if (rel(dpp$sim$y[dpp$sim$sim == k], dk$y) > 1e-8) {
      fail("the simulated densities")
    }
  }
  dps <- plot_data(plot_ppc(statistic = stats::sd, n_sim = 200, seed = 1),
                   f0@spec@distrib, f0, ctx = ctx())
  th0 <- predict(f0, what = "parameter", newdata = data)
  sd_model <- sqrt(mean(th0$sigma^2 + (th0$mu - mean(th0$mu))^2))
  if (abs(dps$observed$value - stats::sd(y)) > 1e-10 ||
      abs(mean(dps$sim$value) / sd_model - 1) > 0.03) {
    fail("the statistic of the simulated responses")
  }
  early <- function(v) stats::sd(v[1:20])
  dpe <- plot_data(plot_ppc(statistic = early, n_sim = 200, seed = 1),
                   f0@spec@distrib, f0, ctx = ctx())
  sims_e <- simulate(f0, nsim = 200, seed = 1)
  ref_e <- vapply(sims_e, early, numeric(1))
  if (rel(dpe$sim$value, ref_e) > 1e-8 ||
      abs(dpe$observed$value - early(y)) > 1e-10 ||
      !(dpe$observed$value < 0.2 * min(dpe$sim$value))) {
    fail("the statistic of the first responses")
  }

  # --- the group of default graphs --------------------------------------------
  if (!inherits(ggplot2::autoplot(f0), "patchwork")) fail("plot()")
  if (!identical(names(diagnostics(f1)),
                 c("ppc", "worm", "calibration", "residuals"))) {
    fail("diagnostics()")
  }
  msg <- tryCatch(ggplot2::autoplot(f0, term = "times"),
                  error = function(e) conditionMessage(e))
  if (!is.character(msg) || !grepl("does not keep its data", msg)) {
    fail("the error of a graph without data")
  }

  # --- the table of the graphs ------------------------------------------------
  ex <- getNamespaceExports("statmodels7")
  types <- setdiff(grep("^plot_", ex, value = TRUE),
                   c("plot_data", "plot_layers", "plot_type"))
  rows <- grep("^\\| `plot_", chapter, value = TRUE)
  listed <- sub("^\\| `(plot_[a-z_]+)\\(\\)`.*$", "\\1", rows)
  if (!setequal(types, listed) || anyDuplicated(listed)) {
    fail("the table of the graphs")
  }
  invisible(TRUE)
}
