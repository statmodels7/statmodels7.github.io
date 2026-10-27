# Coverage of the three prediction intervals of a new observation in a small
# regression (chapter 7, "Prediction"). One call simulates one data set of 12
# observations, fits the regression, and reports for each predictive whether
# the interval contains a new observation at x = 0.5. The bootstrap predictive
# refits the model 100 times, so the 200 data sets cost about 25 minutes in
# one process; chapter 7 reads the result from `predictive-coverage.rds` when
# the file exists and computes it otherwise.

pred_cover_one <- function(i) {
  set.seed(1000 + i)
  d <- data.frame(x = runif(n = 12))
  d$y <- 1 + 2 * d$x + rnorm(n = 12)
  y_new <- 1 + 2 * 0.5 + rnorm(n = 1)
  f <- statmod(y ~ x, distrib = gaussian1_distrib(), data = d)
  nd <- data.frame(x = 0.5)
  predictives <- list(plugin = pred_plugin(), averaged = pred_averaged(),
                      bootstrap = pred_bootstrap(n = 100))
  vapply(predictives, function(pv) {
    set.seed(i)
    p <- predict(f, what = "response", newdata = nd,
                 interval = prediction_interval(predictive = pv))
    p$lower <= y_new && y_new <= p$upper
  }, logical(1))
}
