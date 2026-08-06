# Consistency gates for chapter 8. Each printed formula is transcribed here
# INDEPENDENTLY of the package -- the constants written out, SCAD's rho'
# typed from the displayed equation, the pseudo-determinant recomputed by
# eigen -- so a gate that compares the package with itself cannot happen.
# Injection-checked at development time: a 5% corruption of the ridge
# constant, of the SCAD slope and of the structured theta gradient are each
# caught, and the clean versions pass.

.certify_penalty_definition <- function() {
  out <- character()

  # the ridge value as printed: (1/2s^2) b'b + q log s + (q/2) log 2pi
  q <- 3
  beta <- c(0.8, -0.4, 1.1)
  s <- 1.5
  pen <- ridge_penalty(n_coef = q)
  hand <- sum(beta^2) / (2 * s^2) + q * log(s) + q / 2 * log(2 * pi)
  if (abs(penalty_value(pen, beta, list(sigma = s)) - hand) > 1e-12) {
    out <- c(out, "the printed ridge value disagrees with the package")
  }

  # properness: the value at two scales differs by the constant's difference
  # plus the quadratic part's, both written out here
  s2 <- 2.4
  d_hand <- sum(beta^2) / 2 * (1 / s2^2 - 1 / s^2) + q * log(s2 / s)
  d_pkg <- penalty_value(pen, beta, list(sigma = s2)) -
    penalty_value(pen, beta, list(sigma = s))
  if (abs(d_pkg - d_hand) > 1e-12) {
    out <- c(out, "the normalizing constant does not move as the density's")
  }

  # gradient and mixed block against numDeriv, which shares no code
  g <- penalty_gradient(pen, beta, list(sigma = s))
  gn <- numDeriv::grad(function(b) penalty_value(pen, b, list(sigma = s)), beta)
  if (max(abs(g - gn)) > 1e-8) out <- c(out, "ridge gradient off numDeriv")
  cr <- penalty_cross(pen, beta, list(sigma = s))$sigma
  crn <- as.numeric(numDeriv::jacobian(function(v) {
    penalty_gradient(pen, beta, list(sigma = v))
  }, s))
  if (max(abs(cr - crn)) > 1e-7) out <- c(out, "ridge mixed block off numDeriv")

  out
}

.certify_penalty_branches <- function() {
  out <- character()

  # eq-penalty-quadratic transcribed, with the pseudo-determinant recomputed
  # by eigen here rather than read from the object
  P <- crossprod(diff(diag(6), differences = 2))
  pen <- quadratic_penalty(P)
  beta <- c(0.4, -0.2, 0.9, 0.1, -0.5, 0.3)
  lam <- 1.7
  ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  keep <- ev > 1e-10 * max(ev)
  r <- sum(keep)
  hand <- lam / 2 * sum(beta * as.numeric(P %*% beta)) - r / 2 * log(lam) +
    r / 2 * log(2 * pi) - sum(log(ev[keep])) / 2
  if (abs(penalty_value(pen, beta, list(lambda = lam)) - hand) > 1e-10) {
    out <- c(out, "the printed quadratic value disagrees with the package")
  }
  if (penalty_rank(pen) != r) out <- c(out, "the quadratic rank is not eigen's")

  # the ridge twin: quadratic at lambda = 1/sigma^2 vs the separable
  # gaussian, componentwise, machine precision
  q <- 3
  b2 <- c(0.4, -1.1, 2.2)
  sg <- 1.7
  a <- quadratic_penalty(diag(q))
  b <- ridge_penalty(n_coef = q)
  if (abs(penalty_value(a, b2, list(lambda = 1 / sg^2)) -
          penalty_value(b, b2, list(sigma = sg))) > 1e-12 ||
      max(abs(penalty_gradient(a, b2, list(lambda = 1 / sg^2)) -
              penalty_gradient(b, b2, list(sigma = sg)))) > 1e-12) {
    out <- c(out, "the two constructions of ridge disagree")
  }

  # SCAD's rho' typed from eq-penalty-scad, on points inside each region
  pen <- scad_penalty(n_coef = 1)
  lam <- 1; aa <- 3.7
  th <- list(lambda = lam, a = aa)
  for (t in c(0.3, 0.8, 1.5, 3.0, 4.5, 8.0)) {
    hand <- if (t <= lam) lam else if (t <= aa * lam) {
      (aa * lam - t) / (aa - 1)
    } else 0
    if (abs(penalty_gradient(pen, t, th) - hand) > 1e-12) {
      out <- c(out, sprintf("SCAD rho'(%g) is not the printed one", t))
    }
  }
  if (abs(penalty_value(pen, 0, th)) > 1e-14) {
    out <- c(out, "SCAD is not anchored at rho(0) = 0")
  }
  if (abs(penalty_value(pen, 5, th) - penalty_value(pen, 50, th)) > 1e-14) {
    out <- c(out, "SCAD is not flat beyond the shoulder")
  }
  # MCP's rho' typed from the same display
  pen <- mcp_penalty(n_coef = 1)
  th <- list(lambda = 1, gamma = 3)
  for (t in c(0.4, 1.2, 2.4, 5.0)) {
    hand <- max(1 - t / 3, 0)
    if (abs(penalty_gradient(pen, t, th) - hand) > 1e-12) {
      out <- c(out, sprintf("MCP rho'(%g) is not the printed one", t))
    }
  }

  out
}

.certify_penalty_marginal <- function() {
  out <- character()

  # eq-penalty-logpdet: linear in log lambda with slope exactly r
  P <- crossprod(diff(diag(6), differences = 2))
  pen <- quadratic_penalty(P)
  r <- penalty_rank(pen)
  l1 <- penalty_logpdet(pen, list(lambda = 0.7))$value
  l2 <- penalty_logpdet(pen, list(lambda = 0.7 * exp(2)))$value
  if (abs((l2 - l1) - 2 * r) > 1e-10) {
    out <- c(out, "the log pseudo-determinant is not linear with slope r")
  }
  nb <- penalty_null_basis(pen)
  M <- penalty_matrix(pen, list(lambda = 0.7))
  if (max(abs(M %*% nb)) > 1e-10 * max(abs(M))) {
    out <- c(out, "the null basis does not annihilate the matrix")
  }
  # the refusal: a separable penalty must not answer the marginal pieces
  ref <- tryCatch({
    penalty_matrix(ridge_penalty(n_coef = 2), list(sigma = 1))
    FALSE
  }, error = function(e) TRUE)
  if (!ref) out <- c(out, "a separable penalty answered penalty_matrix")

  # eq-penalty-structured-grad transcribed against parameters7 directly,
  # and independently against numDeriv on the value
  s <- parameters7::ar1(4, role = "precision")
  pen <- structured_penalty(s)
  theta <- list(log_scale = 0.2, z_rho = 0.5)
  beta <- c(0.3, -0.1, 0.4, 0.2)
  eta <- c(0.2, 0.5)
  A <- parameters7::param_d1(s, eta)
  dld <- parameters7::param_dlogdet(s, eta)
  gt <- penalty_grad_theta(pen, beta, theta)
  for (k in 1:2) {
    hand <- sum(beta * as.numeric(A[[k]] %*% beta)) / 2 - dld[[k]] / 2
    if (abs(gt[[k]] - hand) > 1e-12) {
      out <- c(out, "the printed structured theta gradient disagrees")
    }
    nd <- numDeriv::grad(function(v) {
      th <- theta; th[[k]] <- v
      penalty_value(pen, beta, th)
    }, theta[[k]])
    if (abs(gt[[k]] - nd) > 1e-7) {
      out <- c(out, "the structured theta gradient is off numDeriv")
    }
  }
  cr <- penalty_cross(pen, beta, theta)
  if (max(abs(cr[[1]] - as.numeric(A[[1]] %*% beta))) > 1e-12) {
    out <- c(out, "the structured mixed block is not A_k beta")
  }

  # the collapse: log-Cholesky at zero is the plain ridge at lambda = 1
  q <- 3
  sc <- parameters7::log_cholesky(q)
  penS <- structured_penalty(sc)
  th0 <- stats::setNames(as.list(rep(0, sc@n_free)), sc@free_names)
  penR <- quadratic_penalty(diag(q))
  b3 <- c(0.4, -1.1, 2.2)
  if (abs(penalty_value(penS, b3, th0) -
          penalty_value(penR, b3, list(lambda = 1))) > 1e-12) {
    out <- c(out, "the structured prior at zero is not the plain ridge")
  }

  out
}

.penalty_stop <- function(problems, where) {
  if (length(problems)) {
    stop(sprintf("Chapter 8 (%s) disagrees with the package:\n  %s", where,
                 paste(problems, collapse = "\n  ")), call. = FALSE)
  }
  invisible(TRUE)
}

assert_penalty_definition_ok <- function() {
  .penalty_stop(.certify_penalty_definition(), "the definition")
}
assert_penalty_branches_ok <- function() {
  .penalty_stop(.certify_penalty_branches(), "the three constructions")
}
assert_penalty_marginal_ok <- function() {
  .penalty_stop(.certify_penalty_marginal(), "the marginal pieces")
}
