logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Step 1: Predictive Recursion
predictive_recursion <- function(z, grid_u, f0 = dnorm, decay = 0.67, passes = 1) {
  n <- length(z)
  m <- length(grid_u)
  delta <- diff(grid_u)[1]
  G <- rep(1 / m, m)

  for (p in 1:passes) {
    for (i in 1:n) {
      wi <- ((p - 1) * n + i)^(-decay)
      f_kernel <- f0(z[i] - grid_u)
      f_pred <- sum(f_kernel * G) * delta
      G_new <- (1 - wi) * G + wi * f_kernel * G / f_pred
      G <- G_new / sum(G_new)
    }
  }

  list(G = G, grid_u = grid_u, delta = delta)
}

# Step 2–4: Full FDR regression with PR
fdr_regression <- function(z, X, grid_u = seq(-6.5, 6.5, length.out = 300),
                           decay = 0.67, passes = 1) {
  n <- length(z)

  # Estimate f1 using PR
  pr <- predictive_recursion(z, grid_u, decay = decay, passes = passes)

  # Compute f0 and f1 at each z
  f0_z <- dnorm(z, mean = 0, sd = 1)
  f1_z <- sapply(z, function(zi) {
    sum(dnorm(zi - pr$grid_u) * pr$G) * pr$delta
  })

  # Clip for numerical stability
  f0_z <- pmax(f0_z, .Machine$double.eps)
  f1_z <- pmax(f1_z, .Machine$double.eps)

  # Step 3: Estimate c(x) via logistic regression on log-likelihood ratio
  log_ratio <- log(f1_z / f0_z)
  glmfit <- glm(log_ratio ~ X, family = gaussian())

  # Step 4: Compute prior probs and posterior
  eta <- as.numeric(cbind(1, X) %*% coef(glmfit))
  c_x <- inv_logit(eta)
  numer <- c_x * f1_z
  denom <- numer + (1 - c_x) * f0_z
  post_prob <- numer / denom
  local_fdr <- 1 - post_prob

  list(
    z = z,
    X = X,
    f0_z = f0_z,
    f1_z = f1_z,
    prior_prob = c_x,
    post_prob = post_prob,
    local_fdr = local_fdr,
    G = pr$G,
    grid_u = pr$grid_u
  )
}
