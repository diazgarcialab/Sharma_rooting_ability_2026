# ==============================================================================
# 10_prediction_engine.R -- GBLUP in base R. FUNCTIONS ONLY, no side effects.
#
# Sourced by 11_run_predictions.R and by 12_verify_prediction.R. Running it on
# its own does nothing.
#
# ==============================================================================

kernel_G <- function(Z) {
  Z <- as.matrix(Z)
  p <- colMeans(Z, na.rm = TRUE) / 2
  keep <- is.finite(p) & p > 0 & p < 1
  Z <- Z[, keep, drop = FALSE]; p <- p[keep]
  na <- which(is.na(Z), arr.ind = TRUE)
  if (nrow(na)) Z[na] <- 2 * p[na[, "col"]]
  W <- sweep(Z, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  G <- tcrossprod(W) / denom
  # a small ridge keeps G invertible when markers outnumber individuals
  G + diag(1e-6, nrow(G))
}


kernel_S <- function(X) {
  X <- scale(as.matrix(X))
  X <- X[, apply(X, 2, \(v) all(is.finite(v))), drop = FALSE]
  S <- tcrossprod(X) / ncol(X)
  S + diag(1e-6, nrow(S))
}


kernel_T <- function(species) {
  sp <- as.character(species)
  K <- outer(sp, sp, "==") * 1
  K + diag(1e-6, nrow(K))
}

# ---- single-kernel REML ------------------------------------------------------

#' REML for y = 1 mu + u + e via eigendecomposition, profiling over
#' delta = sigma_g^2 / sigma_e^2.
#' @return list(delta, Vg, Ve, mu, h2, logLik)
reml_1k <- function(y, K, log_delta_range = c(-9, 9)) {
  n <- length(y)
  ev <- eigen(K, symmetric = TRUE)
  d <- pmax(ev$values, 0)
  U <- ev$vectors
  ys <- crossprod(U, y)
  Xs <- crossprod(U, matrix(1, n, 1))
  p <- 1L

  negll <- function(log_delta) {
    delta <- exp(log_delta)
    w <- 1 / (delta * d + 1)
    XtWX <- sum(Xs^2 * w)
    if (!is.finite(XtWX) || XtWX <= 0) return(1e10)
    b <- sum(Xs * w * ys) / XtWX
    r <- ys - Xs * b
    s2 <- sum(r^2 * w) / (n - p)
    if (!is.finite(s2) || s2 <= 0) return(1e10)
    0.5 * (sum(log(delta * d + 1)) + (n - p) * log(s2) + log(XtWX))
  }
  opt <- optimize(negll, log_delta_range, tol = 1e-8)
  delta <- exp(opt$minimum)
  w <- 1 / (delta * d + 1)
  XtWX <- sum(Xs^2 * w)
  mu <- sum(Xs * w * ys) / XtWX
  r <- ys - Xs * mu
  Ve <- sum(r^2 * w) / (n - p)
  Vg <- delta * Ve
  list(delta = delta, Vg = Vg, Ve = Ve, mu = mu,
       h2 = Vg / (Vg + Ve), logLik = -opt$objective)
}

#' Fit on a training set and predict a held-out set with one kernel.
#' K must be the FULL kernel over train and test together, in that indexing.
gblup_1k <- function(y_full, K, train, test) {
  y <- y_full[train]
  ok <- is.finite(y)
  train <- train[ok]; y <- y[ok]
  if (length(train) < 10) return(rep(NA_real_, length(test)))
  Kt <- K[train, train, drop = FALSE]
  fit <- reml_1k(y, Kt)
  lambda <- 1 / fit$delta                       # sigma_e^2 / sigma_g^2
  A <- Kt + diag(lambda, nrow(Kt))
  sol <- tryCatch(solve(A, y - fit$mu), error = function(e) NULL)
  if (is.null(sol)) return(rep(NA_real_, length(test)))
  as.vector(fit$mu + K[test, train, drop = FALSE] %*% sol)
}

# ---- multi-kernel REML -------------------------------------------------------

#' REML for y = 1 mu + sum_k u_k + e over several kernels.
#' Parameterised on the log scale so variances stay positive.
reml_mk <- function(y, Ks, tol = 1e-7) {
  n <- length(y); X <- matrix(1, n, 1); p <- 1L
  nk <- length(Ks)
  v0 <- var(y) / (nk + 1)

  negll <- function(theta) {
    s <- exp(theta)
    V <- diag(s[nk + 1], n)
    for (k in seq_len(nk)) V <- V + s[k] * Ks[[k]]
    ch <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    logdetV <- 2 * sum(log(diag(ch)))
    Vi_X <- backsolve(ch, forwardsolve(t(ch), X))
    Vi_y <- backsolve(ch, forwardsolve(t(ch), y))
    XtViX <- crossprod(X, Vi_X)
    ch2 <- tryCatch(chol(XtViX), error = function(e) NULL)
    if (is.null(ch2)) return(1e10)
    beta <- backsolve(ch2, forwardsolve(t(ch2), crossprod(X, Vi_y)))
    r <- y - X %*% beta
    Vi_r <- backsolve(ch, forwardsolve(t(ch), r))
    val <- 0.5 * (logdetV + 2 * sum(log(diag(ch2))) + sum(r * Vi_r))
    if (!is.finite(val)) 1e10 else val
  }
  opt <- optim(rep(log(v0), nk + 1), negll, method = "Nelder-Mead",
               control = list(reltol = tol, maxit = 800))
  s <- exp(opt$par)
  list(Vk = s[seq_len(nk)], Ve = s[nk + 1], logLik = -opt$value,
       converged = opt$convergence == 0)
}

#' Fit on a training set and predict a held-out set with several kernels.
gblup_mk <- function(y_full, Ks, train, test) {
  y <- y_full[train]
  ok <- is.finite(y)
  train <- train[ok]; y <- y[ok]
  if (length(train) < 10) return(rep(NA_real_, length(test)))
  Kt <- lapply(Ks, \(K) K[train, train, drop = FALSE])
  fit <- reml_mk(y, Kt)
  n <- length(y)
  V <- diag(fit$Ve, n)
  for (k in seq_along(Kt)) V <- V + fit$Vk[k] * Kt[[k]]
  ch <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(ch)) return(rep(NA_real_, length(test)))
  Vi_1 <- backsolve(ch, forwardsolve(t(ch), rep(1, n)))
  Vi_y <- backsolve(ch, forwardsolve(t(ch), y))
  mu <- sum(Vi_y) / sum(Vi_1)
  r <- y - mu
  Vi_r <- backsolve(ch, forwardsolve(t(ch), r))
  Cts <- matrix(0, length(test), n)
  for (k in seq_along(Ks)) Cts <- Cts + fit$Vk[k] * Ks[[k]][test, train, drop = FALSE]
  as.vector(mu + Cts %*% Vi_r)
}

# ---- spectral preprocessing --------------------------------------------------
# Reflectance carries instrument and geometry artefacts -- leaf angle, distance,
# stray light -- that scale the whole spectrum without saying anything about the
# plant. These are the two standard corrections.

#' Standard normal variate: centre and scale each SPECTRUM (row), removing
#' multiplicative scatter.
snv <- function(X) t(scale(t(as.matrix(X))))

#' Savitzky-Golay filter along the wavelength axis. Implemented directly so the
#' pipeline needs no extra package; `m` is the derivative order.
sav_gol <- function(X, width = 11, poly = 2, m = 0) {
  X <- as.matrix(X)
  stopifnot(width %% 2 == 1, poly < width)
  half <- (width - 1) %/% 2
  z <- -half:half
  A <- outer(z, 0:poly, "^")
  Cmat <- solve(crossprod(A), t(A))
  coef <- Cmat[m + 1, ] * factorial(m)
  out <- t(apply(X, 1, function(v) {
    padded <- c(rep(v[1], half), v, rep(v[length(v)], half))
    vapply(seq_along(v), \(i) sum(coef * padded[i:(i + width - 1)]), numeric(1))
  }))
  dimnames(out) <- dimnames(X)
  out
}

preprocess_spectra <- function(X, method = c("raw", "snv", "snv_d1", "snv_d2")) {
  method <- match.arg(method)
  switch(method,
    raw    = as.matrix(X),
    snv    = snv(X),
    snv_d1 = sav_gol(snv(X), m = 1),
    snv_d2 = sav_gol(snv(X), m = 2))
}

# ---- accuracy metrics --------------------------------------------------------
#' Correlation is the headline, but a model can correlate well and still be
#' badly scaled. The slope of observed on predicted is 1 for a calibrated
#' prediction; below 1 means the predictions are spread too wide.
pred_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  o <- obs[ok]; p <- pred[ok]
  if (length(o) < 5 || sd(p) < 1e-10 || sd(o) < 1e-10)
    return(c(r = NA_real_, slope = NA_real_, rmse = NA_real_, n = length(o)))
  c(r     = cor(o, p),
    slope = unname(coef(lm(o ~ p))[2]),
    rmse  = sqrt(mean((o - p)^2)),
    n     = length(o))
}

