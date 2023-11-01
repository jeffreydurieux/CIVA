fastIVA_custom <- function (X, nc, source_density = "laplace_diag", student_df = 1, 
          max_iter = 1024, eps = 1e-06, W_init = NA, verbose = FALSE) 
{
  source_density <- match.arg(source_density, c("laplace_diag", 
                                                "student", "entropic"))
  D <- dim(X)[3]
  P <- dim(X)[1]
  
  N <- dim(X)[2]
  whitened <- whiten_custom(X,nc)
  Z <- whitened$Z
  V <- whitened$V
  if (!is.array(W_init)) {
    W <- array(rnorm(nc * nc * D), c(nc, nc, D))
  }
  # else if (all(dim(W_init) == c(P, P, D))) {
  #   W <- W_init
  # }
  # else {
  #   warning("W_init is not in correct shape, identity matrices used as initial W instead.")
  #   W <- array(rep(diag(nrow = P, ncol = P), D), c(P, P, 
  #                                                  D))
  # }
  S <- array(0, c(nc, N, D))
  for (d in 1:D) {
    S[, , d] <- W[, , d] %*% Z[, , d]
  }
  dG <- NULL
  d2G <- NULL
  if (source_density == "student") {
    dG <- function(x) {
      1/(1 + x/student_df)
    }
    d2G <- function(x) {
      x/(student_df + x)^2
    }
  }else if (source_density == "laplace_diag") {
    dG <- function(x) {
      1/(2 * sqrt(x))
    }
    d2G <- function(x) {
      -1/(4 * x^(3/2))
    }
  }else if (source_density == "entropic") {
    dG <- function(x) {
      1/x
    }
    d2G <- function(x) {
      -1/x^2
    }
  }
  "%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
  }
  converged <- FALSE
  for (iter in 1:max_iter) {
    W_old <- W
    for (j in 1:nc) {
      for (d in 1:D) {
        G1 <- sapply(1:N, FUN = function(i) {
          dG(sum(S[j, i, ]^2))
        })
        G2 <- sapply(1:N, FUN = function(i) {
          S[j, i, d]^2 * d2G(sum(S[j, i, ]^2))
        })
        G3 <- sapply(1:N, FUN = function(i) {
          S[j, i, d] * dG(sum(S[j, i, ]^2)) * Z[, i, 
                                                d]
        })
        new_val <- mean(G1 + G2) * W[j, , d] - rowMeans(G3)
        W[j, , d] <- new_val
      }
    }
    for (d in 1:D) {
      W[, , d] <- (W[, , d] %*% t(W[, , d])) %^% (-1/2) %*% 
        W[, , d]
      S[, , d] <- W[, , d] %*% Z[, , d]
    }
    tol <- 0
    for (d in 1:D) {
      tol <- max(tol, 1 - min(diag(tcrossprod(W_old[, , 
                                                    d], W[, , d]))))
    }
    if (verbose) {
      print(paste0("Convergence measure: ", tol))
    }
    if (eps > tol) {
      converged <- TRUE
      break
    }
  }
  #W_nonwhitened <- array(0, dim(W))
  #for (d in 1:D) {
  #  W_nonwhitened[, , d] <- W[, , d] %*% V[, , d]
  #}
  
  dimnames(S)[[1]] <- sapply(1:nc, FUN = function(j) paste0("IC.", 
                                                           j))
  West <- array(data = NA, dim = c(nc,P,D))
  for(i in 1:D){
    West[,,i] <- W[, , i] %*% t(V[,,i])
  }
  
  RES <- list(S = S,W = West, W_whitened = W, V=V,
              X_means = whitened$X_means, niter = iter, converged = converged, 
              source_density = source_density, N = N, D = D, P = P, 
              student_df = student_df, call = deparse(sys.call()), 
              DNAME = paste(deparse(substitute(X))))
  class(RES) <- "iva"
  return(RES)
}
