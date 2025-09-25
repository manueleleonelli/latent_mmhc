label_fun <- function(X) {
  if (any(is.na(X))) {
    stop("Data should not contain NA")
  } else if (!all(is.numeric(as.matrix(X)))) {
    stop("Data should only contain numbers")
  } else {
    y <- c()
    for (i in 1:ncol(X)) {
      temp <- X[, i]
      if (length(table(temp)) == 0) {
        stop("Data contain a variable with no value")
      } else if (length(table(temp)) == 1) {
        stop("Data contain a variable with a constant value")
      } else if (length(table(temp)) == 2) {
        if (all(temp %in% c(0, 1))) {
          y <- c(y, "binary")
        } else {
          y <- c(y, "continuous")
        }
      } else {
        level <- sort(unique(temp))
        if (all(temp %in% 0:max(level))) {
          y <- c(y, "ordinal")
        } else {
          y <- c(y, "continuous")
        }
      }
    }
  }
  return(y)
}

latent_pc <- function(X, label) {
  temp <- transform_fun(X, label)
  Y <- temp[[1]]
  index <- temp[[2]]
  delta <- temp[[3]]
  tau <- tau_fun(Y)
  sig <- matrix(1, ncol(X), ncol(X))
  for (i in 1:(ncol(X) - 1)) {
    for (j in (i + 1):ncol(X)) {
      temp_label <- label[c(i, j)]
      id1 <- which(index == i)
      id2 <- which(index == j)
      temp <- Y[, c(id1, id2)]
      if (all(temp_label == c("continuous", "continuous"))) {
        sig[i, j] <- CC_fun(tau[id1, id2])
      } else if (all(temp_label == c("continuous", "binary"))) {
        sig[i, j] <- optimize(BC_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[j])$minimum
      } else if (all(temp_label == c("binary", "continuous"))) {
        sig[i, j] <- optimize(BC_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[i])$minimum
      } else if (all(temp_label == c("binary", "binary"))) {
        sig[i, j] <- optimize(BB_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[c(i, j)])$minimum
      } else if (all(temp_label == c("continuous", "ordinal"))) {
        sig[i, j] <- CO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[j])
      } else if (all(temp_label == c("ordinal", "continuous"))) {
        sig[i, j] <- CO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[i])
      } else if (all(temp_label == c("binary", "ordinal"))) {
        sig[i, j] <- BO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[c(i, j)])
      } else if (all(temp_label == c("ordinal", "binary"))) {
        sig[i, j] <- BO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[c(i, j)])
      } else if (all(temp_label == c("ordinal", "ordinal"))) {
        sig[i, j] <- OO_fun(tau[c(id1, id2), c(id1, id2)], delta[c(i, j)])
      }
      sig[j, i] <- sig[i, j]
    }
  }
  sig <- lmf::nearPD(sig, corr = TRUE, maxit = 10000)
  return(sig)
}


transform_fun <- function(X, label) {
  Y <- data.frame(rep(NA, nrow(X)))
  index <- c()
  delta <- list()
  for (i in 1:ncol(X)) {
    if (label[i] == "continuous") {
      index <- c(index, i)
      Y <- data.frame(Y, X[, i])
      delta <- append(delta, mean(X[, i]))
    } else if (label[i] == "binary") {
      index <- c(index, i)
      Y <- data.frame(Y, X[, i])
      delta <- append(delta, qnorm(1 - mean(X[, i])))
    } else {
      level <- sort(unique(X[, i]))
      temp <- c()
      for (j in 1:max(level)) {
        a <- X[, i] >= j
        Y <- data.frame(Y, as.numeric(a))
        temp <- c(temp, qnorm(1 - mean(a)))
      }
      index <- c(index, rep(i, max(level)))
      delta <- append(delta, list(temp))
    }
  }
  return(list(Y[, -1], index, delta))
}

tau_fun <- function(X) {
  # Function from the latentcor package
  n_x <- function(x, n) {
    if (length(unique(x) != n)) {
      x.info <- rle(sort(x))
      t_x <- x.info$lengths[x.info$lengths > 1]
      n_x <- sum(t_x * (t_x - 1) / 2)
    } else {
      n_x <- 0
    }
    return(n_x)
  }
  
  n <- nrow(X)
  n0 <- n * (n - 1) / 2
  n_X <- apply(X, 2, function(x) {n_x(x = x, n)})
  n_X_sqd <- sqrt(n0 - n_X)
  K_b <- pcaPP::cor.fk(X)
  K_b.lower <- K_b[lower.tri(K_b)]
  btoa <- n_X_sqd[row(K_b)[lower.tri(K_b)]] * n_X_sqd[col(K_b)[lower.tri(K_b)]] / n0
  K_a.lower <- K_b.lower * btoa
  out <- matrix(0, ncol(X), ncol(X))
  out[lower.tri(out, diag = F)] <- K_a.lower
  out <- out + t(out)
  diag(out) <- 1
  return(out)
}

CC_fun <- function(tau) {
  sin(pi / 2 * tau)
}

BB_fun <- function(x, tau_x, delta_x) {
  delta_x <- unlist(delta_x)
  (2 * pmvnorm(upper = delta_x,
               mean = 0,
               sigma = matrix(c(1, x, x, 1), nrow = 2),
               algorithm = Miwa(steps = 128)) - 
      2 * pnorm(delta_x[1]) * pnorm(delta_x[2]) - tau_x)^2
}


BC_fun <- function(x, tau_x, delta_x) {
  delta_x <- unlist(delta_x)
  (4 * pmvnorm(upper = c(delta_x, 0),
               mean = 0,
               sigma = matrix(c(1, x / sqrt(2), x / sqrt(2), 1), nrow = 2),
               algorithm = Miwa(steps = 128)) -
      2 * pnorm(delta_x) - tau_x)^2
}

CO_fun <- function(tau_o, label_o, delta_o) {
  delta_o <- unlist(delta_o)
  z <- c()
  id <- which(label_o == "ordinal")
  level <- length(delta_o)
  for (h in 1:level) {
    z <- c(z, optimize(BC_fun, c(-0.9999, 0.9999),
                       tau_x = ifelse(id == 1, tau_o[max(level) + 1, h], tau_o[1, h + 1]),
                       delta_x = delta_o[h])$minimum)
  }
  mean(z)
}

BO_fun <- function(tau_o, label_o, delta_o) {
  z <- c()
  id <- which(label_o == "ordinal")
  level <- length(delta_o[[id]])
  for (h in 1:level) {
    z <- c(z, optimize(BB_fun, c(-0.9999, 0.9999),
                       tau_x = ifelse(id == 1, tau_o[max(level) + 1, h], tau_o[1, h + 1]),
                       delta_x = c(delta_o[[-id]], unlist(delta_o[[id]])[h]))$minimum)
  }
  mean(z)
}

OO_fun <- function(tau_o, delta_o) {
  z <- c()
  level1 <- length(delta_o[[1]])
  level2 <- length(delta_o[[2]])
  for (h in 1:max(level1)) { 
    for (k in 1:max(level2)) {
      z <- c(z, optimize(BB_fun, c(-0.9999, 0.9999),
                         tau_x = tau_o[h, level1 + k],
                         delta_x = c(delta_o[[1]][h], delta_o[[2]][k]))$minimum)
    }
  }
  mean(z)
}

make_numeric <- function(df) {
  for (j in seq_along(df)) {
    if (is.factor(df[[j]]) || is.ordered(df[[j]])) {
      df[[j]] <- as.numeric(df[[j]])
    }
  }
  return(df)
}

copula_score <- function(node, parents, data, args) {
  C <- args$corr
  less <- args$less
  gess <- args$gess
  type <- args$n_type
  
  idx_node <- which(colnames(data) == node)
  idx_parents <- which(colnames(data) %in% parents)
  
  if (length(idx_parents) == 0) return(0)
  
  C_pp <- C[idx_parents, idx_parents, drop = FALSE]
  C_np <- C[idx_node, idx_parents, drop = FALSE]
  
  inv_Cpp <- tryCatch(solve(C_pp), error = function(e) return(NULL))
  if (is.null(inv_Cpp)) return(-Inf)
  
  R2 <- as.numeric(C_np %*% inv_Cpp %*% t(C_np))
  R2 <- min(max(R2, 1e-6), 0.999)
  
  if (type == "global") {
    n_eff <- gess
  } else if (type == "average") {
    n_eff <- mean(less[idx_node, idx_parents, drop = FALSE])
  } else if (type == "min") {
    n_eff <- min(less[idx_node, idx_parents, drop = FALSE])
  } else if (type == "min_family") {
    family <- c(idx_node, idx_parents)
    combn_idx <- combn(family, 2)
    n_eff <- min(sapply(1:ncol(combn_idx), function(k) {
      i <- combn_idx[1, k]
      j <- combn_idx[2, k]
      less[i, j]
    }), na.rm = TRUE)
  } else {
    stop("Invalid n_type")
  }
  
  loglik <- - (n_eff / 2) * log(1 - R2)
  penalty <- (length(parents) / 2) * log(n_eff)
  return(loglik - penalty)
}

copula.estimate <- function(Y, n0 = dim(Y)[2] + 1, S0 = diag(dim(Y)[2])/n0,
                            nsamp = 100, odens = max(1, round(nsamp/1000)),
                            impute = any(is.na(Y)),
                            plugin.threshold = 100,
                            plugin.marginal = (apply(Y, 2, function(x) length(unique(x))) > plugin.threshold),
                            seed = NULL, verb = TRUE) {
  require(sbgcop)
  vnames <- colnames(Y)
  Y <- as.matrix(Y)
  colnames(Y) <- vnames
  n <- nrow(Y)
  p <- ncol(Y)
  set.seed(seed)
  R <- matrix(NA, n, p)
  for (j in 1:p) R[, j] <- match(Y[, j], sort(unique(Y[, j])))
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  S <- cov(Z)
  
  LPC <- NULL
  C.psamp <- array(dim = c(p, p, floor(nsamp/odens)))
  dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 1:floor(nsamp/odens))
  
  for (ns in 1:nsamp) {
    for (j in sample(1:p)) {
      Sjc <- S[j, -j] %*% solve(S[-j, -j])
      sdj <- sqrt(S[j, j] - S[j, -j] %*% solve(S[-j, -j]) %*% S[-j, j])
      muj <- Z[, -j] %*% t(Sjc)
      if (!plugin.marginal[j]) {
        for (r in 1:Rlevels[j]) {
          ir <- which(R[, j] == r & !is.na(R[, j]))
          lb <- suppressWarnings(max(Z[R[, j] == r - 1, j], na.rm = TRUE))
          ub <- suppressWarnings(min(Z[R[, j] == r + 1, j], na.rm = TRUE))
          Z[ir, j] <- qnorm(runif(length(ir), pnorm(lb, muj[ir], sdj), pnorm(ub, muj[ir], sdj)), muj[ir], sdj)
        }
      }
      ir <- which(is.na(R[, j]))
      Z[ir, j] <- rnorm(length(ir), muj[ir], sdj)
    }
    
    Z <- t(t(Z) - colMeans(Z))  # recenter
    S <- solve(rwish(solve(S0 * n0 + t(Z) %*% Z), n0 + n))
    
    if (ns %% odens == 0) {
      C <- S / (sqrt(diag(S)) %*% t(sqrt(diag(S))))
      C.psamp[, , ns / odens] <- C
    }
  }
  
  return(list(C.psamp = C.psamp))
}


gaussCItest.local <- function (x, y, S, suffStat) 
{
  matN = suffStat$ESS.mat
  ##
  n = suffStat$n
  if (!(is.null(matN)))
  {
    sub.mat = matN[c(x,y,S),c(x,y,S)]
    n = mean(sub.mat[upper.tri(sub.mat)], na.rm = T)
  }
  ##
  z <- zStat(x, y, S, C = suffStat$C, n = n)
  2 * pnorm(abs(z), lower.tail = FALSE)
}

# --- STRICT 4-ARG SIGNATURE ---
sem_copula_score <- function(node, parents, data, args) {
  # defaults (BIC) + allow override via args$penalty_type
  penalty_type <- if (!is.null(args$penalty_type)) args$penalty_type else "BIC"
  
  Z <- args$Z
  n <- args$n
  
  idx_node    <- which(colnames(Z) == node)
  idx_parents <- which(colnames(Z) %in% parents)
  
  # residual variance (or marginal variance if no parents)
  if (length(idx_parents) == 0) {
    sigma2 <- var(Z[, idx_node])
  } else {
    fit <- lm(Z[, idx_node] ~ Z[, idx_parents])
    sigma2 <- mean(residuals(fit)^2)
  }
  
  sigma2 <- max(sigma2, 1e-6)
  loglik <- - (n / 2) * log(sigma2)
  
  k <- length(parents)  # number of regression coefficients (excluding intercept)
  penalty <- switch(penalty_type,
                    "BIC" = (k / 2) * log(n),
                    "AIC" = k,
                    stop("penalty_type must be 'BIC' or 'AIC'"))
  loglik - penalty
}


estimate_copula_args <- function(data, nsamp = 1000, S0_scale = 100, plugin_threshold = 50, verbose = FALSE) {
  if (!requireNamespace("sbgcop", quietly = TRUE)) {
    stop("The 'sbgcop' package is required but not installed.")
  }
  
  p <- ncol(data)
  
  cop.obj <- copula.estimate(
    Y = data,
    nsamp = nsamp,
    S0 = diag(p) / S0_scale,
    plugin.threshold = plugin_threshold,
    verb = verbose
  )
  
  # Extract posterior samples and retain second half (for convergence)
  C_samples <- cop.obj$C.psamp[, , (nsamp / 2 + 1):nsamp]
  corr.cop <- apply(C_samples, c(1, 2), mean)
  
  # Local and global effective sample size
  less.cop <- ((1 - corr.cop^2)^2) / apply(C_samples, c(1, 2), var)
  gess.cop <- mean(less.cop[upper.tri(less.cop)])
  
  args_list <- list(
    corr = corr.cop,
    less = less.cop,
    gess = gess.cop
  )
  
  return(args_list)
}

estimate_copula_correlation_kendall <- function(data) {
  # Ensure numeric input for copula inversion
  data_num <- make_numeric(data)
  var_names <- colnames(data_num)
  
  # Step 1: Label the variables
  labels <- label_fun(data_num)
  
  # Step 2: Estimate latent correlation matrix from Kendall’s tau
  latent_corr <- latent_pc(data_num, labels)
  
  colnames(latent_corr) <- rownames(latent_corr) <- var_names
  return(latent_corr)
}

