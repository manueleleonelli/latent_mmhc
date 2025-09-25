simulate_copula_data <- function(n = 500, p = 20, s = 0.1, seed = NULL) {
  stopifnot(p %% 2 == 0)  # p must be even
  if (!is.null(seed)) set.seed(seed)
  
  # Step 1: Create weighted adjacency matrix A (lower triangle = DAG)
  A <- matrix(0, p, p)
  for (i in 2:p) {
    for (j in 1:(i - 1)) {
      if (runif(1) < s) {
        A[i, j] <- runif(1, 0.1, 1) * sample(c(-1, 1), 1)
      }
    }
  }
  
  # Step 2: Generate latent Gaussian data Z
  I_A_inv <- solve(diag(p) - A)
  epsilon <- matrix(rnorm(n * p), nrow = n)
  Z <- epsilon %*% t(I_A_inv)
  Z <- scale(Z)
  
  # Step 3: Transform Z into mixed-type data X using Castelletti's copula approach
  X <- matrix(NA, n, p)
  var_types <- character(p)
  
  for (j in 1:(p / 2)) {
    u <- pnorm(Z[, j])  # transform to uniform via standard Gaussian CDF
    
    if (runif(1) < 0.5) {
      # Binary: F_j ~ Bernoulli(θ_j) with θ_j ~ Uniform[0.2, 0.8]
      theta <- runif(1, 0.2, 0.8)
      X[, j] <- as.numeric(u > (1 - theta))  # inverse CDF
      var_types[j] <- "binary"
    } else {
      # Ordinal: F_j ~ Binomial(4, θ_j), then get discrete values 0–4
      theta <- runif(1, 0.2, 0.8)
      qvals <- qbinom(u, size = 4, prob = theta)
      X[, j] <- qvals
      var_types[j] <- "ordinal"
    }
  }
  
  for (j in ((p / 2) + 1):p) {
    X[, j] <- Z[, j]
    var_types[j] <- "continuous"
  }
  
  # Step 4: Construct data.frame with appropriate classes
  data <- as.data.frame(X)
  colnames(data) <- paste0("v", 1:p)
  
  for (j in 1:(p / 2)) {
    if (var_types[j] == "binary") {
      data[[j]] <- factor(data[[j]], levels = c(0, 1))
    } else if (var_types[j] == "ordinal") {
      data[[j]] <- ordered(data[[j]])
    }
  }
  
  for (j in ((p / 2) + 1):p) {
    data[[j]] <- as.numeric(data[[j]])
  }
  
  # Step 5: Construct bnlearn DAG object
  dag <- bnlearn::empty.graph(colnames(data))
  bnlearn::amat(dag) <- t((A != 0) * 1)
  
  return(list(data = data, Z = Z, dag = dag, amat = A, types = var_types))
}
