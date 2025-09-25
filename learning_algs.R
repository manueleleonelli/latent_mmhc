learn_dag_copc <- function(data, corr.cop, gess.cop, alpha = 0.01) {
  require(pcalg)
  
  # Basic checks
  stopifnot(is.matrix(corr.cop) || is.data.frame(corr.cop))
  stopifnot(length(gess.cop) == 1 && is.numeric(gess.cop) && gess.cop > 0)
  
  p <- ncol(data)
  var.names <- colnames(data)
  
  # Construct suffStat with global effective sample size
  suffStat <- list(C = corr.cop, n = gess.cop)
  
  # PC algorithm with Gaussian CI test
  pc_fit <- pc(suffStat = suffStat,
               indepTest = gaussCItest,
               labels = var.names,
               alpha = alpha,
               conservative = TRUE)
  
  # Convert to bnlearn DAG object
  dag <- bnlearn::as.bn(pc_fit, check.cycles = FALSE)
  
  return(dag)
}

learn_dag_latent_pc <- function(data, corr, alpha = 0.01) {
  require(pcalg)
  
  # Step 1: Preprocess
  var_names <- colnames(data)
  n <- nrow(data)
  
  
  # Step 3: PC algorithm with Gaussian CI test
  suffStat <- list(C = corr, n = n)
  pc_fit <- pc(suffStat = suffStat,
               indepTest = gaussCItest,
               labels = var_names,
               alpha = alpha,
               conservative = TRUE)
  
  # Convert to bnlearn DAG object
  dag <- bnlearn::as.bn(pc_fit,check.cycles = FALSE)
  
  return(dag)
}

learn_dag_tau_score <- function(data, corr, blacklist = NULL) {
  n <- nrow(data)
  
  args_list <- list(corr = corr, less = NULL, gess = n, n_type = "global")  # treat ESS = n for now
  
  dag <- bnlearn::tabu(data,
                       score = "custom-score",
                       fun = copula_score,
                       args = args_list,
                       blacklist = blacklist)
  return(dag)
}

learn_dag_copula <- function(data, args_list, blacklist = NULL) {
  # Use minimum ESS between node and each parent (Castelletti-style)
  args_list$n_type <- "min"
  
  dag <- bnlearn::tabu(data,
                       score = "custom-score",
                       fun = copula_score,
                       args = args_list,
                       blacklist = blacklist)
  return(dag)
}



learn_dag_sem <- function(data, blacklist = NULL) {
  # Generate approximate latent Z using rank-based quantiles
  ranks <- apply(data, 2, rank, ties.method = "average") / (nrow(data) + 1)
  Z <- qnorm(ranks)
  
  args_list <- list(Z = Z, n = nrow(Z))
  
  dag <- bnlearn::tabu(data,
                       score = "custom-score",
                       fun = sem_copula_score,
                       args = args_list,
                       blacklist = blacklist)
  return(dag)
}

learn_dag_mmhc_sem <- function(data, corr, alpha = 0.01,
                               use_cholesky = TRUE,
                               verbose = FALSE,
                               penalty_type = "BIC",
                               blacklist = NULL) {
  require(pcalg)
  require(bnlearn)
  
  stopifnot(is.data.frame(data))
  var.names <- colnames(data)
  n <- nrow(data)
  
  # PC skeleton from corr
  suffStat <- list(C = corr, n = n)
  pc_fit <- pc(suffStat = suffStat, indepTest = gaussCItest,
               labels = var.names, alpha = alpha, conservative = TRUE)
  amat_pc <- as(pc_fit@graph, "matrix")
  
  # skeleton -> blacklist
  skeleton <- which((amat_pc + t(amat_pc)) > 0, arr.ind = TRUE)
  skeleton <- skeleton[skeleton[,1] != skeleton[,2], , drop = FALSE]
  skel_names <- data.frame(from = var.names[skeleton[,1]],
                           to   = var.names[skeleton[,2]], 
                           stringsAsFactors = FALSE)
  all_edges <- expand.grid(from = var.names, to = var.names, stringsAsFactors = FALSE)
  all_edges <- all_edges[all_edges$from != all_edges$to, ]
  keep_edges <- rbind(skel_names, skel_names[, c("to","from")])
  edge_str <- paste(keep_edges$from, keep_edges$to, sep = "->")
  all_str  <- paste(all_edges$from, all_edges$to,  sep = "->")
  blacklist <- all_edges[!all_str %in% edge_str, ]
  
  # latent Gaussian Z from ranks; optional Cholesky whiten
  ranks <- apply(data, 2, rank, ties.method = "average") / (n + 1)
  Z_raw <- qnorm(ranks)
  if (use_cholesky) {
    chol_corr <- tryCatch(chol(corr), error = function(e) {
      eps <- 1e-6; chol((1 - eps) * corr + eps * diag(ncol(corr)))
    })
    Z <- scale(Z_raw, center = TRUE, scale = FALSE) %*% chol_corr
  } else {
    Z <- Z_raw
  }
  
  # pass penalty_type here
  args_list <- list(Z = Z, n = n, penalty_type = penalty_type)
  
  bnlearn::tabu(
    data,
    score = "custom-score",
    fun = sem_copula_score,   # <-- 4-arg function above
    args = args_list,
    blacklist = blacklist
  )
}




learn_dag_mmhc_gauss <- function(data, corr, alpha = 0.01, verbose = FALSE) {
  require(pcalg)
  require(bnlearn)
  
  stopifnot(is.data.frame(data))
  var.names <- colnames(data)
  p <- ncol(data)
  n <- nrow(data)
  
  # Step 1: PC Phase using Kendall-based correlation
  suffStat <- list(C = corr, n = n)
  pc_fit <- pc(suffStat = suffStat,
               indepTest = gaussCItest,
               labels = var.names,
               alpha = alpha,
               conservative = TRUE)
  
  amat_pc <- as(pc_fit@graph, "matrix")
  
  # Step 2: Extract skeleton edges
  skeleton <- which((amat_pc + t(amat_pc)) > 0, arr.ind = TRUE)
  skeleton <- skeleton[skeleton[, 1] != skeleton[, 2], , drop = FALSE]
  skeleton_names <- data.frame(
    from = var.names[skeleton[, 1]],
    to   = var.names[skeleton[, 2]],
    stringsAsFactors = FALSE
  )
  
  # Step 3: Construct blacklist (edges not in skeleton)
  all_edges <- expand.grid(from = var.names, to = var.names, stringsAsFactors = FALSE)
  all_edges <- all_edges[all_edges$from != all_edges$to, ]
  keep_edges <- rbind(skeleton_names, skeleton_names[ , c("to", "from")])
  edge_str <- paste(keep_edges$from, keep_edges$to, sep = "→")
  all_str <- paste(all_edges$from, all_edges$to, sep = "→")
  blacklist <- all_edges[!all_str %in% edge_str, ]
  
  # Step 4: Gaussian log-likelihood score based on correlation
  args_list <- list(corr = corr, less = NULL, gess = n, n_type = "global")
  
  if (verbose) cat("Running constrained hill-climbing with Gaussian BIC score...\n")
  dag <- bnlearn::tabu(data,
                       score = "custom-score",
                       fun = copula_score,
                       args = args_list,
                       blacklist = blacklist)
  
  return(dag)
}


