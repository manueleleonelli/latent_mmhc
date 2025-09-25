# Required: bnlearn and your helper functions already loaded:
# - simulate_copula_data(n, p, s)
# - estimate_copula_correlation_kendall(data)
# - learn_dag_mmhc_sem(data, corr, alpha = 0.01, use_cholesky = FALSE, blacklist = NULL)
# - bnlearn::custom.strength, averaged.network, inclusion.threshold, compare, shd

library(bnlearn)

# ----- Bootstrap helpers (as provided, lightly wrapped) -----

bootstrap_sem_copula <- function(data, B = 200, alpha = 0.01, verbose = TRUE, blacklist = NULL) {
  var_names <- colnames(data)
  p <- ncol(data)
  dag_list <- vector("list", B)
  
  for (b in seq_len(B)) {
    if (verbose) cat("Bootstrap iteration", b, "\n")
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    data_boot <- data[idx, , drop = FALSE]
    
    corr_boot <- estimate_copula_correlation_kendall(data_boot)
    
    dag_list[[b]] <- tryCatch(
      learn_dag_mmhc_sem(
        data_boot,
        corr = corr_boot,
        alpha = alpha,
        use_cholesky = FALSE,
        blacklist = blacklist
      ),
      error = function(e) {
        if (verbose) warning(paste("Failed on bootstrap", b, ":", e$message))
        NULL
      }
    )
  }
  
  dag_list
}

learn_dag_mmhc_sem_bootstrap <- function(data, alpha = 0.01, B = 200, blacklist = NULL, verbose = FALSE) {
  dag_boots <- bootstrap_sem_copula(data, B = B, alpha = alpha, verbose = verbose, blacklist = blacklist)
  dag_boots <- Filter(Negate(is.null), dag_boots)
  if (length(dag_boots) == 0) stop("All bootstrap replicates failed.")
  
  strengths <- custom.strength(
    networks = dag_boots,
    nodes = names(data),
    cpdag = TRUE
  )
  averaged.network(strengths, inclusion.threshold(strengths))
}

# ----- Main simulation runner: baseline vs bootstrap -----

run_sim_mmhc_vs_bootstrap <- function(
    sample_sizes     = c(500, 1000),
    p_list           = c(10),
    sparsity_factors = c(2, 5),
    reps             = 2,
    alpha            = 0.01,
    B_boot           = 10,
    out_dir          = "simulationresults_mmvsboot",
    seed             = 1234
) {
  set.seed(seed)
  dir.create(out_dir, showWarnings = FALSE)
  
  fail_log <- data.frame()
  sim_id <- 1
  
  methods <- c("mmhc_sem_raw", "mmhc_sem_boot")
  
  for (p in p_list) {
    for (sf in sparsity_factors) {
      s <- sf / (p - 1)
      
      for (n in sample_sizes) {
        filename <- sprintf("%s/simresults_p%d_s%.4f_n%d.csv", out_dir, p, s, n)
        if (file.exists(filename)) {
          message(sprintf("Skipping existing file: %s", filename))
          next
        }
        
        results <- data.frame()
        
        for (r in seq_len(reps)) {
          message(sprintf("[Run %d] p=%d | s=%.4f | n=%d | rep=%d", sim_id, p, s, n, r))
          
          # ---- simulate mixed data and true DAG ----
          sim <- simulate_copula_data(n = n, p = p, s = s)
          data <- sim$data[, sample(seq_len(ncol(sim$data)))]  # randomize column order
          true_dag <- sim$dag
          
          # Pre-compute Kendall latent corr for baseline (same for all methods on this dataset)
          cop_ken <- tryCatch(
            estimate_copula_correlation_kendall(data),
            error = function(e) {
              fail_log <<- rbind(fail_log, data.frame(
                sim_id = sim_id, p = p, s = s, n = n, rep = r,
                method = "corr_kendall", error = as.character(e$message)
              ))
              return(NULL)
            }
          )
          if (is.null(cop_ken)) next
          
          # ---- Evaluate both methods ----
          for (method in methods) {
            start_time <- Sys.time()
            
            est_dag <- tryCatch({
              if (method == "mmhc_sem_raw") {
                learn_dag_mmhc_sem(
                  data, corr = cop_ken, alpha = alpha,
                  use_cholesky = FALSE, blacklist = NULL
                )
              } else {
                learn_dag_mmhc_sem_bootstrap(
                  data, alpha = alpha, B = B_boot,
                  blacklist = NULL, verbose = FALSE
                )
              }
            }, error = function(e) {
              fail_log <<- rbind(fail_log, data.frame(
                sim_id = sim_id, p = p, s = s, n = n, rep = r,
                method = method, error = as.character(e$message)
              ))
              return(NULL)
            })
            
            runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            
            if (!is.null(est_dag)) {
              comp <- bnlearn::compare(target = true_dag, current = est_dag)
              TP <- comp$tp; FP <- comp$fp; FN <- comp$fn
              
              SEN <- if ((TP + FN) > 0) TP / (TP + FN) else NA
              
              # total possible directed arcs (no self-loops)
              total_arcs <- p * (p - 1)
              TN <- total_arcs - TP - FP - FN
              SPE <- if ((TN + FP) > 0) TN / (TN + FP) else NA
              
              SHD <- bnlearn::shd(est_dag, true_dag)
              
              results <- rbind(results, data.frame(
                sim_id = sim_id, p = p, s = s, n = n, rep = r,
                method = method,
                TP = TP, FP = FP, FN = FN, SEN = SEN, SPE = SPE, SHD = SHD,
                runtime_sec = runtime_sec,
                stringsAsFactors = FALSE
              ))
            }
          }
          
          sim_id <- sim_id + 1
        } # reps
        
        if (nrow(results) > 0) {
          write.csv(results, file = filename, row.names = FALSE)
          message(sprintf("Saved: %s", filename))
        }
      } # n
    } # sf
  } # p
  
  message("=== Failures (if any) ===")
  print(fail_log)
  invisible(list(fail_log = fail_log))
}


res <- run_sim_mmhc_vs_bootstrap()
