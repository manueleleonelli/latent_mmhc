library(bnlearn)
run_full_simulation_study <- function(
    sample_sizes = c(500, 1000, 2500, 5000),
    p_list = c(10, 50),
    sparsity_factors = c(2, 5),
    reps = 100,
    nsamp = 1000,
    seed = 1234
) {
  set.seed(seed)
  dir.create("simulationresults", showWarnings = FALSE)
  fail_log <- data.frame()
  sim_id <- 1
  
  methods <- c("copc", "latent_pc", "tau_score", "copula_score", "sem_score",
               "mmhc_sem_01", "mmhc_sem_05", "mmhc_sem_01_raw", "mmhc_sem_05_raw",
               "mmhc_gauss_01", "mmhc_gauss_05")
  
  for (p in p_list) {
    for (sf in sparsity_factors) {
      s <- sf / (p - 1)
      
      for (n in sample_sizes) {
        filename <- sprintf("simulationresults/simresults_p%d_s%.4f_n%d.csv", p, s, n)
        
        if (file.exists(filename)) {
          message(sprintf("Skipping: %s (already exists)", filename))
          next
        }
        
        results <- data.frame()
        
        for (r in 1:reps) {
          message(sprintf("Running: p = %d, s = %.4f, n = %d, rep = %d", p, s, n, r))
          
          sim <- simulate_copula_data(n = n, p = p, s = s)
          data <- sim$data[, sample(seq_len(ncol(sim$data)))]
          true_dag <- sim$dag
          
          args <- estimate_copula_args(data)
          cop_ken <- estimate_copula_correlation_kendall(data)
          
          for (method in methods) {
            start_time <- Sys.time()
            
            est_dag <- tryCatch({
              switch(method,
                     copc = learn_dag_copc(data, args$corr, args$gess, alpha = 0.01),
                     latent_pc = learn_dag_latent_pc(data, cop_ken, alpha = 0.01),
                     tau_score = learn_dag_tau_score(data, cop_ken),
                     copula_score = learn_dag_copula(data, args),
                     sem_score = learn_dag_sem(data),
                     mmhc_sem_01 = learn_dag_mmhc_sem(data, cop_ken, alpha = 0.01, use_cholesky = TRUE),
                     mmhc_sem_05 = learn_dag_mmhc_sem(data, cop_ken, alpha = 0.05, use_cholesky = TRUE),
                     mmhc_sem_01_raw = learn_dag_mmhc_sem(data, cop_ken, alpha = 0.01, use_cholesky = FALSE),
                     mmhc_sem_05_raw = learn_dag_mmhc_sem(data, cop_ken, alpha = 0.05, use_cholesky = FALSE),
                     mmhc_gauss_01 = learn_dag_mmhc_gauss(data, cop_ken, alpha = 0.01),
                     mmhc_gauss_05 = learn_dag_mmhc_gauss(data, cop_ken, alpha = 0.05),
                     NULL
              )
            }, error = function(e) {
              fail_log <<- rbind(fail_log, data.frame(
                sim_id = sim_id,
                p = p,
                s = s,
                n = n,
                rep = r,
                method = method,
                error = as.character(e$message)
              ))
              return(NULL)
            })
            
            end_time <- Sys.time()
            runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
            
            if (!is.null(est_dag)) {
              comp <- bnlearn::compare(target = true_dag, current = est_dag)
              
              TP <- comp$tp
              FP <- comp$fp
              FN <- comp$fn
              
              SEN <- if ((TP + FN) > 0) TP / (TP + FN) else NA
              total_edges <- choose(p, 2) * 2
              TN <- total_edges - TP - FP - FN
              SPE <- if ((TN + FP) > 0) TN / (TN + FP) else NA
              SHD <- bnlearn::shd(est_dag, true_dag)
              
              results <- rbind(results, data.frame(
                sim_id = sim_id,
                p = p,
                s = s,
                n = n,
                rep = r,
                method = method,
                TP = TP,
                FP = FP,
                FN = FN,
                SEN = SEN,
                SPE = SPE,
                SHD = SHD,
                runtime_sec = runtime_sec
              ))
            }
          }
          
          sim_id <- sim_id + 1
        }
        
        # Save after all reps for this (p, s, n)
        if (nrow(results) > 0) {
          write.csv(results, file = filename, row.names = FALSE)
          message(sprintf("Saved: %s", filename))
        }
      }
    }
  }
  
  message("=== Failures encountered ===")
  print(fail_log)
  
  return(invisible(list(fail_log = fail_log)))
}
