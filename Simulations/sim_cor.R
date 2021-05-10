
## Remove previous objects in the R session
rm(list = ls())



########################
# 1. SIMULATION SETUP
########################

## Load libraries
library(parallel)
library(foreach)
library(doParallel)

## Required packages inside foreach loop
packages <- c("boot", "mvtnorm")

## Register parallel backend
n_cores <- 40
cluster <- makeCluster(n_cores, outfile = "")
registerDoParallel(cluster)

## Simulation design
interval.types <- c("normal", "basic", "student", "percent", "bca")
ns <- c(10, 100, 1000)
alpha <- 0.05
M <- 200  # bootstrap replicates to compute standard error of estimator 
B <- 1.5e4 -1  # bootstrap replicates to compute CIs
R <- 1100  # number of simulations for coverage analysis

## Population parameters
mu <- c(0, 0)
cor <- 0.9
Sigma <- matrix(data = c(1, cor, cor, 1), ncol = 2)



##################################
# 2. PERFORM BOOTSTRAP ANALYSIS
##################################

## Define estimator for second layer bootstrap
estimator_layer2 <- function(x_star, ind){
  
  while (length(unique(x_star[ind, 1])) <= 2) {
    ind <- sample(x = 1:length(ind))  # recompute indices
  }
  
  x_double_star <- x_star[ind, ]  # 2nd layer bootstrap resample
  r_double_star <- cor(x_double_star[, 1], x_double_star[, 2])  # compute correlation coefficient
  
  return(r_double_star)
}

## Define estimator function
estimator_fn <- function(x, indices){
  
  ## Take other indices if the same number is picked all the time 
  #  (could happen for n=10) or only 2 numbers are picked, since
  #  in that case r = 1 or r = -1, and hence z = Inf or z = -Inf
  while (length(unique(indices)) <= 2) {
    indices <- sample(x = 1:length(indices))  # recompute indices
  }
  
  x_star <- x[indices, ]  # bootstrap resample
  r_star <- cor(x_star[, 1], x_star[, 2])  # bootstrap estimate
  
  boot_r_star <- boot(data = x_star, statistic = estimator_layer2, M)
  var_r_star <- var(boot_r_star$t)  # variance of bootstrap estimate
  
  return(c(r_star, var_r_star))
}

## Compute results (R CIs for different methods and ns)
pb <- txtProgressBar(min = 0, max = R, initial = 0, style = 3)

## Run simulations in parallel
results_parallel <- foreach(i=1:R,
                            .packages = packages,
                            .inorder = FALSE,
                            .errorhandling = "pass"
) %dopar% {
  
  ## List of warnings thrown when computing CIs
  boot_warnings <- list()
  boot_warnings_h <- list()
  
  ## Fix seed for reproducible results
  set.seed(i)
  
  ## Generate data X
  # x <- rmvnorm(n = max(ns), mean = mu, sigma = Sigma)
  x <- rmvt(n = max(ns), sigma = Sigma, df = 5)
  
  ## Create array to store bootstrap intervals
  intervals <- array(0, dim = c(2, length(interval.types), length(ns), 2), 
                     dimnames = list(c("untransformed", "transformed"), 
                                     interval.types, 
                                     ns, 
                                     c("lo", "up")))
  
  ## Compute and store bootstrap intervals
  for (j in ns) {
    
    # Compute bootstrap estimates
    boot_rep <- boot(data = x[1:j, ], estimator_fn, B)
    
    
    ## Compute bootstrap CIs for the untransformed correlation coefficient
    tryCatch(expr = {
      cis_boot <- boot.ci(boot_rep, conf = 1-alpha)  # Compute bootstrap CIs
      
    }, warning = function(w) {
      print(w)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
      
      boot_warnings[[length(boot_warnings) + 1]] <<- boot_rep
      
    }, error = function(e) {
      print(e)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
    })
    
    # Store bootstrap confidence intervals for n=ns[j]
    j <- as.character(j)
    intervals["untransformed", "normal", j, ] <- cis_boot$normal[c(2, 3)]
    intervals["untransformed", "basic", j, ] <- cis_boot$basic[c(4, 5)]
    intervals["untransformed", "student", j, ] <- cis_boot$student[c(4, 5)]
    intervals["untransformed", "percent", j, ] <- cis_boot$percent[c(4, 5)]
    intervals["untransformed", "bca", j, ] <- cis_boot$bca[c(4, 5)]
    
    
    ## Compute bootstrap CIs for the transformed correlation coefficient
    tryCatch(expr = {
      cis_boot_h <- boot.ci(boot_rep, conf = 1-alpha, h = atanh, hinv = tanh, 
                          hdot = function(x) {1/(x^2 + 1)})  # Compute bootstrap CIs
      
    }, warning = function(w) {
      print(w)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
      
      boot_warnings_h[[length(boot_warnings_h) + 1]] <<- boot_rep
      
    }, error = function(e) {
      print(e)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
    })
    
    # Store bootstrap confidence intervals for n=ns[j]
    intervals["transformed", "normal", j, ] <- cis_boot_h$normal[c(2, 3)]
    intervals["transformed", "basic", j, ] <- cis_boot_h$basic[c(4, 5)]
    intervals["transformed", "student", j, ] <- cis_boot_h$student[c(4, 5)]
    intervals["transformed", "percent", j, ] <- cis_boot_h$percent[c(4, 5)]
    intervals["transformed", "bca", j, ] <- cis_boot_h$bca[c(4, 5)]
  }
  
  ## Update progress bar
  setTxtProgressBar(pb, i)
  
  ## Return seed, bootstrap intervals and warnings to results_parallel
  return(list(i, intervals, boot_warnings, boot_warnings_h))
}
stopCluster(cluster)  # close parallel backend
close(pb)  # close progressBar



#####################
# STORE RESULTS
#####################

## Store results
name <- "cor0.9_t5"
saveRDS(results_parallel, file = paste0("results_", name, ".rds"))
write(paste(round(proc.time()[3]/60, 2), "mins"), file = paste0("exec_time_", name, ".txt"))  # write processing time in minutes



