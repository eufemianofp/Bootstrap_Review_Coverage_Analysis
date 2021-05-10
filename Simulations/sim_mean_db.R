
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
# packages <- c("boot")
packages <- c("boot", "nor1mix")

## Register parallel backend
n_cores <- 40
cluster <- makeCluster(n_cores, outfile = "")
registerDoParallel(cluster)

## Simulation design
interval.types <- c("percent_db")
ns <- c(10, 100, 1000)

alpha_goal <- 0.05
alphas_step <- 0.02
steps <- -2:2
alphas <- alpha_goal + alphas_step * steps
coverage_goal <- 1 - alpha_goal

M <- 200
N <- 200  # number of CIs to compute empirical coverage for a given alpha
B <- 1.5e4 -1  # bootstrap replicates to compute CIs
R <- 500  # number of simulations for coverage analysis 

## Population parameters
mu <- 0  # normal mean
sd <- 1  # normal sd
# mu <- 1  # exponential and lognormal mean
# sd <- 0.4  # lognormal sd
# dist <- 3*sd  # distance between std normals for bimodal distribution


##################################
# 2. PERFORM BOOTSTRAP ANALYSIS
##################################




## Define estimator function
estimator_fn <- function(x, index, alpha){
  
  theta_hat <- mean(x)  # estimator value
  x_star <- x[index]  # bootstrap resample
  
  boot_star2 <- boot(data = x_star, 
                statistic = function(x_star, ind) { return(mean(x_star[ind])) }, 
                R = B)  # 2nd-layer bootstrap estimates
  
  ci_star <- boot.ci(boot.out = boot_star2, conf = 1-alpha, type = "perc")  # CI for theta_hat
  
  indicator <- ci_star$percent[4] < theta_hat && ci_star$percent[5] > theta_hat  # does CI for "parameter" theta_hat contain theta_hat?
  
  return(indicator)
}

## Compute results (R CIs for different methods and ns)
pb <- txtProgressBar(min = 0, max = R, initial = 0, style = 3)

## Run simulations in parallel
results_parallel <- foreach(i=1:R,
                            .packages = packages,
                            .inorder = FALSE, 
                            .errorhandling = "pass"
) %dopar% {
  
  boots_warnings <- list()
  
  ## Fix seed for reproducible results
  set.seed(i)
  
  ## X drawn from normal population with mean mu and standard deviation sd
  x <- rnorm(n = max(ns), mu, sd)  # normal
  # x <- rlnorm(n = max(ns), mu, sd)  # lognormal
  # x <- rexp(n = max(ns), rate = mu)  # exponential
  # x <- rt(n = max(ns), df = 5)  # student t
  # mix <- norMix(mu = c(mu - dist/2, mu + dist/2), sigma = c(sd, sd))
  # x <- rnorMix(n = max(ns), mix)  # normal mixture
  
  ## Create array to store bootstrap intervals
  intervals <- array(0, dim = c(length(interval.types), length(ns), 2), 
                     dimnames = list(interval.types, ns, c("lo", "up")))
  
  # Vector of empirical coverages
  emp_coverages <- numeric(length = length(alphas))
  
  ## Compute and store bootstrap intervals
  for (j in ns) {
    
    ## Compute empirical coverages for different values of alpha
    for (k in 1:length(alphas)) {
      
      # Compute bootstrap estimates
      boot_rep <- boot(data = x[1:j], statistic = estimator_fn, R = N, alpha = alphas[k])
      emp_coverages[k] <- mean(boot_rep$t)
    }
    
    ## Choose best value of alpha from vector of alphas
    pos_best_alpha <- which.min(abs(emp_coverages - coverage_goal))
    best_alpha <- alphas[pos_best_alpha]
    
    ## Let's find a better value of best_alpha by linear interpolation
    if (pos_best_alpha != 1 && pos_best_alpha != length(alphas)) {
      
      # Avoid interpolating outside the range of obtained coverages
      if (max(emp_coverages) < coverage_goal || min(emp_coverages) > coverage_goal) {
        alpha_star <- best_alpha
      } else {
        
        diff_to_goal <- emp_coverages[pos_best_alpha] - coverage_goal
        
        if (diff_to_goal <= 0) {
          m <- (emp_coverages[pos_best_alpha - 1] - emp_coverages[pos_best_alpha]) / (-alphas_step)
        } else {
          m <- (emp_coverages[pos_best_alpha + 1] - emp_coverages[pos_best_alpha]) / alphas_step
        }
        
        # Avoiding dividing by zero
        if (m != 0) {
          alpha_star <- best_alpha - 1/m * diff_to_goal
        } else {
          alpha_star <- best_alpha
        }
      }
      
    } else {
      alpha_star <- best_alpha
    }
    
    ## If the interpolation procedure went wrong
    if (alpha_star < min(alphas) || alpha_star > max(alphas)) {
      alpha_star <- best_alpha
    }
    
    ## Compute bootstrap percentile interval
    boot_rep <- boot(data = x[1:j], function(x, ind) {return(mean(x[ind]))}, R = B)
    # ci_boot <- boot.ci(boot.out = boot_rep, conf = 1-alpha_star, type = "perc")
    
    tryCatch(expr = {
      # Compute bootstrap confidence intervals
      ci_boot <- boot.ci(boot.out = boot_rep, conf = 1-alpha_star, type = "perc")
    }, warning = function(w) {
      print(w)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
      
      boots_warnings[[length(boots_warnings) + 1]] <<- boot_rep
    }, error = function(e) {
      print(e)
      print(paste0("Seed: ", i))
      print(paste0("Sample size: ", j))
      print(paste0("Alpha star: ", alpha_star))
    })
    
    
    ## Store bootstrap confidence intervals for n=ns[j]
    j <- as.character(j)
    intervals["percent_db", j, ] <- ci_boot$percent[c(4, 5)]
  }
  
  ## Update progress bar
  setTxtProgressBar(pb, i)
  
  ## Return seed used and bootstrap intervals to results_parallel
  return(list(i, intervals, boots_warnings))
}
stopCluster(cluster)  # close parallel backend
close(pb)  # close progressBar


name <- "mean_bimodal_db"
saveRDS(results_parallel, file = paste0("results_parallel_", name, ".rds"))
write(paste(round(proc.time()[3]/60, 2), "mins"), file = paste0("exec_time_", name, ".txt"))  # write processing time in minutes


## Create array for bootstrap intervals
bootstrap_intervals <- array(dim = c(length(interval.types), length(ns), 2, R),
                             dimnames = list(interval.types,
                                             ns,
                                             c("lo", "up"),
                                             paste0("R", 1:R)))

## Save bootstrap intervals
for (i in 1:R) {
  j <- results_parallel[[i]][[1]]  # seed used (range from 1:R)
  bootstrap_intervals[, , , j] <- results_parallel[[i]][[2]]
}


#####################
# STORE RESULTS
#####################

## Store results
# name <- "mean_bimodal_db"
saveRDS(bootstrap_intervals, file = paste0("boot_CIs_", name, ".rds"))



