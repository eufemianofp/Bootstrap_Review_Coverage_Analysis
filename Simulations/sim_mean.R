
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
packages <- c("boot")
# packages <- c("boot", "nor1mix")

## Register parallel backend
n_cores <- 20
cluster <- makeCluster(n_cores, outfile = "")
registerDoParallel(cluster)

## Simulation design
interval.types <- c("normal", "basic", "student", "percent", "bca")
ns <- c(10, 100, 1000)
alpha <- 0.05
B <- 1.5e4 -1  # bootstrap replicates to compute CIs
R <- 1e4  # number of simulations for coverage analysis 

## Population parameters
mu <- 0  # normal and bimodal mean
# mu <- 1  # exponential and lognormal mean
sd <- 1  # normal and bimodal sd
# sd <- 0.4  # lognormal sd
# mu <- exp(mu + (sigma^2)/2)  # lognormal mean (log(N(mu, sd)))
# dist <- 3*sd  # distance between std normals for bimodal distribution


##################################
# 2. PERFORM BOOTSTRAP ANALYSIS
##################################

## Define estimator function
estimator_fn <- function(x, index){
  
  x <- x[index]  # bootstrap resample
  
  mu_hat <- mean(x)  # bootstrap estimate
  var_mu_hat <- var(x)/length(x)  # variance of bootstrap estimate
  
  return(c(mu_hat, var_mu_hat))
}

## Compute results (R CIs for different methods and ns)
pb <- txtProgressBar(min = 0, max = R, initial = 0, style = 3)

## Run simulations in parallel
results_parallel <- foreach(i=1:R,
                            .packages = packages,
                            .inorder = FALSE
) %dopar% {
  
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
  
  ## Compute and store bootstrap intervals
  for (j in ns) {
    
    # Compute bootstrap estimates
    boot_rep <- boot(data = x[1:j], estimator_fn, B)
    
    # Compute bootstrap confidence intervals
    ci.boot <- boot.ci(boot_rep, 1-alpha)
    
    # Store bootstrap confidence intervals for n=ns[j]
    j <- as.character(j)
    intervals["normal", j, ] <- ci.boot$normal[c(2, 3)]
    intervals["basic", j, ] <- ci.boot$basic[c(4, 5)]
    intervals["student", j, ] <- ci.boot$student[c(4, 5)]
    intervals["percent", j, ] <- ci.boot$percent[c(4, 5)]
    intervals["bca", j, ] <- ci.boot$bca[c(4, 5)]
  }
  
  ## Update progress bar
  setTxtProgressBar(pb, i)
  
  ## Return seed used and bootstrap intervals to results_parallel
  return(list(i, intervals))
}
stopCluster(cluster)  # close parallel backend
close(pb)  # close progressBar

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
name <- "mean_normal"
saveRDS(bootstrap_intervals, file = paste0("boot_CIs_", name, ".rds"))
write(paste(round(proc.time()[3]/60, 2), "mins"), file = paste0("exec_time_", name, ".txt"))  # write processing time in minutes



