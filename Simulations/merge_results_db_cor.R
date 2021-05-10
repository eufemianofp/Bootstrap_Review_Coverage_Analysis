
## Remove previous objects in the R session
rm(list = ls())



####################################
##  Setup
####################################

path <- file.path("/home", "eufe", "master_thesis")
setwd(file.path(path, "results_plots_cor"))

estimators <- c("cor0.5", "cor0.9")
populations <- c("normal", "t5")

interval_types <- c("normal", "basic", "student", "percent", "bca", "percent_db")
ns <- c(10, 100, 1000)
R <- 1000



##########################################
##  Function for plot of noncoverages
##########################################

plot_cis <- function (noncoverage, median_lengths, right, y_lims) {
  
  ## Basic plot
  ns_plot <- c(2,4,6)
  
  ylab = "One-sided noncoverage"
  toplab <- "Left"
  par(mar=c(4, 4, 2, 0))
  if (right == TRUE) {
    ylab = ""
    toplab <- "Right"
    par(mar=c(4, 3, 2, 0.5))
  }
  
  plot.new()
  plot.window(xlim = c(1.25, 6.75), ylim = y_lims)
  axis(1, at = ns_plot, labels = ns)  # set x axis values
  ticks <- sort(c(0.025, seq(0, 0.5, by = 0.01)))
  axis(2, at = ticks)
  box()
  title(xlab = "Sample size n", ylab = ylab, main = toplab, cex.main = 0.9)  # font.main = 1 to remove bold main title
  abline(h = 0.025, col = "black", lty = 2)
  
  ## Interpolation function to get proportional widths according to value of median lengths
  min_medians <- apply(median_lengths, 2, min)
  max_medians <- apply(median_lengths, 2, max)
  min_max_medians <- matrix(c(min_medians, max_medians), ncol = 3, byrow = TRUE)
  min_max_proportions <- matrix(rep(c(0,1), each = 3), ncol = 3, byrow = TRUE)
  interpolate <- approxfun(x = min_max_medians, y = min_max_proportions)
  
  ## Get widths (deltas)
  prop <- t(apply(median_lengths, 1, interpolate))
  min_width <- 0.2
  max_width <- 0.6
  deltas <- min_width + (max_width - min_width) * prop
  
  ## Create x and y coordinates for segments function, to plot horizontal intervals
  ns_deltas <- matrix(rep(ns_plot, times = 6), ncol = 3, byrow = TRUE)
  from.x <- ns_deltas - deltas
  to.x   <- ns_deltas + deltas
  to.y   <- from.y <- noncoverage  # from and to y coords the same
  
  ## Plot horizontal lines
  colours <- c("red3", "green4", "blue3", "cyan", "magenta", "darkorange3")
  segments(x0 = from.x, y0 = from.y, x1 = to.x, y1 = to.y, 
           col = rep(colours, times = 3), lwd = 1.5)
  
  interval_codes <- c("zB-BC", "R", "BT", "P", "BCa", "DB")
  matplot(x = t(ns_deltas), 
          y = t(noncoverage), 
          type = rep("b", 6), lty = 3, pch = 4, col = colours, add = TRUE, cex = 1.4)
  
  ## Add legend
  if (right == TRUE) {
    legend(x = 4.9, y = y_lims[2], legend = interval_codes,
           col = colours, lty = 1, lwd = 3, cex = 0.8)
  }
}



###################################################
##  Loop through all estimators and populations, 
##  and create the plots
###################################################

for (estimator in estimators) {
  for (population in populations) {
    
    # estimator <- estimators[1]
    # population <- populations[2]
    estpop <- paste0(estimator, "_", population)
    
    
    ## Set true parameter value and skip variable to skip warnings
    if (estimator == "cor0.5") {
      true_parameter <- 0.5
      if (population == "normal") {
        skip <- NULL
      } else {  # t5
        skip <- c(301, 623, 784)
      }
    } else {  # cor0.9
      true_parameter <- 0.9
      if (population == "normal") {
        skip <- c(301)
      } else {  # t5
        skip <- c(194, 239, 301, 605, 623, 836, 960, 1034)
      }
    }
    
    
    ## Get R seeds with no warnings
    seeds <- 1:R
    if (!is.null(skip)) {
      seeds <- 1:(R + length(skip))
      seeds <- seeds[-skip]
    }
    
    ## Read intervals
    res_cis_path <- file.path(path, "sim_cor", estpop, paste0("results_", estpop, ".rds"))
    res_cis_db_path <- file.path(path, "sim_cor", estpop, paste0("results_parallel_", estpop, "_db.rds"))
    res_cis_db_atanh_path <- file.path(path, "sim_cor", estpop, paste0("results_parallel_", estpop, "_atanh_db.rds"))
    
    res_cis <- readRDS(res_cis_path)
    res_cis_db <- readRDS(res_cis_db_path)
    res_cis_db_atanh <- readRDS(res_cis_db_atanh_path)
    
    
    ## If there were warnings, add the 10 extra simulated for double bootstrap
    if (!is.null(skip)) {
      res_cis_db_extra_path <- file.path(path, "sim_cor", estpop, paste0("results_parallel_", estpop, "_db_extra.rds"))
      res_cis_db_atanh_extra_path <- file.path(path, "sim_cor", estpop, paste0("results_parallel_", estpop, "_atanh_db_extra.rds"))
      
      res_cis_db_extra <- readRDS(res_cis_db_extra_path)
      res_cis_db_atanh_extra <- readRDS(res_cis_db_atanh_extra_path)
      
      res_cis_db <- c(res_cis_db, res_cis_db_extra)
      res_cis_db_atanh <- c(res_cis_db_atanh, res_cis_db_atanh_extra)
    }
    
    
    ## Create array for bootstrap intervals
    boot_cis <- array(data = NA, 
                      dim = c(2, length(interval_types), length(ns), 2, R),
                      dimnames = list(c("untransformed", "transformed"), 
                                      interval_types,
                                      ns,
                                      c("lo", "up"),
                                      paste0("R", 1:R)))
    
    
    ## Rearrange untransformed and transformed standard bootstrap intervals
    for (i in 1:R) {
      j <- seeds[i]
      boot_cis[, interval_types[1:5], , , i] <- res_cis[[j]][[2]]
    }
    
    # Add untransformed double bootstrap intervals
    for (i in 1:(R/2)) {
      j <- seeds[i]
      boot_cis["untransformed", interval_types[6], , , i] <- res_cis_db[[j]][[2]]
    }
    
    # Add transformed double bootstrap intervals, and transform it back using tanh()
    for (i in 1:(R/2)) {
      j <- seeds[i]
      boot_cis["transformed", interval_types[6], , , i] <- tanh(res_cis_db_atanh[[j]][[2]])
    }
    
    
    ## Create matrices to store noncoverage and length results
    empty_matrix <- matrix(nrow = length(interval_types), 
                           ncol = length(ns), 
                           dimnames = list(interval_types, ns))
    noncoverage_left <- noncoverage_right <- median_lengths <- empty_matrix
    noncoverage_left_atanh <- noncoverage_right_atanh <- median_lengths_atanh <- empty_matrix
    
    
    
    ## -------------------------------
    ## Untransformed correlation
    ## -------------------------------
    
    ## Compute noncoverage to the left and to the right sides of CI endpoints and median length of CIs
    for (i in interval_types) {
      for (j in ns) {
        j <- as.character(j)
        
        lower_endpoints <- boot_cis["untransformed", i, j, "lo", ]
        upper_endpoints <- boot_cis["untransformed", i, j, "up", ]
        
        if (i == "percent_db") {
          lower_endpoints <- lower_endpoints[!is.na(lower_endpoints)]
          upper_endpoints <- upper_endpoints[!is.na(upper_endpoints)]
        }
        
        ci_lengths <- abs(upper_endpoints - lower_endpoints)
        median_lengths[i, j] <- median(ci_lengths)
        
        noncvg_left <- lower_endpoints > true_parameter
        noncoverage_left[i, j] <- mean(noncvg_left)
        
        noncvg_right <- upper_endpoints < true_parameter
        noncoverage_right[i, j] <- mean(noncvg_right)
      }
    }
    
    
    ## Save plots as png and pdf
    y_min <- floor(min(noncoverage_left, noncoverage_right)*1000)/1000
    y_max <- ceiling(max(noncoverage_left, noncoverage_right)*1000)/1000
    y_lims <- c(y_min, y_max)
    
    # png
    png(filename = paste0(estpop, ".png"), 
        width = 180*3.7, height = 135*3.7, pointsize = 15)
    par(mfrow = c(1, 2))
    plot_cis(noncoverage_left, median_lengths, right = FALSE, y_lims)
    plot_cis(noncoverage_right, median_lengths, right = TRUE, y_lims)
    par(mfrow = c(1, 1))
    dev.off()
    
    # pdf
    pdf(file = paste0(estpop, ".pdf"), 
        width = 18/2.54, height = 13.5/2.54, pointsize = 10)
    par(mfrow = c(1, 2))
    plot_cis(noncoverage_left, median_lengths, right = FALSE, y_lims)
    plot_cis(noncoverage_right, median_lengths, right = TRUE, y_lims)
    par(mfrow = c(1, 1))
    dev.off()
    
    
    ## Save table of noncoverages and median lengths
    merged <- round(cbind(noncoverage_left, noncoverage_right, median_lengths), digits = 4)
    rownames(merged) <- c("zB-BC", "R", "BT", "P", "BCa", "DB")
    merged <- as.data.frame(merged)
    merged <- cbind(merged, rep("\\\\[0.3ex]", 6))
    
    write.table(format(merged, digits=4, nsmall=4), file = paste0(estpop, "_table.csv"), 
                col.names = FALSE, quote = FALSE, sep = " & ")
    
    
    ## Get ranks for given estimator and population
    errors_left <- round(abs(noncoverage_left - 0.025), digits = 4)
    errors_right <- round(abs(noncoverage_right - 0.025), digits = 4)
    
    ranks_left <- apply(errors_left, 2, rank)
    ranks_right <- apply(errors_right, 2, rank)
    
    sum_ranks <- apply(cbind(ranks_left, ranks_right), 1, sum)
    ranks_total <- rank(sum_ranks)
    
    write.csv(ranks_total, file = paste0(estpop, "_ranks.csv"),
              quote = FALSE)
    
    
    
    ## -------------------------------
    ## Transformed correlation
    ## -------------------------------
    
    ## Compute noncoverage to the left and to the right sides of CI endpoints and median length of CIs
    for (i in interval_types) {
      for (j in ns) {
        j <- as.character(j)
        
        lower_endpoints <- boot_cis["transformed", i, j, "lo", ]
        upper_endpoints <- boot_cis["transformed", i, j, "up", ]
        
        if (i == "percent_db") {
          lower_endpoints <- lower_endpoints[!is.na(lower_endpoints)]
          upper_endpoints <- upper_endpoints[!is.na(upper_endpoints)]
        }
        
        ci_lengths <- abs(upper_endpoints - lower_endpoints)
        median_lengths_atanh[i, j] <- median(ci_lengths)
        
        noncvg_left <- lower_endpoints > true_parameter
        noncoverage_left_atanh[i, j] <- mean(noncvg_left)
        
        noncvg_right <- upper_endpoints < true_parameter
        noncoverage_right_atanh[i, j] <- mean(noncvg_right)
      }
    }
    
    
    ## Save plots as png and pdf
    y_min_atanh <- floor(min(noncoverage_left_atanh, noncoverage_right_atanh)*1000)/1000
    y_max_atanh <- ceiling(max(noncoverage_left_atanh, noncoverage_right_atanh)*1000)/1000
    y_lims_atanh <- c(y_min_atanh, y_max_atanh)
    
    # png
    png(filename = paste0(estpop, "_atanh.png"), 
        width = 180*3.7, height = 135*3.7, pointsize = 15)
    par(mfrow = c(1, 2))
    plot_cis(noncoverage_left_atanh, median_lengths_atanh, right = FALSE, y_lims_atanh)
    plot_cis(noncoverage_right_atanh, median_lengths_atanh, right = TRUE, y_lims_atanh)
    par(mfrow = c(1, 1))
    dev.off()
    
    # pdf
    pdf(file = paste0(estpop, "_atanh.pdf"), 
        width = 18/2.54, height = 13.5/2.54, pointsize = 10)
    par(mfrow = c(1, 2))
    plot_cis(noncoverage_left_atanh, median_lengths_atanh, right = FALSE, y_lims_atanh)
    plot_cis(noncoverage_right_atanh, median_lengths_atanh, right = TRUE, y_lims_atanh)
    par(mfrow = c(1, 1))
    dev.off()
    
    
    ## Save table of noncoverages and median lengths
    merged_atanh <- round(cbind(noncoverage_left_atanh, noncoverage_right_atanh, median_lengths_atanh), digits = 4)
    rownames(merged_atanh) <- c("zB-BC", "R", "BT", "P", "BCa", "DB")
    merged_atanh <- as.data.frame(merged_atanh)
    merged_atanh <- cbind(merged_atanh, rep("\\\\[0.3ex]", 6))
    
    write.table(format(merged_atanh, digits=4, nsmall=4), file = paste0(estpop, "_atanh_table.csv"), 
                col.names = FALSE, quote = FALSE, sep = " & ")
    
    
    ## Get ranks for given estimator and population
    errors_left_atanh <- round(abs(noncoverage_left_atanh - 0.025), digits = 4)
    errors_right_atanh <- round(abs(noncoverage_right_atanh - 0.025), digits = 4)
    
    ranks_left_atanh <- apply(errors_left_atanh, 2, rank)
    ranks_right_atanh <- apply(errors_right_atanh, 2, rank)
    
    sum_ranks_atanh <- apply(cbind(ranks_left_atanh, ranks_right_atanh), 1, sum)
    ranks_total_atanh <- rank(sum_ranks_atanh)
    
    write.csv(ranks_total_atanh, file = paste0(estpop, "_atanh_ranks.csv"),
              quote = FALSE)
  }
}






