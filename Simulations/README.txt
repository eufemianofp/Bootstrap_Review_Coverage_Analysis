The files to run the simulations are:
	- sim_mean.R
	- sim_mean_db.R (double bootstrap)
	- sim_median.R 
	- sim_median_db.R (double bootstrap)
	- sim_cor.R
	- sim_cor_db.R (double bootstrap)

The results of the first two files are stored in folder "sim_mean", the following two in folder "sim_median" and the last two in folder "sim_cor".

The files to generate the plots and tables in Chapter 4 are:
	- merge_results_db.R
	- merge_results_db_cor.R

The results are stored in folders "results_plots" and "results_plots_cor".

Folder "miscellaneous" contains two R files: 
	- plot_montecarlo_cis.R
	- sample_cor_distr.R

The first R file computes the first two images in Chapter 4 with a 100 confidence intervals with a 95% confidence level, one with symmetric noncoverage and another with asymmetric noncoverage. The second R file computes the plot in Chapter 4 showing how the distribution of the sample correlation coefficient varies with the true correlation coefficient and the sample size.

