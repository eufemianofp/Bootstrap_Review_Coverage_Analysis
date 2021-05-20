# Review of bootstrap principles and coverage analysis of bootstrap confidence intervals for common estimators
This project was developed in Autumn 2019 as part of my master thesis for the MSc in Statistics at ETH Zurich.

## Abstract
The bootstrap is a statistical technique that has been around for 40 years since it was
introduced by Efron (1979). Its use in practice is widespread, but many practitioners
do not fully understand its limits and under which circumstances it works or does not
work. This thesis tries to address this issue, first by exploring the theoretical underpinnings
of the bootstrap and then by analysing its performance in some practical scenarios for common estimators.

* __Chapter 1__ starts with an introduction to the bootstrap, its fundamental principles and how it works.
* __Chapter 2__ gets into when the bootstrap works (consistency) and when it does not, and if it works at which rate it does so (accuracy).
* __Chapter 3__ explores a number of first- and second-order accurate bootstrap confidence intervals, introduces a general technique
  to improve bootstrap intervals called double bootstrap and presents the two most well-known R packages to implement the bootstrap.
* __Chapter 4__ illustrates the results from a coverage analysis of bootstrap intervals in different scenarios for the
  sample mean, sample median and sample (Pearson) correlation coefficient, computed via simulations.
* __Chapter 5__ makes a summary of the thesis and presents a list of conclusions.
* __Chapter 6__ outlines interesting points and topics to explore further.

## Code
The code available in this repository has been used to run the simulations for Chapter 4: Coverage
Analysis of Bootstrap CIs for Common Estimators.

## License
MIT License.
