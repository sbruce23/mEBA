# mEBA (multivariate Empirical Band Analysis)
Simulated examples and electroencephalography analysis accompanying manuscript "Frequency Band Analysis of Nonstationary Multivariate Time Series" by Raanju R. Sundararajan and Scott A. Bruce

R code developed and tested using R version 4.2.2 on both macOS Monterey and Windows 11 operating systems.

## Simulated examples
Separate R scripts are provided to generate simulated time series following the six settings introduced in the manuscript (WN1B, TD1B, L3B, S3B, M3B-1, and M3B-2).  These files use the naming convention `mEBA_simulatedexamples_<abbreviation>.R` where abbreviation is the corresponding abbreviation noted above.  These files also include a description of the output from the method and visualizations of output.  To get started, follow these steps:
1. Download the files into a directory on your local.
2. Open any one of the simulated example R scripts.
3. Update the file path in the `setwd()` function to point to the directory containing the files on your local.
4. Run the entire R script to produce visualizations of the simulated time series, test statistics, bootstrap p-values, and estimated frequency band structure.

## EEG analysis
R script provided to reproduce elements of the EEG data analysis presented in the manuscript (specifically Figures 1, 2b, 2c, and 3).  As in the simulated examples, all that is needed is to update the file path in the `setwd()` function in the R script to point to the directory containing the files on your local.  Then the entire script can be run to reproduce visualizations and results.
