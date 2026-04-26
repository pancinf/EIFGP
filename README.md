#  General overview
Extended isolation forest based gene prioritization (EIFGP) is a flexible approach for gene prioritization using unsupervised machine learning. It is designed to integrate association statistics from cohort studies, tissue/cell type specificity values, and random walk with restart probabilities. However additional datatypes can be provided aswell.
The underlying framework is an extended isolation forest (Hariri et al. 2018). Furthermore the script computes Approximate Shapley Values (https://cran.r-project.org/web/packages/fastshap/index.html) for the top n predictions to obtain feature explainability.

##  Required R-packages (R v. 4.3.0)
*  data.table v. 1.17.0
*  optparse v. 1.7.5
*  isotree v. 0.6.1.4
*  fastshap v. 0.1.1
*  ggplot2 v. 3.5.2

##  Required inputs
*  inData.txt. This is a file providing information about the locations of the input files, which columns to consider, their names, if gene symbol/and or ensemblID are given, and if log transformation should be applied for the feature prior to EIFGP. An example can be found in the inLists folder.

##  How to run:
Rscript --vanilla EIFGP.R --dataList ../data/inLists/inData.txt --shapNo 20 --outName myResults

##  Output
*  Extended isolation forest results
*  Approximate Shapley values for top n predictions
