#  General overview
EIFGPI is a flexible approach for gene prioritization using unsupervised machine learning. It is designed to integrate association statistics from cohort studies, tissue/cell type specificity values and random walk with restart probabilities.
The underlying framework is an extended isolation forest (Hariri et al. 2018). Furthermore the script computes Approximate Shapley Values (https://cran.r-project.org/web/packages/fastshap/index.html) for the top n predictions to obtain feature explainability.

##  Required inputs
*  inData.txt. This is a file providing information about the locations of the input files, which columns to consider, their names, if gene symbol/and or ensemblID are given, and if log transformation should be applied for the feature prior to EIFGPI. An example can be found in the inLists folder.

##  How to run:
Rscript --vanilla eifgpi.R --dataList ../data/inLists/inData.txt --shapNo 20 --outName myResults

##  Output
*  Extended isolaiton forest results
*  Approximate Shapley values for top n predictions
