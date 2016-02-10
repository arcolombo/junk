#' qusage single run
#' @param eset a matrix of log2(expression values), with rows of features and columns of samples. OR an object of class ExpressionSet 
#' @param labels vector of labels representing each column of eset.
#' @param contrast a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'. 
#' @param geneSets a list of pathways to be compared. Each item in the list is a vector of names that correspond to the row names of eset.
#' @param pairVector A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
#' @param filter.genes a boolean indicating whether the genes in eset should be filtered to remove genes with low mean and sd.
#' @param n.points The number of points to sample the convolution at. Passed to aggregateGeneSet
#' @export
#' @return qusage results 

qusage.single = function(eset,       ##a matrix of log2(expression values), with rows of features and columns of samples. OR an object of class ExpressionSet 
                         labels,            ##vector of labels representing each column of eset.
                         contrast,          ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'. 
                         geneSets,          ##a list of pathways to be compared. Each item in the list is a vector of names that correspond to the row names of eset.
                         pairVector=NULL,   ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                         var.equal=FALSE,   ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
                         filter.genes=FALSE,##a boolean indicating whether the genes in eset should be filtered to remove genes with low mean and sd.
                         n.points=2^12      ##The number of points to sample the convolution at. Passed to aggregateGeneSet
                 ){
  cat("Calculating gene-by-gene comparisons...")
  results = makeComparison(eset, labels, contrast, pairVector=pairVector,var.equal=var.equal)
  if(filter.genes){
 results = filterGenes(results)
  }

  cat("Done.\nAggregating gene data for gene sets.")
  nu = floor(min(results$dof,na.rm=T))
  if(nu<5){cat("\nLow sample size detected. Increasing n.points in aggregateGeneSet.")}
  results = aggregateGeneSet(results, geneSets, silent=F, n.points=n.points)
  cat("Done.\nCalculating variance inflation factors...")
  results = calcVIF(eset, results)
  #cat("Done.\nCalculating homogeneity scores...")
  #results = calcHomogeneity(results, silent=TRUE, addVIF=FALSE)
  #cat("Done.\nCalculating correlation matrix...")  
  #results = calcPCor(eset, results)
  cat("Done.\n")
  results
}



