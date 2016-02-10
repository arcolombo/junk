#' homogenity score
#' @param Qsarray   qusage object container
#' @export
#' @return score of homogeneity
homogeneityScore = function(QSarray){
  cat("Calculating Individual Gene Pvals.")
  gene.pvals = absoluteTest.genePvals(QSarray, compareTo="p")
  cat("Done.\n")

  ##aggregate p-values 
  PVALS<-(sapply(1:length(gene.pvals),function(i){
    Chi2<-(-2)*sum(log(abs(gene.pvals[[i]])))
    k = length(gene.pvals[[i]])
    P<-1-pchisq(Chi2,2*k)
  }))

  homogeneity = 1/(1-log10(PVALS))

  return(newQSarray(QSarray, homogeneity = homogeneity))
}

