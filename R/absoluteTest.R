#' absolute test homogeneity score 
#' @param eset expression set matrix
#' @param QSarray qusage object
#' @param p.adjust.method method for correcting falses
#' @param silent verbose
#' @export
#' @return homogeneity score
absoluteTest = function(eset, QSarray, p.adjust.method="fdr", silent=F){
  cat("Calculating Individual Gene Pvals.")
  gene.pvals = absoluteTest.genePvals(QSarray, compareTo="z",silent=silent)
  cat("Done.\nCalculating VIF...")
  corMats = calcPCor(eset, QSarray)
  cat("Done.\n")

  ##aggregate p-values 
  PVALS<-(sapply(1:length(gene.pvals),function(i){
    Chi2<-(-2)*sum(log(abs(gene.pvals[[i]])))
    k = length(gene.pvals[[i]])
    corMat = corMats[[i]]
    COR = 4*k+sum(ifelse(corMat > 0,corMat*(3.25 + 0.75*corMat),corMat*(3.27 + 0.71*corMat) ))
    f = 8*(k)^2/COR
    c = COR/4/k
    Chi2<-Chi2/c
    P<-1-pchisq(Chi2,f)
  }))

  return(newQSarray(QSarray, absolute.p = PVALS))
}

