#' computes the 95% CI for a pdf
#' @param QSarray qusage object
#' @param low  percentile
#' @param up  upper percentile
#' @param addVIF vif to add
#' @export 
#' @return CI for pdf
calcBayesCI <- function(QSarray,low=0.025,up=1-low,addVIF=!is.null(QSarray$vif)){
  cis = sapply(1:ncol(QSarray$path.PDF), function(i){
    if( (!is.null(QSarray$pathways) && length(QSarray$pathways[[i]])==0 ) ||
        any(is.na(QSarray$path.PDF[,i]))){return(c(NA,NA))}
    x = getXcoords(QSarray,i,addVIF=addVIF)
    cdf = cumsum(QSarray$path.PDF[,i])
    cdf = cdf/cdf[length(cdf)]
    INDEX_LOW<-findInterval(low,cdf)
    INDEX_UP<-findInterval(up,cdf)
    return( c(  x[INDEX_LOW]+ ((low-cdf[INDEX_LOW])/(cdf[INDEX_LOW+1]-cdf[INDEX_LOW]))*(x[INDEX_LOW+1]-x[INDEX_LOW]) ,
                x[INDEX_UP] + ((up-cdf[INDEX_UP])/(cdf[INDEX_UP+1]-cdf[INDEX_UP]))*(x[INDEX_UP+1]-x[INDEX_UP])
            ) )
#     return( c(x[findInterval(low,cdf)-1] , x[findInterval(up,cdf)]) )
  })
  colnames(cis) = colnames(QSarray$path.PDF)
  rownames(cis) = c("low","up")
  return(cis)
}

