#' Compute p-value of two distributions not centered at zero
#' @param dens1 density 1
#' @param dens2 density 2 
#' @param alternative alternative pdf
#' @export
#' @return compared pdf
compareTwoDistsFaster <-function(dens1=runif(256*8,0,1), dens2=runif(256*8,0,1),alternative="two.sided"){
  if(length(dens1)>1 & length(dens2)>1 ){
    dens1<-dens1/sum(dens1)
    dens2<-dens2/sum(dens2)
    cum2 <- cumsum(dens2)-dens2/2
    tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
    if(alternative=="two.sided"){
      if(tmp>0.5)tmp<-tmp-1
      return( tmp*2 )
    }
    if(alternative=="less"){return(1-tmp)}
    return(tmp)
  }
  else {
    return(NA)
  }
}


