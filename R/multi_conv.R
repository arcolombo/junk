#' fast fourier transform method
#' @param x input to fft
#' @export
#' @return fft object
multi_conv<-function(x){
  x_fft<-apply(x,2,function(x)fft(x, inverse = FALSE))
  M<-max(Mod(x_fft))
  x_fft<-x_fft/M
  Prod_fft<-apply(x_fft,1,prod)
  p1<-Re(fft(Prod_fft,inverse=TRUE))
  N<-nrow(x)
  #     if(ncol(x)%%2==0)p1<-c(p1[(N/2+1):N],p1[1:(N/2)])
  Mp1<-which.max(p1)
  Delta<-N/2-Mp1
  if(Delta>0){
    p1<-c(p1[(N-Delta):N],p1[1:(N-Delta-1)])
  }
  if(Delta<0){
    p1<-c(p1[(1-Delta):N],p1[1:(1-Delta-1)])
  }
  p2 = abs(p1)
  return(p2/sum(p2))
}

