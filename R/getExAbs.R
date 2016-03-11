#' approximates Nu values
#' @param Nu  values for Nu
#' @export 
#' @return values for Nu
getExAbs<-function(Nu){
  approx(c(seq(1,10,0.1),seq(11,100,1),seq(110,1000,10)),approximatedNu,Nu)$y
}

