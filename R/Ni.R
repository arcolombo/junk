#' helper funciton to calculate DOF
#' @param SigmaOld first sigma to compare
#' @param SigmaYoung second sigma to compare
#' @param NOld Non-nas of old, ncolumns
#' @param NYoung non-nas of columns of young
#' @export
#' @return degrees of freedom
Ni<-function(SigmaOld,SigmaYoung,NOld,NYoung){
 # (SigmaYoung^2/NYoung+SigmaOld^2/NOld)^2/(((SigmaYoung^4/(NYoung^2*(NYoung-1))))+((SigmaOld^4/(NOld^2*(NOld-1)))))
  (SigmaYoung/NYoung+SigmaOld/NOld)^2/(((SigmaYoung^2/(NYoung^2*(NYoung-1))))+((SigmaOld^2/(NOld^2*(NOld-1)))))
}

