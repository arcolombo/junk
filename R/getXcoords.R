#'  Calculates the x-coordinates for the PDF of a given pathway. 
#' @param QSarray qusage object
#' @param path.index path.index can either be an integer between 1 and length(path.means), or the name of the pathway.
#' @param addVIF is a boolean determining whether information on the VIF should be used in the calculation of the x-coordinates.
#' @export 
#' @return coordinates for pdf in pathway
getXcoords = function(QSarray,path.index=1, addVIF=!is.null(QSarray$vif)){ #,absolute=FALSE){
  if(length(path.index)>1){stop("path.index must be of length 1")}
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}

  sif = ifelse(addVIF,sqrt(QSarray$vif[path.index]),1)
  if(is.na(sif)){sif=1}

#   if(!absolute){
  seq(-1,1,length.out=QSarray$n.points)* QSarray$ranges[path.index]* sif + QSarray$path.mean[path.index]
#   }
#    else {
#    ###First calculate the new mean of the pathway based on the absolute values of the means
#    MeanAbs<-mean(abs(QSarray$mean[QSarray$pathways[[path.index]]]))
#    seq(-1,1,length.out=QSarray$n.points)* QSarray$ranges[path.index]* sif / QSarray$path.size[path.index] + MeanAbs
#   }
}

