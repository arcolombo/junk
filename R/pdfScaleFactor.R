#' calculate a scaling factor to multiply the PDF by (so that it actually integrates to 1)
#' @param QSarray qusage object
#' @param addVIF vif value
#' @export 
#' @return  pdf plot
pdfScaleFactor = function(QSarray, addVIF=!is.null(QSarray$vif)){
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = sapply(1:numPathways(QSarray),function(i){ifelse(addVIF,sqrt(QSarray$vif[i]),1)})
  sif[is.na(sif)] = 1
  pdfSum = colSums(QSarray$path.PDF)

  ##get the pdf range
  ranges = QSarray$ranges * 2

  ##the scale factor is essentially the distance between points in the x coordinates times the (current) sum of the pdf
  scaleFactor = (ranges*sif) / (QSarray$n.points-1) * pdfSum
  scaleFactor
}

