#' enforced lm.series authored by Gordon Smyth, a master class worker.  this enforces ndupes =1 and spacing = 1 with weigths as NULL, skipping alot of checks, and going into a direct computation over C++ using lm.fit
#' @param M  expression set class object 
#' @param design  a model.matrix formed using formula
#' @import stats
#' @import limma
#' @return a differential expression list
#' @export
enforced.lm.series <- function(M,design=NULL)
        #member variables assumed in an enforced case.
        ndups=1
        spacing=1
        weights=NULL
        NoProbeWts<-TRUE
#       Fit linear model for each gene to a series of arrays
#       Gordon Smyth
#       18 Apr 2002. Revised 26 June 2015.
{
#       Check expression matrix
        M <- as.matrix(M)
        narrays <- ncol(M)

#       Check design
        if(is.null(design))
                design <- matrix(1,narrays,1)
        else
                design <- as.matrix(design)
        nbeta <- ncol(design)
        coef.names <- colnames(design)
        if(is.null(coef.names)) coef.names <- paste("x",1:nbeta,sep="")

#       Check weights : assumed NULL is.na(weights) is TRUE
      #  if(!is.null(weights)) {
       #         weights <- asMatrixWeights(weights,dim(M))
        #        weights[weights <= 0] <- NA
         #       M[!is.finite(weights)] <- NA
       # }

#       Reform duplicated rows into columns : dupes not tolerated
       # if(ndups>1) {
        #        M <- unwrapdups(M,ndups=ndups,spacing=spacing)
         #       design <- design %x% rep(1,ndups)
  # s.matrix(M)
        narrays <- ncol(M)     if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
# }
 # Initialize standard errors
        ngenes <- nrow(M)
        stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))

#       Check whether QR-decomposition is constant for all genes
#       If so, fit all genes in one sweep
        #NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights,"arrayweights")))
       # if(NoProbeWts) {
               # if(is.null(weights)) weights is defaulted to NULL
                        fit <- stats::lm.fit(design, t(M))
                #else {
                 #       fit <- lm.wfit(design, t(M), weights[1,])
                  #      fit$weights <- NULL
               # }
                if(fit$df.residual>0) {
                        if(is.matrix(fit$effects))
                                fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays,,drop=FALSE]^2))
                        else
                                fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
                } else
                        fit$sigma <- rep(NA,ngenes)
                


                #FIX ME : this must be conducted in C++ 
                fit$fitted.values <- fit$residuals <- fit$effects <- NULL
                fit$coefficients <- t(fit$coefficients)
                fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
                est <- fit$qr$pivot[1:fit$qr$rank]
                dimnames(fit$cov.coefficients) <- list(coef.names[est],coef.names[est])
 stdev.unscaled[,est] <- matrix(sqrt(diag(fit$cov.coefficients)),ngenes,fit$qr$rank,byrow = TRUE)
                fit$stdev.unscaled <- stdev.unscaled
                fit$df.residual <- rep.int(fit$df.residual,ngenes)
                dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
                fit$pivot <- fit$qr$pivot
                return(fit)
        }

