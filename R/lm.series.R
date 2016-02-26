lm.Series <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL)
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

#       Check weights
        if(!is.null(weights)) {
                weights <- asMatrixWeights(weights,dim(M))
                weights[weights <= 0] <- NA
                M[!is.finite(weights)] <- NA
        }
 #Reform duplicated rows into columns
        if(ndups>1) {
                M <- unwrapdups(M,ndups=ndups,spacing=spacing)
                design <- design %x% rep(1,ndups)
                if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
        }

#       Initialize standard errors
        ngenes <- nrow(M)
        stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))

#       Check whether QR-decomposition is constant for all genes
#       If so, fit all genes in one sweep
        NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights,"arrayweights")))
        if(NoProbeWts) {
                if(is.null(weights))
                        fit <- stats::lm.fit(design, t(M))
                else {
                        fit <- lm.wfit(design, t(M), weights[1,])
                        fit$weights <- NULL
                }
                if(fit$df.residual>0) {
                        if(is.matrix(fit$effects))
                                fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays,,drop=FALSE]^2))
                        else
                                fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
                } else
  fit$sigma <- rep(NA,ngenes)
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

#       Genewise QR-decompositions are required, so iterate through genes
        beta <- stdev.unscaled
        sigma <- rep(NA,ngenes)
        df.residual <- rep(0,ngenes)
        for (i in 1:ngenes) {
                y <- as.vector(M[i,])
                obs <- is.finite(y)
                if(sum(obs) > 0) {
                        X <- design[obs,,drop=FALSE]
                        y <- y[obs]
                        if(is.null(weights))
                                out <- lm.fit(X,y)
                        else {
  w <- as.vector(weights[i,obs])
                                out <- lm.wfit(X,y,w)
                        }
                        est <- !is.na(out$coef)
                        beta[i,] <- out$coef
                        stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
                        df.residual[i] <- out$df.residual
                        if(df.residual[i] > 0) sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
                }
        }

#       Correlation matrix of coefficients
        QR <- qr(design)
        cov.coef <- chol2inv(QR$qr,size=QR$rank)
        est <- QR$pivot[1:QR$rank]
        dimnames(cov.coef) <- list(coef.names[est],coef.names[est])

        list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot,rank=QR$rank)
}

