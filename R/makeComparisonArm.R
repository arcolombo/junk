#' A wrapper function to calculate comparisons between groups in a dataset. A note on var.equal: LIMMA's linear model function can only be run when assuming equal  variances. If var.equal==TRUE, then a linear model will be created on the entire dataset at once. Otherwise, calcIndividualExpressions will be called for the comparison of interest, and means and SDs will be calculated directly.One benefit of using LIMMA's pooled variance calculation is that the linear models allow for more complicated comparisons (e.g. "(A+B)-C" or similar). This may be of interest to some users, but in order to do this, you must assume equal variances between all groups. One caveat regarding paired samples: LIMMA can not fit a linear model when the paired samples are convoluted with the groups (e.g. one set of paired (trt vs mock) samples in patients with  disease, combined with a set of paired samples from healthy controls). If var.equal==TRUE, these groups must be run separately to correctly fit the model (e.g. run disease first, then healthy controls).
#' @param eset a matrix of log2(expression values), with rows of features and columns of samples 
#' @param labels vector of labels representing each column of eset.
#' @param contrast a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'.
#' @param pairVector A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
#' @param bayesEstimation if true, use a bayesian framework to estimate the standard deviation (via limma's eBayes function). 
#' @param min.variance.factor a factor to add to the SDs to ensure that none are equal to 0. Only used if var.equal==FALSE or bayesEstimation==FALSE.
#' @param implicitContrast boolean to explicitly contrast the design fit using an intercept constant, this increases for performance by minimizing calls to limma
#' @import limma
#' @export
#' @return qusage data object
makeComparisonArm <- function(eset,       ##a matrix of log2(expression values), with rows of features and columns of samples 
                           labels,     ##vector of labels representing each column of eset.
                           contrast,   ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'. 
                           pairVector=NULL,  ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                           var.equal = FALSE, ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used. 
                           bayesEstimation = TRUE, ##if true, use a bayesian framework to estimate the standard deviation (via limma's eBayes function). 
                           min.variance.factor=10^-8  ##a factor to add to the SDs to ensure that none are equal to 0. Only used if var.equal==FALSE or bayesEstimation==FALSE. 
){
  if(is(eset, "ExpressionSet")){eset = exprs(eset)}
  ##check that input is formatted correctly
  if(length(labels)!=ncol(eset)){stop("labels length does not match columns of eset")}
  labels = as.factor(as.vector(labels))
  if(!is.character(contrast)){stop("Contrast must be a character vector of length 1.")}
  if(length(contrast)!=1){
    warning("Multiple contrasts provided. Using first contrast only.")
    contrast = contrast[1]
  }
  if(is.null(rownames(eset))){
    stop("Rownames for eset not found")
  }
  if(length(unique(rownames(eset)))!=nrow(eset) | any(rownames(eset)=="")){
    stop("The rownames of eset are invalid. Rownames must be unique and must not contain any empty values")
  }


  params = list(labels=labels, contrast = contrast)

  if((paired = !is.null(pairVector))){
    if(length(pairVector)!=ncol(eset)){stop("PairVector length does not match columns of eset")}
    pairVector = as.factor(as.vector(pairVector))
    params[["pairVector"]] = pairVector
  }


  if(var.equal){
    ###################################
    ## Pooled Variance (Linear Model) method

    ##create design matrix
     f = "~0+labels"

    if( class(labels)!="factor" ){
    designNames = levels(factor(labels))
      }
    if( class(labels)=="factor") {
     designNames<-levels(labels) 
     }

    if(paired){
      f = paste(f,"+pairVector",sep="")
  designNames = c(designNames, paste("P",levels(pairVector)[-1],sep=""))
    }
    design <- model.matrix(formula(f))
    colnames(design) <- designNames

    ## Fit the linear model with the given deign matrix
    ## 'fit' contains info on each coefficient (i.e. column of design matrix) in the model.
  
    fit<-lmFit(eset,design=design)
    contrast.matrix <- makeContrasts( contrasts=contrast, levels=design)
    fit2 <- contrasts.fit(fit,contrast.matrix)

    ##calculate number of samples use in the contrast
    n.samples = sum(labels %in% rownames(contrast.matrix)[contrast.matrix!=0])
    
    ##if using Bayes estimation, calculate the moderated t-statistics for each comparison
   #FIX ME: Armadillo can perform these comps 
   if(bayesEstimation){
      fit2b <- eBayes(fit2)
     # SD = (fit2b$coefficients/(fit2b$t))[,1]
  #sd.alpha = SD/(fit2b$sigma*fit2b$stdev.unscaled)
  #    sd.alpha[is.infinite(sd.alpha)] = 1
   #   dof = fit2b$df.total
    arrayResult<-bayesEstimation(fit2b$coefficients,fit2b$t,fit2b$sigma,fit2b$stdev.unscaled,fit2b$df.total)
    }else{
     # SD = sqrt((fit2$sigma*fit2b$stdev.unscaled)^2 + min.variance.factor)
     # sd.alpha = SD/(fit2$sigma*fit2$stdev.unscaled)
      #sd.alpha[is.infinite(sd.alpha)] = 1
     # dof = fit2$df.residual
     arrayResult<-notbayesEstimation(fit2$sigma,fit2b$stdev.unscaled,min.variance.factor,as.matrix(fit2$df.residual))
    }

    #FIX ME : convert to vectors, and name everything
    ##format 
    results = newQSarray(params,
                      mean = fit2$coefficients[,1],
                      SD = SD,
                      sd.alpha = sd.alpha,
                      dof = dof,
                      var.method="Pooled",
                      n.samples=n.samples
                     )
  }
 if(!var.equal){
    ###################################
    ## Welch's method
    ##parse the contrast
    grps = strsplit(contrast,"-")[[1]]
    grps = sub("\\s","",grps)          ##remove whitespace
    if(length(grps)!=2){stop("Only contrasts of the form 'A-B' are allowed when var.equal is FALSE.")}
    grp.1 = labels==grps[1]  ##PostTreatment
    grp.2 = labels==grps[2]  ##Baseline
    if(sum(grp.1)==0 | sum(grp.2)==0){stop("Contrast groups do not match labels")}
    params$n.samples = sum(grp.1) + sum(grp.2)

    eset.1 = eset[,grp.1]
    eset.2 = eset[,grp.2]

    if(paired){
      colnames(eset.1) = pairVector[grp.1]
      colnames(eset.2) = pairVector[grp.2]

      #if(ncol(eset.1)!=ncol(eset.2)){
  eset.1 = eset.1[,colnames(eset.1) %in% colnames(eset.2)]
        eset.2 = eset.2[,colnames(eset.2) %in% colnames(eset.1)]
        eset.1 = eset.1[,match(colnames(eset.2),colnames(eset.1))]
      #}
    }

    results = newQSarray(c(params, calcIndividualExpressionsArm(eset.2,eset.1,paired=paired,min.variance.factor=min.variance.factor)))
  }
  return(results)
}

