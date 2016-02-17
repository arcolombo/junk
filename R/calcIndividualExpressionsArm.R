#' modificed expression calculation that enforces clean data without NAs for speeding up qusage
#' @param Baseline               expression values before treatment
#' @param PostTreatment          expression values after treatment
#' @param paired                 boolean, paired reads
#' @param min.variance.factor    numeric, variance factor for pooled variance
#' @import parallel
#' @export 
calcIndividualExpressionsArm<-function(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6){
  ###Baseline is the matix of gene expressions at baseline, row names are gene names
  ###PostTreatment is the matix of gene expressions after treatment, row names are gene names
  ###paired: logical, whether the data is paired or not
#enforce no NAs  not a flexible solution
  ##########Some error checks
  if(length(dim(Baseline))!=2 | length(dim(PostTreatment))!=2){stop("Input Matrices need to be matrices... \n")}
  if(nrow(Baseline)!=nrow(PostTreatment)){stop("Input Matrices need to have the same number of genes \n")}
  if(sum(!(rownames(Baseline)%in%rownames(PostTreatment)))){stop("Input Matrices need to have the same list of genes. Gene names should be the row names \n")}
  if(ncol(Baseline)<2 | ncol(PostTreatment)<2 ){stop("Input Matrices need to have at least two columns \n")}
 
 #the design could be improved: by checking if there are any NA values and ensureing there do not exist any NA values, then we can say that the row sums of the NA values is equal to the number of columns. and do not have to check that. assuming this
 if(any(is.na(Baseline)) || any(is.na(PostTreatment)) ){
 stop("NA values are present in the expression matrix, please pluck out ...")
 }


 #########Reorder PostTreatment
  PostTreatment<-PostTreatment[rownames(Baseline),]

  ###########Paired
    if(paired){
      if(ncol(Baseline)!=ncol(PostTreatment)){
        stop("Input Matrices need to have the same number of columns when paired flag is on \n")
      }
  if(sum(!colnames(Baseline)%in%colnames(PostTreatment))){
        stop("Input Matrices need to have the same list of samples when paired flag is on \n")
      }
      PostTreatment = PostTreatment[,colnames(Baseline)]
      ##########First calculate the differential expression for individual genes
#      Sums_Base<-rowSums(Baseline)
      # Sums_Post<-rowSums(PostTreatment)
     # Ns<-ncol(Baseline-PostTreatment)#the ncol will be equal to all of the remove NA
      #if(min(Ns)!=ncol(Baseline)){warning("Some NA's in data")} ignore by assumption
     Ns<-ncol(Baseline) #ncol identical(ncol(Baseline,ncol(PostTreatment)) if NA enforced 
  #  Sigmas_Base<-rowSums((Baseline-PostTreatment-(Sums_Base-Sums_Post)/Ns)^2)/(Ns-1)
    qsList<-sigmaArm(Baseline,PostTreatment,Ns,min.variance.factor)
    qv<-lapply(qsList,function(x) as.vector(x))
    names(qv$Mean)<-attr(qsList$Mean,"names")
    names(qv$SD)<-attr(qsList$SD,"names")
     names(qv$SDAlpha)<-attr(qsList$SDAlpha,"names")

      if(any(qv$DOF<3)){warning("Some degrees of freedom are below minimum. They have been set to 3.\nPlease refer to section 3.4 of the vignette for information on running qusage with small sample sizes.")}
      qv$DOF[qv$DOF<3]<-3
   
     
     # Mean=(Sums_Post-Sums_Base)/Ns
    #  SD1=sqrt(Sigmas_Base/Ns)
   #  SD2=sqrt(SD1^2+min.variance.factor)
#FIX ME: assign from qsLIst dont' need to save twice
  # sd.alpha = sqrt(SD1^2+min.variance.factor)/SD1
  #sd.alpha[is.infinite(sd.alpha)] = 1

   qv$SDAlpha[is.infinite(qv$SDAlpha)]=1
      dat = newQSarray(mean=qv$Mean,
                SD=qv$SD,
                sd.alpha = qv$SDAlpha,
                dof=qv$DOF,
                var.method="Welch's"
  )



    }
    ###########Non Paired
    if(!paired){
      ##########First calculate the differential expression for individual genes
     # Sums_Base<-rowSums(Baseline)
     # Sums_Post<-rowSums(PostTreatment)
     # Ns_Base<-ncol(Baseline) #no NAs numeric, not vector
     # Ns_Post<-ncol(PostTreatment) #no NAs , numeric not vector
      sumsList<-list(Baseline,PostTreatment)
      names(sumsList)<-c("Sums_Base","Sums_Post")
      out<-lapply(sumsList,function(x) rowSums(x))
      Ns<-list(Baseline,PostTreatment)
      Ns<-lapply(Ns,function(x) ncol(x))
      names(Ns)<-c("Ns_Base","Ns_Post")
      #if(min(Ns_Base)!=ncol(Baseline) | min(Ns_Post)!=ncol(PostTreatment)){warning("Some NA's in data")}  we assume this: because we enforce no existence of NA values, then the sum of each row will have the ncol.

   

#FIX ME: slow for eset.1 eset.2  
       Sigmas_Base<-sigmasCpp(Baseline,out$Sums_Base/Ns$Ns_Base,Ns$Ns_Base)
       Sigmas_Post<-sigmasCpp(PostTreatment,out$Sums_Post/Ns$Ns_Post,Ns$Ns_Post)

     #Sigmas_Base<-rowSums((Baseline-(Sums_Base)/Ns_Base)^2)/(Ns_Base-1)
      #Sigmas_Post<-rowSums((PostTreatment-(Sums_Post)/Ns_Post)^2)/(Ns_Post-1)
      ROWS<-rownames(Baseline)
      DOF<-Ni(Sigmas_Post+min.variance.factor,Sigmas_Base+min.variance.factor,Ns_Post,Ns_Base)
      #calculate degrees of freedom


  if(any(DOF<3)){warning("Some degrees of freedom are below minimum. They have been set to 3.")}
      DOF[DOF<3]<-3
      Mean=(out$Sums_Post/Ns$Ns_Post-out$Sums_Base/Ns$Ns_Base)
      SD=sqrt(Sigmas_Base/Ns$Ns_Base+Sigmas_Post/Ns$Ns_Post)
    }
  sd.alpha = sqrt(SD^2+min.variance.factor)/SD
  sd.alpha[is.infinite(sd.alpha)] = 1

  dat = newQSarray(mean=Mean,
                SD=sqrt(SD^2+min.variance.factor),
                sd.alpha = sd.alpha,
                dof=DOF,
                var.method="Welch's"
  )
  dat
}


