calcIndividualExpressionsC<-function(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6){
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
      Sums_Base<-rowSums(Baseline)
      Sums_Post<-rowSums(PostTreatment)
      #Ns<-rowSumsC(!is.na(Baseline-PostTreatment))
      #if(min(Ns)!=ncol(Baseline)){warning("Some NA's in data")}
      Sigmas_Base<-rowSums((Baseline-PostTreatment-(Sums_Base-Sums_Post)/(Baseline-PostTreatment))^2)/((Baseline-PostTreatment)-1)
      DOF<-(Baseline-PostTreatment)
      if(any(DOF<3, na.rm=T)){warning("Some degrees of freedom are below minimum. They have been set to 3.\nPlease refer to section 3.4 of the vignette for information on running qusage with small sample sizes.")}
      DOF[DOF<3]<-3
      Mean=(Sums_Post-Sums_Base)/(Baseline-PostTreatment)
      SD=sqrt(Sigmas_Base/(Baseline-PostTreatment))
    }
    ###########Non Paired
    if(!paired){
      ##########First calculate the differential expression for individual genes
      Sums_Base<-rowSums(Baseline)
      Sums_Post<-rowSums(PostTreatment)
      #Ns_Base<-rowSums(!is.na(Baseline))
      #Ns_Post<-rowSums(!is.na(PostTreatment))
      #if(min(Ns_Base)!=ncol(Baseline) | min(Ns_Post)!=ncol(PostTreatment)){warning("Some NA's in data")}  we assume this: because we enforce no existence of NA values, then the sum of each row will have the ncol.
      
      Sigmas_Base<-rowSums((Baseline-(Sums_Base)/ncol(Baseline))^2)/(ncol(Baseline)-1)
      Sigmas_Post<-rowSums((PostTreatment-(Sums_Post)/ncol(PostTreatment))^2)/(ncol(PostTreatment)-1)
      ROWS<-rownames(Baseline)
      DOF<-Ni(Sigmas_Post+min.variance.factor,Sigmas_Base+min.variance.factor,ncol(Baseline),ncol(PostTreatment))
      #calculate degrees of freedom


  if(any(DOF<3, na.rm=T)){warning("Some degrees of freedom are below minimum. They have been set to 3.")}
      DOF[DOF<3]<-3
      Mean=(Sums_Post/ncol(PostTreatment)-Sums_Base/ncol(Baseline))
      SD=sqrt(Sigmas_Base/ncol(Baseline)+Sigmas_Post/ncol(PostTreatment))
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


