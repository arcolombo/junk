library(qusage)
#tests C version which enforces no NA in baseline or PT
#baseline and PT is created in qusage 
# enforcing no NAs, not a flexible solution
results1<-calcIndividualExpressionsC(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6)
results2<-calcIndividualExpressions(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6)
identical(results1,results2)

microbenchmark(
results1<-calcIndividualExpressionsC(Baseline,PostTreatment,min.variance.factor=10^-6),
results2<-calcIndividualExpressions(Baseline,PostTreatment,min.variance.factor=10^-6) 
)

#add NAs and test
testPT<-PostTreatment[1:20,]
testPT<-cbind(rbind(testPT,NaN),NA)
testB<-Baseline[1:20,]
testB<-cbind(rbind(testB,NaN),NA)


