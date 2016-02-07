#SpeedSage Intro
qusage is published software that is slow for large runs, SpeedSage corrects for speed and efficiency at large orders
#Bottlenecking of Functions
Qusage can improve the speed of its algorithm by minimizing the cost of computaiton.

##changes calcIndividualExpressionsC
trading NA flexibility slows down qusage runs, but having the user input no NAs enforcing good input, this speeds up calcIndividualExpressionsC 2X

![qusage profile](/demo/qusageSingleBottleNeck.pdf "Plot of Qusage Profile")


#Individual Expression Function
This test the local version which enforces no NA in Baseline or PostTreatment object, this reduces the flexibility.


```{r} 
library(speedSage)
library(qusage)
eset<-system.file("extdata","eset.RData",package="speedSage")
load(eset)
labels<-c(rep("t0",134),rep("t1",134))
contrast<-"t1-t0"
fileISG<-system.file("extdata","c2.cgp.v5.1.symbols.gmt",package="speedSage")
ISG.geneSet<-read.gmt(fileISG)
ISG.geneSet<-ISG.geneSet[grepl("DER_IFN_GAMMA_RESPONSE_UP",names(ISG.geneSet))]




result<-calcIndividualExpressionsC(Baseline,PostTreatment,paired=FALSE,min.variance.factor=10^-6)
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
```


