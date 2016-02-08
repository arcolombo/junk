#SpeedSage Intro
qusage is published software that is slow for large runs, SpeedSage corrects for speed and efficiency at large orders
#Bottlenecking of Functions
Qusage can improve the speed of its algorithm by minimizing the cost of computaiton.

##changes calcIndividualExpressionsC
trading NA flexibility slows down qusage runs, but having the user input no NAs enforcing good input, this speeds up calcIndividualExpressionsC 2X

![qusage profile](/demo/qusageSingleBottleNeck.pdf "Plot of Qusage Profile")


#calculate Individual Expression Function
This test the local version which enforces no NA in Baseline or PostTreatment object, this reduces the flexibility.
for non-paired end results:
  min       lq     mean   median       uq      max neval cld
  176.5508 224.7822 223.4475 226.2532 229.9082 242.5602   100   b
  112.5428 116.8076 142.3476 163.1434 165.7608 170.5896   100  a

for paired end results:
 min       lq     mean   median       uq      max neval cld
 151.4904 175.8198 204.5262 191.4459 229.1922 332.8214   100   b
  129.6539 143.8261 166.9900 155.2507 191.1199 248.1528   100  a

for input expression matrix with NA enforced:
calcIndividualExpressionsC(testB,testPT)
Error in calcIndividualExpressionsC(testB, testPT) : 
  NA values are present in the expression matrix, please pluck out ...
#makeComparison changes

