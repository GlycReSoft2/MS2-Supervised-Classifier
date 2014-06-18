

## Set Parameters

```r
N = 150
K = .30
```


## Generate training and testing data sets for leave K% out cross-validation.


```r
dataPartitions <- lapply(LETTERS[1:N], function(i){
  USSRLabel$randu <- runif(nrow(USSRLabel), 0, 1)
  trainData <- USSRLabel[USSRLabel$randu > K,]
  testData <- USSRLabel[USSRLabel$randu <= K,]
  return(list(trainData = trainData, testData = testData))
})
names(dataPartitions) <- paste(LETTERS,1:N,sep='')
```


```r
fitAndCheck <- lapply(dataPartitions, function(dataSet){
  fitForest <- randomForest(complexModelFormulaWithInteractions, dataSet$trainData, mtry=6, ntree=5000, importance=T)
  pred <- predict(fitForest, dataSet$testData, type='prob')
  check <- checkModel(dataSet$testData, pred)
  return(list(check = check, fit = fitForest))
})

totalFit <- randomForest(complexModelFormulaWithInteractions, USSRLabel, mtry=6, ntree=5000, importance=T)
totalPred <- predict(totalFit, USSRLabel, type='prob')
totalCheck <- as.data.frame(checkModel(USSRLabel, totalPred))
totalCheck$type <- "Total"

fitTable <- as.data.frame(do.call(rbind, lapply(fitAndCheck, function(fc)unlist(fc$check))))
fitTable$type <- "K-out"
fitTable <- rbind(fitTable, totalCheck)
ggplot(rbind(fitTable, totalCheck), aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, shape = type)) + geom_point() + scale_size(range = c(4,10)) + labs(title = "USSR label prediction by leave K out fits")
```

![plot of chunk fitAndCheck](figure/fitAndCheck.png) 


```r
bestFitIDs <- row.names(fitTable[fitTable$ErrorRate < 0.15, ])
print(bestFitIDs)
```

```
##  [1] "B2"   "C3"   "K11"  "M13"  "N14"  "O15"  "V22"  "X24"  "Z26"  "H34" 
## [11] "I35"  "K37"  "P42"  "T46"  "V48"  "X50"  "B54"  "E57"  "I61"  "J62" 
## [21] "P68"  "R70"  "T72"  "X76"  "Y77"  "Z78"  "G85"  "J88"  "K89"  "A105"
## [31] "D108" "G111" "H112" "J114" "L116" "T124" "Y129" "Z130" "F136" "G137"
## [41] "M143" "N144" "O145" "Q147" "151"
```

```r
bestFits <- lapply(bestFitIDs, function(i){
  if(i %in% names(fitAndCheck)){fitAndCheck[[i]]$fit}
  else {totalFit}
})

names(bestFits) <- bestFitIDs

SolomonTests <- lapply(bestFits, function(f){
  checkModel(SolomonLabel, predict(f, SolomonLabel, type="prob"))
})

SolomonTestTable <- as.data.frame(do.call(rbind, lapply(SolomonTests, unlist)))
row.names(SolomonTestTable) <- bestFitIDs
SolomonTestTable$type <- fitTable[fitTable$ErrorRate < 0.15, "type"]
ggplot(SolomonTestTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, shape = type)) + geom_point() + scale_size(range = c(4,10)) + labs(title = "Solomon Islands label prediction by best performing USSR predictors")
```

![plot of chunk crossSolomonTests](figure/crossSolomonTests.png) 


```r
STT <- as.data.table(name_rows(SolomonTestTable))
bestAcrossFitIDs <- STT[(FalseNegative < 12) & (TruePositive/.numYes > 0.84) & type != "Total", .rownames]
ggplot(STT[(FalseNegative < 12) & (TruePositive/.numYes > 0.84) & type != "Total",], aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, shape = type)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```r
bestData <- dataPartitions[bestAcrossFitIDs]
bestForests <- bestFits[bestAcrossFitIDs]
saveRDS(list(data =  bestData, models = bestForests), paste("../SavedData/leave-k-out-crossvalidate-USSR-", format(now(), format="%Y-%m-%d_%H-%M-%S"),".RDS", sep = ''))
```

