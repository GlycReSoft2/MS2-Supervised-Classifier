

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
##  [1] "B2"   "F6"   "H8"   "I9"   "J10"  "K11"  "L12"  "M13"  "N14"  "O15" 
## [11] "P16"  "Q17"  "S19"  "U21"  "V22"  "X24"  "Y25"  "Z26"  "A27"  "B28" 
## [21] "E31"  "G33"  "H34"  "I35"  "K37"  "M39"  "N40"  "P42"  "R44"  "T46" 
## [31] "V48"  "W49"  "X50"  "Y51"  "Z52"  "B54"  "E57"  "H60"  "I61"  "J62" 
## [41] "N66"  "P68"  "R70"  "S71"  "T72"  "U73"  "X76"  "Z78"  "D82"  "E83" 
## [51] "J88"  "K89"  "L90"  "M91"  "S97"  "U99"  "A105" "B106" "D108" "E109"
## [61] "H112" "J114" "K115" "M117" "N118" "P120" "T124" "X128" "Y129" "Z130"
## [71] "C133" "D134" "E135" "F136" "G137" "J140" "M143" "N144" "O145" "Q147"
## [81] "S149" "T150" "151"
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

