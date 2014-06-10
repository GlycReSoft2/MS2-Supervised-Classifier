

## Data Partitioning and Bootstrapping

Bootstrap 10 random partitionings of the naively labeled data


```r
dataPartitions <- lapply(1:10, function(i){
  combinedData$randu <- runif(nrow(combinedData), 0,1)
  trainData <- combinedData[combinedData$randu >= 0.5,]
  testData <- combinedData[combinedData$randu < 0.5,]
  list(trainData = trainData, testData = testData)
})
```

Fit an RPART CART model and a Random Forest Model on each training set, and predict the test set labels with them. This first set of data will use a simple predictive formula/ 


```r
fitsAndTest <- lapply(dataPartitions, function(data){
  suppressWarnings(forest <- randomForest(formula = modelFormula, data = data$trainData, mtry = 15, ntree=5000))
  results <- list()
  results$data = data
  results$forest = suppressWarnings(randomForest(formula = modelFormula, data = data$trainData, mtry = 20, ntree=5000))
  results$tree = rpart(formula = modelFormula, data = data$trainData)
  results$forestPred = predict(results$forest, data$testData, type='prob')
  results$treePred = predict(results$tree, data$testData)
  results$forestPredLabAGP = predict(results$forest, AGPLabel, type='prob')
  results$forestPredLabUSSR = predict(results$forest, USSRLabel, type='prob')
  results$treePredLabAGP = predict(results$tree, AGPLabel)
  results$treePredLabUSSR = predict(results$tree, USSRLabel)
  
  return(results)
})
```


```r
checkFits <- lapply(fitsAndTest, function(component){
  list(forestError <- checkModel(component$data$testData, component$forestPred),
       treeError <- checkModel(component$data$testData, component$treePred),
       forestAGPError <- checkModel(AGPLabel, component$forestPredLabAGP),
       forestUSSRError <- checkModel(USSRLabel, component$forestPredLabUSSR),
       treeAGPError <- checkModel(AGPLabel, component$treePredLabAGP),
       treeUSSRError <- checkModel(USSRLabel, component$treePredLabUSSR)
       )
})
```


```r
fitTable <- as.data.frame(do.call(rbind, lapply(do.call(c, checkFits), unlist)))
lenNames <- length(names(fitsAndTest[[1]])[c(-1:-3)])
fitTable$type <- rep(names(fitsAndTest[[1]])[c(-1:-3)], nrow(fitTable)/lenNames)
#kable(fitTable)
ggplot(fitTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, shape = type, color = type)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


```r
fitsAndTestComplex <- lapply(dataPartitions, function(data){
  suppressWarnings(forest <- randomForest(formula = complexModelFormula, data = data$trainData, mtry = 15, ntree=5000))
  results <- list()
  results$data = data
  results$forest = suppressWarnings(randomForest(formula = complexModelFormula, data = data$trainData, mtry = 20, ntree=5000))
  results$tree = rpart(formula = modelFormula, data = data$trainData)
  results$forestPred = predict(results$forest, data$testData, type='prob')
  results$treePred = predict(results$tree, data$testData)
  results$forestPredLabAGP = predict(results$forest, AGPLabel, type='prob')
  results$forestPredLabUSSR = predict(results$forest, USSRLabel, type='prob')
  results$treePredLabAGP = predict(results$tree, AGPLabel)
  results$treePredLabUSSR = predict(results$tree, USSRLabel)
  
  return(results)
})
```


```r
checkFitsComplex <- lapply(fitsAndTestComplex, function(component){
  list(forestError <- checkModel(component$data$testData, component$forestPred),
       treeError <- checkModel(component$data$testData, component$treePred),
       forestAGPError <- checkModel(AGPLabel, component$forestPredLabAGP),
       forestUSSRError <- checkModel(USSRLabel, component$forestPredLabUSSR),
       treeAGPError <- checkModel(AGPLabel, component$treePredLabAGP),
       treeUSSRError <- checkModel(USSRLabel, component$treePredLabUSSR)
       )
})
```


```r
fitTableComplex <- as.data.frame(do.call(rbind, lapply(do.call(c, checkFitsComplex), unlist)))
lenNames <- length(names(fitsAndTest[[1]])[c(-1:-3)])
fitTableComplex$type <- rep(names(fitsAndTest[[1]])[c(-1:-3)], nrow(fitTable)/lenNames)
#kable(fitTableComplex)
ggplot(fitTableComplex, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, shape = type, color = type)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


```r
# complexModelFormula <- call ~ numOxIons + (X.percent_b.ion_coverage * 
#   X.percent_y.ion_coverage * X.b_ion_with_HexNAc_coverage * 
#   X.y_ion_with_HexNAc_coverage *  numStubs * peptideLens)

USSRForest <- randomForest(formula=complexModelFormula, data = USSRLabel, mtry=6, ntree=5000, importance=T)

USSRNaivePred <- predict(USSRForest, combinedData, type="prob")
USSRAGPPred <- predict(USSRForest, AGPLabel, type="prob")
USSRUSSRPred <- predict(USSRForest, USSRLabel, type="prob")

USSRNaiveCheck <- checkModel(combinedData, USSRNaivePred)
USSRAGPCheck <- checkModel(AGPLabel, USSRAGPPred)
USSRUSSRCheck <- checkModel(USSRLabel, USSRUSSRPred)

USSRFitTable <- as.data.frame(do.call(rbind, lapply(list(USSRNaiveCheck, USSRAGPCheck, USSRUSSRCheck), unlist)))
USSRFitTable$target <- c("Naive", "AGP", "USSR")
kable(USSRFitTable)
```

```
## 
## 
## | TruePositive| FalsePositive| FalseNegative| TrueNegative|    FPR|   FDR| .numYes| .numNo|target |
## |------------:|-------------:|-------------:|------------:|------:|-----:|-------:|------:|:------|
## |          958|          1452|           530|         5731| 0.2021| 2e-04|    1488|   7183|Naive  |
## |           56|            10|             6|          165| 0.0571| 1e-03|      62|    175|AGP    |
## |          150|             4|             4|          383| 0.0103| 1e-04|     154|    387|USSR   |
```

```r
ggplot(USSRFitTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, shape = target, color = target)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk fitUSSRComplex](figure/fitUSSRComplex.png) 



