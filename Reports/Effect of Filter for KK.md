

## Filtered Data v.s. Unfiltered Technical Replicate

```r
filtUSSR <- prepareModel("/Users/jaklein/Dropbox/Glycresoft 2014/For ASMS using glycomics for hypo/20131219_005_results-scored.fix.csv", call =F)
repUSSR <- prepareModel("../USSR//20131220_006_results-scored.fix.csv", call =F)
labelUSSR <- prepareAnnotatedModel("../USSR//20131219_005_results-scored-annotated.fix.csv")
matchedUSSR <- shimPercentCalcs(prepareModel("/Users/jaklein/Dropbox/Glycomics Sandbox/Test Data/Sample outputs for MS2 data/20131219_005_results.fix.csv", call =F))

labelForest <- randomForest(complexModelFormulaWithInteractions, labelUSSR, mtry=6, ntree=5000, importance=T)

filtUSSR$call <- predict(labelForest, filtUSSR)
repUSSR$call <- predict(labelForest, repUSSR)
matchedUSSR$call <- predict(labelForest, matchedUSSR)

filtForest <- randomForest(complexModelFormulaWithInteractions, filtUSSR, mtry=6, ntree=5000, importance=T)
repForest <- randomForest(complexModelFormulaWithInteractions, repUSSR, mtry=6, ntree=5000, importance=T)
matchForest <- randomForest(complexModelFormulaWithInteractions, matchedUSSR, mtry=6, ntree=5000, importance=T)

checkSelf <- unlist(checkModel(labelUSSR, predict(labelForest, labelUSSR, 'prob')))
checkFilt <- unlist(checkModel(labelUSSR, predict(filtForest, labelUSSR, 'prob')))
checkRep <- unlist(checkModel(labelUSSR, predict(repForest, labelUSSR, 'prob')))
checkMatch <- unlist(checkModel(labelUSSR, predict(matchForest, labelUSSR, 'prob')))

errorTable <- name_rows(as.data.frame(rbind(checkSelf, checkMatch,checkFilt, checkRep)))
ggplot(errorTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, shape = .rownames)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

```r
print(errorTable)
```

```
##   TruePositive FalsePositive FalseNegative TrueNegative      FPR       FDR
## 1          149             1             5          386 0.002584 1.734e-05
## 2          149             1             5          386 0.002584 1.734e-05
## 3          147             1             7          386 0.002584 1.758e-05
## 4          125            31            29          356 0.080103 6.404e-04
##   ErrorRate .numYes .numNo  .rownames
## 1   0.01109     154    387  checkSelf
## 2   0.01109     154    387 checkMatch
## 3   0.01479     154    387  checkFilt
## 4   0.11091     154    387   checkRep
```

```r
sum(filtUSSR$call == "Yes")
```

```
## [1] 166
```


