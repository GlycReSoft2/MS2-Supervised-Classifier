


```r
USSRLabel <- prepareAnnotatedModel("../USSR//20131219_005_results-scored-annotated.fix.csv")
USSRData <- loadDir("../USSR/", includeAnno=F)

labelForest <- randomForest(complexModelFormulaWithInteractions, USSRLabel, mtry=6, ntree = 5000, importance = T)

newLabels <- lapply(USSRData, function(dataSet){
  dataSet$call <- predict(labelForest, dataSet)
  fit <- randomForest(complexModelFormulaWithInteractions, dataSet, mtry=6, ntree = 5000, importance = T)
  list(check = checkModel(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fNeg = getFalseNegatives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fPos = getFalsePositives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")))
  
})
```

```
## Warning: row names were found from a short variable and have been
## discarded
```

```r
errorTable <- as.data.frame(do.call(rbind, lapply(newLabels, function(x){unlist(x$check)})))

ggplot(errorTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


```r
falseNegatives <- lapply(newLabels, function(x){
  mod <- model.frame(complexModelFormulaWithInteractions,x[["fNeg"]])
  x[["fNeg"]] <- name_rows(x[["fNeg"]])
  print(x[["fNeg"]][,c('.rownames','Yes', "Notes")])
  x[["fNeg"]]
})
```

```
##   .rownames    Yes
## 1       336 0.0564
## 2       351 0.0458
## 3       389 0.0162
## 4       473 0.0130
## 5       498 0.0452
##                                                                                                             Notes
## 1 can't say unless we lookat the fragments. Also other matches show that this site is glycosylated not the other.
## 2                        this is the site found to be correct from previous matches. Check the neighbor's score!!
## 3                                                                                           based on neighbors!!!
## 4                                                                        based on the site found in other matches
## 5                                                                                                                
##    .rownames    Yes
## 1        144 0.1020
## 2        204 0.4892
## 3        301 0.1348
## 4        336 0.2744
## 5        337 0.4802
## 6        351 0.2864
## 7        367 0.4288
## 8        368 0.0048
## 9        389 0.0498
## 10       407 0.1572
## 11       421 0.0002
## 12       442 0.1990
## 13       444 0.1078
## 14       469 0.4094
## 15       473 0.1120
## 16       485 0.4548
## 17       493 0.3618
## 18       498 0.0518
##                                                                                                              Notes
## 1                                                                                           very little difference
## 2                                                                                                                 
## 3                                                                                         little backbone coverage
## 4  can't say unless we lookat the fragments. Also other matches show that this site is glycosylated not the other.
## 5                                                                                                        low score
## 6                         this is the site found to be correct from previous matches. Check the neighbor's score!!
## 7                                                                                                                 
## 8                                                  probably true, glycan composition has high probability + 1 stub
## 9                                                                                            based on neighbors!!!
## 10                                                                                                     most likely
## 11                                                                                                                
## 12                                                                                                                
## 13                                                                                                                
## 14                                                                                                                
## 15                                                                        based on the site found in other matches
## 16                                                                                                                
## 17                                                                                                                
## 18                                                                                                                
##    .rownames    Yes
## 1        144 0.1978
## 2        169 0.4170
## 3        301 0.2104
## 4        336 0.2558
## 5        337 0.4852
## 6        340 0.3780
## 7        351 0.1926
## 8        368 0.0142
## 9        372 0.4730
## 10       389 0.1694
## 11       398 0.4594
## 12       407 0.0754
## 13       421 0.0006
## 14       423 0.2688
## 15       442 0.3092
## 16       444 0.0672
## 17       469 0.3958
## 18       473 0.1426
## 19       485 0.0544
## 20       488 0.4106
## 21       498 0.4550
## 22       515 0.0682
## 23       516 0.1262
## 24       534 0.4162
##                                                                                                              Notes
## 1                                                                                           very little difference
## 2                                                                                               most probably true
## 3                                                                                         little backbone coverage
## 4  can't say unless we lookat the fragments. Also other matches show that this site is glycosylated not the other.
## 5                                                                                                        low score
## 6                                                                                   most likely glycan composition
## 7                         this is the site found to be correct from previous matches. Check the neighbor's score!!
## 8                                                  probably true, glycan composition has high probability + 1 stub
## 9                                                                                                                 
## 10                                                                                           based on neighbors!!!
## 11                                                                                                                
## 12                                                                                                     most likely
## 13                                                                                                                
## 14                                                                                                                
## 15                                                                                                                
## 16                                                                                                                
## 17                                                                                                                
## 18                                                                        based on the site found in other matches
## 19                                                                                                                
## 20                                                                                               again second site
## 21                                                                                                                
## 22                                                                                                                
## 23                                                            most likely true, I never see stubs for this peptide
## 24
```

```r
falsePositives <- lapply(newLabels, function(x){
  mod <- model.frame(complexModelFormulaWithInteractions,x[["fPos"]])
  x[["fPos"]] <- name_rows(x[["fPos"]])
  print(x[["fPos"]][,c('.rownames','Yes', "Notes")])
  x[["fPos"]]
})
```

```
## Error: undefined columns selected
```

