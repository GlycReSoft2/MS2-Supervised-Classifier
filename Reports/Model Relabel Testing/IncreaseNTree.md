---
title: "Relabel Testing"
author: "Joshua Klein"
date: "Tuesday, June 17, 2014"
output:
  html_document:
    fig_height: 8
    fig_width: 8
    number_sections: yes
    theme: journal
    toc: yes
---
# Model Fit Test



```r
testFormula <- call ~ meanCoverage + percentUncovered + abs_ppm_error + X.b_ion_with_HexNAc_coverage + 
  X.y_ion_with_HexNAc_coverage + numStubs + peptideLens

mtry <- 5

ntree <- 10000
```


```
## [1] "testFormula ~"                                                                                                                                     
## [2] "testFormula call"                                                                                                                                  
## [3] "testFormula meanCoverage + percentUncovered + abs_ppm_error + X.b_ion_with_HexNAc_coverage + X.y_ion_with_HexNAc_coverage + numStubs + peptideLens"
```

```
## [1] "mtry 5"
```

```
## [1] "ntree 10000"
```

## Load Data

```r
USSRLabel <- prepareAnnotatedModel("../../USSR//20131219_005_results-scored-annotated.fix.csv")
SolomonLabel <- prepareAnnotatedModel("../../Solomon Islands//20131222_004_results-scored-annotated.fix.csv")
SouthCarolinaLabel <- prepareAnnotatedModel("../../South Carolina//20131220_002_results-scored-annotated.fix.csv")

labelData <- list(USSRLabel = USSRLabel, SolomonLabel = SolomonLabel, 
                  SouthCarolinaLabel = SouthCarolinaLabel, Merged = 
                    rbind(model.frame(testFormula, USSRLabel), model.frame(testFormula, SolomonLabel), model.frame(testFormula, SouthCarolinaLabel)))

#USSRData <- loadDir("../../USSR/", includeAnno = F)
#SolomonData <- loadDir("../../Solomon Islands//", includeAnno = F)
#SouthCarolina <- loadDir("../../South Carolina/", includeAnno = F)
```

## Fit Forests on Labeled Data

```r
print(names(labelData))
```

```
## [1] "USSRLabel"          "SolomonLabel"       "SouthCarolinaLabel"
## [4] "Merged"
```

```r
labelForests <- lapply(names(labelData), function(d){
  rf <- randomForest(formula = testFormula,data = labelData[[d]], mtry = mtry, ntree = ntree, importance = T)
  varImpPlot(rf, main = d)
  return(rf)
})
```

<img src="figure/IncreaseNTree-Rmdfit-forests1.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/IncreaseNTree-Rmdfit-forests2.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/IncreaseNTree-Rmdfit-forests3.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/IncreaseNTree-Rmdfit-forests4.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" />

```r
names(labelForests) <- names(labelData)
```


```r
newLabels <- lapply(labelData, function(dataSet){
  lapply(labelForests, function(fit){
      list(check = checkModel(dataSet=dataSet, modelPred=predict(fit, dataSet, type="prob")),
       fNeg = getFalseNegatives(dataSet=dataSet, modelPred=predict(fit, dataSet, type="prob")),
       fPos = getFalsePositives(dataSet=dataSet, modelPred=predict(fit, dataSet, type="prob")),
       tPos = getTruePositives(dataSet=dataSet, modelPred=predict(fit, dataSet, type="prob")),
       tNeg = getTrueNegatives(dataSet=dataSet, modelPred=predict(fit, dataSet, type="prob"))
       )
  })  
})

newLabels <- unlist(newLabels, recursive = F)
```



```r
errorTable <- as.data.frame(do.call(rbind, lapply(newLabels, function(x){unlist(x$check)})))
errorTable <- name_rows(errorTable)
errorTable$TPR <- with(errorTable, TruePositive/(TruePositive + FalseNegative))
errorTable <- arrange(errorTable, TPR, -FPR)

kable(errorTable)
```



| TruePositive| FalsePositive| FalseNegative| TrueNegative|    FPR|   FDR| ErrorRate| .numYes| .numNo|.rownames                             |    TPR|
|------------:|-------------:|-------------:|------------:|------:|-----:|---------:|-------:|------:|:-------------------------------------|------:|
|          215|            28|            28|          574| 0.0465| 2e-04|    0.0663|     243|    602|USSRLabel.SouthCarolinaLabel          | 0.8848|
|          107|             0|            13|          581| 0.0000| 0e+00|    0.0185|     120|    581|SolomonLabel.SouthCarolinaLabel       | 0.8917|
|          387|            28|            41|         2188| 0.0126| 0e+00|    0.0261|     428|   2216|Merged.SouthCarolinaLabel             | 0.9042|
|           60|            32|             5|         1001| 0.0310| 5e-04|    0.0337|      65|   1033|SouthCarolinaLabel.USSRLabel          | 0.9231|
|          114|             4|             6|          577| 0.0069| 1e-04|    0.0143|     120|    581|SolomonLabel.USSRLabel                | 0.9500|
|          231|            45|            12|          557| 0.0748| 3e-04|    0.0675|     243|    602|USSRLabel.SolomonLabel                | 0.9506|
|          414|            66|            14|         2150| 0.0298| 1e-04|    0.0303|     428|   2216|Merged.SolomonLabel                   | 0.9673|
|           63|            21|             2|         1012| 0.0203| 3e-04|    0.0209|      65|   1033|SouthCarolinaLabel.SolomonLabel       | 0.9692|
|          417|            36|            11|         2180| 0.0162| 0e+00|    0.0178|     428|   2216|Merged.USSRLabel                      | 0.9743|
|          243|             0|             0|          602| 0.0000| 0e+00|    0.0000|     243|    602|USSRLabel.USSRLabel                   | 1.0000|
|          243|             0|             0|          602| 0.0000| 0e+00|    0.0000|     243|    602|USSRLabel.Merged                      | 1.0000|
|          120|             0|             0|          581| 0.0000| 0e+00|    0.0000|     120|    581|SolomonLabel.SolomonLabel             | 1.0000|
|          120|             0|             0|          581| 0.0000| 0e+00|    0.0000|     120|    581|SolomonLabel.Merged                   | 1.0000|
|           65|             0|             0|         1033| 0.0000| 0e+00|    0.0000|      65|   1033|SouthCarolinaLabel.SouthCarolinaLabel | 1.0000|
|           65|             0|             0|         1033| 0.0000| 0e+00|    0.0000|      65|   1033|SouthCarolinaLabel.Merged             | 1.0000|
|          428|             0|             0|         2216| 0.0000| 0e+00|    0.0000|     428|   2216|Merged.Merged                         | 1.0000|


```r
ggplot(errorTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, label = zapsmall(TruePositive/.numYes, 3))) + geom_point() + scale_size(range = c(4,10))
```

<img src="figure/IncreaseNTree-Rmderror-plot1.png" title="plot of chunk error-plot" alt="plot of chunk error-plot" style="display: block; margin: auto;" />

```r
ggplot(errorTable, aes(y = TPR, x = FPR, size = ErrorRate, color = .rownames)) + geom_point(position="jitter") + scale_size(range = c(4,10)) 
```

<img src="figure/IncreaseNTree-Rmderror-plot2.png" title="plot of chunk error-plot" alt="plot of chunk error-plot" style="display: block; margin: auto;" />

