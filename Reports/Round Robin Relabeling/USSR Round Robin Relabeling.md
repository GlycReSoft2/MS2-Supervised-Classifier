


```r
USSRLabel <- prepareAnnotatedModel("../USSR//20131219_005_results-scored-annotated.fix.csv")
USSRData <- loadDir("../USSR/", includeAnno=F)
SolomonLabel <- prepareAnnotatedModel("../Solomon Islands//20131222_004_results-scored-annotated.fix.csv")
SouthCarolinaLabel <- prepareAnnotatedModel("../South Carolina//20131220_002_results-scored-annotated.fix.csv")

labelForest <- randomForest(testFormula, USSRLabel, mtry=mtry, ntree = 5000, importance = T)

varImpPlot(labelForest)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

```r
newLabels <- lapply(USSRData, function(dataSet){
  dataSet$call <- predict(labelForest, dataSet)
  fit <- randomForest(testFormula, dataSet, mtry=mtry, ntree = 5000, importance = T)
  list(check = checkModel(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fNeg = getFalseNegatives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fPos = getFalsePositives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       tPos = getTruePositives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       tNeg = getTrueNegatives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob"))
       )
  
})
```


```r
USSRLabelsOther <- lapply(list(SolomonLabelByUSSR = SolomonLabel, SouthCarolinaLabelByUSSR = SouthCarolinaLabel), function(dataSet){
  list(check = checkModel(dataSet=dataSet, modelPred=predict(labelForest, dataSet, type="prob")),
       fNeg = getFalseNegatives(dataSet=dataSet, modelPred=predict(labelForest, dataSet, type="prob")),
       fPos = getFalsePositives(dataSet=dataSet, modelPred=predict(labelForest, dataSet, type="prob")))  
})
```


```r
otherLabelsUSSR <- lapply(list(SolomonLabeledUSSR = SolomonLabel, SouthCarolinaLabeledUSSR = SouthCarolinaLabel), function(dataSet){
  fit <- randomForest(testFormula, dataSet, mtry=mtry, ntree = 5000, importance = T)
  list(check = checkModel(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fNeg = getFalseNegatives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")),
       fPos = getFalsePositives(dataSet=USSRLabel, modelPred=predict(fit, USSRLabel, type="prob")))  
})
```


```r
errorTableRep <- as.data.frame(do.call(rbind, lapply(newLabels, function(x){unlist(x$check)})))

errorTableUSSRLabelOther <- as.data.frame(do.call(rbind, lapply(USSRLabelsOther, function(x){unlist(x$check)})))

errorTableOtherLabelUSSR <- as.data.frame(do.call(rbind, lapply(otherLabelsUSSR, function(x){unlist(x$check)})))

errorTable <- name_rows(rbind(errorTableRep, errorTableUSSRLabelOther, errorTableOtherLabelUSSR))

errorTable <- arrange(errorTable, TruePositive/.numYes, -ErrorRate)

kable(errorTable)
```

```
## 
## 
## | TruePositive| FalsePositive| FalseNegative| TrueNegative|    FPR|   FDR| ErrorRate| .numYes| .numNo|.rownames                           |
## |------------:|-------------:|-------------:|------------:|------:|-----:|---------:|-------:|------:|:-----------------------------------|
## |          215|            28|            28|          574| 0.0465| 2e-04|    0.0663|     243|    602|SouthCarolinaLabeledUSSR            |
## |           60|            31|             5|         1002| 0.0300| 5e-04|    0.0328|      65|   1033|SouthCarolinaLabelByUSSR            |
## |          114|             4|             6|          577| 0.0069| 1e-04|    0.0143|     120|    581|SolomonLabelByUSSR                  |
## |          231|            45|            12|          557| 0.0748| 3e-04|    0.0675|     243|    602|SolomonLabeledUSSR                  |
## |          234|             4|             9|          598| 0.0066| 0e+00|    0.0154|     243|    602|20131220_006_results-scored.fix.csv |
## |          236|            14|             7|          588| 0.0233| 1e-04|    0.0249|     243|    602|20131219_006_results-scored.fix.csv |
## |          243|             0|             0|          602| 0.0000| 0e+00|    0.0000|     243|    602|20131219_005_results-scored.fix.csv |
```

```r
ggplot(errorTable, aes(size = TruePositive/.numYes,  y = FalsePositive, x = FalseNegative, color = ErrorRate, label = zapsmall(TruePositive/.numYes, 3))) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) 

```r
ggplot(errorTable, aes(y = TruePositive/(TruePositive + FalseNegative), x = FalsePositive/(FalsePositive + TrueNegative), size = ErrorRate, color =  .rownames)) + geom_point() + scale_size(range = c(4,10))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 


```r
falseNegatives <- lapply(newLabels, function(x){
  kable(x[["fNeg"]][,c('Glycopeptide_identifier', "ppm_error", 'peptideLens',"meanCoverage", "percentUncovered", "numStubs", 'fNeg')])
})
```

```
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
```

```
## 
## 
## |Glycopeptide_identifier | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs| fNeg|
## |:-----------------------|---------:|-----------:|------------:|----------------:|--------:|----:|
## 
## 
## |    |Glycopeptide_identifier                                                                      | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   fNeg|
## |:---|:--------------------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |143 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;5;4;1;0]                                      |    -2.598|          19|        3.947|           0.1053|        0| 0.3704|
## |150 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;4;1;0]                                      |    -3.296|          19|        4.895|           0.0000|        1| 0.3564|
## |163 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;0;0]                                      |     1.653|          19|        3.105|           0.1053|        1| 0.1718|
## |176 |N(HexNAc)GSYPN(Deamidated)LSK[0;6;6;1;0]                                                     |     9.431|           9|        1.778|           0.0000|        1| 0.0584|
## |187 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;0;0]                                      |    -5.907|          19|        4.316|           0.0000|        1| 0.0128|
## |213 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQN(Deamidated)IHPVTIGEC(Carbamidomethyl)PK[1;6;5;3;0] |     5.575|          27|        3.074|           0.1481|        0| 0.1160|
## |282 |N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                     |    -2.425|           9|        2.000|           0.1111|        1| 0.1790|
## 
## 
## |    |Glycopeptide_identifier                                                                      | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   fNeg|
## |:---|:--------------------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |110 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[3;7;6;1;0]                           |    -8.104|          22|        6.682|           0.0000|        1| 0.3332|
## |143 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;5;4;1;0]                                      |    -2.598|          19|        3.947|           0.1053|        0| 0.3574|
## |148 |N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                                     |     8.121|           9|        3.889|           0.1111|        0| 0.4898|
## |150 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;4;1;0]                                      |    -3.296|          19|        4.895|           0.0000|        1| 0.4778|
## |163 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;0;0]                                      |     1.653|          19|        3.105|           0.1053|        1| 0.3446|
## |176 |N(HexNAc)GSYPN(Deamidated)LSK[0;6;6;1;0]                                                     |     9.431|           9|        1.778|           0.0000|        1| 0.1256|
## |187 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;0;0]                                      |    -5.907|          19|        4.316|           0.0000|        1| 0.0446|
## |213 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQN(Deamidated)IHPVTIGEC(Carbamidomethyl)PK[1;6;5;3;0] |     5.575|          27|        3.074|           0.1481|        0| 0.2684|
## |282 |N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                     |    -2.425|           9|        2.000|           0.1111|        1| 0.0520|
```

```r
falsePositives <- lapply(newLabels, function(x){
  kable(x[["fPos"]][,c('Glycopeptide_identifier', "ppm_error",'peptideLens',"meanCoverage", "percentUncovered", "numStubs", 'fPos')])
})
```

```
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
```

```
## 
## 
## |Glycopeptide_identifier | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs| fPos|
## |:-----------------------|---------:|-----------:|------------:|----------------:|--------:|----:|
## 
## 
## |    |Glycopeptide_identifier                                                        | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   fPos|
## |:---|:------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |175 |N(Deamidated)LLWLTEKNGSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.957|           0.0000|        0| 0.9552|
## |194 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)NK[1;5;4;2;0] |    -3.987|          23|        5.304|           0.0000|        0| 0.7548|
## |195 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVNN(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.304|           0.0000|        0| 0.7548|
## |196 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0] |    -3.987|          23|        5.783|           0.0000|        0| 0.8544|
## |197 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.783|           0.0000|        0| 0.8544|
## |198 |NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.217|           0.0000|        0| 0.9696|
## |199 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVNN(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.696|           0.0000|        0| 0.9504|
## |210 |NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.000|           0.0000|        0| 0.6472|
## |218 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)NK[1;7;6;4;0] |    -3.710|          23|        5.783|           0.0000|        0| 0.8532|
## |249 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        5.087|           0.0000|        0| 0.7796|
## |252 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                 |    -7.327|          23|        5.913|           0.0000|        1| 0.6518|
## |295 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.7496|
## |296 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.7496|
## |314 |N(HexNAc)GSYPN(HexNAc)LSK[3;14;13;0;0]                                         |     7.494|           9|        1.222|           0.5556|        2| 0.9880|
## 
## 
## |    |Glycopeptide_identifier                                                        | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   fPos|
## |:---|:------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |175 |N(Deamidated)LLWLTEKNGSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.957|           0.0000|        0| 0.6490|
## |295 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.8664|
## |296 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.8664|
## |314 |N(HexNAc)GSYPN(HexNAc)LSK[3;14;13;0;0]                                         |     7.494|           9|        1.222|           0.5556|        2| 0.9888|
```

```r
truePositives<-lapply(newLabels, function(x){
  kable(x[["tPos"]][,c('Glycopeptide_identifier', "ppm_error", 'peptideLens',"meanCoverage", "percentUncovered", "numStubs", 'tPos')])
  })
```

```
## 
## 
## |    |Glycopeptide_identifier                                                                      | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   tPos|
## |:---|:--------------------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |1   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                            |   -5.9097|          18|      22.5556|           0.0000|      8.0| 1.0000|
## |2   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                   |    1.3485|          22|      26.5909|           0.0000|      8.0| 1.0000|
## |3   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;2;0]                                   |    1.3485|          22|      26.0455|           0.0000|      8.0| 1.0000|
## |4   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                            |   -3.9412|          18|      20.2222|           0.0000|      9.0| 1.0000|
## |5   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                            |   -6.5906|          18|      18.0000|           0.0000|      7.0| 1.0000|
## |6   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                            |   -8.5683|          18|      21.6667|           0.0000|      8.0| 0.9996|
## |7   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                            |   -6.7784|          18|      18.5000|           0.0000|      6.0| 1.0000|
## |8   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                   |    2.8744|          22|      25.0455|           0.0000|      9.0| 1.0000|
## |9   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;2;0]                                   |    2.8744|          22|      24.5000|           0.0000|      9.0| 1.0000|
## |10  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;1;0]                           |    2.9786|          22|      22.5000|           0.0000|      9.0| 1.0000|
## |11  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;1;0]                                      |   -3.9610|          19|      22.8947|           0.0000|      6.0| 1.0000|
## |12  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                                   |   -2.7807|          22|      20.6364|           0.0000|      8.0| 1.0000|
## |13  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                   |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |14  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                                   |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |15  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;1;0]                                      |   -0.8840|          19|      22.7895|           0.0000|      9.0| 1.0000|
## |16  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;14;4;0;0]                          |    9.6761|          22|      29.1364|           0.0000|      7.0| 0.9884|
## |17  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;0;0]                                      |   -2.7961|          19|      20.9474|           0.0000|      5.0| 1.0000|
## |18  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;0;0]                                   |   -2.7807|          22|      20.0909|           0.0000|      8.0| 1.0000|
## |19  |NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                                 |    0.1308|           9|       9.8889|           0.0000|     10.0| 1.0000|
## |20  |N(HexNAc)GSYPNLSK[1;5;4;2;0]                                                                 |    0.1308|           9|       9.0000|           0.0000|     10.0| 1.0000|
## |21  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;3;0]                                      |   -2.8410|          19|      21.0526|           0.0000|      8.0| 1.0000|
## |22  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                   |   -0.1015|          22|      22.5000|           0.0000|      7.0| 1.0000|
## |23  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;1;0]             |   -4.3482|          27|      21.9630|           0.0000|      6.0| 0.9988|
## |24  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;2;0]                             |   -1.5373|          28|      22.3929|           0.0000|      1.0| 0.9976|
## |25  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                   |   -0.1015|          22|      21.9545|           0.0000|      7.0| 1.0000|
## |26  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;1;0]                                      |   -3.5609|          19|      18.5789|           0.0000|      8.0| 1.0000|
## |27  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;5;4;2;0]             |   -0.1807|          27|      22.2222|           0.0000|      5.0| 0.9996|
## |28  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;2;0]                                      |    0.7980|          19|      18.5789|           0.0000|      5.0| 1.0000|
## |29  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]                             |   -0.4458|          28|      22.9643|           0.0000|      2.0| 0.9950|
## |30  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]                             |    0.4326|          28|      22.5000|           0.0000|      1.0| 0.9970|
## |31  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                            |   -3.0106|          18|      16.5556|           0.0000|      6.0| 1.0000|
## |32  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                                   |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |33  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;0;0]                                   |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |34  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;2;0]                                      |   -0.0219|          19|      19.6842|           0.0000|      8.0| 0.9996|
## |35  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]                             |   -2.8434|          28|      20.0000|           0.0000|      2.0| 0.9932|
## |36  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;5;4;1;0]                                                        |   -7.4680|          18|      12.3889|           0.0000|      4.0| 0.9804|
## |37  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;2;0]                                      |    0.7444|          19|      18.9474|           0.0000|      5.0| 1.0000|
## |38  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;2;0]                                      |   -1.2420|          19|      17.6842|           0.0000|      5.0| 1.0000|
## |39  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]                       |   -1.5816|          22|      15.9545|           0.0000|      6.0| 1.0000|
## |40  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]                       |   -1.5416|          22|      17.0000|           0.0000|      6.0| 1.0000|
## |41  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;3;0]                       |   -1.5816|          22|      15.5000|           0.0000|      6.0| 1.0000|
## |42  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;2;0]                                      |    5.9064|          19|      15.6842|           0.0000|      5.0| 1.0000|
## |43  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]                             |   -0.7608|          28|      19.1429|           0.0000|      2.0| 0.9950|
## |44  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                            |   -2.7114|          18|      14.1111|           0.0000|      6.0| 1.0000|
## |45  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;4;0]                                      |   -3.7451|          19|      17.9474|           0.0000|      6.0| 1.0000|
## |46  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;2;0]                       |   -1.5416|          22|      15.9091|           0.0000|      6.0| 1.0000|
## |47  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;3;0]                                      |   -2.2429|          19|      17.1053|           0.0000|      8.0| 1.0000|
## |48  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                   |   -4.0561|          22|      14.7727|           0.0000|      6.0| 1.0000|
## |49  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;2;0]                                   |   -4.0561|          22|      14.7727|           0.0000|      6.0| 1.0000|
## |50  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                            |   -8.2882|          18|      14.7222|           0.0000|      5.0| 0.9994|
## |51  |NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                                 |   -1.4933|           9|       8.3333|           0.0000|     10.0| 1.0000|
## |52  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;2;0]                                      |   -0.5044|          19|      15.4737|           0.0000|      4.0| 1.0000|
## |53  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                            |   -7.0563|          18|      12.6667|           0.0000|      8.0| 0.9988|
## |54  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;1;0]                             |   -4.1984|          28|      18.5714|           0.0000|      2.0| 0.9948|
## |55  |NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                                 |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |56  |HN(HexNAc)VTR[1;5;4;1;0]                                                                     |   -1.3984|           5|       2.8000|           0.0000|      3.0| 0.9934|
## |57  |N(HexNAc)GSYPNLSK[1;4;5;1;0]                                                                 |   -1.4933|           9|       6.8889|           0.0000|     10.0| 1.0000|
## |58  |N(HexNAc)GSYPNLSK[1;6;5;0;0]                                                                 |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |59  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                                            |   -2.2699|          18|      14.2222|           0.0000|      2.0| 0.9940|
## |60  |N(HexNAc)GSYPN(HexNAc)LSK[2;6;6;0;0]                                                         |   -1.0568|           9|       9.6667|           0.0000|      9.0| 1.0000|
## |61  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                            |   -4.3064|          18|      15.0556|           0.0000|      5.0| 1.0000|
## |62  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                   |   -5.7919|          22|      14.2273|           0.0000|      7.0| 1.0000|
## |63  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;1;0]                                   |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |64  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                   |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |65  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;1;0]                                      |   -2.0419|          19|      13.8947|           0.0000|      3.0| 1.0000|
## |66  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;5;2;0]             |   -3.7703|          27|      14.1852|           0.0000|      3.0| 0.9780|
## |67  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;3;0]                                      |    2.5239|          19|      11.0000|           0.0000|      7.0| 1.0000|
## |68  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;0;0]                                      |    0.4413|          19|      15.4211|           0.0000|      4.0| 1.0000|
## |69  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;3;0]                                   |   -5.7919|          22|      13.6818|           0.0000|      7.0| 1.0000|
## |70  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                                   |   -9.2195|          22|      13.2727|           0.0000|      3.0| 0.9918|
## |71  |NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                                 |   -1.8672|           9|       7.8889|           0.0000|      9.0| 1.0000|
## |72  |N(HexNAc)GSYPNLSK[1;5;4;1;0]                                                                 |   -1.8672|           9|       5.7778|           0.0000|      9.0| 0.9996|
## |73  |N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                                 |   -1.5766|           9|       6.4444|           0.0000|     10.0| 0.9998|
## |74  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                            |   -5.9097|          18|      14.1111|           0.0000|      8.0| 1.0000|
## |75  |N(HexNAc)GSYPN(HexNAc)LSK[1;7;6;1;0]                                                         |   -0.9177|           9|       6.3333|           0.0000|      8.0| 0.9998|
## |76  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;6;6;0;0]                           |    0.1096|          22|      13.4545|           0.0000|      8.0| 1.0000|
## |77  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                            |   -3.3875|          18|      15.5000|           0.0000|      4.0| 1.0000|
## |78  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;1;0]                             |   -5.4029|          28|      15.2500|           0.0000|      1.0| 0.9802|
## |79  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;0;0]                                   |   -9.2195|          22|      12.7273|           0.0000|      3.0| 0.9908|
## |80  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                            |   -4.4787|          18|      11.5000|           0.0000|      2.0| 0.9994|
## |81  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;1;0]                                      |   -4.0059|          19|      11.0000|           0.0000|      4.0| 1.0000|
## |82  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;0;0]                                      |   -0.4635|          19|      12.6842|           0.0000|      3.0| 1.0000|
## |83  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;0;0]                           |   -0.8716|          22|      12.5000|           0.0000|      8.0| 1.0000|
## |84  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                            |   -2.7455|          18|      13.2778|           0.0000|      6.0| 1.0000|
## |85  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;1;0]                                            |   -3.9412|          18|      14.5556|           0.0000|      9.0| 1.0000|
## |86  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;2;0]                                            |   -6.5906|          18|      12.8333|           0.0000|      7.0| 1.0000|
## |87  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                            |   -6.7784|          18|      14.5000|           0.0000|      6.0| 1.0000|
## |88  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                            |   -8.5683|          18|      14.2778|           0.0000|      8.0| 0.9996|
## |89  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;1;0]                                                        |   -6.4433|          18|       9.7778|           0.0000|      5.0| 1.0000|
## |90  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;1;0]                                      |    0.2467|          19|      10.6316|           0.0000|      6.0| 1.0000|
## |91  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;1;0]                                      |   -3.4625|          19|      10.2632|           0.0000|      3.0| 1.0000|
## |92  |NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                                 |   -4.2508|           9|       6.1111|           0.0000|      7.0| 0.9998|
## |93  |N(HexNAc)GSYPNLSK[2;4;5;0;0]                                                                 |   -4.2508|           9|       5.4444|           0.0000|      7.0| 1.0000|
## |94  |NGSYPN(HexNAc)LSK[1;6;5;1;0]                                                                 |   -1.5766|           9|       6.2222|           0.0000|     10.0| 0.9998|
## |95  |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]                 |   -0.0952|          28|      13.0000|           0.0000|      1.0| 0.9894|
## |96  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;0;0]                                                        |   -4.0895|          18|      10.4444|           0.0000|      5.0| 1.0000|
## |97  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;0;0]                                      |   -3.0749|          19|      11.8947|           0.0000|      3.0| 1.0000|
## |98  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;0;0]                                            |   -3.0106|          18|      12.0000|           0.0000|      6.0| 1.0000|
## |99  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;3;0]             |    6.3070|          27|      10.2222|           0.0000|      2.0| 0.9320|
## |100 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                            |   -4.7866|          18|      10.0000|           0.0000|      2.0| 0.9996|
## |101 |NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                                 |   -1.3818|           9|       6.7778|           0.0000|      6.0| 1.0000|
## |102 |NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                                 |   -1.0611|           9|       6.1111|           0.0000|      9.0| 0.9998|
## |103 |N(HexNAc)GSYPNLSK[1;6;5;2;0]                                                                 |   -1.0611|           9|       5.2222|           0.0000|      9.0| 1.0000|
## |104 |N(HexNAc)GSYPNLSK[2;5;4;0;0]                                                                 |   -1.3818|           9|       5.3333|           0.0000|      6.0| 1.0000|
## |105 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;2;0]             |    4.7947|          27|      16.3333|           0.0000|      0.0| 0.9964|
## |106 |N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                                 |   -1.0866|           9|       4.1111|           0.0000|      8.0| 0.9984|
## |107 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;0;0]                                      |   -6.0673|          19|       7.6316|           0.0000|      2.0| 0.9990|
## |108 |N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                                     |    9.4307|           9|       5.5556|           0.0000|      1.0| 0.8534|
## |109 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;6;6;1;0]             |   -6.6991|          27|      11.8148|           0.0000|      1.0| 0.9966|
## |110 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[3;7;6;1;0]                           |   -8.1037|          22|       6.6818|           0.0000|      1.0| 0.7948|
## |111 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;1;0]                                            |   -4.3064|          18|       9.2778|           0.0000|      5.0| 1.0000|
## |112 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;6;0;0]                                                        |   -6.7592|          18|       8.3333|           0.0000|      4.0| 1.0000|
## |113 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;1;0]                                            |   -8.2882|          18|       9.6667|           0.0000|      5.0| 0.9992|
## |114 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]                 |   -1.6746|          28|       8.5000|           0.0000|      1.0| 0.9910|
## |115 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;2;0]                                      |    1.0948|          19|       7.0526|           0.0000|      2.0| 0.9968|
## |116 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;1;0]                                            |   -7.0563|          18|      10.3333|           0.0000|      8.0| 0.9988|
## |117 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                                   |   -1.6205|          22|       6.8182|           0.0000|      6.0| 1.0000|
## |118 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;1;0]                                   |   -1.6205|          22|       6.8182|           0.0000|      6.0| 1.0000|
## |119 |NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                                 |   -1.8675|           9|       5.1111|           0.0000|      8.0| 1.0000|
## |120 |N(HexNAc)GSYPNLSK[1;5;4;0;0]                                                                 |   -1.8675|           9|       5.1111|           0.0000|      8.0| 1.0000|
## |121 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;0;0]                                                        |   -4.6859|          18|       8.1667|           0.0000|      5.0| 1.0000|
## |122 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;0;0]                                                         |   -3.5034|           9|       5.3333|           0.0000|      8.0| 1.0000|
## |123 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                                     |   -1.3280|           9|       5.5556|           0.0000|      6.0| 1.0000|
## |124 |NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                                 |   -0.0550|           9|       4.8889|           0.0000|      6.0| 1.0000|
## |125 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]                       |   -9.7324|          22|       7.6364|           0.0000|      3.0| 0.9662|
## |126 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]                       |   -9.7324|          22|       7.5455|           0.0000|      3.0| 0.9654|
## |127 |N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                                 |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |128 |NGSYPN(HexNAc)LSK[1;5;6;1;0]                                                                 |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |129 |NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                                 |   -1.0866|           9|       3.8889|           0.0000|      8.0| 1.0000|
## |130 |N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                                 |   -1.5617|           9|       4.5556|           0.0000|      6.0| 0.9998|
## |131 |NGSYPN(HexNAc)LSK[1;7;6;2;0]                                                                 |   -1.5617|           9|       4.5556|           0.0000|      6.0| 0.9998|
## |132 |N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                                 |   -0.0550|           9|       4.0000|           0.0000|      6.0| 0.9998|
## |133 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                                            |   -2.2699|          18|       8.0000|           0.0000|      2.0| 0.9990|
## |134 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;1;0]             |    0.1401|          27|       9.7778|           0.0000|      0.0| 0.9942|
## |135 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;0;0]                                            |   -2.7114|          18|       8.4444|           0.0000|      6.0| 1.0000|
## |136 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;0;0]                                            |   -2.7455|          18|       8.8889|           0.0000|      6.0| 1.0000|
## |137 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;9;0;0]                                                         |   -4.4694|           9|       6.5556|           0.0000|      7.0| 1.0000|
## |138 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;3;0]                                                        |   -5.3203|          18|       6.9444|           0.0000|      4.0| 0.9998|
## |139 |N(HexNAc)GSYPN(HexNAc)LSK[0;14;4;0;0]                                                        |    7.7308|           9|       4.5556|           0.0000|      5.0| 0.9854|
## |140 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;1;0]                                      |    1.6958|          19|       5.0526|           0.0000|      0.0| 0.8108|
## |141 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;8;0;0]                                                         |   -3.0549|           9|       3.2222|           0.0000|     10.0| 0.9994|
## |142 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[2;7;6;2;0]                 |    1.2589|          28|       5.8571|           0.0000|      0.0| 0.9004|
## |143 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;5;4;1;0]                                      |   -2.5981|          19|       3.9474|           0.1053|      0.0| 0.7054|
## |144 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                                            |   -4.4787|          18|       7.5000|           0.0000|      2.0| 0.9996|
## |145 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]                 |    4.3525|          28|       6.5000|           0.0000|      0.0| 0.9420|
## |146 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]                 |   -3.4230|          28|       5.1786|           0.0000|      0.0| 0.8418|
## |147 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]                       |   -5.5353|          22|       5.2727|           0.0000|      5.0| 0.9998|
## |148 |N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                                     |    8.1210|           9|       3.8889|           0.1111|      0.0| 0.8118|
## |149 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;1;0]                                      |   -4.0430|          19|       6.8947|           0.0000|      2.0| 0.9996|
## |150 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;4;1;0]                                      |   -3.2963|          19|       4.8947|           0.0000|      1.0| 0.7406|
## |151 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;3;0]                                      |   -0.6209|          19|       5.9474|           0.0000|      3.0| 0.9992|
## |152 |NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                                 |   -7.2148|           9|       3.3333|           0.1111|      4.0| 0.9974|
## |153 |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]                             |   -3.9020|          28|       7.1786|           0.0000|      0.0| 0.9736|
## |154 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;0;0]                                      |   -4.3502|          19|       5.5789|           0.0000|      2.0| 0.9986|
## |155 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;0;0]                                            |   -4.7866|          18|       7.7222|           0.0000|      2.0| 0.9996|
## |156 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]                 |   -3.9614|          28|       6.3214|           0.0000|      0.0| 0.9454|
## |157 |NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]               |   -3.7097|          23|       8.6957|           0.0000|      0.0| 0.7336|
## |158 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                                   |    5.7869|          22|       5.5909|           0.0000|      4.0| 1.0000|
## |159 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;6;1;0]                                   |    5.7869|          22|       5.5909|           0.0000|      4.0| 1.0000|
## |160 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;2;0]                       |   -5.5353|          22|       4.7273|           0.0000|      5.0| 1.0000|
## |161 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;5;1;0]             |   -8.0917|          27|       6.0741|           0.0000|      0.0| 0.9320|
## |162 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;0;0]                                            |   -3.3875|          18|       9.1111|           0.0000|      4.0| 1.0000|
## |163 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;0;0]                                      |    1.6531|          19|       3.1053|           0.1053|      1.0| 0.7062|
## |164 |N(Deamidated)LLWLTEKN(HexNAc)GSYPNLSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]               |   -3.7097|          23|       6.9565|           0.0000|      0.0| 0.6428|
## |165 |NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                                 |    3.0849|           9|       3.6667|           0.0000|      6.0| 0.9996|
## |166 |N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                                 |    1.8418|           9|       2.5556|           0.0000|      8.0| 1.0000|
## |167 |NGSYPN(HexNAc)LSK[1;6;5;3;0]                                                                 |    1.8418|           9|       2.5556|           0.0000|      8.0| 1.0000|
## |168 |N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                                 |    3.0849|           9|       3.0000|           0.0000|      6.0| 1.0000|
## |169 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;6;1;0]             |   -9.7402|          27|       5.0741|           0.0000|      0.0| 0.8824|
## |170 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;1;0]                                                     |   -1.3280|           9|       2.6667|           0.0000|      6.0| 1.0000|
## |171 |HN(HexNAc)VTR[1;6;6;1;0]                                                                     |    3.5881|           5|       0.8000|           0.4000|      1.8| 0.9984|
## |172 |N(HexNAc)GSYPN(HexNAc)LSK[2;12;6;0;0]                                                        |    2.9929|           9|       2.6667|           0.0000|      7.0| 0.9994|
## |173 |HN(HexNAc)VTR[1;6;5;0;0]                                                                     |    0.6254|           5|       1.0000|           0.0000|      2.7| 0.9998|
## |174 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;6;0;0]                                      |   -3.5028|          19|       7.2105|           0.0000|      2.0| 1.0000|
## |176 |N(HexNAc)GSYPN(Deamidated)LSK[0;6;6;1;0]                                                     |    9.4307|           9|       1.7778|           0.0000|      1.0| 0.6900|
## |177 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;7;6;0;0]                           |   -6.6968|          22|       4.6818|           0.0000|      2.0| 0.9204|
## |178 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]                       |   -2.1214|          22|       4.5909|           0.0000|      4.0| 0.9954|
## |179 |N(HexNAc)GSYPNLSK[2;3;6;0;0]                                                                 |   -7.2148|           9|       2.4444|           0.2222|      4.0| 0.9972|
## |180 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                                   |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9998|
## |181 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;0;0]                                   |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9998|
## |184 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;6;1;0]             |   -0.2344|          27|       4.2222|           0.1111|      0.0| 0.8190|
## |186 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHN(Deamidated)GK[1;7;6;3;0]                                |   -5.4866|          18|       2.5556|           0.0000|      2.0| 0.9176|
## |187 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;0;0]                                      |   -5.9074|          19|       4.3158|           0.0000|      1.0| 0.6740|
## |189 |N(HexNAc)GSYPN(HexNAc)LSK[2;8;6;1;0]                                                         |   -2.9545|           9|       2.4444|           0.1111|      8.0| 0.9996|
## |190 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;0;0]                                   |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9998|
## |191 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                                   |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9996|
## |193 |HN(HexNAc)VTR[0;4;2;0;0]                                                                     |   -8.8012|           5|       1.6000|           0.0000|      1.8| 0.8206|
## |201 |N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                                     |    2.5274|           9|       3.0000|           0.1111|      4.0| 0.9990|
## |202 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;7;0;0]                                                         |   -0.6787|           9|       2.5556|           0.0000|      5.0| 0.9998|
## |203 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]                       |   -2.1214|          22|       3.5909|           0.0000|      4.0| 0.9974|
## |204 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                                   |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9996|
## |205 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;1;0]                                   |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9992|
## |206 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                            |   -4.9256|          18|       5.7778|           0.0000|      3.0| 0.9994|
## |208 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)KEK[4;8;11;2;0]                                    |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9974|
## |209 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                                    |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9992|
## |212 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                                     |    2.8345|           9|       2.1111|           0.1111|      6.0| 0.9998|
## |213 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQN(Deamidated)IHPVTIGEC(Carbamidomethyl)PK[1;6;5;3;0] |    5.5752|          27|       3.0741|           0.1481|      0.0| 0.6370|
## |214 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;9;0;0]                                                         |   -0.7798|           9|       2.2222|           0.0000|      5.0| 0.9998|
## |215 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[3;8;8;1;0]                           |   -2.2628|          15|       4.7333|           0.0000|      2.0| 0.9838|
## |216 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                                     |   -0.2389|           9|       2.2222|           0.0000|      5.0| 0.9998|
## |217 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;0;0]                                                         |   -1.8379|           9|       3.5556|           0.0000|      8.0| 0.9994|
## |233 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;1;0]                                                        |   -6.3796|          18|       2.8333|           0.0000|      3.0| 0.9978|
## |234 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;1;0]                                                         |   -1.2417|           9|       2.0000|           0.0000|      9.0| 0.9964|
## |245 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;1;0]                                                         |    6.9927|           9|       1.7778|           0.3333|      6.0| 0.9990|
## |248 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;7;1;0]                                                         |   -0.4717|           9|       1.3333|           0.0000|      5.0| 0.9998|
## |250 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                                     |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9974|
## |251 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[0;10;10;0;0]                                     |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9974|
## |253 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;7;0;0]                                            |   -4.9256|          18|       5.2222|           0.0000|      3.0| 1.0000|
## |255 |N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                                     |    0.9561|           9|       1.4444|           0.1111|      2.0| 0.9880|
## |264 |N(HexNAc)GSYPN(HexNAc)LSK[1;6;6;0;0]                                                         |   -1.6219|           9|       2.0000|           0.1111|      5.0| 0.9998|
## |265 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                                   |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9868|
## |266 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;5;1;0]                                   |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9868|
## |282 |N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                     |   -2.4248|           9|       2.0000|           0.1111|      1.0| 0.6638|
## |291 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;8;6;0;0]                           |    5.8987|          22|       5.4545|           0.0000|      3.0| 0.9998|
## |305 |HN(HexNAc)VTR[1;5;4;0;0]                                                                     |   -2.0406|           5|       0.2000|           0.8000|      2.7| 0.9968|
## |311 |HN(HexNAc)VTR[0;7;2;0;0]                                                                     |   -3.4059|           5|       0.6000|           0.4000|      1.8| 0.9986|
## |312 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                         |   -2.1511|           9|       0.8889|           0.1111|      5.0| 1.0000|
## |313 |HN(HexNAc)VTR[1;5;5;0;0]                                                                     |   -3.9621|           5|       0.6000|           0.4000|      2.1| 0.9984|
## |319 |N(HexNAc)GSYPN(HexNAc)LSK[4;9;11;2;0]                                                        |    7.9440|           9|       0.8889|           0.1111|      3.0| 0.9570|
## |320 |N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                                 |   -1.6596|           9|       1.2222|           0.3333|      4.0| 1.0000|
## |321 |NGSYPN(HexNAc)LSK[1;5;5;0;0]                                                                 |   -1.6596|           9|       1.2222|           0.3333|      4.0| 1.0000|
## |327 |N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                                 |   -1.4113|           9|       0.5556|           0.4444|      3.0| 1.0000|
## |328 |NGSYPN(HexNAc)LSK[1;5;6;2;0]                                                                 |   -1.4113|           9|       0.5556|           0.4444|      3.0| 1.0000|
## |329 |N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                                 |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9902|
## |330 |NGSYPN(HexNAc)LSK[1;4;5;0;0]                                                                 |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9902|
## |340 |N(HexNAc)GSYPN(HexNAc)LSK[1;9;6;1;0]                                                         |   -2.0313|           9|       0.5556|           0.4444|      4.0| 0.9998|
## |345 |N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                                 |   -3.0949|           9|       0.5556|           0.4444|      4.0| 1.0000|
## |346 |NGSYPN(HexNAc)LSK[1;6;6;2;0]                                                                 |   -3.0949|           9|       0.5556|           0.4444|      4.0| 1.0000|
## |350 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;0;0]                                                     |   -0.2389|           9|       0.7778|           0.2222|      5.0| 0.9992|
## |351 |N(HexNAc)GSYPN(HexNAc)LSK[0;9;8;1;0]                                                         |    7.0308|           9|       0.5556|           0.4444|      3.0| 0.9992|
## |454 |N(HexNAc)GSYPN(HexNAc)LSK[0;12;4;0;0]                                                        |    8.6872|           9|       0.5556|           0.4444|      3.0| 0.9336|
## |471 |N(HexNAc)GSYPN(HexNAc)LSK[0;13;4;0;0]                                                        |    6.5838|           9|       0.5556|           0.4444|      5.0| 0.9996|
## |485 |N(HexNAc)GSYPN(Deamidated)LSK[1;7;6;2;0]                                                     |    2.5274|           9|       0.0000|           1.0000|      4.0| 0.9998|
## |487 |N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                                     |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9994|
## |488 |N(HexNAc)GSYPN(Deamidated)LSK[1;4;5;1;0]                                                     |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9994|
## |491 |N(HexNAc)GSYPN(Deamidated)LSK[1;6;5;1;0]                                                     |    0.9561|           9|       0.0000|           1.0000|      2.0| 0.9828|
## |492 |HN(HexNAc)VTR[2;4;5;0;0]                                                                     |   -5.3622|           5|       0.0000|           1.0000|      3.0| 0.9926|
## |516 |HN(HexNAc)VTR[0;5;2;0;0]                                                                     |   -6.1880|           5|       0.0000|           1.0000|      1.8| 0.9678|
## |525 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;2;0]                                                     |    2.8345|           9|       0.0000|           1.0000|      6.0| 0.9994|
## |554 |HN(HexNAc)VTR[1;4;5;1;0]                                                                     |   -3.1177|           5|       0.0000|           1.0000|      1.8| 0.9944|
## |562 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;0;0]                                      |    1.0918|          19|       0.7895|           0.5263|      2.0| 0.9236|
## |627 |N(HexNAc)GSYPN(HexNAc)LSK[0;6;6;0;0]                                                         |   -2.6125|           9|       0.0000|           1.0000|      5.0| 0.9996|
## |687 |HN(HexNAc)VTR[2;5;4;0;0]                                                                     |   -2.7862|           5|       0.0000|           1.0000|      1.8| 0.9952|
## 
## 
## |    |Glycopeptide_identifier                                                          | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   tPos|
## |:---|:--------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |1   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                |   -5.9097|          18|      22.5556|           0.0000|      8.0| 1.0000|
## |2   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                       |    1.3485|          22|      26.5909|           0.0000|      8.0| 1.0000|
## |3   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;2;0]                       |    1.3485|          22|      26.0455|           0.0000|      8.0| 1.0000|
## |4   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                |   -3.9412|          18|      20.2222|           0.0000|      9.0| 1.0000|
## |5   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                |   -6.5906|          18|      18.0000|           0.0000|      7.0| 0.9982|
## |6   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                |   -8.5683|          18|      21.6667|           0.0000|      8.0| 0.9994|
## |7   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                |   -6.7784|          18|      18.5000|           0.0000|      6.0| 0.9950|
## |8   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                       |    2.8744|          22|      25.0455|           0.0000|      9.0| 1.0000|
## |9   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;2;0]                       |    2.8744|          22|      24.5000|           0.0000|      9.0| 1.0000|
## |10  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;1;0]               |    2.9786|          22|      22.5000|           0.0000|      9.0| 1.0000|
## |11  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;1;0]                          |   -3.9610|          19|      22.8947|           0.0000|      6.0| 1.0000|
## |12  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                       |   -2.7807|          22|      20.6364|           0.0000|      8.0| 1.0000|
## |13  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                       |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |14  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                       |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |15  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;1;0]                          |   -0.8840|          19|      22.7895|           0.0000|      9.0| 1.0000|
## |16  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;14;4;0;0]              |    9.6761|          22|      29.1364|           0.0000|      7.0| 0.8816|
## |17  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;0;0]                          |   -2.7961|          19|      20.9474|           0.0000|      5.0| 1.0000|
## |18  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;0;0]                       |   -2.7807|          22|      20.0909|           0.0000|      8.0| 1.0000|
## |19  |NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                     |    0.1308|           9|       9.8889|           0.0000|     10.0| 0.9994|
## |20  |N(HexNAc)GSYPNLSK[1;5;4;2;0]                                                     |    0.1308|           9|       9.0000|           0.0000|     10.0| 0.9994|
## |21  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;3;0]                          |   -2.8410|          19|      21.0526|           0.0000|      8.0| 1.0000|
## |22  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                       |   -0.1015|          22|      22.5000|           0.0000|      7.0| 0.9994|
## |23  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;1;0] |   -4.3482|          27|      21.9630|           0.0000|      6.0| 0.9956|
## |24  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;2;0]                 |   -1.5373|          28|      22.3929|           0.0000|      1.0| 0.9722|
## |25  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                       |   -0.1015|          22|      21.9545|           0.0000|      7.0| 0.9994|
## |26  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;1;0]                          |   -3.5609|          19|      18.5789|           0.0000|      8.0| 1.0000|
## |27  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;5;4;2;0] |   -0.1807|          27|      22.2222|           0.0000|      5.0| 0.9976|
## |28  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;2;0]                          |    0.7980|          19|      18.5789|           0.0000|      5.0| 1.0000|
## |29  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]                 |   -0.4458|          28|      22.9643|           0.0000|      2.0| 0.9980|
## |30  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]                 |    0.4326|          28|      22.5000|           0.0000|      1.0| 0.9596|
## |31  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                |   -3.0106|          18|      16.5556|           0.0000|      6.0| 1.0000|
## |32  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                       |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |33  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;0;0]                       |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |34  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;2;0]                          |   -0.0219|          19|      19.6842|           0.0000|      8.0| 0.9994|
## |35  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]                 |   -2.8434|          28|      20.0000|           0.0000|      2.0| 0.9898|
## |36  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;5;4;1;0]                                            |   -7.4680|          18|      12.3889|           0.0000|      4.0| 0.9996|
## |37  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;2;0]                          |    0.7444|          19|      18.9474|           0.0000|      5.0| 1.0000|
## |38  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;2;0]                          |   -1.2420|          19|      17.6842|           0.0000|      5.0| 1.0000|
## |39  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]           |   -1.5816|          22|      15.9545|           0.0000|      6.0| 1.0000|
## |40  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]           |   -1.5416|          22|      17.0000|           0.0000|      6.0| 1.0000|
## |41  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;3;0]           |   -1.5816|          22|      15.5000|           0.0000|      6.0| 1.0000|
## |42  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;2;0]                          |    5.9064|          19|      15.6842|           0.0000|      5.0| 1.0000|
## |43  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]                 |   -0.7608|          28|      19.1429|           0.0000|      2.0| 0.9980|
## |44  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                |   -2.7114|          18|      14.1111|           0.0000|      6.0| 1.0000|
## |45  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;4;0]                          |   -3.7451|          19|      17.9474|           0.0000|      6.0| 1.0000|
## |46  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;2;0]           |   -1.5416|          22|      15.9091|           0.0000|      6.0| 1.0000|
## |47  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;3;0]                          |   -2.2429|          19|      17.1053|           0.0000|      8.0| 1.0000|
## |48  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                       |   -4.0561|          22|      14.7727|           0.0000|      6.0| 1.0000|
## |49  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;2;0]                       |   -4.0561|          22|      14.7727|           0.0000|      6.0| 1.0000|
## |50  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                |   -8.2882|          18|      14.7222|           0.0000|      5.0| 0.9956|
## |51  |NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                     |   -1.4933|           9|       8.3333|           0.0000|     10.0| 1.0000|
## |52  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;2;0]                          |   -0.5044|          19|      15.4737|           0.0000|      4.0| 1.0000|
## |53  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                |   -7.0563|          18|      12.6667|           0.0000|      8.0| 0.9856|
## |54  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;1;0]                 |   -4.1984|          28|      18.5714|           0.0000|      2.0| 0.9894|
## |55  |NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                     |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |56  |HN(HexNAc)VTR[1;5;4;1;0]                                                         |   -1.3984|           5|       2.8000|           0.0000|      3.0| 0.9968|
## |57  |N(HexNAc)GSYPNLSK[1;4;5;1;0]                                                     |   -1.4933|           9|       6.8889|           0.0000|     10.0| 0.9784|
## |58  |N(HexNAc)GSYPNLSK[1;6;5;0;0]                                                     |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |59  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                                |   -2.2699|          18|      14.2222|           0.0000|      2.0| 1.0000|
## |60  |N(HexNAc)GSYPN(HexNAc)LSK[2;6;6;0;0]                                             |   -1.0568|           9|       9.6667|           0.0000|      9.0| 1.0000|
## |61  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                |   -4.3064|          18|      15.0556|           0.0000|      5.0| 1.0000|
## |62  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                       |   -5.7919|          22|      14.2273|           0.0000|      7.0| 1.0000|
## |63  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;1;0]                       |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |64  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                       |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |65  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;1;0]                          |   -2.0419|          19|      13.8947|           0.0000|      3.0| 1.0000|
## |66  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;5;2;0] |   -3.7703|          27|      14.1852|           0.0000|      3.0| 0.9094|
## |67  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;3;0]                          |    2.5239|          19|      11.0000|           0.0000|      7.0| 1.0000|
## |68  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;0;0]                          |    0.4413|          19|      15.4211|           0.0000|      4.0| 1.0000|
## |69  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;3;0]                       |   -5.7919|          22|      13.6818|           0.0000|      7.0| 1.0000|
## |70  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                       |   -9.2195|          22|      13.2727|           0.0000|      3.0| 0.8270|
## |71  |NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                     |   -1.8672|           9|       7.8889|           0.0000|      9.0| 0.9996|
## |72  |N(HexNAc)GSYPNLSK[1;5;4;1;0]                                                     |   -1.8672|           9|       5.7778|           0.0000|      9.0| 1.0000|
## |73  |N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                     |   -1.5766|           9|       6.4444|           0.0000|     10.0| 0.9986|
## |74  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                |   -5.9097|          18|      14.1111|           0.0000|      8.0| 1.0000|
## |75  |N(HexNAc)GSYPN(HexNAc)LSK[1;7;6;1;0]                                             |   -0.9177|           9|       6.3333|           0.0000|      8.0| 0.9992|
## |76  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;6;6;0;0]               |    0.1096|          22|      13.4545|           0.0000|      8.0| 0.9990|
## |77  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                |   -3.3875|          18|      15.5000|           0.0000|      4.0| 1.0000|
## |78  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;1;0]                 |   -5.4029|          28|      15.2500|           0.0000|      1.0| 0.9978|
## |79  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;0;0]                       |   -9.2195|          22|      12.7273|           0.0000|      3.0| 0.8268|
## |80  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                |   -4.4787|          18|      11.5000|           0.0000|      2.0| 1.0000|
## |81  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;1;0]                          |   -4.0059|          19|      11.0000|           0.0000|      4.0| 0.9998|
## |82  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;0;0]                          |   -0.4635|          19|      12.6842|           0.0000|      3.0| 1.0000|
## |83  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;0;0]               |   -0.8716|          22|      12.5000|           0.0000|      8.0| 1.0000|
## |84  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                |   -2.7455|          18|      13.2778|           0.0000|      6.0| 1.0000|
## |85  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;1;0]                                |   -3.9412|          18|      14.5556|           0.0000|      9.0| 1.0000|
## |86  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;2;0]                                |   -6.5906|          18|      12.8333|           0.0000|      7.0| 0.9982|
## |87  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                |   -6.7784|          18|      14.5000|           0.0000|      6.0| 0.9950|
## |88  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                |   -8.5683|          18|      14.2778|           0.0000|      8.0| 0.9994|
## |89  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;1;0]                                            |   -6.4433|          18|       9.7778|           0.0000|      5.0| 0.9944|
## |90  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;1;0]                          |    0.2467|          19|      10.6316|           0.0000|      6.0| 1.0000|
## |91  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;1;0]                          |   -3.4625|          19|      10.2632|           0.0000|      3.0| 1.0000|
## |92  |NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                     |   -4.2508|           9|       6.1111|           0.0000|      7.0| 0.9994|
## |93  |N(HexNAc)GSYPNLSK[2;4;5;0;0]                                                     |   -4.2508|           9|       5.4444|           0.0000|      7.0| 0.9994|
## |94  |NGSYPN(HexNAc)LSK[1;6;5;1;0]                                                     |   -1.5766|           9|       6.2222|           0.0000|     10.0| 0.9992|
## |95  |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]     |   -0.0952|          28|      13.0000|           0.0000|      1.0| 0.9478|
## |96  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;0;0]                                            |   -4.0895|          18|      10.4444|           0.0000|      5.0| 1.0000|
## |97  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;0;0]                          |   -3.0749|          19|      11.8947|           0.0000|      3.0| 0.9994|
## |98  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;0;0]                                |   -3.0106|          18|      12.0000|           0.0000|      6.0| 1.0000|
## |99  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;3;0] |    6.3070|          27|      10.2222|           0.0000|      2.0| 0.8130|
## |100 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                |   -4.7866|          18|      10.0000|           0.0000|      2.0| 1.0000|
## |101 |NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                     |   -1.3818|           9|       6.7778|           0.0000|      6.0| 0.9978|
## |102 |NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                     |   -1.0611|           9|       6.1111|           0.0000|      9.0| 0.9984|
## |103 |N(HexNAc)GSYPNLSK[1;6;5;2;0]                                                     |   -1.0611|           9|       5.2222|           0.0000|      9.0| 0.9980|
## |104 |N(HexNAc)GSYPNLSK[2;5;4;0;0]                                                     |   -1.3818|           9|       5.3333|           0.0000|      6.0| 0.9974|
## |105 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;2;0] |    4.7947|          27|      16.3333|           0.0000|      0.0| 0.9988|
## |106 |N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                     |   -1.0866|           9|       4.1111|           0.0000|      8.0| 1.0000|
## |107 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;0;0]                          |   -6.0673|          19|       7.6316|           0.0000|      2.0| 0.9982|
## |108 |N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                         |    9.4307|           9|       5.5556|           0.0000|      1.0| 0.6746|
## |109 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;6;6;1;0] |   -6.6991|          27|      11.8148|           0.0000|      1.0| 0.9960|
## |110 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[3;7;6;1;0]               |   -8.1037|          22|       6.6818|           0.0000|      1.0| 0.8650|
## |111 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;1;0]                                |   -4.3064|          18|       9.2778|           0.0000|      5.0| 1.0000|
## |112 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;6;0;0]                                            |   -6.7592|          18|       8.3333|           0.0000|      4.0| 0.9900|
## |113 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;1;0]                                |   -8.2882|          18|       9.6667|           0.0000|      5.0| 0.9956|
## |114 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]     |   -1.6746|          28|       8.5000|           0.0000|      1.0| 0.9880|
## |115 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;2;0]                          |    1.0948|          19|       7.0526|           0.0000|      2.0| 0.9974|
## |116 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;1;0]                                |   -7.0563|          18|      10.3333|           0.0000|      8.0| 0.9856|
## |117 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                       |   -1.6205|          22|       6.8182|           0.0000|      6.0| 0.9918|
## |118 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;1;0]                       |   -1.6205|          22|       6.8182|           0.0000|      6.0| 0.9896|
## |119 |NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                     |   -1.8675|           9|       5.1111|           0.0000|      8.0| 0.9996|
## |120 |N(HexNAc)GSYPNLSK[1;5;4;0;0]                                                     |   -1.8675|           9|       5.1111|           0.0000|      8.0| 0.9998|
## |121 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;0;0]                                            |   -4.6859|          18|       8.1667|           0.0000|      5.0| 0.9992|
## |122 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;0;0]                                             |   -3.5034|           9|       5.3333|           0.0000|      8.0| 0.9976|
## |123 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                         |   -1.3280|           9|       5.5556|           0.0000|      6.0| 0.9996|
## |124 |NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                     |   -0.0550|           9|       4.8889|           0.0000|      6.0| 0.9972|
## |125 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]           |   -9.7324|          22|       7.6364|           0.0000|      3.0| 0.8244|
## |126 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]           |   -9.7324|          22|       7.5455|           0.0000|      3.0| 0.8244|
## |127 |N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                     |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |128 |NGSYPN(HexNAc)LSK[1;5;6;1;0]                                                     |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |129 |NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                     |   -1.0866|           9|       3.8889|           0.0000|      8.0| 1.0000|
## |130 |N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                     |   -1.5617|           9|       4.5556|           0.0000|      6.0| 0.9968|
## |131 |NGSYPN(HexNAc)LSK[1;7;6;2;0]                                                     |   -1.5617|           9|       4.5556|           0.0000|      6.0| 0.9968|
## |132 |N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                     |   -0.0550|           9|       4.0000|           0.0000|      6.0| 0.9990|
## |133 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                                |   -2.2699|          18|       8.0000|           0.0000|      2.0| 1.0000|
## |134 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;1;0] |    0.1401|          27|       9.7778|           0.0000|      0.0| 0.9828|
## |135 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;0;0]                                |   -2.7114|          18|       8.4444|           0.0000|      6.0| 0.9984|
## |136 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;0;0]                                |   -2.7455|          18|       8.8889|           0.0000|      6.0| 1.0000|
## |137 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;9;0;0]                                             |   -4.4694|           9|       6.5556|           0.0000|      7.0| 0.9998|
## |138 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;3;0]                                            |   -5.3203|          18|       6.9444|           0.0000|      4.0| 0.9940|
## |139 |N(HexNAc)GSYPN(HexNAc)LSK[0;14;4;0;0]                                            |    7.7308|           9|       4.5556|           0.0000|      5.0| 0.9696|
## |140 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;1;0]                          |    1.6958|          19|       5.0526|           0.0000|      0.0| 0.8486|
## |141 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;8;0;0]                                             |   -3.0549|           9|       3.2222|           0.0000|     10.0| 0.9996|
## |142 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[2;7;6;2;0]     |    1.2589|          28|       5.8571|           0.0000|      0.0| 0.9236|
## |144 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                                |   -4.4787|          18|       7.5000|           0.0000|      2.0| 0.9992|
## |145 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]     |    4.3525|          28|       6.5000|           0.0000|      0.0| 0.9518|
## |146 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]     |   -3.4230|          28|       5.1786|           0.0000|      0.0| 0.9564|
## |147 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]           |   -5.5353|          22|       5.2727|           0.0000|      5.0| 0.9960|
## |148 |N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                         |    8.1210|           9|       3.8889|           0.1111|      0.0| 0.5156|
## |149 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;1;0]                          |   -4.0430|          19|       6.8947|           0.0000|      2.0| 0.9772|
## |151 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;3;0]                          |   -0.6209|          19|       5.9474|           0.0000|      3.0| 0.9986|
## |152 |NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                     |   -7.2148|           9|       3.3333|           0.1111|      4.0| 0.9446|
## |153 |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]                 |   -3.9020|          28|       7.1786|           0.0000|      0.0| 0.9892|
## |154 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;0;0]                          |   -4.3502|          19|       5.5789|           0.0000|      2.0| 0.9986|
## |155 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;0;0]                                |   -4.7866|          18|       7.7222|           0.0000|      2.0| 1.0000|
## |156 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]     |   -3.9614|          28|       6.3214|           0.0000|      0.0| 0.9976|
## |157 |NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]   |   -3.7097|          23|       8.6957|           0.0000|      0.0| 0.9864|
## |158 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                       |    5.7869|          22|       5.5909|           0.0000|      4.0| 0.9974|
## |159 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;6;1;0]                       |    5.7869|          22|       5.5909|           0.0000|      4.0| 0.9964|
## |160 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;2;0]           |   -5.5353|          22|       4.7273|           0.0000|      5.0| 0.9940|
## |161 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;5;1;0] |   -8.0917|          27|       6.0741|           0.0000|      0.0| 0.8334|
## |162 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;0;0]                                |   -3.3875|          18|       9.1111|           0.0000|      4.0| 1.0000|
## |164 |N(Deamidated)LLWLTEKN(HexNAc)GSYPNLSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]   |   -3.7097|          23|       6.9565|           0.0000|      0.0| 0.9554|
## |165 |NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                     |    3.0849|           9|       3.6667|           0.0000|      6.0| 1.0000|
## |166 |N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                     |    1.8418|           9|       2.5556|           0.0000|      8.0| 0.9984|
## |167 |NGSYPN(HexNAc)LSK[1;6;5;3;0]                                                     |    1.8418|           9|       2.5556|           0.0000|      8.0| 0.9984|
## |168 |N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                     |    3.0849|           9|       3.0000|           0.0000|      6.0| 1.0000|
## |169 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;6;1;0] |   -9.7402|          27|       5.0741|           0.0000|      0.0| 0.7152|
## |170 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;1;0]                                         |   -1.3280|           9|       2.6667|           0.0000|      6.0| 0.9982|
## |171 |HN(HexNAc)VTR[1;6;6;1;0]                                                         |    3.5881|           5|       0.8000|           0.4000|      1.8| 1.0000|
## |172 |N(HexNAc)GSYPN(HexNAc)LSK[2;12;6;0;0]                                            |    2.9929|           9|       2.6667|           0.0000|      7.0| 0.9980|
## |173 |HN(HexNAc)VTR[1;6;5;0;0]                                                         |    0.6254|           5|       1.0000|           0.0000|      2.7| 1.0000|
## |174 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;6;0;0]                          |   -3.5028|          19|       7.2105|           0.0000|      2.0| 0.9958|
## |177 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;7;6;0;0]               |   -6.6968|          22|       4.6818|           0.0000|      2.0| 0.9116|
## |178 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]           |   -2.1214|          22|       4.5909|           0.0000|      4.0| 0.9874|
## |179 |N(HexNAc)GSYPNLSK[2;3;6;0;0]                                                     |   -7.2148|           9|       2.4444|           0.2222|      4.0| 0.9554|
## |180 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                       |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9878|
## |181 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;0;0]                       |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9874|
## |184 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;6;1;0] |   -0.2344|          27|       4.2222|           0.1111|      0.0| 0.5684|
## |186 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHN(Deamidated)GK[1;7;6;3;0]                    |   -5.4866|          18|       2.5556|           0.0000|      2.0| 0.9932|
## |189 |N(HexNAc)GSYPN(HexNAc)LSK[2;8;6;1;0]                                             |   -2.9545|           9|       2.4444|           0.1111|      8.0| 0.9992|
## |190 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;0;0]                       |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9882|
## |191 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                       |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9898|
## |193 |HN(HexNAc)VTR[0;4;2;0;0]                                                         |   -8.8012|           5|       1.6000|           0.0000|      1.8| 0.7210|
## |201 |N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                         |    2.5274|           9|       3.0000|           0.1111|      4.0| 0.9942|
## |202 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;7;0;0]                                             |   -0.6787|           9|       2.5556|           0.0000|      5.0| 0.9980|
## |203 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]           |   -2.1214|          22|       3.5909|           0.0000|      4.0| 0.9918|
## |204 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                       |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9782|
## |205 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;1;0]                       |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9802|
## |206 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                |   -4.9256|          18|       5.7778|           0.0000|      3.0| 1.0000|
## |208 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)KEK[4;8;11;2;0]                        |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9988|
## |209 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                        |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9990|
## |212 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                         |    2.8345|           9|       2.1111|           0.1111|      6.0| 1.0000|
## |214 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;9;0;0]                                             |   -0.7798|           9|       2.2222|           0.0000|      5.0| 0.9998|
## |215 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[3;8;8;1;0]               |   -2.2628|          15|       4.7333|           0.0000|      2.0| 0.9946|
## |216 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                         |   -0.2389|           9|       2.2222|           0.0000|      5.0| 1.0000|
## |217 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;0;0]                                             |   -1.8379|           9|       3.5556|           0.0000|      8.0| 1.0000|
## |233 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;1;0]                                            |   -6.3796|          18|       2.8333|           0.0000|      3.0| 0.9840|
## |234 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;1;0]                                             |   -1.2417|           9|       2.0000|           0.0000|      9.0| 0.9996|
## |245 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;1;0]                                             |    6.9927|           9|       1.7778|           0.3333|      6.0| 0.9628|
## |248 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;7;1;0]                                             |   -0.4717|           9|       1.3333|           0.0000|      5.0| 1.0000|
## |250 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                         |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9890|
## |251 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[0;10;10;0;0]                         |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9890|
## |253 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;7;0;0]                                |   -4.9256|          18|       5.2222|           0.0000|      3.0| 0.9976|
## |255 |N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                         |    0.9561|           9|       1.4444|           0.1111|      2.0| 0.9986|
## |264 |N(HexNAc)GSYPN(HexNAc)LSK[1;6;6;0;0]                                             |   -1.6219|           9|       2.0000|           0.1111|      5.0| 0.9998|
## |265 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                       |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9504|
## |266 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;5;1;0]                       |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9504|
## |291 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;8;6;0;0]               |    5.8987|          22|       5.4545|           0.0000|      3.0| 0.9880|
## |305 |HN(HexNAc)VTR[1;5;4;0;0]                                                         |   -2.0406|           5|       0.2000|           0.8000|      2.7| 0.9996|
## |311 |HN(HexNAc)VTR[0;7;2;0;0]                                                         |   -3.4059|           5|       0.6000|           0.4000|      1.8| 1.0000|
## |312 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;5;0;0]                                             |   -2.1511|           9|       0.8889|           0.1111|      5.0| 1.0000|
## |313 |HN(HexNAc)VTR[1;5;5;0;0]                                                         |   -3.9621|           5|       0.6000|           0.4000|      2.1| 1.0000|
## |319 |N(HexNAc)GSYPN(HexNAc)LSK[4;9;11;2;0]                                            |    7.9440|           9|       0.8889|           0.1111|      3.0| 0.9836|
## |320 |N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                     |   -1.6596|           9|       1.2222|           0.3333|      4.0| 0.9998|
## |321 |NGSYPN(HexNAc)LSK[1;5;5;0;0]                                                     |   -1.6596|           9|       1.2222|           0.3333|      4.0| 0.9998|
## |327 |N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                     |   -1.4113|           9|       0.5556|           0.4444|      3.0| 0.9994|
## |328 |NGSYPN(HexNAc)LSK[1;5;6;2;0]                                                     |   -1.4113|           9|       0.5556|           0.4444|      3.0| 0.9994|
## |329 |N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                     |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9996|
## |330 |NGSYPN(HexNAc)LSK[1;4;5;0;0]                                                     |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9996|
## |340 |N(HexNAc)GSYPN(HexNAc)LSK[1;9;6;1;0]                                             |   -2.0313|           9|       0.5556|           0.4444|      4.0| 0.9992|
## |345 |N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                     |   -3.0949|           9|       0.5556|           0.4444|      4.0| 0.9998|
## |346 |NGSYPN(HexNAc)LSK[1;6;6;2;0]                                                     |   -3.0949|           9|       0.5556|           0.4444|      4.0| 0.9998|
## |350 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;0;0]                                         |   -0.2389|           9|       0.7778|           0.2222|      5.0| 0.9998|
## |351 |N(HexNAc)GSYPN(HexNAc)LSK[0;9;8;1;0]                                             |    7.0308|           9|       0.5556|           0.4444|      3.0| 0.9676|
## |454 |N(HexNAc)GSYPN(HexNAc)LSK[0;12;4;0;0]                                            |    8.6872|           9|       0.5556|           0.4444|      3.0| 0.9642|
## |471 |N(HexNAc)GSYPN(HexNAc)LSK[0;13;4;0;0]                                            |    6.5838|           9|       0.5556|           0.4444|      5.0| 0.9984|
## |485 |N(HexNAc)GSYPN(Deamidated)LSK[1;7;6;2;0]                                         |    2.5274|           9|       0.0000|           1.0000|      4.0| 0.9998|
## |487 |N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                         |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9992|
## |488 |N(HexNAc)GSYPN(Deamidated)LSK[1;4;5;1;0]                                         |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9992|
## |491 |N(HexNAc)GSYPN(Deamidated)LSK[1;6;5;1;0]                                         |    0.9561|           9|       0.0000|           1.0000|      2.0| 0.9984|
## |492 |HN(HexNAc)VTR[2;4;5;0;0]                                                         |   -5.3622|           5|       0.0000|           1.0000|      3.0| 0.9986|
## |516 |HN(HexNAc)VTR[0;5;2;0;0]                                                         |   -6.1880|           5|       0.0000|           1.0000|      1.8| 0.9964|
## |525 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;2;0]                                         |    2.8345|           9|       0.0000|           1.0000|      6.0| 0.9998|
## |554 |HN(HexNAc)VTR[1;4;5;1;0]                                                         |   -3.1177|           5|       0.0000|           1.0000|      1.8| 0.9992|
## |562 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;0;0]                          |    1.0918|          19|       0.7895|           0.5263|      2.0| 0.9960|
## |627 |N(HexNAc)GSYPN(HexNAc)LSK[0;6;6;0;0]                                             |   -2.6125|           9|       0.0000|           1.0000|      5.0| 0.9998|
## |687 |HN(HexNAc)VTR[2;5;4;0;0]                                                         |   -2.7862|           5|       0.0000|           1.0000|      1.8| 0.9996|
## 
## 
## |    |Glycopeptide_identifier                                                          | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   tPos|
## |:---|:--------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |1   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                |   -5.9097|          18|      22.5556|           0.0000|      8.0| 0.9998|
## |2   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                       |    1.3485|          22|      26.5909|           0.0000|      8.0| 1.0000|
## |3   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;2;0]                       |    1.3485|          22|      26.0455|           0.0000|      8.0| 1.0000|
## |4   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                |   -3.9412|          18|      20.2222|           0.0000|      9.0| 0.9454|
## |5   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                |   -6.5906|          18|      18.0000|           0.0000|      7.0| 0.9988|
## |6   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                |   -8.5683|          18|      21.6667|           0.0000|      8.0| 1.0000|
## |7   |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                |   -6.7784|          18|      18.5000|           0.0000|      6.0| 0.9988|
## |8   |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                       |    2.8744|          22|      25.0455|           0.0000|      9.0| 1.0000|
## |9   |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;2;0]                       |    2.8744|          22|      24.5000|           0.0000|      9.0| 1.0000|
## |10  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;1;0]               |    2.9786|          22|      22.5000|           0.0000|      9.0| 1.0000|
## |11  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;1;0]                          |   -3.9610|          19|      22.8947|           0.0000|      6.0| 0.9442|
## |12  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                       |   -2.7807|          22|      20.6364|           0.0000|      8.0| 1.0000|
## |13  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                       |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |14  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                       |   -4.9353|          22|      17.5909|           0.0000|      7.0| 1.0000|
## |15  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;1;0]                          |   -0.8840|          19|      22.7895|           0.0000|      9.0| 1.0000|
## |16  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;14;4;0;0]              |    9.6761|          22|      29.1364|           0.0000|      7.0| 0.8270|
## |17  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;0;0]                          |   -2.7961|          19|      20.9474|           0.0000|      5.0| 1.0000|
## |18  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;0;0]                       |   -2.7807|          22|      20.0909|           0.0000|      8.0| 1.0000|
## |19  |NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                     |    0.1308|           9|       9.8889|           0.0000|     10.0| 0.9998|
## |20  |N(HexNAc)GSYPNLSK[1;5;4;2;0]                                                     |    0.1308|           9|       9.0000|           0.0000|     10.0| 0.9998|
## |21  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;3;0]                          |   -2.8410|          19|      21.0526|           0.0000|      8.0| 1.0000|
## |22  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                       |   -0.1015|          22|      22.5000|           0.0000|      7.0| 1.0000|
## |23  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;1;0] |   -4.3482|          27|      21.9630|           0.0000|      6.0| 0.9942|
## |24  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;2;0]                 |   -1.5373|          28|      22.3929|           0.0000|      1.0| 0.9888|
## |25  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                       |   -0.1015|          22|      21.9545|           0.0000|      7.0| 1.0000|
## |26  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;1;0]                          |   -3.5609|          19|      18.5789|           0.0000|      8.0| 1.0000|
## |27  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;5;4;2;0] |   -0.1807|          27|      22.2222|           0.0000|      5.0| 0.9970|
## |28  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;2;0]                          |    0.7980|          19|      18.5789|           0.0000|      5.0| 1.0000|
## |29  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]                 |   -0.4458|          28|      22.9643|           0.0000|      2.0| 0.9992|
## |30  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]                 |    0.4326|          28|      22.5000|           0.0000|      1.0| 0.9882|
## |31  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                |   -3.0106|          18|      16.5556|           0.0000|      6.0| 1.0000|
## |32  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                       |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |33  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;0;0]                       |   -3.5683|          22|      19.4091|           0.0000|      6.0| 1.0000|
## |34  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;2;0]                          |   -0.0219|          19|      19.6842|           0.0000|      8.0| 1.0000|
## |35  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]                 |   -2.8434|          28|      20.0000|           0.0000|      2.0| 0.9992|
## |36  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;5;4;1;0]                                            |   -7.4680|          18|      12.3889|           0.0000|      4.0| 0.9998|
## |37  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;2;0]                          |    0.7444|          19|      18.9474|           0.0000|      5.0| 1.0000|
## |38  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;2;0]                          |   -1.2420|          19|      17.6842|           0.0000|      5.0| 1.0000|
## |39  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]           |   -1.5816|          22|      15.9545|           0.0000|      6.0| 1.0000|
## |40  |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]           |   -1.5416|          22|      17.0000|           0.0000|      6.0| 1.0000|
## |41  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;3;0]           |   -1.5816|          22|      15.5000|           0.0000|      6.0| 1.0000|
## |42  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;2;0]                          |    5.9064|          19|      15.6842|           0.0000|      5.0| 0.9996|
## |43  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]                 |   -0.7608|          28|      19.1429|           0.0000|      2.0| 0.9992|
## |44  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                |   -2.7114|          18|      14.1111|           0.0000|      6.0| 1.0000|
## |45  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;4;0]                          |   -3.7451|          19|      17.9474|           0.0000|      6.0| 1.0000|
## |46  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;7;6;2;0]           |   -1.5416|          22|      15.9091|           0.0000|      6.0| 1.0000|
## |47  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;5;3;0]                          |   -2.2429|          19|      17.1053|           0.0000|      8.0| 1.0000|
## |48  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                       |   -4.0561|          22|      14.7727|           0.0000|      6.0| 0.9992|
## |49  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;2;0]                       |   -4.0561|          22|      14.7727|           0.0000|      6.0| 0.9992|
## |50  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                |   -8.2882|          18|      14.7222|           0.0000|      5.0| 1.0000|
## |51  |NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                     |   -1.4933|           9|       8.3333|           0.0000|     10.0| 1.0000|
## |52  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;2;0]                          |   -0.5044|          19|      15.4737|           0.0000|      4.0| 1.0000|
## |53  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                |   -7.0563|          18|      12.6667|           0.0000|      8.0| 1.0000|
## |54  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;1;0]                 |   -4.1984|          28|      18.5714|           0.0000|      2.0| 0.9930|
## |55  |NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                     |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |56  |HN(HexNAc)VTR[1;5;4;1;0]                                                         |   -1.3984|           5|       2.8000|           0.0000|      3.0| 1.0000|
## |57  |N(HexNAc)GSYPNLSK[1;4;5;1;0]                                                     |   -1.4933|           9|       6.8889|           0.0000|     10.0| 0.9996|
## |58  |N(HexNAc)GSYPNLSK[1;6;5;0;0]                                                     |   -2.2131|           9|       7.2222|           0.0000|      9.0| 1.0000|
## |59  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                                |   -2.2699|          18|      14.2222|           0.0000|      2.0| 1.0000|
## |60  |N(HexNAc)GSYPN(HexNAc)LSK[2;6;6;0;0]                                             |   -1.0568|           9|       9.6667|           0.0000|      9.0| 1.0000|
## |61  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                |   -4.3064|          18|      15.0556|           0.0000|      5.0| 1.0000|
## |62  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                       |   -5.7919|          22|      14.2273|           0.0000|      7.0| 0.9864|
## |63  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;1;0]                       |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |64  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                       |    3.6927|          22|      15.3182|           0.0000|      7.0| 1.0000|
## |65  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;1;0]                          |   -2.0419|          19|      13.8947|           0.0000|      3.0| 1.0000|
## |66  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;5;2;0] |   -3.7703|          27|      14.1852|           0.0000|      3.0| 0.9950|
## |67  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;3;0]                          |    2.5239|          19|      11.0000|           0.0000|      7.0| 1.0000|
## |68  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;7;6;0;0]                          |    0.4413|          19|      15.4211|           0.0000|      4.0| 1.0000|
## |69  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;7;6;3;0]                       |   -5.7919|          22|      13.6818|           0.0000|      7.0| 0.9864|
## |70  |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                       |   -9.2195|          22|      13.2727|           0.0000|      3.0| 0.9974|
## |71  |NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                     |   -1.8672|           9|       7.8889|           0.0000|      9.0| 1.0000|
## |72  |N(HexNAc)GSYPNLSK[1;5;4;1;0]                                                     |   -1.8672|           9|       5.7778|           0.0000|      9.0| 0.9994|
## |73  |N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                     |   -1.5766|           9|       6.4444|           0.0000|     10.0| 1.0000|
## |74  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                |   -5.9097|          18|      14.1111|           0.0000|      8.0| 0.9998|
## |75  |N(HexNAc)GSYPN(HexNAc)LSK[1;7;6;1;0]                                             |   -0.9177|           9|       6.3333|           0.0000|      8.0| 0.9900|
## |76  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;6;6;0;0]               |    0.1096|          22|      13.4545|           0.0000|      8.0| 0.9994|
## |77  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                |   -3.3875|          18|      15.5000|           0.0000|      4.0| 1.0000|
## |78  |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;1;0]                 |   -5.4029|          28|      15.2500|           0.0000|      1.0| 0.9944|
## |79  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;0;0]                       |   -9.2195|          22|      12.7273|           0.0000|      3.0| 0.9972|
## |80  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                |   -4.4787|          18|      11.5000|           0.0000|      2.0| 1.0000|
## |81  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;5;1;0]                          |   -4.0059|          19|      11.0000|           0.0000|      4.0| 0.9774|
## |82  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;0;0]                          |   -0.4635|          19|      12.6842|           0.0000|      3.0| 1.0000|
## |83  |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;7;6;0;0]               |   -0.8716|          22|      12.5000|           0.0000|      8.0| 1.0000|
## |84  |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                |   -2.7455|          18|      13.2778|           0.0000|      6.0| 1.0000|
## |85  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;1;0]                                |   -3.9412|          18|      14.5556|           0.0000|      9.0| 0.9454|
## |86  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;2;0]                                |   -6.5906|          18|      12.8333|           0.0000|      7.0| 0.9988|
## |87  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                |   -6.7784|          18|      14.5000|           0.0000|      6.0| 0.9988|
## |88  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                |   -8.5683|          18|      14.2778|           0.0000|      8.0| 1.0000|
## |89  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;1;0]                                            |   -6.4433|          18|       9.7778|           0.0000|      5.0| 0.9730|
## |90  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;1;0]                          |    0.2467|          19|      10.6316|           0.0000|      6.0| 1.0000|
## |91  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;1;0]                          |   -3.4625|          19|      10.2632|           0.0000|      3.0| 1.0000|
## |92  |NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                     |   -4.2508|           9|       6.1111|           0.0000|      7.0| 0.9924|
## |93  |N(HexNAc)GSYPNLSK[2;4;5;0;0]                                                     |   -4.2508|           9|       5.4444|           0.0000|      7.0| 0.9978|
## |94  |NGSYPN(HexNAc)LSK[1;6;5;1;0]                                                     |   -1.5766|           9|       6.2222|           0.0000|     10.0| 0.9902|
## |95  |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;3;0]     |   -0.0952|          28|      13.0000|           0.0000|      1.0| 0.9882|
## |96  |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;0;0]                                            |   -4.0895|          18|      10.4444|           0.0000|      5.0| 0.9998|
## |97  |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;0;0]                          |   -3.0749|          19|      11.8947|           0.0000|      3.0| 1.0000|
## |98  |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;0;0]                                |   -3.0106|          18|      12.0000|           0.0000|      6.0| 1.0000|
## |99  |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;7;6;3;0] |    6.3070|          27|      10.2222|           0.0000|      2.0| 0.9612|
## |100 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                |   -4.7866|          18|      10.0000|           0.0000|      2.0| 1.0000|
## |101 |NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                     |   -1.3818|           9|       6.7778|           0.0000|      6.0| 0.9960|
## |102 |NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                     |   -1.0611|           9|       6.1111|           0.0000|      9.0| 0.9958|
## |103 |N(HexNAc)GSYPNLSK[1;6;5;2;0]                                                     |   -1.0611|           9|       5.2222|           0.0000|      9.0| 1.0000|
## |104 |N(HexNAc)GSYPNLSK[2;5;4;0;0]                                                     |   -1.3818|           9|       5.3333|           0.0000|      6.0| 1.0000|
## |105 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;2;0] |    4.7947|          27|      16.3333|           0.0000|      0.0| 0.9768|
## |106 |N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                     |   -1.0866|           9|       4.1111|           0.0000|      8.0| 0.9996|
## |107 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;4;0;0]                          |   -6.0673|          19|       7.6316|           0.0000|      2.0| 0.9988|
## |108 |N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                         |    9.4307|           9|       5.5556|           0.0000|      1.0| 0.6376|
## |109 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;6;6;1;0] |   -6.6991|          27|      11.8148|           0.0000|      1.0| 0.9844|
## |111 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;1;0]                                |   -4.3064|          18|       9.2778|           0.0000|      5.0| 1.0000|
## |112 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;6;0;0]                                            |   -6.7592|          18|       8.3333|           0.0000|      4.0| 0.9982|
## |113 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;1;0]                                |   -8.2882|          18|       9.6667|           0.0000|      5.0| 1.0000|
## |114 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;7;6;2;0]     |   -1.6746|          28|       8.5000|           0.0000|      1.0| 0.9606|
## |115 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;2;0]                          |    1.0948|          19|       7.0526|           0.0000|      2.0| 0.9998|
## |116 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;1;0]                                |   -7.0563|          18|      10.3333|           0.0000|      8.0| 1.0000|
## |117 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                       |   -1.6205|          22|       6.8182|           0.0000|      6.0| 0.9958|
## |118 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;1;0]                       |   -1.6205|          22|       6.8182|           0.0000|      6.0| 0.9958|
## |119 |NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                     |   -1.8675|           9|       5.1111|           0.0000|      8.0| 1.0000|
## |120 |N(HexNAc)GSYPNLSK[1;5;4;0;0]                                                     |   -1.8675|           9|       5.1111|           0.0000|      8.0| 1.0000|
## |121 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;6;5;0;0]                                            |   -4.6859|          18|       8.1667|           0.0000|      5.0| 0.9992|
## |122 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;0;0]                                             |   -3.5034|           9|       5.3333|           0.0000|      8.0| 0.9978|
## |123 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                         |   -1.3280|           9|       5.5556|           0.0000|      6.0| 0.9994|
## |124 |NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                     |   -0.0550|           9|       4.8889|           0.0000|      6.0| 0.9990|
## |125 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]           |   -9.7324|          22|       7.6364|           0.0000|      3.0| 0.8258|
## |126 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]           |   -9.7324|          22|       7.5455|           0.0000|      3.0| 0.8256|
## |127 |N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                     |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |128 |NGSYPN(HexNAc)LSK[1;5;6;1;0]                                                     |   -2.0741|           9|       3.5556|           0.0000|      6.0| 1.0000|
## |129 |NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                     |   -1.0866|           9|       3.8889|           0.0000|      8.0| 0.9956|
## |130 |N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                     |   -1.5617|           9|       4.5556|           0.0000|      6.0| 1.0000|
## |131 |NGSYPN(HexNAc)LSK[1;7;6;2;0]                                                     |   -1.5617|           9|       4.5556|           0.0000|      6.0| 1.0000|
## |132 |N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                     |   -0.0550|           9|       4.0000|           0.0000|      6.0| 0.9996|
## |133 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                                |   -2.2699|          18|       8.0000|           0.0000|      2.0| 1.0000|
## |134 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;5;5;1;0] |    0.1401|          27|       9.7778|           0.0000|      0.0| 0.9646|
## |135 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;0;0]                                |   -2.7114|          18|       8.4444|           0.0000|      6.0| 1.0000|
## |136 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;6;0;0]                                |   -2.7455|          18|       8.8889|           0.0000|      6.0| 1.0000|
## |137 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;9;0;0]                                             |   -4.4694|           9|       6.5556|           0.0000|      7.0| 0.9922|
## |138 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;3;0]                                            |   -5.3203|          18|       6.9444|           0.0000|      4.0| 0.9992|
## |139 |N(HexNAc)GSYPN(HexNAc)LSK[0;14;4;0;0]                                            |    7.7308|           9|       4.5556|           0.0000|      5.0| 0.9998|
## |140 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[0;6;6;1;0]                          |    1.6958|          19|       5.0526|           0.0000|      0.0| 0.5720|
## |141 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;8;0;0]                                             |   -3.0549|           9|       3.2222|           0.0000|     10.0| 0.9998|
## |142 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[2;7;6;2;0]     |    1.2589|          28|       5.8571|           0.0000|      0.0| 0.6692|
## |144 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                                |   -4.4787|          18|       7.5000|           0.0000|      2.0| 0.9992|
## |145 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]     |    4.3525|          28|       6.5000|           0.0000|      0.0| 0.8340|
## |146 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;3;0]     |   -3.4230|          28|       5.1786|           0.0000|      0.0| 0.8204|
## |147 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]           |   -5.5353|          22|       5.2727|           0.0000|      5.0| 0.9952|
## |149 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;6;1;0]                          |   -4.0430|          19|       6.8947|           0.0000|      2.0| 0.9978|
## |151 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;6;6;3;0]                          |   -0.6209|          19|       5.9474|           0.0000|      3.0| 0.9986|
## |152 |NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                     |   -7.2148|           9|       3.3333|           0.1111|      4.0| 0.9994|
## |153 |SWSYIAETPNSEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;6;5;0;0]                 |   -3.9020|          28|       7.1786|           0.0000|      0.0| 0.8940|
## |154 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;5;5;0;0]                          |   -4.3502|          19|       5.5789|           0.0000|      2.0| 0.9938|
## |155 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;4;0;0]                                |   -4.7866|          18|       7.7222|           0.0000|      2.0| 1.0000|
## |156 |SWSYIAETPN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGYFADYEELR[1;5;4;2;0]     |   -3.9614|          28|       6.3214|           0.0000|      0.0| 0.7024|
## |157 |NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]   |   -3.7097|          23|       8.6957|           0.0000|      0.0| 0.7236|
## |158 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                       |    5.7869|          22|       5.5909|           0.0000|      4.0| 0.9860|
## |159 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;6;1;0]                       |    5.7869|          22|       5.5909|           0.0000|      4.0| 0.9856|
## |160 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;2;0]           |   -5.5353|          22|       4.7273|           0.0000|      5.0| 0.9986|
## |161 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;5;1;0] |   -8.0917|          27|       6.0741|           0.0000|      0.0| 0.8384|
## |162 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;6;0;0]                                |   -3.3875|          18|       9.1111|           0.0000|      4.0| 1.0000|
## |164 |N(Deamidated)LLWLTEKN(HexNAc)GSYPNLSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]   |   -3.7097|          23|       6.9565|           0.0000|      0.0| 0.7304|
## |165 |NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                     |    3.0849|           9|       3.6667|           0.0000|      6.0| 0.9992|
## |166 |N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                     |    1.8418|           9|       2.5556|           0.0000|      8.0| 1.0000|
## |167 |NGSYPN(HexNAc)LSK[1;6;5;3;0]                                                     |    1.8418|           9|       2.5556|           0.0000|      8.0| 1.0000|
## |168 |N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                     |    3.0849|           9|       3.0000|           0.0000|      6.0| 1.0000|
## |169 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[1;6;6;1;0] |   -9.7402|          27|       5.0741|           0.0000|      0.0| 0.8240|
## |170 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;1;0]                                         |   -1.3280|           9|       2.6667|           0.0000|      6.0| 1.0000|
## |171 |HN(HexNAc)VTR[1;6;6;1;0]                                                         |    3.5881|           5|       0.8000|           0.4000|      1.8| 0.9992|
## |172 |N(HexNAc)GSYPN(HexNAc)LSK[2;12;6;0;0]                                            |    2.9929|           9|       2.6667|           0.0000|      7.0| 0.9998|
## |173 |HN(HexNAc)VTR[1;6;5;0;0]                                                         |    0.6254|           5|       1.0000|           0.0000|      2.7| 1.0000|
## |174 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[1;4;6;0;0]                          |   -3.5028|          19|       7.2105|           0.0000|      2.0| 1.0000|
## |177 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[0;7;6;0;0]               |   -6.6968|          22|       4.6818|           0.0000|      2.0| 0.9980|
## |178 |DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]           |   -2.1214|          22|       4.5909|           0.0000|      4.0| 0.9996|
## |179 |N(HexNAc)GSYPNLSK[2;3;6;0;0]                                                     |   -7.2148|           9|       2.4444|           0.2222|      4.0| 0.9984|
## |180 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                       |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9968|
## |181 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;0;0]                       |   -5.5389|          22|       3.9545|           0.0000|      5.0| 0.9968|
## |184 |C(Carbamidomethyl)QTPQGAIN(HexNAc)SSLPFQNIHPVTIGEC(Carbamidomethyl)PK[0;4;6;1;0] |   -0.2344|          27|       4.2222|           0.1111|      0.0| 0.5954|
## |186 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHN(Deamidated)GK[1;7;6;3;0]                    |   -5.4866|          18|       2.5556|           0.0000|      2.0| 0.9994|
## |189 |N(HexNAc)GSYPN(HexNAc)LSK[2;8;6;1;0]                                             |   -2.9545|           9|       2.4444|           0.1111|      8.0| 1.0000|
## |190 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;0;0]                       |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9988|
## |191 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                       |   -2.2514|          22|       4.0455|           0.0000|      4.0| 0.9992|
## |193 |HN(HexNAc)VTR[0;4;2;0;0]                                                         |   -8.8012|           5|       1.6000|           0.0000|      1.8| 0.9944|
## |201 |N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                         |    2.5274|           9|       3.0000|           0.1111|      4.0| 1.0000|
## |202 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;7;0;0]                                             |   -0.6787|           9|       2.5556|           0.0000|      5.0| 1.0000|
## |203 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]           |   -2.1214|          22|       3.5909|           0.0000|      4.0| 0.9996|
## |204 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                       |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9992|
## |205 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;1;0]                       |   -2.7934|          22|       4.2727|           0.0000|      3.0| 0.9994|
## |206 |N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                |   -4.9256|          18|       5.7778|           0.0000|      3.0| 0.9988|
## |208 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)KEK[4;8;11;2;0]                        |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9998|
## |209 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                        |    2.5723|          17|       3.0000|           0.0000|      3.0| 0.9996|
## |212 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                         |    2.8345|           9|       2.1111|           0.1111|      6.0| 1.0000|
## |214 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;9;0;0]                                             |   -0.7798|           9|       2.2222|           0.0000|      5.0| 1.0000|
## |215 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[3;8;8;1;0]               |   -2.2628|          15|       4.7333|           0.0000|      2.0| 0.9978|
## |216 |N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                         |   -0.2389|           9|       2.2222|           0.0000|      5.0| 1.0000|
## |217 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;0;0]                                             |   -1.8379|           9|       3.5556|           0.0000|      8.0| 0.9998|
## |233 |N(HexNAc)VTVTHSVNLLEDSHNGK[1;7;6;1;0]                                            |   -6.3796|          18|       2.8333|           0.0000|      3.0| 0.9922|
## |234 |N(HexNAc)GSYPN(HexNAc)LSK[2;9;6;1;0]                                             |   -1.2417|           9|       2.0000|           0.0000|      9.0| 0.9990|
## |245 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;8;1;0]                                             |    6.9927|           9|       1.7778|           0.3333|      6.0| 0.9986|
## |248 |N(HexNAc)GSYPN(HexNAc)LSK[2;7;7;1;0]                                             |   -0.4717|           9|       1.3333|           0.0000|      5.0| 1.0000|
## |250 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                         |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9976|
## |251 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[0;10;10;0;0]                         |   -4.3522|          15|       3.2000|           0.0000|      3.0| 0.9976|
## |253 |N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;7;0;0]                                |   -4.9256|          18|       5.2222|           0.0000|      3.0| 0.9994|
## |255 |N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                         |    0.9561|           9|       1.4444|           0.1111|      2.0| 0.9986|
## |264 |N(HexNAc)GSYPN(HexNAc)LSK[1;6;6;0;0]                                             |   -1.6219|           9|       2.0000|           0.1111|      5.0| 0.9992|
## |265 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                       |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9930|
## |266 |DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;5;1;0]                       |   -3.8433|          22|       1.3636|           0.2727|      3.0| 0.9930|
## |291 |DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEK[1;8;6;0;0]               |    5.8987|          22|       5.4545|           0.0000|      3.0| 0.9882|
## |305 |HN(HexNAc)VTR[1;5;4;0;0]                                                         |   -2.0406|           5|       0.2000|           0.8000|      2.7| 0.9616|
## |311 |HN(HexNAc)VTR[0;7;2;0;0]                                                         |   -3.4059|           5|       0.6000|           0.4000|      1.8| 0.9990|
## |312 |N(HexNAc)GSYPN(HexNAc)LSK[2;6;5;0;0]                                             |   -2.1511|           9|       0.8889|           0.1111|      5.0| 1.0000|
## |313 |HN(HexNAc)VTR[1;5;5;0;0]                                                         |   -3.9621|           5|       0.6000|           0.4000|      2.1| 0.9558|
## |319 |N(HexNAc)GSYPN(HexNAc)LSK[4;9;11;2;0]                                            |    7.9440|           9|       0.8889|           0.1111|      3.0| 0.9988|
## |320 |N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                     |   -1.6596|           9|       1.2222|           0.3333|      4.0| 1.0000|
## |321 |NGSYPN(HexNAc)LSK[1;5;5;0;0]                                                     |   -1.6596|           9|       1.2222|           0.3333|      4.0| 1.0000|
## |327 |N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                     |   -1.4113|           9|       0.5556|           0.4444|      3.0| 1.0000|
## |328 |NGSYPN(HexNAc)LSK[1;5;6;2;0]                                                     |   -1.4113|           9|       0.5556|           0.4444|      3.0| 1.0000|
## |329 |N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                     |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9992|
## |330 |NGSYPN(HexNAc)LSK[1;4;5;0;0]                                                     |   -2.3321|           9|       1.2222|           0.3333|      2.0| 0.9992|
## |340 |N(HexNAc)GSYPN(HexNAc)LSK[1;9;6;1;0]                                             |   -2.0313|           9|       0.5556|           0.4444|      4.0| 1.0000|
## |345 |N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                     |   -3.0949|           9|       0.5556|           0.4444|      4.0| 1.0000|
## |346 |NGSYPN(HexNAc)LSK[1;6;6;2;0]                                                     |   -3.0949|           9|       0.5556|           0.4444|      4.0| 1.0000|
## |350 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;0;0]                                         |   -0.2389|           9|       0.7778|           0.2222|      5.0| 1.0000|
## |351 |N(HexNAc)GSYPN(HexNAc)LSK[0;9;8;1;0]                                             |    7.0308|           9|       0.5556|           0.4444|      3.0| 0.9910|
## |454 |N(HexNAc)GSYPN(HexNAc)LSK[0;12;4;0;0]                                            |    8.6872|           9|       0.5556|           0.4444|      3.0| 0.9852|
## |471 |N(HexNAc)GSYPN(HexNAc)LSK[0;13;4;0;0]                                            |    6.5838|           9|       0.5556|           0.4444|      5.0| 0.9908|
## |485 |N(HexNAc)GSYPN(Deamidated)LSK[1;7;6;2;0]                                         |    2.5274|           9|       0.0000|           1.0000|      4.0| 0.9992|
## |487 |N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                         |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9992|
## |488 |N(HexNAc)GSYPN(Deamidated)LSK[1;4;5;1;0]                                         |    2.3905|           9|       0.0000|           1.0000|      3.0| 0.9992|
## |491 |N(HexNAc)GSYPN(Deamidated)LSK[1;6;5;1;0]                                         |    0.9561|           9|       0.0000|           1.0000|      2.0| 0.9964|
## |492 |HN(HexNAc)VTR[2;4;5;0;0]                                                         |   -5.3622|           5|       0.0000|           1.0000|      3.0| 0.9992|
## |516 |HN(HexNAc)VTR[0;5;2;0;0]                                                         |   -6.1880|           5|       0.0000|           1.0000|      1.8| 0.9876|
## |525 |N(HexNAc)GSYPN(Deamidated)LSK[1;5;4;2;0]                                         |    2.8345|           9|       0.0000|           1.0000|      6.0| 0.9992|
## |554 |HN(HexNAc)VTR[1;4;5;1;0]                                                         |   -3.1177|           5|       0.0000|           1.0000|      1.8| 0.9970|
## |562 |GFGSGIITSN(HexNAc)ASMDEC(Carbamidomethyl)DTK[2;7;6;0;0]                          |    1.0918|          19|       0.7895|           0.5263|      2.0| 0.9810|
## |627 |N(HexNAc)GSYPN(HexNAc)LSK[0;6;6;0;0]                                             |   -2.6125|           9|       0.0000|           1.0000|      5.0| 0.9992|
## |687 |HN(HexNAc)VTR[2;5;4;0;0]                                                         |   -2.7862|           5|       0.0000|           1.0000|      1.8| 0.9976|
```

```r
trueNegatives<-lapply(newLabels, function(x){
  kable(x[["tNeg"]][,c('Glycopeptide_identifier', "ppm_error", 'peptideLens',"meanCoverage", "percentUncovered", "numStubs", 'tNeg')])
  })
```

```
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: no non-missing arguments to max; returning -Inf
```

```
## 
## 
## |Glycopeptide_identifier | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs| tNeg|
## |:-----------------------|---------:|-----------:|------------:|----------------:|--------:|----:|
## 
## 
## |    |Glycopeptide_identifier                                                        | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   tNeg|
## |:---|:------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |175 |N(Deamidated)LLWLTEKNGSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.957|           0.0000|        0| 0.9552|
## |194 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)NK[1;5;4;2;0] |    -3.987|          23|        5.304|           0.0000|        0| 0.7548|
## |195 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVNN(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.304|           0.0000|        0| 0.7548|
## |196 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0] |    -3.987|          23|        5.783|           0.0000|        0| 0.8544|
## |197 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.783|           0.0000|        0| 0.8544|
## |198 |NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.217|           0.0000|        0| 0.9696|
## |199 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVNN(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.696|           0.0000|        0| 0.9504|
## |210 |NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;5;4;2;0] |    -3.987|          23|        5.000|           0.0000|        0| 0.6472|
## |218 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)NK[1;7;6;4;0] |    -3.710|          23|        5.783|           0.0000|        0| 0.8532|
## |249 |N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        5.087|           0.0000|        0| 0.7796|
## |252 |N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                 |    -7.327|          23|        5.913|           0.0000|        1| 0.6518|
## |295 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.7496|
## |296 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.7496|
## |314 |N(HexNAc)GSYPN(HexNAc)LSK[3;14;13;0;0]                                         |     7.494|           9|        1.222|           0.5556|        2| 0.9880|
## 
## 
## |    |Glycopeptide_identifier                                                        | ppm_error| peptideLens| meanCoverage| percentUncovered| numStubs|   tNeg|
## |:---|:------------------------------------------------------------------------------|---------:|-----------:|------------:|----------------:|--------:|------:|
## |175 |N(Deamidated)LLWLTEKNGSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |    -3.710|          23|        6.957|           0.0000|        0| 0.6490|
## |295 |N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.8664|
## |296 |N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)K[2;11;7;0;0]                        |     9.544|          15|        2.133|           0.2000|        2| 0.8664|
## |314 |N(HexNAc)GSYPN(HexNAc)LSK[3;14;13;0;0]                                         |     7.494|           9|        1.222|           0.5556|        2| 0.9888|
```

