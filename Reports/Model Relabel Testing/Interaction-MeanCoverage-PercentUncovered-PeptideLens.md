---
title: "Relabel Testing"
author: "Joshua Klein"
date: "Tuesday, June 17, 2014"
output:
  html_document:
    fig_height: 8
    fig_width: 15
    number_sections: yes
    theme: journal
    toc: yes
---
<style>body{max-width:100% !important;}</style>
# Model Fit Test



```r
testFormula <- call ~ (meanCoverage * percentUncovered  *  peptideLens) + abs_ppm_error + X.b_ion_with_HexNAc_coverage + 
  X.y_ion_with_HexNAc_coverage + numStubs

mtry <- 5

ntree <- 5000

RUN_ID <- "Interaction-MeanCoverage-PercentUncovered-PeptideLens"
```


```
## [1] "testFormula ~"                                                                                                                                       
## [2] "testFormula call"                                                                                                                                    
## [3] "testFormula (meanCoverage * percentUncovered * peptideLens) + abs_ppm_error + X.b_ion_with_HexNAc_coverage + X.y_ion_with_HexNAc_coverage + numStubs"
```

```
## [1] "mtry 5"
```

```
## [1] "ntree 5000"
```

## Load Data


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

<img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmdfit-forests1.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmdfit-forests2.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmdfit-forests3.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" /><img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmdfit-forests4.png" title="plot of chunk fit-forests" alt="plot of chunk fit-forests" style="display: block; margin: auto;" />

```r
names(labelForests) <- names(labelData)
```




### Error Table of Labeled Components

| TruePositive| FalsePositive| FalseNegative| TrueNegative|    FPR|   FDR| ErrorRate| .numYes| .numNo|.rownames                             |    TPR|
|------------:|-------------:|-------------:|------------:|------:|-----:|---------:|-------:|------:|:-------------------------------------|------:|
|          215|            28|            28|          574| 0.0465| 2e-04|    0.0663|     243|    602|USSRLabel.SouthCarolinaLabel          | 0.8848|
|          107|             0|            13|          581| 0.0000| 0e+00|    0.0185|     120|    581|SolomonLabel.SouthCarolinaLabel       | 0.8917|
|          387|            28|            41|         2188| 0.0126| 0e+00|    0.0261|     428|   2216|Merged.SouthCarolinaLabel             | 0.9042|
|           60|            29|             5|         1004| 0.0281| 5e-04|    0.0310|      65|   1033|SouthCarolinaLabel.USSRLabel          | 0.9231|
|          114|             4|             6|          577| 0.0069| 1e-04|    0.0143|     120|    581|SolomonLabel.USSRLabel                | 0.9500|
|          231|            45|            12|          557| 0.0748| 3e-04|    0.0675|     243|    602|USSRLabel.SolomonLabel                | 0.9506|
|          414|            66|            14|         2150| 0.0298| 1e-04|    0.0303|     428|   2216|Merged.SolomonLabel                   | 0.9673|
|           63|            21|             2|         1012| 0.0203| 3e-04|    0.0209|      65|   1033|SouthCarolinaLabel.SolomonLabel       | 0.9692|
|          417|            33|            11|         2183| 0.0149| 0e+00|    0.0166|     428|   2216|Merged.USSRLabel                      | 0.9743|
|          243|             0|             0|          602| 0.0000| 0e+00|    0.0000|     243|    602|USSRLabel.USSRLabel                   | 1.0000|
|          243|             0|             0|          602| 0.0000| 0e+00|    0.0000|     243|    602|USSRLabel.Merged                      | 1.0000|
|          120|             0|             0|          581| 0.0000| 0e+00|    0.0000|     120|    581|SolomonLabel.SolomonLabel             | 1.0000|
|          120|             0|             0|          581| 0.0000| 0e+00|    0.0000|     120|    581|SolomonLabel.Merged                   | 1.0000|
|           65|             0|             0|         1033| 0.0000| 0e+00|    0.0000|      65|   1033|SouthCarolinaLabel.SouthCarolinaLabel | 1.0000|
|           65|             0|             0|         1033| 0.0000| 0e+00|    0.0000|      65|   1033|SouthCarolinaLabel.Merged             | 1.0000|
|          428|             0|             0|         2216| 0.0000| 0e+00|    0.0000|     428|   2216|Merged.Merged                         | 1.0000|

## Summary Plots

```r
ggplot(errorTable, aes(size = TPR,  y = FPR, x = FalseNegative/(FalseNegative + TruePositive), color = .rownames, shape = sapply(sapply(.rownames, function(rn){strsplit(rn, ".", T)}), `[`,1), label = zapsmall(TruePositive/.numYes, 3))) + geom_point() + scale_size(range = c(4,10)) + guides(col = guide_legend(ncol = 3, byrow = TRUE))
```

<img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmderror-plot1.png" title="plot of chunk error-plot" alt="plot of chunk error-plot" style="display: block; margin: auto;" />

```r
ggplot(errorTable, aes(y = TPR, x = FPR, size = ErrorRate, color = .rownames)) + geom_point(position="jitter") + scale_size(range = c(4,10)) + guides(col = guide_legend(ncol = 3, byrow = TRUE))
```

<img src="figure/Interaction-MeanCoverage-PercentUncovered-PeptideLens-Rmderror-plot2.png" title="plot of chunk error-plot" alt="plot of chunk error-plot" style="display: block; margin: auto;" />

## Positive Entities
###  USSRLabel.USSRLabel 
Entities Positive:  69 


| MS1_Score| Calc_mass| count|bestSeq                                                                        | worstScore| bestScore|
|---------:|---------:|-----:|:------------------------------------------------------------------------------|----------:|---------:|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                   |     0.9970|    0.9988|
|    0.3989|      4061|     2|N(HexNAc)GSYPN(Deamidated)LSK[1;7;6;2;0]                                       |     0.9992|    1.0000|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                   |     1.0000|    1.0000|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                   |     0.9922|    0.9922|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                   |     1.0000|    1.0000|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                       |     0.9994|    0.9996|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                       |     1.0000|    1.0000|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                              |     0.9998|    1.0000|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                       |     0.9852|    0.9876|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                   |     1.0000|    1.0000|
|    0.5380|      3259|     2|N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                       |     0.0086|    0.6712|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                   |     1.0000|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                   |     1.0000|    1.0000|
|    0.5614|      3183|     2|N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                   |     0.9996|    1.0000|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                   |     1.0000|    1.0000|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                       |     1.0000|    1.0000|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                       |     1.0000|    1.0000|
|    0.5762|      4351|     2|N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                   |     0.9998|    1.0000|
|    0.5799|      2876|     2|NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                   |     0.9992|    1.0000|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                     |     1.0000|    1.0000|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                              |     0.9932|    0.9980|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                   |     0.9998|    0.9998|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                   |     0.9998|    0.9998|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                   |     0.9998|    0.9998|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                       |     0.6766|    0.8540|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                     |     0.9878|    0.9882|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                   |     1.0000|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                     |     0.9848|    0.9848|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                              |     0.9992|    0.9992|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]         |     0.9956|    0.9980|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                     |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                   |     1.0000|    1.0000|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                     |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                              |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]         |     0.9998|    0.9998|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]         |     0.9584|    0.9588|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                     |     0.9998|    1.0000|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                     |     1.0000|    1.0000|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                              |     0.9994|    0.9998|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                       |     0.0066|    0.8038|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                   |     1.0000|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                              |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                     |     1.0000|    1.0000|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                     |     1.0000|    1.0000|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                              |     0.9996|    0.9996|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                              |     1.0000|    1.0000|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]         |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                              |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                     |     1.0000|    1.0000|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                     |     1.0000|    1.0000|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                   |     0.9998|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                   |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                              |     1.0000|    1.0000|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                   |     1.0000|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                     |     1.0000|    1.0000|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                     |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                              |     1.0000|    1.0000|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                              |     0.9996|    0.9996|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                     |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                              |     1.0000|    1.0000|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                     |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]         |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                              |     1.0000|    1.0000|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                     |     1.0000|    1.0000|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                       |     0.9974|    0.9974|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                              |     0.9996|    0.9996|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                      |     0.9990|    0.9996|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                              |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |     0.0100|    0.7384|
###  USSRLabel.SolomonLabel 
Entities Positive:  76 


| MS1_Score| Calc_mass| count|bestSeq                                                                        | worstScore| bestScore|
|---------:|---------:|-----:|:------------------------------------------------------------------------------|----------:|---------:|
|    0.3777|      2975|     2|N(HexNAc)GSYPNLSK[2;3;6;0;0]                                                   |     0.9810|    0.9956|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                       |     0.9366|    0.9918|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                   |     0.9532|    0.9532|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                   |     0.5096|    0.5096|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                   |     0.9468|    0.9468|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                       |     0.9520|    0.9944|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                       |     0.9302|    0.9302|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                              |     0.9998|    0.9998|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                       |     0.1596|    0.7266|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                   |     0.9856|    0.9856|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                   |     0.9990|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                   |     1.0000|    1.0000|
|    0.5614|      3183|     2|NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                   |     0.9860|    0.9978|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                   |     0.9850|    0.9850|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                       |     0.9892|    1.0000|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                       |     0.9366|    0.9958|
|    0.5762|      4351|     2|N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                   |     0.9110|    0.9964|
|    0.5799|      2876|     2|N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                   |     0.9074|    0.9948|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;0;0]                     |     0.9880|    0.9886|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                              |     0.9596|    0.9596|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                   |     1.0000|    1.0000|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                   |     0.9996|    0.9996|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                   |     0.9850|    0.9850|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                       |     0.1312|    0.8752|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                     |     0.9748|    0.9750|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                   |     0.9990|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                     |     0.9816|    0.9816|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                              |     0.9512|    0.9512|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]         |     0.9358|    0.9870|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                     |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                   |     0.9932|    0.9932|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                     |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                              |     0.9964|    0.9964|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]         |     0.9754|    0.9994|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]         |     0.9746|    0.9746|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                     |     0.9928|    0.9934|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                     |     0.9406|    0.9416|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                              |     0.9506|    0.9656|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                       |     0.0102|    0.5952|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                   |     0.9982|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                              |     0.9996|    0.9996|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                     |     1.0000|    1.0000|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                     |     1.0000|    1.0000|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                              |     0.9746|    0.9748|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                              |     1.0000|    1.0000|
|    0.6264|      5400|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;6;5;2;0] |     0.0584|    0.7408|
|    0.6264|      5400|     3|NLLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;8;7;0;0]                 |     0.6598|    0.7288|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]         |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                              |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                     |     1.0000|    1.0000|
|    0.6299|      5181|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                        |     0.6770|    0.6770|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                     |     1.0000|    1.0000|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                   |     1.0000|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                   |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;9;9;0;0]                         |     0.5522|    0.5522|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                              |     0.9952|    0.9952|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                   |     0.9990|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                     |     0.9982|    0.9982|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                     |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                              |     0.9956|    0.9956|
|    0.6452|      5343|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                 |     0.2804|    0.9304|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                              |     0.9668|    0.9668|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                     |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                              |     0.9990|    1.0000|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                     |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]         |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                              |     0.9982|    0.9982|
|    0.6635|      5035|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0] |     0.0494|    0.9642|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                     |     1.0000|    1.0000|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                       |     0.9622|    0.9622|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                              |     0.9766|    0.9766|
|    0.7180|      6365|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[4;8;6;2;0]     |     0.5288|    0.7836|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                      |     0.9704|    0.9894|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                              |     0.9986|    0.9986|
|    0.7364|      6348|     8|NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |     0.6570|    0.9706|
|    0.7409|      5782|     3|NLLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[4;8;6;0;0]     |     0.6822|    0.8120|
###  USSRLabel.SouthCarolinaLabel 
Entities Positive:  68 


| MS1_Score| Calc_mass| count|bestSeq                                                                        | worstScore| bestScore|
|---------:|---------:|-----:|:------------------------------------------------------------------------------|----------:|---------:|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                   |     0.6066|    0.6390|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                       |     0.4448|    0.7708|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                   |     0.5876|    0.5876|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                   |     0.5876|    0.5876|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                       |     0.7330|    0.7980|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                              |     0.9670|    0.9908|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                   |     0.7452|    0.7452|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                   |     0.9150|    0.9664|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                   |     0.9966|    0.9970|
|    0.5614|      3183|     2|N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                   |     0.8364|    0.8490|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                   |     0.9518|    0.9518|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                       |     0.7978|    0.8992|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                       |     0.4448|    0.7538|
|    0.5762|      4351|     2|NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                   |     0.8234|    0.8254|
|    0.5799|      2876|     2|N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                   |     0.8572|    0.9830|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                     |     0.7630|    0.7638|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                              |     0.9420|    0.9476|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                   |     0.9954|    0.9954|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                   |     0.9856|    0.9856|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                   |     0.8256|    0.8256|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                       |     0.2182|    0.6602|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                     |     0.9120|    0.9120|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                   |     0.9168|    0.9706|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                     |     0.6506|    0.6506|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                              |     0.9840|    0.9842|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]         |     0.8154|    0.9460|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                     |     0.9628|    0.9628|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                   |     0.7980|    0.7980|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                     |     0.9958|    0.9958|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                              |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]         |     0.9048|    0.9186|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]         |     0.8920|    0.9100|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;1;0]                     |     0.8744|    0.8834|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;0;0]                     |     0.6208|    0.6212|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                              |     0.9626|    0.9820|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                   |     0.9544|    0.9792|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                              |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                     |     0.9966|    0.9966|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                     |     0.9966|    0.9966|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                              |     0.9570|    0.9570|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                              |     1.0000|    1.0000|
|    0.6264|      5400|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;6;5;2;0] |     0.0048|    0.7484|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]         |     0.9944|    0.9944|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                              |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                     |     0.9946|    0.9946|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                     |     0.9962|    0.9962|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                   |     0.9798|    0.9846|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                   |     0.9946|    0.9954|
|    0.6383|      4971|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;9;9;0;0]                         |     0.5536|    0.5536|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                              |     0.9960|    0.9960|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                   |     0.9544|    0.9764|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                     |     0.9530|    0.9530|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                     |     0.9964|    0.9964|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                              |     0.9992|    0.9992|
|    0.6452|      5343|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                 |     0.0008|    0.9290|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                              |     0.9960|    0.9960|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                     |     0.9966|    0.9966|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                              |     0.9768|    0.9994|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                     |     0.9858|    0.9858|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]         |     0.9944|    0.9944|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                              |     0.9996|    0.9996|
|    0.6635|      5035|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0] |     0.0084|    0.9382|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                     |     0.9946|    0.9946|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                       |     0.9538|    0.9538|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                              |     0.9060|    0.9060|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)KEK[4;8;11;2;0]                      |     0.9532|    0.9600|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                              |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |     0.7120|    0.9720|
###  USSRLabel.Merged 
Entities Positive:  69 


| MS1_Score| Calc_mass| count|bestSeq                                                                        | worstScore| bestScore|
|---------:|---------:|-----:|:------------------------------------------------------------------------------|----------:|---------:|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                   |     0.9988|    0.9990|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                       |     0.9968|    1.0000|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                   |     1.0000|    1.0000|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                   |     0.9704|    0.9704|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                   |     1.0000|    1.0000|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                       |     1.0000|    1.0000|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                       |     0.9898|    0.9898|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                              |     1.0000|    1.0000|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                       |     0.9630|    0.9684|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                   |     1.0000|    1.0000|
|    0.5380|      3259|     2|N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                       |     0.0134|    0.6540|
|    0.5384|      2893|     2|N(HexNAc)GSYPNLSK[2;5;4;0;0]                                                   |     0.9996|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                   |     1.0000|    1.0000|
|    0.5614|      3183|     2|NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                   |     0.9996|    1.0000|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                   |     0.9988|    0.9988|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                       |     0.9996|    1.0000|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                       |     0.9982|    1.0000|
|    0.5762|      4351|     2|NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                   |     1.0000|    1.0000|
|    0.5799|      2876|     2|NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                   |     0.9992|    0.9996|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                     |     0.9996|    0.9996|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                              |     0.9934|    0.9940|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                   |     0.9994|    1.0000|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                   |     1.0000|    1.0000|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                   |     1.0000|    1.0000|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                       |     0.6578|    0.8326|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                     |     0.9998|    0.9998|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                   |     1.0000|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                     |     0.9988|    0.9988|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                              |     0.9990|    0.9990|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]         |     0.9958|    0.9992|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                     |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                   |     1.0000|    1.0000|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                     |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                              |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]         |     0.9986|    0.9992|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]         |     0.9970|    0.9974|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                     |     0.9996|    0.9998|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                     |     1.0000|    1.0000|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                              |     0.9954|    0.9978|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                       |     0.0190|    0.7334|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                   |     1.0000|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                              |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                     |     1.0000|    1.0000|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                     |     1.0000|    1.0000|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                              |     0.9996|    0.9996|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                              |     1.0000|    1.0000|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]         |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                              |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                     |     1.0000|    1.0000|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                     |     1.0000|    1.0000|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                   |     1.0000|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                   |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                              |     1.0000|    1.0000|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                   |     1.0000|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;6;1;0]                     |     0.9996|    0.9998|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                     |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                              |     1.0000|    1.0000|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                              |     0.9996|    0.9996|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                     |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                              |     1.0000|    1.0000|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                     |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]         |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                              |     1.0000|    1.0000|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                     |     1.0000|    1.0000|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                       |     0.9992|    0.9992|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                              |     0.9994|    0.9994|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                      |     0.9996|    1.0000|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                              |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0] |     0.0236|    0.8252|
###  SolomonLabel.USSRLabel 
Entities Positive:  14 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                       | worstScore| bestScore|
|---------:|---------:|-----:|:-------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                              |     0.5498|    0.6322|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                             |     0.9884|    0.9934|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                              |     0.9984|    0.9996|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0] |     0.5306|    0.6246|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                              |     0.6440|    0.6806|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                              |     0.9996|    0.9996|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                    |     0.9916|    0.9916|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;1;0]                                                    |     0.9908|    0.9974|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                                                    |     0.9914|    0.9918|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                             |     0.6578|    0.6632|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                             |     0.9858|    0.9864|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                              |     0.9866|    0.9882|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                             |     0.9990|    0.9996|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                             |     1.0000|    1.0000|
###  SolomonLabel.SolomonLabel 
Entities Positive:  14 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                       | worstScore| bestScore|
|---------:|---------:|-----:|:-------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                              |     0.0738|    0.6832|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                                             |     0.9950|    0.9974|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                              |     0.9952|    1.0000|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEKN(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0] |     0.9882|    0.9998|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                              |     0.1676|    0.8314|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                              |     1.0000|    1.0000|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                    |     0.9978|    0.9978|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                    |     0.9964|    0.9964|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                                                    |     0.9948|    0.9952|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                             |     0.1594|    0.7952|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                             |     0.9894|    0.9948|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                              |     0.9966|    0.9996|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                                             |     0.9550|    0.9978|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                             |     0.9998|    0.9998|
###  SolomonLabel.SouthCarolinaLabel 
Entities Positive:  11 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                       | worstScore| bestScore|
|---------:|---------:|-----:|:-------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                             |     0.7500|    0.7586|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                              |     0.7018|    0.8892|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0] |     0.7488|    0.8448|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                              |     0.9968|    0.9968|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                    |     0.5914|    0.5914|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                    |     0.8656|    0.8658|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                    |     0.9196|    0.9196|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                             |     0.9280|    0.9280|
|    0.7328|      4850|     2|N(Deamidated)LLWLTGKNGLYPN(HexNAc)LSK[1;6;6;2;0]                                                              |     0.6566|    0.6570|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                             |     0.8328|    0.8640|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                             |     1.0000|    1.0000|
###  SolomonLabel.Merged 
Entities Positive:  14 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                       | worstScore| bestScore|
|---------:|---------:|-----:|:-------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                              |     0.1512|    0.6768|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                                             |     0.9992|    0.9998|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                              |     0.9998|    1.0000|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0] |     0.9894|    1.0000|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                              |     0.2248|    0.8444|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                              |     1.0000|    1.0000|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                    |     0.9994|    0.9994|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                    |     0.9992|    0.9992|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                    |     0.9982|    0.9994|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                             |     0.1280|    0.7372|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                             |     0.9998|    0.9998|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                              |     0.9992|    0.9998|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                                             |     0.9998|    1.0000|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                             |     1.0000|    1.0000|
###  SouthCarolinaLabel.USSRLabel 
Entities Positive:  33 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9860|    0.9870|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9998|    0.9998|
|    0.2810|      8521|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;5;4;0;0] |     0.1338|    0.5198|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     0.9990|    0.9996|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9356|    0.9358|
|    0.3273|      8504|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;4;4;0;0]             |     0.8636|    0.8936|
|    0.3359|      5191|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;2;0]                                                                                                         |     0.6514|    0.6524|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.5022|    0.9158|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.4804|    0.8884|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;5;5;1;0] |     0.8532|    0.8690|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;5;5;0;0]             |     0.8102|    0.8138|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9326|    0.9518|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.3480|    0.5364|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.3062|    0.6914|
|    0.4645|      7705|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;6;0;0]            |     0.4728|    0.5534|
|    0.4687|      8319|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;7;6;1;0]                        |     0.2158|    0.5228|
|    0.4795|      7283|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;4;4;0;0]            |     0.1308|    0.6072|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8046|    0.9016|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.8780|    0.9478|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0110|    0.5180|
|    0.5856|      8029|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;7;6;0;0]            |     0.4882|    0.5730|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.5272|    0.8582|
|    0.6159|      7298|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]                        |     0.0100|    0.5858|
|    0.6314|      7663|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]                        |     0.0450|    0.5438|
|    0.6490|      7443|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;4;1;0]                        |     0.1630|    0.8192|
###  SouthCarolinaLabel.SolomonLabel 
Entities Positive:  29 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     0.9964|    0.9964|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;8;7;0;0]                                                                                                         |     0.9578|    0.9596|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9996|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     0.9982|    0.9982|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.8656|    0.8656|
|    0.3273|      8504|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;4;4;0;0]             |     0.8728|    0.8738|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.4486|    0.8694|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.6060|    0.8960|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;5;5;1;0] |     0.8414|    0.8500|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;5;5;0;0]             |     0.8432|    0.8458|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.8414|    0.8772|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.4350|    0.6282|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.4570|    0.8148|
|    0.4795|      7283|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;4;4;0;0]            |     0.0190|    0.6202|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8604|    0.9212|
|    0.5163|      8357|     3|TIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKNVTVTHSVNLLEDSHNGKLC(Carbamidomethyl)K[2;7;11;0;0]                                                  |     0.5618|    0.6598|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.8496|    0.8912|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.4418|    0.6268|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0244|    0.5518|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.5406|    0.8456|
|    0.6314|      7663|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]                        |     0.0206|    0.5752|
|    0.6490|      7443|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;4;1;0]                        |     0.0650|    0.7962|
###  SouthCarolinaLabel.SouthCarolinaLabel 
Entities Positive:  22 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     0.9996|    0.9998|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     0.9984|    0.9984|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9104|    0.9400|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     0.9966|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     0.9994|    0.9994|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9960|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     0.9970|    1.0000|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9694|    0.9734|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.0070|    0.9872|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.2282|    0.8602|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9570|    0.9800|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.6702|    0.8550|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.0338|    0.8432|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8296|    0.9570|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.9836|    0.9984|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.7630|    0.9064|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0000|    0.7432|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.1990|    0.9612|
###  SouthCarolinaLabel.Merged 
Entities Positive:  22 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9986|    0.9994|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9830|    0.9854|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.0744|    0.9914|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.2686|    0.9144|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9604|    0.9850|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.7316|    0.8786|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.0650|    0.8436|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8248|    0.9642|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.9852|    0.9998|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.7646|    0.8864|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0000|    0.8016|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.0880|    0.9616|
###  Merged.USSRLabel 
Entities Positive:  116 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9860|    0.9870|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9998|    0.9998|
|    0.2810|      8521|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;5;4;0;0] |     0.1338|    0.5198|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     0.9990|    0.9996|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9356|    0.9358|
|    0.3273|      8504|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;4;4;0;0]             |     0.8636|    0.8936|
|    0.3359|      5191|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;2;0]                                                                                                         |     0.6514|    0.6524|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                                                                                              |     0.9970|    0.9988|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.5022|    0.9158|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.4804|    0.8884|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;5;5;1;0] |     0.8532|    0.8690|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;5;5;0;0]             |     0.8102|    0.8138|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9326|    0.9518|
|    0.3989|      4061|     2|N(HexNAc)GSYPN(Deamidated)LSK[1;7;6;2;0]                                                                                                                  |     0.9992|    1.0000|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.3480|    0.5364|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.3062|    0.6914|
|    0.4645|      7705|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;6;0;0]            |     0.4728|    0.5534|
|    0.4687|      8319|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;7;6;1;0]                        |     0.2158|    0.5228|
|    0.4795|      7283|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;4;4;0;0]            |     0.1308|    0.6072|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                                                                                              |     0.9922|    0.9922|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8046|    0.9016|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                  |     0.9994|    0.9996|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                  |     1.0000|    1.0000|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.8780|    0.9478|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                                                                                         |     0.9998|    1.0000|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                                                                                                  |     0.9852|    0.9876|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5380|      3259|     2|N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                                                                                  |     0.0086|    0.6712|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                                                                          |     0.5498|    0.6322|
|    0.5614|      3183|     2|N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                                                                                              |     0.9996|    1.0000|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                  |     1.0000|    1.0000|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0110|    0.5180|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                  |     1.0000|    1.0000|
|    0.5762|      4351|     2|N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                                                                                              |     0.9998|    1.0000|
|    0.5799|      2876|     2|NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                                                                                              |     0.9992|    1.0000|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                                                                                                |     1.0000|    1.0000|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                                                                                                         |     0.9932|    0.9980|
|    0.5856|      8029|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;7;6;0;0]            |     0.4882|    0.5730|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                                                                                              |     0.9998|    0.9998|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.5272|    0.8582|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                              |     0.9998|    0.9998|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                                                                                              |     0.9998|    0.9998|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                                                                                                  |     0.6766|    0.8540|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                                                                                                |     0.9878|    0.9882|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                                                                                                |     0.9848|    0.9848|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     0.9992|    0.9992|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]                                                                                    |     0.9956|    0.9980|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                    |     0.9998|    0.9998|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     0.9884|    0.9934|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]                                                                                    |     0.9584|    0.9588|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                                                                                                |     0.9998|    1.0000|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                                                                                                         |     0.9994|    0.9998|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                                                                                                  |     0.0066|    0.8038|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6159|      7298|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]                        |     0.0100|    0.5858|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                                                                          |     0.9984|    0.9996|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                    |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6314|      7663|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]                        |     0.0450|    0.5438|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                              |     0.9998|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     1.0000|    1.0000|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                                                                                         |     1.0000|    1.0000|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6490|      7443|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;4;1;0]                        |     0.1630|    0.8192|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                                |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0]                                             |     0.5306|    0.6246|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                    |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                                                                          |     0.6440|    0.6806|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                          |     0.9996|    0.9996|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                                                                |     0.9916|    0.9916|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;6;1;0]                                                                                                |     0.9908|    0.9974|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                                                                                                |     0.9914|    0.9918|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                                                                                                  |     0.9974|    0.9974|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.9996|    0.9996|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                                                                         |     0.6578|    0.6632|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                                                                         |     0.9858|    0.9864|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                                                                                                 |     0.9990|    0.9996|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                                                                          |     0.9866|    0.9882|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]                                                                            |     0.0100|    0.7384|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.9990|    0.9996|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
###  Merged.SolomonLabel 
Entities Positive:  119 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     0.9964|    0.9964|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;8;7;0;0]                                                                                                         |     0.9578|    0.9596|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9996|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     0.9982|    0.9982|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.8656|    0.8656|
|    0.3273|      8504|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;4;4;0;0]             |     0.8728|    0.8738|
|    0.3777|      2975|     2|N(HexNAc)GSYPNLSK[2;3;6;0;0]                                                                                                                              |     0.9810|    0.9956|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.4486|    0.8694|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.6060|    0.8960|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;5;5;1;0] |     0.8414|    0.8500|
|    0.3963|      8869|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[2;5;5;0;0]             |     0.8432|    0.8458|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.8414|    0.8772|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                                                                                                  |     0.9366|    0.9918|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                                                                                              |     0.9532|    0.9532|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.4350|    0.6282|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.4570|    0.8148|
|    0.4795|      7283|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;4;4;0;0]            |     0.0190|    0.6202|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                                                                                              |     0.5096|    0.5096|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8604|    0.9212|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                                                                                              |     0.9468|    0.9468|
|    0.5163|      8357|     3|TIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKNVTVTHSVNLLEDSHNGKLC(Carbamidomethyl)K[2;7;11;0;0]                                                  |     0.5618|    0.6598|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                  |     0.9520|    0.9944|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                  |     0.9302|    0.9302|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.8496|    0.8912|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                                                                                         |     0.9998|    0.9998|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.4418|    0.6268|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                                                                                                  |     0.1596|    0.7266|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                                                                                              |     0.9856|    0.9856|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                                                                                              |     0.9990|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                                                                          |     0.0738|    0.6832|
|    0.5614|      3183|     2|NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                                                                                              |     0.9860|    0.9978|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                                                                                              |     0.9850|    0.9850|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                  |     0.9892|    1.0000|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0244|    0.5518|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                  |     0.9366|    0.9958|
|    0.5762|      4351|     2|N(HexNAc)GSYPNLSK[1;7;6;3;0]                                                                                                                              |     0.9110|    0.9964|
|    0.5799|      2876|     2|N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                                                                                              |     0.9074|    0.9948|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;0;0]                                                                                                |     0.9880|    0.9886|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                                                                                                         |     0.9596|    0.9596|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.5406|    0.8456|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                              |     0.9996|    0.9996|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                                                                                              |     0.9850|    0.9850|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                                                                                                  |     0.1312|    0.8752|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                                                                                                |     0.9748|    0.9750|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                                              |     0.9990|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                                                                                                |     0.9816|    0.9816|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     0.9512|    0.9512|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                    |     0.9358|    0.9870|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                                                                                              |     0.9932|    0.9932|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     0.9964|    0.9964|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                    |     0.9754|    0.9994|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                                                                                         |     0.9950|    0.9974|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]                                                                                    |     0.9746|    0.9746|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                                                                                                |     0.9928|    0.9934|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                                                                                                |     0.9406|    0.9416|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;5;5;0;0]                                                                                                         |     0.9506|    0.9656|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                                                                                                  |     0.0102|    0.5952|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                                                                                              |     0.9982|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                                                                          |     0.9952|    1.0000|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                                                                                         |     0.9746|    0.9748|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6264|      5400|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;6;5;2;0]                                                                            |     0.0584|    0.7408|
|    0.6264|      5400|     3|NLLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;8;7;0;0]                                                                                            |     0.6598|    0.7288|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                    |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6299|      5181|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[2;11;7;0;0]                                                                                                   |     0.6770|    0.6770|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6314|      7663|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]                        |     0.0206|    0.5752|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;9;9;0;0]                                                                                                    |     0.5522|    0.5522|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     0.9952|    0.9952|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                              |     0.9990|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                                                                                                |     0.9982|    0.9982|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                                                                                         |     0.9956|    0.9956|
|    0.6452|      5343|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                                                                                            |     0.2804|    0.9304|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                                                                                         |     0.9668|    0.9668|
|    0.6490|      7443|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;4;1;0]                        |     0.0650|    0.7962|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                                |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     0.9990|    1.0000|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEKN(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0]                                             |     0.9882|    0.9998|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                    |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     0.9982|    0.9982|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                                                                          |     0.1676|    0.8314|
|    0.6635|      5035|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0]                                                                            |     0.0494|    0.9642|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                          |     1.0000|    1.0000|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                                                                |     0.9978|    0.9978|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     0.9964|    0.9964|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;4;1;0]                                                                                                |     0.9948|    0.9952|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                                                                                                  |     0.9622|    0.9622|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.9766|    0.9766|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                                                                         |     0.1594|    0.7952|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                                                                         |     0.9894|    0.9948|
|    0.7180|      6365|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[4;8;6;2;0]                                                                                |     0.5288|    0.7836|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                                                                                                 |     0.9704|    0.9894|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                                                                          |     0.9966|    0.9996|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     0.9986|    0.9986|
|    0.7364|      6348|     8|NLLWLTEKN(HexNAc)GSYPN(Deamidated)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]                                                                            |     0.6570|    0.9706|
|    0.7409|      5782|     3|NLLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[4;8;6;0;0]                                                                                |     0.6822|    0.8120|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                                                                                         |     0.9550|    0.9978|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     0.9998|    0.9998|
###  Merged.SouthCarolinaLabel 
Entities Positive:  101 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     0.9996|    0.9998|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     0.9984|    0.9984|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9104|    0.9400|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     0.9966|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     0.9994|    0.9994|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     0.9960|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     0.9970|    1.0000|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9694|    0.9734|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                                                                                              |     0.6066|    0.6390|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.0070|    0.9872|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.2282|    0.8602|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9570|    0.9800|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                                                                                                  |     0.4448|    0.7708|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                                                                                              |     0.5876|    0.5876|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.6702|    0.8550|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.0338|    0.8432|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8296|    0.9570|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                                                                                              |     0.5876|    0.5876|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                  |     0.7330|    0.7980|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.9836|    0.9984|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                                                                                         |     0.9670|    0.9908|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.7630|    0.9064|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                                                                                              |     0.7452|    0.7452|
|    0.5384|      2893|     2|NGSYPN(HexNAc)LSK[2;5;4;0;0]                                                                                                                              |     0.9150|    0.9664|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                                                                                              |     0.9966|    0.9970|
|    0.5614|      3183|     2|N(HexNAc)GSYPNLSK[0;5;4;2;0]                                                                                                                              |     0.8364|    0.8490|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                                                                                              |     0.9518|    0.9518|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                  |     0.7978|    0.8992|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0000|    0.7432|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                  |     0.4448|    0.7538|
|    0.5762|      4351|     2|NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                                                                                              |     0.8234|    0.8254|
|    0.5799|      2876|     2|N(HexNAc)GSYPNLSK[1;4;4;1;0]                                                                                                                              |     0.8572|    0.9830|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                                                                                                |     0.7630|    0.7638|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;2;0;0]                                                                                                         |     0.9420|    0.9476|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                                                                                              |     0.9954|    0.9954|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.1990|    0.9612|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                              |     0.9856|    0.9856|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                                                                                              |     0.8256|    0.8256|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                                                                                                  |     0.2182|    0.6602|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                                                                                                |     0.9120|    0.9120|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                                              |     0.9168|    0.9706|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                                                                                                |     0.6506|    0.6506|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     0.9840|    0.9842|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                    |     0.8154|    0.9460|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                                                                                                |     0.9628|    0.9628|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                                                                                              |     0.7980|    0.7980|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                                                                                                |     0.9958|    0.9958|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                    |     0.9048|    0.9186|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     0.7500|    0.7586|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;5;4;2;0]                                                                                    |     0.8920|    0.9100|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;1;0]                                                                                                |     0.8744|    0.8834|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;4;6;0;0]                                                                                                |     0.6208|    0.6212|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     0.9626|    0.9820|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                                                                                              |     0.9544|    0.9792|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                                                                                |     0.9966|    0.9966|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     0.9966|    0.9966|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                                                                          |     0.7018|    0.8892|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                                                                                         |     0.9570|    0.9570|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6264|      5400|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;6;5;2;0]                                                                            |     0.0048|    0.7484|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                    |     0.9944|    0.9944|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                                                                                                |     0.9946|    0.9946|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                                |     0.9962|    0.9962|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                              |     0.9798|    0.9846|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                              |     0.9946|    0.9954|
|    0.6383|      4971|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;9;9;0;0]                                                                                                    |     0.5536|    0.5536|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     0.9960|    0.9960|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                              |     0.9544|    0.9764|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;6;1;0]                                                                                                |     0.9530|    0.9530|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     0.9964|    0.9964|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                                                                                         |     0.9992|    0.9992|
|    0.6452|      5343|     3|N(Deamidated)LLWLTEKN(HexNAc)GSYPN(HexNAc)LSKSYVNNK[1;8;6;0;0]                                                                                            |     0.0008|    0.9290|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                                                                                         |     0.9960|    0.9960|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                                |     0.9966|    0.9966|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     0.9768|    0.9994|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0]                                             |     0.7488|    0.8448|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                                |     0.9858|    0.9858|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                    |     0.9944|    0.9944|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6635|      5035|     8|N(Deamidated)LLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[1;5;4;2;0]                                                                            |     0.0084|    0.9382|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                                |     0.9946|    0.9946|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                          |     0.9968|    0.9968|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                                                                |     0.5914|    0.5914|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     0.8656|    0.8658|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     0.9196|    0.9196|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                                                                                                  |     0.9538|    0.9538|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.9060|    0.9060|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                                                                         |     0.9280|    0.9280|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVNN(Deamidated)KEK[4;8;11;2;0]                                                                                                 |     0.9532|    0.9600|
|    0.7328|      4850|     2|N(Deamidated)LLWLTGKNGLYPN(HexNAc)LSK[1;6;6;2;0]                                                                                                          |     0.6566|    0.6570|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]                                                                            |     0.7120|    0.9720|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.8328|    0.8640|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
###  Merged.Merged 
Entities Positive:  105 


| MS1_Score| Calc_mass| count|bestSeq                                                                                                                                                   | worstScore| bestScore|
|---------:|---------:|-----:|:---------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|    0.0000|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0001|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0155|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0217|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0275|      4609|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.0458|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1028|      4828|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;8;7;0;0]                                                                                                         |     0.9986|    0.9994|
|    0.1418|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.1772|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.2137|      3774|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;4;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.3033|      4389|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.3273|      8504|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[0;4;4;1;0] |     0.9830|    0.9854|
|    0.3777|      2975|     2|NGSYPN(HexNAc)LSK[2;3;6;0;0]                                                                                                                              |     0.9988|    0.9990|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;2;0]            |     0.0744|    0.9914|
|    0.3920|      7938|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]                        |     0.2686|    0.9144|
|    0.3969|      8886|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELREQLSSVSSFEK[1;6;5;0;0] |     0.9604|    0.9850|
|    0.3989|      4061|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;7;6;2;0]                                                                                                                  |     0.9968|    1.0000|
|    0.4327|      3898|     2|N(HexNAc)GSYPNLSK[1;6;6;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.4346|      7647|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[0;5;5;1;0]            |     0.7316|    0.8786|
|    0.4346|      7647|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;0;0]                        |     0.0650|    0.8436|
|    0.4914|      2788|     2|N(HexNAc)GSYPNLSK[1;4;5;0;0]                                                                                                                              |     0.9704|    0.9704|
|    0.5060|      7664|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSNSEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;0;0]            |     0.8248|    0.9642|
|    0.5064|      3735|     2|N(HexNAc)GSYPNLSK[1;5;6;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5168|      2748|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                  |     1.0000|    1.0000|
|    0.5173|      3080|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                  |     0.9898|    0.9898|
|    0.5207|      7955|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;6;5;1;0]            |     0.9852|    0.9998|
|    0.5290|      4682|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[0;8;7;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.5311|      7939|     3|C(Carbamidomethyl)N(Deamidated)IAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[2;5;5;1;0]            |     0.7646|    0.8864|
|    0.5328|      3404|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;6;5;1;0]                                                                                                                  |     0.9630|    0.9684|
|    0.5333|      2950|     2|N(HexNAc)GSYPNLSK[1;5;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5380|      3259|     2|N(Deamidated)GSYPN(HexNAc)LSK[2;6;5;0;0]                                                                                                                  |     0.0134|    0.6540|
|    0.5384|      2893|     2|N(HexNAc)GSYPNLSK[2;5;4;0;0]                                                                                                                              |     0.9996|    1.0000|
|    0.5494|      3112|     2|NGSYPN(HexNAc)LSK[1;6;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5525|      5594|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;4;0]                                                                                                          |     0.1512|    0.6768|
|    0.5614|      3183|     2|NGSYPN(HexNAc)LSK[0;5;4;2;0]                                                                                                                              |     0.9996|    1.0000|
|    0.5647|      4060|     2|N(HexNAc)GSYPNLSK[1;7;6;2;0]                                                                                                                              |     0.9988|    0.9988|
|    0.5655|      3039|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                  |     0.9996|    1.0000|
|    0.5677|      7589|     3|C(Carbamidomethyl)NIAGWLLGNPEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;1;0]                        |     0.0000|    0.8016|
|    0.5710|      3330|     2|N(Deamidated)GSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                  |     0.9982|    1.0000|
|    0.5762|      4351|     2|NGSYPN(HexNAc)LSK[1;7;6;3;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5799|      2876|     2|NGSYPN(HexNAc)LSK[1;4;4;1;0]                                                                                                                              |     0.9992|    0.9996|
|    0.5809|      4234|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;0;0]                                                                                                |     0.9996|    0.9996|
|    0.5834|      3667|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[0;8;2;0;0]                                                                                                         |     0.9934|    0.9940|
|    0.5902|      3403|     2|N(HexNAc)GSYPNLSK[1;6;5;1;0]                                                                                                                              |     0.9994|    1.0000|
|    0.5918|      7299|     3|C(Carbamidomethyl)NIAGWLLGN(Deamidated)PEC(Carbamidomethyl)DLLLTASSWSYIVETSN(Deamidated)SEN(HexNAc)GTC(Carbamidomethyl)YPGDFIDYEELR[1;5;4;0;0]            |     0.0880|    0.9616|
|    0.5922|      2747|     2|NGSYPN(HexNAc)LSK[1;5;4;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5946|      3444|     2|N(HexNAc)GSYPNLSK[1;5;6;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.5966|      3461|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;6;6;1;0]                                                                                                                  |     0.6578|    0.8326|
|    0.5982|      4437|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;5;0;0]                                                                                                |     0.9998|    0.9998|
|    0.5984|      3694|     2|NGSYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6017|      4566|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;5;1;0]                                                                                                |     0.9988|    0.9988|
|    0.6034|      3733|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;0;0]                                                                                                         |     0.9990|    0.9990|
|    0.6035|      4891|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(Deamidated)STDTVDTVLEK[1;6;5;1;0]                                                                                    |     0.9958|    0.9992|
|    0.6040|      4728|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;5;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6043|      3986|     2|N(HexNAc)GSYPNLSK[1;6;5;3;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6045|      4640|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6055|      4301|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6065|      5182|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                    |     0.9986|    0.9992|
|    0.6073|      4971|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;6;5;3;0]                                                                                                         |     0.9992|    0.9998|
|    0.6082|      4817|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                                                                    |     0.9970|    0.9974|
|    0.6087|      4769|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;1;0]                                                                                                |     0.9996|    0.9998|
|    0.6095|      4478|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;4;6;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6118|      3936|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;0;0]                                                                                                         |     0.9954|    0.9978|
|    0.6133|      2934|     2|N(Deamidated)GSYPN(HexNAc)LSK[0;4;5;1;0]                                                                                                                  |     0.0190|    0.7334|
|    0.6133|      2934|     2|NGSYPN(HexNAc)LSK[2;4;5;0;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6147|      4592|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.6148|      4816|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6167|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6180|      4688|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;5;6;2;0]                                                                                                          |     0.9998|    1.0000|
|    0.6196|      4227|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;5;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6253|      4139|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6271|      5547|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                    |     1.0000|    1.0000|
|    0.6291|      4098|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6297|      4599|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;0;0]                                                                                                |     1.0000|    1.0000|
|    0.6299|      5181|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6340|      3329|     2|NGSYPN(HexNAc)LSK[1;5;4;2;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6381|      3079|     2|NGSYPN(HexNAc)LSK[1;4;5;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6383|      4971|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;3;0]                                                                                                         |     1.0000|    1.0000|
|    0.6389|      3038|     2|NGSYPN(HexNAc)LSK[1;5;4;1;0]                                                                                                                              |     1.0000|    1.0000|
|    0.6428|      5093|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;6;1;0]                                                                                                |     0.9996|    0.9998|
|    0.6433|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6451|      4315|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;2;0]                                                                                                         |     1.0000|    1.0000|
|    0.6453|      4430|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;6;1;0]                                                                                                         |     0.9996|    0.9996|
|    0.6506|      5837|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                                |     1.0000|    1.0000|
|    0.6512|      4463|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;0;0]                                                                                                         |     1.0000|    1.0000|
|    0.6541|      7869|     9|DTIC(Carbamidomethyl)IGYHAN(HexNAc)N(HexNAc)STDTVDTVLEKN(Deamidated)VTVTHSVNLLEDSHN(Deamidated)GK[2;7;10;0;0]                                             |     0.9894|    1.0000|
|    0.6549|      4890|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;6;5;1;0]                                                                                                |     1.0000|    1.0000|
|    0.6552|      5838|     2|DTIC(Carbamidomethyl)IGYHAN(Deamidated)N(HexNAc)STDTVDTVLEK[1;7;6;3;0]                                                                                    |     1.0000|    1.0000|
|    0.6589|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.6621|      5303|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;7;6;3;0]                                                                                                          |     0.2248|    0.8444|
|    0.6642|      5546|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;7;6;2;0]                                                                                                |     1.0000|    1.0000|
|    0.6691|      4647|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;5;2;0]                                                                                                          |     1.0000|    1.0000|
|    0.6754|      4890|     2|DTIC(Carbamidomethyl)IGYHAN(HexNAc)NSTDTVDTVLEK[1;6;5;1;0]                                                                                                |     0.9994|    0.9994|
|    0.6817|      4931|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;6;1;0]                                                                                                |     0.9992|    0.9992|
|    0.6843|      4525|     2|DTIC(Carbamidomethyl)IGYHANN(HexNAc)STDTVDTVLEK[1;5;4;1;0]                                                                                                |     0.9982|    0.9994|
|    0.6850|      5336|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NK[0;10;10;0;0]                                                                                                  |     0.9992|    0.9992|
|    0.6850|      5336|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;3;0]                                                                                                         |     0.9994|    0.9994|
|    0.7079|      4024|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;5;4;1;0]                                                                                                         |     0.1280|    0.7372|
|    0.7125|      5627|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;4;0]                                                                                                         |     0.9998|    0.9998|
|    0.7219|      6639|     2|N(HexNAc)GSYPN(HexNAc)LSKSYVN(Deamidated)NKEK[4;8;11;2;0]                                                                                                 |     0.9996|    1.0000|
|    0.7328|      4850|     2|NLLWLTGKN(Deamidated)GLYPN(HexNAc)LSK[1;6;6;2;0]                                                                                                          |     0.9992|    0.9998|
|    0.7359|      4754|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;7;6;1;0]                                                                                                         |     1.0000|    1.0000|
|    0.7364|      6348|     8|NLLWLTEKN(Deamidated)GSYPN(HexNAc)LSKSYVN(Deamidated)N(Deamidated)K[1;7;6;4;0]                                                                            |     0.0236|    0.8252|
|    0.8216|      5336|     2|N(HexNAc)VTVTHSVN(Deamidated)LLEDSHNGK[1;7;6;3;0]                                                                                                         |     0.9998|    1.0000|
|    0.8268|      4389|     2|N(HexNAc)VTVTHSVNLLEDSHN(Deamidated)GK[1;6;5;1;0]                                                                                                         |     1.0000|    1.0000|

