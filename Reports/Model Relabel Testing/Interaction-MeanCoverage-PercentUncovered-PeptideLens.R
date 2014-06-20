testFormula <- call ~ (meanCoverage * percentUncovered  *  peptideLens) + abs_ppm_error + X.b_ion_with_HexNAc_coverage + 
  X.y_ion_with_HexNAc_coverage + numStubs

mtry <- 5

ntree <- 5000

RUN_ID <- "Interaction-MeanCoverage-PercentUncovered-PeptideLens"