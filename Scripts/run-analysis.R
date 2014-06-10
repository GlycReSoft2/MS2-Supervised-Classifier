require(reshape2)
require(ggplot2)


#' Returns a structure of the form
#' - list
#'  - list
#'   - rpart length(fit x data) (redundant)
#'  - list 
#'   - data.frame length(fit x data) (redundant)
#'  - list
#'   - list (accuracy data)
withinSetAnalyze = function(preparedDataSets){
  
  treeModels <- lapply(preparedDataSets, fitModel)
  testPairs <- expand.grid(model = treeModels, testData = preparedDataSets)
  testNames <- apply(expand.grid(names(preparedDataSets), names(preparedDataSets)), 
                     1, function(row){
    paste(row[[1]], row[[2]], sep = "_on_")
  })
  row.names(testPairs) <- testNames
  predictions <- apply(testPairs, 1, function(row){
    model <- row[[1]]
    data <- row[[2]]
    predict(model, data)
  })
  
  testPairs$predictions <-  predictions
  
  resultsData <- as.data.frame(t(t(testPairs)))
  
  accuracy <- apply(resultsData, 1, function(row){
    checkModel(row[[2]], row[[3]])
  })
  
  resultsData$accuracy <- accuracy
  
  return(resultsData)
}

getSetError <- function(modelSet){
  errorTable <- as.data.frame(
      do.call(rbind, lapply(modelSet$accuracy, unlist))
    )
  nameData <- strsplit(x=row.names(errorTable),split="_on_")
  errorTable$model <- sapply(nameData, function(x)x[[1]])
  errorTable$data <- sapply(nameData, function(x)x[[2]])
  
  return(errorTable)
}

plotSetError <- function(errorTable){
  ggplot(errorTable, aes(x = FalseNegative, y = FalsePositive,
                        size = TruePositive/.numYes, color = data, 
                        shape = model)) + 
    geom_point() + scale_size(range = c(4,10)) + 
    labs(title = paste(substitute(errorTable), " Results"))
}

loadDir <- function(dir, includeAnno = F){
  dataFiles = list.files(dir, pattern = ".csv", full.names=T)
  preparedDataSets <- lapply(dataFiles, prepareModelHandler, includeAnno)
  names(preparedDataSets) <- basename(dataFiles)
  check <- sapply(preparedDataSets, function(x)class(x)=="data.frame")
  return(preparedDataSets[check])  
}

MAIN <- function(){
  parent <- parent.frame()
  dataDirs <- c("USSR","AGP","Transferrin","Solomon Islands","South Carolina")
  
  # All within-class pair-wise tests
  results <- sapply(dataDirs, function(dir){
    print(paste("On ", dir))
    parent[[dir]] <- loadDir(dir)
    parent[[paste(dir, "Within", sep = '')]] <- withinSetAnalyze(parent[[dir]])
    parent[[paste(dir, "Error", sep = '')]] <- getSetError(parent[[paste(dir, "Within", sep = '')]])
  })
  
  transformResults <- apply(results, 2, as.data.frame)
  
  return(transformResults)
}