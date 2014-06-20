runAll <- function(){
  scripts <- list.files(pattern = "^[^_]*.R$")
  sapply(scripts, function(src){
    stitch(src, template = "Model_Relabel_Testing.Rmd", envir = new.env())
  })
}

loadAllErrorTables <- function(){
  files <- list.files(pattern = "errorTable.*.rds", recursive = T)
  tableNames <- make.names(list.files(pattern = "errorTable.*.rds", recursive = T))
  tableNames <- gsub("SavedData.|.rds|\\.", "", tableNames)
  errorTables <- lapply(files, function(f){
    tab <- readRDS(f)
    row.names(tab) <- tab$.rownames
    tab
  })
  names(errorTables) <- tableNames  
  errorTables
}

orderAllErrorTables <- function(tables){
  ord <- row.names(tables[[1]])
  lapply(tables, function(tab){
    tab[ord,]
  })
}


MAIN <- function(){
  apply(cbind(errorTables[[1]]$ErrorRate, errorTables[[2]]$ErrorRate, errorTables[[3]]$ErrorRate),1, function(i)all(i[1] == i))
}