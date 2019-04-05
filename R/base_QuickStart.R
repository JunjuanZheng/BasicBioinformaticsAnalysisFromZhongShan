

## QuickStart help run series R script automatically
# script.path# the path of script
# data.path#the path of data
# result.path#the path of result
#' @export
QuickStart <- function(script.path,
                       data.path,
                       result.path){
  ## create result path
  dir.create(result.path,recursive = T)
  result.path <<- result.path
  data.path <<- data.path

  ## run scripts
  s <- list.files(path = script.path,pattern = ".R",full.names = T)
  for(i in s) source(i,encoding = "UTF-8")

}









