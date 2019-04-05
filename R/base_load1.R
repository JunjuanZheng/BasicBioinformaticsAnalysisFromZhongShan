
#' @export
## load1 help load local .rda files by a pattern
load1 <- function(pattern,path=".",envir=parent.frame()){
  for(i in pattern){
    file = list.files(path = path,pattern = i,full.names = T)
    if(length(file) > 1){
      print("Attention!Multiple files and the first was selected.")
    }
    load(file[1],envir = envir)
  }
}
