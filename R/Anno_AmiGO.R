




AmiGO <- function(GO.term){

  ## 加载必要的包
  nd <- c("XML","plyr")
  Plus.library(nd)

  ## unix时间戳
  unix.time <- as.numeric(as.POSIXct(Sys.Date(), format="%Y-%m-%d"))
  us <- sample(0:9,3,replace = T)
  us <- paste0(us,collapse = "")
  us <- paste(us,us,sep = ".")

  ## find GO term
  url1 <- "http://amigo.geneontology.org/amigo/search/bioentity?q="
  url <- paste0(url1,GO.term,";time=",unix.time,us)


  tables <- readHTMLTable(url)







}




