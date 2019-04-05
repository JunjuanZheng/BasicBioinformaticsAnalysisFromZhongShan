

## a fast way to extra part of elements in a vetor.
#vt # a vector with common split
#split # the common split
#n # which splited string should be extrated
#' @export
Fastextra <- function(vt,split,n=NULL){
  vt <- as.character(vt)
  get1 <- function(i,split,n=NULL){
    if(is.null(n)){
      vt1.i <- unlist(strsplit(i,split)) #全部输出
    } else {
      vt1.i <- unlist(strsplit(i,split))[n]
    }
    return(vt1.i)
  }
  vt1 <- apply(as.matrix(vt),1,function(z)get1(z,split = split,n = n))
  vt1 <- as.vector(vt1)
  return(vt1)
}

#vt <- c("abc_1","abc_2","abc_3","abc_4")
#split = "_";n =NULL
#Fastextra(a,"_")
#get1(a[1],"_")

# Fastextra <- function(vt,split,n=NULL){
#  vt <- as.character(vt)
#  vt1 <- NULL
#  for(i in vt){
#    if(is.null(n)){
#      vt1.i <- unlist(strsplit(i,split)) #全部输出
#    } else {
#      vt1.i <- unlist(strsplit(i,split))[n]
#    }
#    vt1 <- c(vt1,vt1.i)
#  }
#  return(vt1)
#}









