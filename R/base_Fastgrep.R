


# A Fast grep when lots of pattern used.return an unique position of all the pattern.Fastgrep() is based on base::grep().

## arguement
# pattern #a string
# x # a vector that may contain patterns.

####=====================grep======================####

#Fastgrep <- function(pattern,x){
#   p <- NULL
#   for(i in pattern){
#     p.i <- grep(i,x)
#     p <- c(p,p.i)
#   }
#   p <- unique(p)
#   return(p)
#}
#' @export
Fastgrep <- function(pattern,x){
  p <- apply(as.matrix(pattern),1,function(z)grep(z,x))
  p <- unlist(p)
  p <- as.vector(as.matrix(p))
  return(p)
}

#Fastgrep(c("a","c"),c("a","a","a","a","a","b","b","c","d","e"))

#' @export
grep2 <- function(pattern,x){
  p <- grep(pattern,x)
  ifelse(length(p)==0,{p <- NA},{p <- p})
  return(p)
}

#Fastgrep的升级版。利用了apply函数。支持并行运算
#' @export
Fastgrep2 <- function(pattern,x,parrallel = F,logical=F){
  if(parrallel == F){
    p <- apply(as.matrix(pattern),1,function(i)grep2(i,x))
  } else {
    #使用并行运行以加快速度
    parrallel.apply <- function(median.pattern,median.x){
      require(parallel)
      ncore <- detectCores(logical = logical)
      cl <- makeCluster(mc <- getOption("cl.cores",ncore))
      clusterExport(cl,c("median.pattern","median.x","grep2","pattern","x"))
      i1 <- parApply(cl = cl,as.matrix(median.pattern),1,function(i)grep2(i,median.x))
      return(i1)
    }
    p <-  parrallel.apply(pattern,x)
  }
  return(p)
}
#Fastgrep2(c("a","c"),c("b","c","d"),parrallel = T)



####===================match======================####
#Fastmatch <- function(pattern,x){
#  p <- NULL
#  for(i in pattern){
#    p.i <- match(i,x)
#    p <- c(p,p.i)
#  }
#  return(p)
#}
#' @export
Fastmatch <- function(pattern,x){
  p <- apply(as.matrix(pattern),1,function(z)match(z,x))
  return(p)
}
#Fastmatch(c("a","c"),c("a","a","a","a","a","b","b","c","d","e"))


#' @export
match2 <- function(pattern,x){
  p <- match(pattern,x)
  ifelse(length(p)==0,{p <- NA},{p <- p})
  return(p)
}

#Fastmatch的升级版。利用了apply函数。
#' @export
Fastmatch2 <- function(pattern,x){
  p <- apply(as.matrix(pattern),1,function(i)match2(i,x))
  return(p)
}



## compared with grep and Fastgrep
# grep("a",c("b","c","d"))
# Fastgrep("a",c("b","c","d"))
# grep("b",c("b","c","d"))
# Fastgrep("b",c("b","c","d"))
# Fastgrep(c("a","c"),c("b","c","d"))
# Fastgrep(c("b","c"),c("b","c","d"))
