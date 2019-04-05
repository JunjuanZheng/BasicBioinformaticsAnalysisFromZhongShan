

# stringposition help find the position of a specified string pattern
#' @export
stringposition <- function(pattern,x,unlist=T){
  p2 <- NULL
  for(i in 1:length(pattern)){
    p <- gregexpr(pattern[i],x)[[1]]#-中的位置
    p1 <- p[1:length(p)]
    if(length(p)==1){
      if(p1 == -1){
        #说明并无此字符串
        p1 <- NA
      }
    }
    p2 <- c(p2,list(p1))
    names(p2)[i] <- pattern[i]
  }
  if(unlist==T){
    p3 <- as.numeric(unlist(p2))
    return(p3)
  } else {
    return(p2)
  }
}

