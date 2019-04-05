

# usage: run order() for lots of cols in a matrix
# data # a matrix or data frame
# col # colnames of the data
#' @export
Plus.order <- function(data,col,decreasing=T){

  x2 <- as.matrix(data)
  order.col <- NULL;
  for(a in col){
    order.col.i <- match(a,colnames(x2))
    order.col <- c(order.col,order.col.i)
  }
  library(dplyr)
  query <- paste('x2[,',order.col,']',sep='') %>%
    paste(.,collapse=',') %>%
    paste('order(',.,paste0(',decreasing = ',decreasing,')',collapse = ""),sep='')
  class(query);x2 <- x2[eval(parse(text = query)),]
  x2 <- as.data.frame(x2)
  return(x2)
}



