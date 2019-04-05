
## Usage:filter.low() help quickly filtering lots of different dataset via lots of condition based on lower cut off number.It returns,by every dataset,a targets vetor and the related sub-data via list style.
# list.data # a list of data frame
# list.target.col # a list of character vectors containing target marker colnames of related data.
# list.cut.off# a list of numeric vectors containing target marker cut-off value.
# range #the range you want to do the filter.
#' @export
filter.low <- function(list.data,
                       list.cutoff.col,
                       list.cut.off,
                       list.target.col,
                       range){
  gene <- NULL
  for(i in 1:length(list.data)){
    ## 提取数据
    data.i <- list.data[[i]]
    cutoff.col.i <- list.cutoff.col[[i]]
    cut.off.i <- list.cutoff[[i]]
    target.col.i <- list.target.col[[i]]

    ## 根据range对数据进行初步的筛选
    data.i<- data.i[data.i[,target.col.i] %in% range,]

    ## 对多个列进行排序
    library(dplyr)
    data.i1 <- as.matrix(data.i)
    position1 <- Fastmatch(cutoff.col.i,colnames(data.i1))
    query <- paste('data.i1[,',position1,']',sep='') %>%
      paste(.,collapse=',') %>%
      paste('order(',.,',decreasing =T)',sep='')
    class(query);data.i1 <- data.i1[eval(parse(text = query)),]

    ## 按cut.off值进行筛选
    query1 <- paste('data.i1[,',position1,']',sep='') %>%
      paste("abs(as.numeric(as.character(",.,")))",sep = "") %>%
      paste(.,">=",cut.off.i,collapse=' & ')
    class(query1);data.i2 <- data.i1[eval(parse(text = query1)),]

    ## 最后输出候选列的id
    a <- as.character(data.i2[,target.col.i])
    a1 <- list(select.marker = a,select.data = data.i2)
    gene <- c(gene,list(a1))
  }
  names(gene) <- names(list.data)
  return(gene)
}






