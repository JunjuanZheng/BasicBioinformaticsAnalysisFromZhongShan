

## Usage: translation between cormatirx and related long data frame
# data # a cormatrix or a long data frame from a cormatirx
# source.col #If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of source columns.
# target.col #If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of target columns.
# value #If data is a cormatrix,it isn't needed.If data is a long data frame,it is the colnames of value columns.
#' @export
Make.cormatrix <- function(data,
                           source.col=NULL,
                           target.col=NULL,
                           value=NULL,
                           report = T){
  ## 判断是相关矩阵还是长型数据
  if(is.null(source.col)|is.null(target.col)|is.null(value)){
    if(report == T) {print("data是相关矩阵，将相关矩阵转化为长型矩阵。")}
    l1 <- list()
    for(i in 1:ncol(data)){
      l.i <- data[,i]
      l1 <- c(l1,list(l.i))
    }
    names(l1) <- colnames(data)
    cm1 <- stack(l1)
    cm1$target <- rep(rownames(data),ncol(data))
    colnames(cm1)[1:2] <- c("value","source")
    # 将值为1和重复的结果去除
    cm2 <- cm1[cm1$value !=1,]
    cm2 <- cm2[order(cm2$value,decreasing = T),]
    select <- seq(2,nrow(cm2),by=2) #选择偶数行
    cm3 <- cm2[select,];rownames(cm3) <- 1:nrow(cm3)
    cm3 <- subset(cm3,select = c("source","target","value"))
    logic1 <- length(unique(as.character(cm3$value))) == length(as.character(cm3$value))
    if(logic1){print("value值并无重复,结果可信。")} else {print("value值有重复，结果可疑，请检查算法。")}
    return(cm3)
  }      else {
    if(report == T){print("data是长型矩阵，将长型矩阵转化为相关矩阵。")}
    colnames(data)[Fastmatch(c(source.col,target.col,value),colnames(data))] <- c("source","target","value")
    data1=data.frame(source = data[,"target"],
                     target = data[,"source"],
                     value = data[,"value"])
    co.genes <- unique(c(as.character(data[,"target"]),as.character(data[,"source"])))
    data2 = data.frame(source = co.genes,
                       target = co.genes,
                       value = 1)
    data3 <- rbind(data,data1,data2)
    Plus.library("tidyr")
    options(warn = -1)
    a <- spread(data3,
                key = source,
                value = value)
    options(warn = 0)
    rownames(a) <- a$target;a1 <- a[,-1];
    a1 <- a1[colnames(a1),]
    a1 <- as.matrix(a1)
    a1[is.na(a1)] <- 1
    return(a1)
  }
}

