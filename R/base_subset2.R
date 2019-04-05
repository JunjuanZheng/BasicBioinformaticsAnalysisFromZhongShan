
#' @export
## Usage: subset2 can get a subset of a data frame,and keep data frame style regardless of the nrow or ncol of this subset.

# data # a data frame
# colnames #the colnames vector of a subset
# rownames #the rownames vector of a subset

subset2 <- function(data,colnames = NULL,rownames = NULL){
  ## 设置默认值
  default <- c(rep(list(NULL),2))
  input <- list(colnames = colnames,rownames =rownames)
  do <- list(colnames = colnames(data),rownames = rownames(data))
  library(lucky)
  output <- set.default(input,default,do)
  colnames = output[["colnames"]]
  rownames = output[["rownames"]]

  ## 选择子集:列
  data1 <- cbind(data,sos.c=0)
  data1 <- rbind(data1,sos.r=0)
  data2 <- subset(data1,select = colnames)

  ## 选择子集：行
  data3 <- as.data.frame(t(data2))
  data3 <- subset(data3,select = rownames)
  data3 <- as.data.frame(t(data3))

  ## 输出
  return(data3)
}


