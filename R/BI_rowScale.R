

####===========================rowScale===========================####
#对矩阵按行进行归一化处理。
# matrix
# log.convert是否对matrix进行log2的转换
#' @export
rowScale <- function(matrix,log.convert = T){
  if(log.convert ==F){m1 <- matrix} else {
    m1 <- log2(matrix+1)
  }
  x1 <- apply(m1,1,scale)
  rownames(x1) <- colnames(matrix)
  return(t(x1))
}
