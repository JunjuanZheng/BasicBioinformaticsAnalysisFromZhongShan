

## boxplot.MA is based on boxplot.genedata2 and a quick way to draw a boxplot for specified contrast via MAList object.boxplot.MA2() is based on boxplot.MA and appropriate for multiple significant markers plotting.

# MA1#MA object after lucky::MApipeline2()
# design #a design matrix from model.matrix().
# matrix.type # M or A matrix
# contrast#the group in every facet
# select#marker names
# p.val.matrix#a result object containing p.val
# p.val.col#colnames of p value in the p.val.matrix
# p.val.position#the position of significant symbol in the plot
#' @export
boxplot.MA <- function(MA1,
                       design,
                       matrix.type="M",
                       contrast,
                       select,
                       p.val.matrix,
                       p.val.col="P.Value",
                       p.val.position=NULL){

  ##矩阵
  mt <- MA1[[matrix.type]]

  ##构建condition对象
  getcondition1 <- function(design1){
    c.i <- names(design1)[grep(1,design1)]
    return(c.i)
  }
  condition <- apply(design,1,getcondition1)
  condition <- data.frame(condition = condition,
                          row.names = colnames(mt))

  ##构建矩阵
  ss <- condition$condition %in% contrast
  condition1 <- subset(condition,subset = ss)
  mt1 <- mt[,rownames(condition1)]

  ##绘制图
  p <- boxplot.genedata2(select.genes=select,
                         design=condition1,
                         expr.matrix=mt1,
                         p.val.matrix=p.val.matrix,
                         p.val.col=p.val.col,
                         p.val.position=NULL,
                         contrast = "condition",
                         contrast.list = contrast)

  ##输出图片
  return(p)

}


boxplot.MA2 <- function(MA1,
                        matrix.type="M",
                        cut.n,
                        design,
                        contrast,
                        p.val.matrix,
                        p.val.col="P.Value",
                        p.val.position=NULL,
                        names="test1"){

  ## 箱式图观察差异基因
  NO.list <- cut.vector(1:nrow(p.val.matrix),cut.n)
  pdf(paste0("boxplot.sig.genes_",names,".pdf"),width = width,height = height)
  for(i in 1:length(NO.list)){
    p.i <- boxplot.MA(MA1=MA1,
                      design=design,
                      matrix.type="M",
                      contrast=contrast,
                      select=rownames(p.val.matrix)[NO.list[[i]]],
                      p.val.matrix=p.val.matrix,
                      p.val.col=p.val.col,
                      p.val.position=p.val.position)
    print(p.i)
  }
  dev.off()
  print("完成boxplot绘制并保存于当前工作空间.")
}
















