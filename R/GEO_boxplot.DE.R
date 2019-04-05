
# Usage:quick way to draw a boxplot for specified contrast via DesignEset object

#boxplot.DE() is based on boxplot.genedata2.

# DesignEset1# a DesignEset object.It is always be one after lucky::AnnotateDesignEset()
# contrast# the group in every facet
# select# marker names
# p.val.matrix#  a result object containing p.val
# p.val.col# colnames of p value in the p.val.matrix
# p.val.position# the position of significant symbol in the plot

#' @export
boxplot.DE <- function(DesignEset1,
                       contrast,
                       select,
                       p.val.matrix=sig.exprs,
                       p.val.col="P.Value",
                       p.val.position=NULL){
  ## 加载必要的包
  library(Biobase)
  library(limma)

  ## 普通项
  eset <- DesignEset1[["ESet"]]
  expr <- exprs(eset)
  condition <- DesignEset1[["DesignList"]][["condition"]]

  ## 构建condition object
  condition1 <- c(paste(contrast[1],contrast[2],sep = "-"),paste(contrast[2],contrast[1],sep = "-"))
  l1 <- names(condition) %in% condition1
  condition1 <- condition[[match(T,l1)]]
  expr1 <- expr[,rownames(condition1)]

  ## 画图
  p <- boxplot.genedata2(select.genes=select,
                         design=condition1,
                         expr.matrix=expr1,
                         p.val.matrix=p.val.matrix,
                         p.val.col=p.val.col,
                         p.val.position=p.val.position,
                         contrast="condition",
                         contrast.list=contrast)
  return(p)
}


boxplot.DE2 <- function(DesignEset1,
                        cut.n=11,
                        contrast,
                        p.val.matrix=sig.exprs,
                        p.val.col="P.Value",
                        p.val.position=NULL,
                        names="test1",
                        width = 12,height = 8){
  ## 箱式图观察差异基因
  NO.list <- cut.vector(1:nrow(p.val.matrix),cut.n)
  pdf(paste0("boxplot.sig.genes_",names,".pdf"),width = width,height = height)
  for(i in 1:length(NO.list)){
    NO.i <- NO.list[[i]]
    contrast=c("Np","N0")
    p.val.matrix=sig.exprs
    select <- rownames(sig.exprs)[NO.i]
    p.val.col="P.Value"
    library(ggplot2)
    p <- boxplot.DE(DesignEset1,
                    contrast,
                    select,
                    p.val.matrix=sig.exprs,
                    p.val.col,
                    p.val.position=NULL)
    p <- p +
      ggtitle(paste0("Significant Genes_",names(NO.list)[[i]])) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
  dev.off()
  print("完成boxplot绘制并保存于当前工作空间.")
}







