
## Usage: quick way to draw boxplot for multiple markers.
# boxplot.genedata2 is a plus version of boxplot.genedata.It use the exogenous p values from p.val.matrix from limma package or other differentially analysis tools.

# select.genes# marker names
# design# a design object of a expression matrix.
# expr.matrix# a expression matrix.
# p.val.matrix# a result object containing p.val
# p.val.col#colnames of p value in the p.val.matrix
# p.val.position#the position of significant symbol in the plot
# contrast #the colname of the group in every facet.
# contrast.list#the group in every facet
# returndata #whether output the data of boxplot based.
#' @export
boxplot.genedata2 <- function(select.genes,
                              design,
                              expr.matrix,
                              p.val.matrix,
                              p.val.col="P.Value",
                              p.val.position=NULL,
                              contrast = "N.status",
                              contrast.list = c("N0","Np")){
  ## 加载必要包
  library(ggplot2);library(ggpubr)

  ## 生成boxplot.data
  expr.matrix <- expr.matrix[,rownames(design)]
  x <- data.frame(expr=NULL,genes =NULL);
  for(i in select.genes){
    x.i <- expr.matrix[i,]
    x.i <- as.data.frame(x.i)
    colnames(x.i) <- "expr"
    x.i$genes <- i
    x.i$symbols <- convert(x.i$genes)
    p <- match(contrast,colnames(design))
    x.i$col1 <- design[,p]
    colnames(x.i)[match("col1",colnames(x.i))] <- contrast
    #外源性p值
    x.i$P.val <- p.val.matrix[i,p.val.col]
    #expr最大值
    x.i$max.expr <- max(x.i$expr)
    x <- rbind(x,x.i)
  }
  x <- x[order(x$symbols),]

  ## significant label
  p.symbol <- unique(x$symbols) #基因
  p.ensembl <- convert(p.symbol,
                       fromtype = "SYMBOL",
                       totype = "ENSEMBL")
  p.label <- as.numeric(as.character(p.val.matrix[p.ensembl,p.val.col]))
  p.label <- cut(p.label,breaks = c(0,0.0001,0.001,0.01,0.05,Inf),
                 labels = c("****","***","**","*","ns"),
                 right = T)
  p.label <- as.character(p.label) #显著性标记

  ## 构建颜色向量
  library(lucky)
  palette = rep(mycolor,100)
  palette1 = palette[1:length(select.genes)]

  ## 画图
  p <- ggboxplot(x,
                 x = contrast, y = "expr",
                 fill = "symbols",#箱体填充颜色
                 color = "black",
                 palette = palette1,
                 add = "jitter") + #加上jitter点
    labs(x = contrast,y = 'The expression level of genes') +
    theme(axis.text.x = element_text(size = 12,colour = "black",face = "bold")) +
    facet_grid(. ~ symbols) +
    theme(axis.title.x = element_text(size = 15,colour = "black",face = "bold"),
          axis.title.y = element_text(size = 15,colour = "black",face = "bold"),
          legend.position='none',
          strip.background = element_rect(fill="white")) +
    rotate_x_text(angle = 45) #x轴的text以逆时针旋转45度角

  #添加显著性标记
  if(is.null(p.val.position)){
    pos <- max(x$expr)
  } else {
    pos <- p.val.position
  }
  f_labels <- data.frame(symbols = p.symbol, label = p.label);
  p1 <- p + geom_text(x=1.5, y=pos, aes(label=label), data=f_labels)
  print(p1)

  ## 输出数据框或图片
  l <- list(
    plot = p1,
    metadata = x
  )
  return(l)
}






