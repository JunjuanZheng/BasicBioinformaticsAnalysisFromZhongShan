

###======QC.density
## QC.density可以观察不同分组（treat vs. control）的芯片表达量密度。基于ggpubr的ggdensity函数。
#' @export
QC.density <- function(DesignEset,
                       contrast,
                       savefile=F,
                       names="test1",
                       limit.array = 30,
                       limit.genes=3000,
                       seed=2018){
  ## 包
  library(ggpubr)

  ## contrast的设置
  if(is.null(contrast)){
    #如果不提供contrast，就直接用全部的contrast
    contrast.matrix <- DesignEset[["DesignList"]][["contrast.matrix"]]
    c1 <- colnames(contrast.matrix)
    c2 <- list()
    for(i in c1){
      c2.i <- unlist(strsplit(i,"-"))
      c2 <- c(c2,list(c2.i))
    }
  } else {
    #提供了contrast，直接使用定义的contrast
    c2 <- contrast
  }

  ## 提取某个DesignList
  p1 <- NULL
  for(i in 1:length(c2)){
    ## 提取某个DesignList和contrast
    contrast.a <- c2[[i]]
    DesignEset.i <- getOneDesignEset(DesignEset,contrast.a)
    names.i <- paste0(contrast.a,collapse = "-")

    ## 赋值
    library(Biobase)
    condition <- DesignEset.i[["DesignList"]][["condition"]]
    expr1 <- exprs(DesignEset.i[["ESet"]])

    ## 计算个数
    n.condition <- length(unique(condition$condition))
    n.array <- table(condition$condition)
    if(sum(n.array) > limit.array){
      print(str_c("按limit.array=",limit.array,"的标准,属于巨大矩阵。随机抽样",limit.genes,"个基因进行分析。"))
      set.seed(seed);select.row <- sample(1:nrow(expr1),limit.genes,replace=F)
      expr <- expr1[select.row,]
    } else {
      expr <- expr1
    }
    n.expr <- dim(expr)[1]

    ## 密度图data
    n.condition1 <- NULL;n.array1 <- NULL;n.expr1 <- NULL;
    for(i in 1:length(contrast.a)){
      ##condition向量
      contrast.i <- contrast.a[i]
      c.contrast.i <- n.array[match(contrast.i,names(n.array))]*n.expr
      n.condition1i <- rep(contrast.i,c.contrast.i)
      n.condition1 <- c(n.condition1, n.condition1i)

      ##array向量
      array1.i <- rownames(condition)[condition$condition %in% contrast.i]
      n.array.i <- NULL
      for(a in 1:length(array1.i)){
        array.i.a <- array1.i[a]
        n.array.i.a <- rep(array.i.a,n.expr)
        n.array.i <- c(n.array.i,n.array.i.a)
      }
      n.array1 <- c(n.array1,n.array.i)

      ##expr向量
      n.expr.i <- NULL
      for(a in 1:length(array1.i)){
        array.i.a <- array1.i[a]
        n.expr.i.a <- expr[,array.i.a]
        #n.expr.i.a <- expr[,a]
        n.expr.i <- c(n.expr.i,n.expr.i.a)
      }
      n.expr1 <- c(n.expr1,n.expr.i)
    }
    df1 <- data.frame(
      condition = n.condition1,
      array = n.array1,
      expr = n.expr1
    )
    df1$array <- factor(df1$array,levels = rownames(condition))

    ## 画密度图
    if(length(rownames(condition))<=10){
      legend.position = "top"
    } else {
      legend.position = "none"
    }
    library(ggpubr)
    p <- ggdensity(data = df1,
                   x = "expr",rug = F,size = 1,
                   color = "array",
                   facet.by = "condition")+
      labs(x = "Expression Intensity", y= "Density")+
      ggtitle(str_c("DensityPlot:",paste0(contrast.a,collapse = "-"))) +
      theme(panel.grid =element_blank(),
            plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
            axis.title = element_text(face = "bold",size = 12),
            axis.text = element_text(face = "bold",size = 10),
            legend.title=element_blank(),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = legend.position
      )
    p1 <- c(p1,list(p))

    ##预览图片
    win.graph(width = 9,height = 7);print(p)
  }

  #保存文件
  if(savefile==T){
    pdf(str_c(names,"_QualityDensityPlot.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_QualityDensityPlot.pdf已保存在本地"))
    #End
    return(p1)
  } else {
    #End
    return(p1)
  }
}

# data(DesignEset)
# contrast=list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")) # a list that contain lots of contrasts.
# QC.density(DesignEset,
#            contrast,
#            savefile=F,
#            names="test1",
#            limit.array = 30,
#            limit.genes=3000,
#            seed=2018)

