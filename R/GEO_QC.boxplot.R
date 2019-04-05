


###=====QC.boxplot
# \item{contrast}{list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")).a list that contain lots of contrasts}

# limit.array = 30,limit.genes=3000指的是如果array个数>30，就不进行全基因分析，只随机选取3000个基因进行质量分析。
# contrast=list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")) # a list that contain lots of contrasts.
#' @export
QC.boxplot <- function(DesignEset,
                       contrast,
                       savefile=F,
                       names="test1",
                       limit.array = 30,
                       limit.genes=3000,
                       seed=2018){
  ## 包
  library(ggpubr);library(stringr)

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

  ## 找到某个contrast.i的信息
  p1 <- NULL
  for(i in 1:length(c2)){
    ## 提取某一个contrast信息
    contrast.a <- c2[[i]]
    DesignEset.i <- getOneDesignEset(DesignEset,contrast.a)
    names.i <- paste0(contrast.a,collapse ="-")

    ## 赋值
    condition <- DesignEset.i[["DesignList"]][["condition"]]
    library(Biobase)
    expr1 <- exprs(DesignEset.i[["ESet"]])

    ## 计算个数
    n.condition <- length(unique(condition$condition))
    n.array <- table(condition$condition)
    if(sum(n.array) > limit.array){
      print(str_c("按limit.array=",limit.array,"的标准,属于巨大矩阵。随机抽样",limit.genes,"个基因进行分析。"))
      set.seed(seed);select.row <- sample(1:nrow(expr1),limit.genes,replace=F)
      expr <- expr1[select.row,]
      x.names=scale_x_discrete(breaks=NULL)#不显示下标
    } else {
      expr <- expr1
      x.names=scale_x_discrete()#默认显示下标
    }
    n.expr <- dim(expr)[1]

    ## data构建
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

    ## 画箱式图
    # c <- c(4,1,6,10,2,3,5,7,8,9,11,12)
    #if(length(unique(n.array1))<=12){
    #  palette <- brewer.pal(12,"Set3")[c]
    #} else {palette <- NULL}
    library(ggpubr)
    p.i <- ggboxplot(data=df1,
                   x="array",y="expr",
                   color = "black",
                   fill = "condition",
                   palette = "nature") +
      labs(x = "Arrays",y = 'The expression level of genes') +
      #ggtitle(str_c("BoxPlot:",paste0(contrast,collapse = "-"))) +
      theme(plot.title = element_text(face = "bold",size = 15,hjust = 0.5),axis.text.x = element_text(size = 12,colour = "black",face = "bold"),axis.title.x = element_text(size = 15,colour = "black",face = "bold"),axis.title.y = element_text(size = 15,colour = "black",face = "bold"),axis.text = element_text(face = "bold",size = 12),
            axis.ticks.x = element_blank(),
            legend.position='top',
            legend.title=element_blank(),
            legend.text =element_text(face = "bold",size = 12),
            panel.border=element_rect(fill='transparent', color='black')) +
      rotate_x_text(angle = 45) +
      x.names

    # 新窗口预览图片
    win.graph(width = 9,height = 7);print(p.i)
    p1 <- c(p1,list(p.i))

    }

  ##保存文件
  if(savefile==T){
    pdf(str_c(names,"_QualityBoxPlot.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_QualityBoxPlot.pdf已保存在本地"))
    return(p1)
  } else {
    return(p1)
  }
}

## example
# QC.boxplot(DesignEset,
#            contrast,
#            savefile=F,
#            names="test1",
#            limit.array = 30,
#            limit.genes=3000,
#            seed=2018)



