


###======QC.PCA
# QC.PCA用于绘制PCA图以观察不同组别的区分情况
#' @export
QC.PCA <- function(DesignEset,
                   contrast,
                   savefile=F,
                   names="test1"){
  ## 加载包
  library(ggplot2)

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

  p1 <- NULL
  for(i in 1:length(c2)){
    ## 提取某个DesignList和contrast
    contrast.a <- c2[[i]]
    DesignEset.i <- getOneDesignEset(DesignEset,contrast.a)
    names.i <- paste0(contrast.a,collapse = "-")

    ## 赋值
    condition <- DesignEset.i[["DesignList"]][["condition"]]
    library(Biobase)
    mt <- exprs(DesignEset.i[["ESet"]])
    mt1 <- as.matrix(mt)

    ##PCA分析
    print("正在进行主成分分析...")
    library(stats)
    pca.mt <- prcomp(mt1)
    print("完成主成分分析!")
    pcadata <- pca.mt[["rotation"]]
    percentVar <- pca.mt$sdev^2/sum(pca.mt$sdev^2)
    pcadata <- pcadata[,1:2]

    ##合并pca数据与分组信息
    pcadata1 <- rbind(t(pcadata),t(condition))
    pcadata1 <- as.data.frame(t(pcadata1))
    pcadata1$PC1 <- as.numeric(as.character(pcadata1$PC1))
    pcadata1$PC2 <- as.numeric(as.character(pcadata1$PC2))

    ##PCA作图
    library(ggplot2)
    p <- ggplot(pcadata1, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", round(percentVar[1]*100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) +
      ggtitle(str_c("PCA:",paste0(contrast.a,collapse = "-"))) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
            axis.title = element_text(face = "bold",size = 12),
            axis.text = element_text(face = "bold",size = 10),
            legend.title=element_blank(),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = "right")
    p1 <- c(p1,list(p))

    ##预览图片
    win.graph(width = 9,height = 7);print(p)
  }

  ## 保存文件
  if(savefile==T){
    pdf(str_c(names,"_Quality.PCAPlot.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_Quality.PCAPlot.pdf已保存在本地"))
    #End
    return(p1)
  } else {
    #End
    return(p1)
  }


}



