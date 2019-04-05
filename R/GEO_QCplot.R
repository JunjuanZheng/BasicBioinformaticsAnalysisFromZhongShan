



##======QCplot
# 5种QC图像函数的快速通道
#' @export
QCplot <- function(DesignEset,
                   contrast=NULL,
                   names,
                   cluster_cols=T,
                   log.convert = F,
                   savefile=T){
  ## 加载包
  library(Biobase);library(limma);library(stringr)

  ##预览
  ## MA QC
  ma <- QC.plotMA(DesignEset,
                  contrast,
                  savefile=F,
                  names = names)
  print("完成QC.plotMA")

  ##Heatmap QC:
  hm <- QC.cluster(DesignEset,
                   contrast,
                   savefile=F,
                   names=names,
                   log.convert = log.convert,
                   cluster_cols =cluster_cols)
  print("完成QC.cluster")

  ##Boxplot QC
  bo <- QC.boxplot(DesignEset,
                   contrast,
                   savefile=F,
                   names=names)
  print("完成QC.boxplot")

  ##Density QC:
  ds <- QC.density(DesignEset,
                   contrast,
                   savefile=F,
                   names=names)
  print("完成QC.density")

  ##PCA QC
  pca <- QC.PCA(DesignEset,
                contrast,
                savefile=F,
                names=names)
  print("完成QC.PCA")
  if(savefile==F){
    print("完成QCplots展示。")
  } else {
    ##保存为PDF格式
    p1 <- c(ma,hm,bo,ds,pca)
    pdf(str_c(names,"_QCplot.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_QCPlot.pdf已保存在本地"))
  }

}

# QCplot(DesignEset,
#       contrast,
#       names="test1",
#        savefile=F)




