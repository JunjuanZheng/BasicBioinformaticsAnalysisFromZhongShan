

###======QC.plotMA()
# contrast=list(c("N1.plus","N1"),c("N1.plus","N0"),c("N1","N0")) # a list that contain lots of contrasts.
## QC.plotMA()对于给定的DesignEset和contrast可以进行MA图的绘制。
#' @export
QC.plotMA <- function(DesignEset,
                      contrast=NULL,
                      savefile=F,
                      names = "test1"){
  ## 加载包
  library(stringr);
  library(limma);library(Biobase)
  library(ggplot2);library(RColorBrewer)

  ## 赋值
  design <- DesignEset[["DesignList"]][["design"]]
  eset <- exprs(DesignEset[["ESet"]])
  anotation0 <- fData(DesignEset[["ESet"]])
  contrast.matrix <- DesignEset[["DesignList"]][["contrast.matrix"]]

  ## 计算logFC、mean值和P值
  fit <- lmFit(eset, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  results <- decideTests(fit);
  print(summary(results)) #查看up和down的信息
  results.matrix <- results@.Data #包含探针名与表达差异的对应关系矩阵

  ## contrast的设置
  if(is.null(contrast)){
    #提取DesignEset中的contrast信息
    c1 <- colnames(contrast.matrix)
  } else {
    c1 <- NULL
    for(i in 1:length(contrast)){
      contrast.a <- contrast[[i]]
      contrast.a1 <- c(paste(contrast.a[1],contrast.a[2],sep = "-"),paste(contrast.a[2],contrast.a[1],sep = "-"))
      l1 <- contrast.a1 %in% colnames(contrast.matrix)
      contrast.a1 <- contrast.a1[l1]
      c1 <- c(c1,contrast.a1)
    }
    c1
  }

  p1 <- NULL
  for(i in 1:length(c1)){
    contrast.i=colnames(contrast.matrix)[i]
    sig.exprs <- topTable(fit, coef=contrast.i, n=Inf) #输出
    select.cols <- c("logFC","AveExpr","t","P.Value","adj.P.Val")
    sig.exprs1 <- subset(sig.exprs,select = select.cols)
    sig.exprs1$sig <- ifelse(sig.exprs1$adj.P.Val < 0.1 & abs(sig.exprs1$logFC) > 1,"sig","nonsig");table(sig.exprs1$sig)
    print(paste0(contrast.i,":完成logFC与表达水平的计算。",collapse = ""))

    ### 对探针颜色进行标注
    control.probes = DesignEset[["DesignList"]][["control.probes"]]
    ## 标记显著性表达的基因
    all.position <- 1:nrow(sig.exprs1) #探针位置
    sig.position <- ifelse(sig.exprs1$sig == "sig",T,F)
    sig.position <- all.position[sig.position]
    #results.matrix.i <- results.matrix[,contrast.i] #带探针名的向量
    #sig.position <- ifelse(abs(results.matrix.i) == 1,T,F)
    ## 分有/无control探针两种情况
    if(!is.null(control.probes)){
      ## 按要求特殊标记control.probes
      con.position <- NULL
      for(i in 1:length(control.probes)){
        con.i <- control.probes[[i]]
        #找到control.probes的位置
        con.position.i <- NULL
        # s="control"
        for(s in con.i$control.symbol){
          con.position.s <- grep(s,anotation0[,con.i$control.col])
          con.position.i <- c(con.position.i,con.position.s)
        }
        con.position <- c(con.position,con.position.i)
      }
      con.position <- unique(con.position) #control在全探针里的位置
      ## 设置颜色向量:
      probe.type = ifelse(all.position %in% con.position,"control",ifelse(all.position %in% sig.position,"sig","non-sig"))#table(probe.type)
    } else {
      ##没有control.probes
      probe.type = ifelse(all.position %in% sig.position,"sig","non-sig")#table(probe.type)
    }
    ## control探针和显著探针标亮色，非显著基因为浅色。
    if(!is.null(control.probes)){
      c <- c("control"="yellow","sig"="#FB8072","non-sig"="#BEBADA")
    } else {
      c <- c("sig"="#FB8072","non-sig"="#BEBADA")
    }
    data.MA <- cbind(probe.type,sig.exprs1)
    data.MA$Significant <- 1/data.MA$adj.P.Val
    print(paste0(contrast.i,":完成探针的颜色标注。"))

    ## 绘制MA图：y=logFC,x=AveExpr。
    p <- ggplot(data.MA,aes(x=AveExpr,y=logFC,colour=probe.type,size = Significant)) +
      geom_point(alpha = 1) +
      scale_colour_manual(values = c,aesthetics = "colour") +
      xlab("Average Expression") + ylab("log Fold Change") + ggtitle(str_c("MAplot:",contrast.i))+
      theme_bw() +
      theme(panel.grid =element_blank(),
            plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
            axis.title = element_text(face = "bold",size = 12),
            axis.text = element_text(face = "bold",size = 10),
            legend.title = element_text(face = "bold",size = 12),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = "right")

    print(paste0(contrast.i,":完成MA图绘制。",collapse = ""))
    ##预览图像
    win.graph(width=9, height=8);print(p)
    p1 <- c(p1,list(p))
  }

  ##保存文件
  if(savefile==T){
    pdf(str_c(names,"_QualityMAplot.pdf"),width=9,height=8)
    for(i in 1:length(p1)){
      p.i <- p1[[i]]
      print(p.i)
    }
    dev.off()
    print(str_c(names,"_QualityMAplot.pdf已保存在本地"))
    #End
    return(p1)
  }  else {
    #End
    return(p1)
  }
}

## example:
# DesignEset <- AddDesignList(eset,factor,levels,control.probes)
# contrast = c("Vemurafenib","Control")#对照放后面
# QC.plotMA(DesignEset,contrast,savefile=F)

# DesignEset # a DesignEset
# contrast=c("Vemurafenib","Control") #对比
# savefile=F #是否将图片保存为pdf
# names="test1" #部分文件名



