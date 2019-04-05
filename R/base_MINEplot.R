

###=======MINEplot
# Usage # MINE.corrplot() help to draw some plot base on the result from FastMINE.
#' @export
MINEplot <- function(data,
                     transposition = F,
                     select,
                     order = T,#是否对相关矩阵排序
                     QC=T,outer.ratio = 0.1,
                     lower.col = NULL,
                     upper.col =NULL,
                     upper = NULL,
                     tl.pos = NULL,
                     tl.col=NULL,
                     tl.srt=NULL,
                     diag = NULL,
                     savefile=T,
                     names="love"){
  ## 加载相应的包
  library(corrplot);library(rJava)

  ## 矩阵设置
  expr1 <- as.matrix(data)
  if(transposition==T){
    # 矩阵要转置
    expr2 <- t(expr1)
  } else {
    expr2 <- expr1
  }
  expr2 <- expr2[,select]

  ## 计算矩阵MIC系数
  mycor <- FastMINE(data=expr2,
                    transposition = F,
                    method = "all.pairs",
                    control.markers=NULL,
                    target.markers=NULL)
  mycor1 <- mycor[["MIC.matirx"]]

  ## 计算矩阵相关系数并决定是否自动排序
  if(order == F){
    mycor2 <- mycor1
  } else {
    ord <- corrMatOrder(mycor1, order = "AOE")
    mycor2 <- mycor1[ord,ord]
  }
  mycor[["MIC.matirx"]] <- mycor2

  ## 加载默认或者自定义信息
  library(grDevices)
  #library(RColorBrewer);display.brewer.pal(12, "Set3")
  # "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072"
  # "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5"
  # "#D9D9D9" "#BC80BD" "#CCEBC5" "#FFED6F"
  col <- colorRampPalette(c("#D9D9D9","#80B1D3","#FB8072"))(100)
  default <- c(rep(list(NULL),7)) #
  input <- list(lower.col = lower.col,
                upper.col = upper.col,
                upper = upper,
                tl.pos = tl.pos,
                tl.col=tl.col,
                tl.srt=tl.srt,
                diag = diag)
  do <- list(lower.col = col,
             upper.col =col,
             upper = "pie",
             tl.pos = "lt",
             tl.col="black",
             tl.srt=45,
             diag = "l")
  output <- set.default(input,default,do)#base_set.default


  ## 质量检测：观察MINE与其它方法的表现差别；
  if(QC == F){
    print("Quality Control plots would not be printed.")
    out <- "Quality Control plots would not be printed."
  } else {
    ### 进行质检
    print("Here is some Quality Control plots to compare MINE with other correlation methods...")
    ## MIC与Pearson系数的散点图
    df.cor1 <- mycor[["MINE.result"]] #colnames(df.cor1)
    newcol <- c("x1","y1","MIC","mic-p2","MAS","MEV","MCN","Pearson")
    colnames(df.cor1) <- newcol
    library(stringr)
    df.cor2 <- df.cor1
    #构建颜色向量
    color1 <- ifelse(df.cor2$MIC <= outer.ratio &abs(df.cor2$Pearson)>=(1-outer.ratio),NA,1)
    color2 <- ifelse(df.cor2$MIC >= (1-outer.ratio) & abs(df.cor2$Pearson)<=outer.ratio,NA,1)
    co <- cbind(color1,color2)
    color <- apply(co,1,function(x)is.one.na(x))
    color <- ifelse(color==F,"normal","outer")
    l1 <- "outer" %in% color
    if(!l1){
      #没有异常值
      c = c("normal"="#8DD3C7")
    } else {
      #有异常值。高亮
      c = c("normal"="#8DD3C7","outer"="#FB8072")
    }
    df.cor2$status <-  color #colnames(df.cor2)
    out <- subset(df.cor2,status=="outer");
    colnames(out)[1:2] <- c("pair1","pair2")

    ## 绘制散点图：
    library(ggplot2)
    p <- ggplot(df.cor2,aes(x=MIC,y=Pearson,colour=status,size = `mic-p2`)) +
      geom_point(alpha = 1) +
      scale_colour_manual(values = c,aesthetics = "colour") +
      xlab("Average Expression") + ylab("log Fold Change") + ggtitle("MIC-Pearson")+ xlim(0,1) + ylim(-1,1) +
      theme_bw() +
      theme(panel.grid =element_blank(),
            plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
            axis.title = element_text(face = "bold",size = 12),
            axis.text = element_text(face = "bold",size = 10),
            legend.title = element_text(face = "bold",size = 12),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = "right") +
      geom_hline(aes(yintercept=0),#加上直线
                 colour="#80B1D3" ,
                 linetype="dashed");
  if(savefile==F){
    win.graph(width = 10,height = 10);print(p)
    ggsave(paste0(names,"_","MINE.QCplot.pdf"),p,width = 10,height = 10)
    print(paste0(names,"_","PointPlot.pdf已经成功保存!"))
  } else {
    win.graph(width = 10,height = 10);print(p)
  }
}

  ## 绘制 corrplot
  mycor.x <- c(mycor,list(out));names(mycor.x)[3] <- "outer"
  if(savefile==T){
    win.graph(width = 10,height = 10);
    corrplot.mixed(mycor2,
                   cl.lim=c(0,1),#设置颜色的范围
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    file.names <- paste0(names,"_MINE.Corrplot.pdf",collapse = "")
    pdf(file.names,width = 10,height = 10)
    corrplot.mixed(mycor2,
                   cl.lim=c(0,1),
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    dev.off()
    print(paste0(file.names,"已经保存至当前工作目录",collapse = ""))
    return(mycor.x)
  } else {
    win.graph(width = 10,height = 10)
    corrplot.mixed(mycor2,
                   cl.lim=c(0,1),
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    return(mycor.x)
  }
  #End
}



