

###=====compareROC2():用于普通二分类变量的多条ROC曲线绘制；支持glm的多marker融合变量的ROC计算。
#' @title get ROC Curve of one or more variates for a binary status
#' @description compareROC2 helps get ROC curves of one or more variates for a binary status.It supports the merge of lots of variates via glm function.
#' @param data a data frame
#' @param markers  target markers
#' @param status the dependent variable in ROC.like "N.status"
#' @param merge.markers  A list of markers you want to merge.Default is NULL.The merge stragegy is based on Generalized Linear Models(glm)
#' @param roc.type one of "ggplot" and "pROC"."ggplot" is recommanded and the default setting
#' @param title  the plot title
#' @param color Default is NULL.You can set other colors like "#8DD3C7"
#' @param half.border whether to show half border style
#' @param reference whether to show a reference line in ROC plot
#' @param legend.position the position of legend in ggplot
#' @param show.auc  whether to show related auc in legend labels
#' @param auc.digits the digits of auc in legend labels
#' @param width,height the size of saved PDF plot
#' @param names part of saved PDF plot names
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and NOT RUN
#'train.LM.ROC <- compareROC2(
#'  data = design.x,
#'  markers = mgenes.s,
#'  status ="lymphatic.metastasis",
#'  merge.markers=list(merge = mgenes.s),
#'  roc.type = "ggplot",
#'  title="ROC Curves of Lymphatic Metastasis",
#'  color = mycolor[c(1,20,3,5,6,7,4)],
#'  legend.position=c(0.8,0.2),
#'  show.auc = T,auc.digits = 2,
#'  width = 10,height = 10,
#'  names =paste0("seed",seed,"_train.mgenes")
#')
#'## Quick Start
#'data(fat)
#'markers=c("BMI","Waist","WHR","belly fat thickness") #marker #'colnames
#'status =  "outcomes" #survival status
#'merge.markers = list(c("BMI","Waist"))# the markers list that you want to merge as a co-prognostic factor
#'return.plot = T # if T，return plot;if F，return data
#'roc.type = c("ggplot","pROC")[1] # ggplot style is recommanded
#'title="this is a title" # plot title
#'color = brewer.pal(12, "Set3")[1:5] #curve colors.The lenght of colors must >= the number of curves.
#'output.name = "ROC of something" # part of PDF name
#'compareROC2(data=fat,
#'            markers,
#'            status,
#'            merge.markers=NULL,
#'            return.plot=T,
#'            roc.type,
#'            title,
#'            color,
#'            output.name="test1")
#' @export
compareROC2 <- function(data,
                        markers,
                        status,
                        merge.markers=NULL,
                        roc.type = "ggplot",
                        title="ROC Curves",
                        color=NULL,
                        half.border = F,reference = T,
                        legend.position=c(0.8,0.2),
                        show.auc = T,auc.digits = 2,
                        width = 10,height = 10,
                        names ="love"){
  ##加载必要的包
  need <- c("pROC","ggplot2","RColorBrewer","stringr")
  Plus.library(need)

  ##df的初步处理
  df <- as.data.frame(data)
  df1 <- df[!is.na(df[,match(status,colnames(df))]),]


  ##单marker计算每个marker的roc对象
  if(!is.null(markers)){
    roc.single <- list() #i=1
    for(i in 1:length(markers)){
      outcome.i <- df1[,match(status,colnames(df1))]
      marker.i <- df1[,Fastmatch(markers[i],colnames(df1))]
      roc.i <- roc(outcome.i,marker.i)
      roc.single <- c(roc.single,list(roc.i))
      names(roc.single)[i] <- markers[i]
    }
  } else {
    #没有markers值输入
    roc.single <- NULL
  }

  ##merge.marker，计算roc对象 j=1
  if(is.null(merge.markers)){roc.models <- roc.single} else{
    #有联合marker的需求
    multi.models <- list()
    m.names <- names(merge.markers)
    for(j in 1:length(merge.markers)){
      data1 <- df1[,c(status,merge.markers[[j]])]
      colnames(data1)[1] <- "status"
      logic1 <- apply(data1,1,is.one.na)
      data1 <- data1[!logic1,]
      model.j <- glm(status ~ ., data=data1, family='binomial')
      pre.j <- predict(model.j)
      roc.j <- roc(data1$status,pre.j)
      multi.models <- c(multi.models,list(roc.j))
      if(is.null(m.names[j])){
        #未指定merge.marker的名字
        names(multi.models)[j] <- paste0(merge.markers[[j]],collapse = "_")
      } else {
        names(multi.models)[j] <- m.names[j]
      }

    }

    ###合并单/多变量的model
    roc.models <- c(roc.single,multi.models)
  }

  ## 颜色
  if(is.null(color)){
    c <- rep(mycolor,100)
    color1 <- c[1:length(roc.models)]
  } else {
    color1 <- color
  }

  ## 半边框
  if(half.border == T){
    panel.border=element_rect(fill='transparent',color='transparent')
    axis.line = element_line(colour = "black")
  } else {
    panel.border=element_rect()
    axis.line = element_line()
  }

  ## 是否添加斜线
  if(reference == T){
    r <- geom_abline(intercept=0,slope=1,linetype = 2) #斜率为1，截距为0的直线
  } else {
    r <- NULL
  }

  ## 是否在legend上添加AUC相关信息
  if(show.auc == T){
    # 展示auc曲线面积
    auc <- NULL
    for(i in 1:length(roc.models)){ # i=1
      auc.i <- roc.models[[i]]
      auc.i <- as.numeric(auc.i$auc)
      auc.i <- round(auc.i,auc.digits)
      auc <- c(auc,auc.i)
    }
    auc.labels <- str_c(names(roc.models),"(AUC=",auc,")")
    gg.auc <- scale_color_manual(
      values = color1,
      aesthetics = "colour",
      breaks = names(roc.models),
      labels = auc.labels)
  } else {
    gg.auc <- scale_color_manual(
      values = color1,
      aesthetics = "colour")
  }

  ##绘制制ROC图
  if(roc.type == "ggplot"){
    #ggplot风格    roc.list = roc.models; i=1
    extraroc <- function(roc.list){
      df <- NULL
      for(i in 1:length(roc.list)){
        sensitivity <- as.numeric(as.character(roc.list[[i]][["sensitivities"]]))
        specificity <- as.numeric(as.character(roc.list[[i]][["specificities"]]))
        specificity1 <- 1-specificity
        Groups <- names(roc.list)[i]
        df.i <- cbind(specificity1,sensitivity,Groups,specificity)
        # data frame
        df.i <- as.data.frame(df.i)
        df.i$specificity1 <- as.numeric(as.character(df.i$specificity1))
        df.i$sensitivity <- as.numeric(as.character(df.i$sensitivity))
        df.i$specificity <- as.numeric(as.character(df.i$specificity))
        df.i <- df.i[order(df.i$sensitivity,df.i$specificity),]
        df <- rbind(df,df.i)
      }
      df <- as.data.frame(df)
      df$Groups <- factor(df$Groups,levels = names(roc.list))
      return(df)
    }
    df2 <- extraroc(roc.models);#str(df2)
    g <- ggplot(df2, aes(x=specificity1,y=sensitivity,color = Groups)) +
      scale_x_continuous(limits=c(0,1), breaks=seq(1,0,-0.25)) +
      geom_line(size = 1.2) +
      ggtitle(title) +
      r + gg.auc +
      labs(x = "1 - Specificity",y = "Sensitivity") +
      theme_bw() +
      theme(panel.grid =element_blank(),
            plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
            axis.title = element_text(face = "bold",size = 15),
            axis.text = element_text(face = "bold",size = 12),
            legend.title = element_text(face = "bold",size = 12),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = legend.position,
            panel.border=panel.border,
            axis.line = axis.line);
    print(g)
    ggsave(str_c(names,"_compareROC_gg.pdf"),g,width = 10,height = 8)
    print(str_c("完成",str_c(names,"_compareROC_gg.pdf"),"的绘制。"))
  } else {
    if(roc.type == "pROC"){
      #pROC的风格~研究不深哈，就随便画一个了~
      g <- NULL
      pdf(str_c(output.name,"_compareROC.pdf"),width = width,height = width)
      plot(roc.models[[1]], print.auc=TRUE,
           grid=c(0.1, 0.2),
           grid.col=c("green", "red"), max.auc.polygon=TRUE,
           auc.polygon=F,
           #auc.polygon.col="skyblue",
           print.thres=TRUE)
      if(length(roc.models) >= 2){
        for(j in 2:length(roc.models)){
          plot.roc(roc.models[[j]], add=TRUE, col=color[j])#在上图中继续添加ROC曲线
        }
      }
      dev.off()
      print(str_c("完成",str_c(names,"_compareROC.pdf"),"的绘制。"))
    }
  }

  ##输出数据
  n = 1:length(roc.models);
  if(length(n)==1){
    print("仅单个变量，无法比较");
    d <- list(
      roc.plot = g,
      roc.models = roc.models
    )
  } else {
    #可进行两两比较
    n.mt <- combn(n,2)
    x <- NULL;
    for(j in 1:ncol(n.mt)){
      s1 =n.mt[1,j];s2 = n.mt[2,j]
      x.j <- roc.test(roc.models[[s1]],roc.models[[s2]],method = "delong")#ROC曲线下面积的比较
      x <- c(x,list(x.j));
      names(x)[j] <- paste(names(roc.models)[s1],names(roc.models)[s2],sep = " vs ")
    }
    d <- NULL
    for(j in 1:length(x)){
      di <- data.frame(contrast = names(x)[j],
                       Z_val = x[[j]]$statistic,
                       p_val = x[[j]]$p.value)
      d <- rbind(d,di)
    }
    d <- list(compare = d,
              roc.plot = g,
              roc.models = roc.models)

  }
  return(d)

  #End
}






