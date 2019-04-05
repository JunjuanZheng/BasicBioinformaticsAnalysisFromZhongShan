

#' @title Fast way to deal with results of qRT-PCR
#' @description Fast way to deal with results of qRT-PCR,including compared plot,statistics and quality control.
#' @param data the result of qRT-PCR.
#' @param sample the colnames of sample.Note that different samples should use different sample id.
#' @param marker the colnames of marker
#' @param bioRepeat the colnames of bioRepeat.The data in "bioRepeat" would be considered as available repeat in afterward statistics.
#' @param parallelRepeat the colnames of parallelRepeat.The data in "parallelRepeat" would be treated via mean strategy.
#' @param group the colnames of group.Like "treatment" and "control".
#' @param group.control the names of control group.Default is "control".
#' @param internal the name of internal reference gene.Default is "GAPDH".
#' @param value the colnames of Ct value
#' @param plot.type one of 1 and 2.1 is used in multiple markers,but 2 is more suitable for one marker.
#' @param palette the color of groups.The number of color must be equal to the number of groups.Default is NULL,which mean that the strategy of lucky package would be use.
#' @param size ggplot parameters.the size of plot
#' @param label.position ggplot parameters.the position of p significance
#' @param label.type character string specifying label type. Allowed values include "p.signif" (shows the significance levels), "p.format" (shows the formatted p value).
#' @param method method use in statistics."t.test" is always recommended.You can also use "wilcox.test" for a try.Alternative choices are "anova" and "kruskal.test".
#' @param x.title ggplot parameters.the x axis title
#' @param y.title ggplot parameters.the y axis title
#' @param legend.position ggplot parameters.legend position
#' @importFrom plyr ddply summarise
#' @importFrom ggpubr ggbarplot stat_compare_means
#' @importFrom ggplot2 labs theme
#' @return a Lucky Objects
#' @seealso \code{\link[ggpubr]{stat_compare_means}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' data("qpcr",package = "lucky") # see the example of qpcr result
#' result.pcr <- FastQPCR(qpcr)
#' result.pcr <- FastQPCR(data=qpcr,
#'                        sample = "samples",
#'                        marker = "markers",
#'                        bioRepeat = "biorepeat",
#'                        parallelRepeat = "parepeat",
#'                        group = "groups",
#'                        group.control = "control",
#'                        internal = "GAPDH",
#'                        value = "Ct",
#'                        size=20,
#'                        label.y = 650,
#'                        method = "t.test",
#'                        x.title = "Genes",
#'                        y.title = "The Relative Expression of Genes",
#'                        legend.position = "top")
#' ## Use t.test.However,it's not recommand.
#' result.pcr <- FastQPCR(data=qpcr,
#'                        sample = "samples",
#'                        marker = "markers",
#'                        bioRepeat = "biorepeat",
#'                        parallelRepeat = "parepeat",
#'                        group = "groups",
#'                        group.control = "control",
#'                        internal = "GAPDH",
#'                        value = "Ct",
#'                        size=20,
#'                        label.y = 650,
#'                        method = "t.test", #t.test
#'                        x.title = "Genes",
#'                        y.title = "The Relative Expression of Genes",
#'                        legend.position = "top")
#'
#' ## View result
#' View(result_qpcr$Data$statistc$whole)
#' View(result_qpcr$Data$statistc$pair)
#'
#' ## test
#' a = result.pcr$Data$metadata
#' x <- a$fc[a$markers %in% "marker1" & a$groups %in% "control"]
#' y <- a$fc[a$markers %in% "marker1" & a$groups %in% "treat"]
#' t.test(x,y)
#' @export
FastQPCR <- function(data,
                     sample = "samples",
                     marker = "markers",
                     bioRepeat = "biorepeat",
                     parallelRepeat = "parepeat",
                     group = "groups",
                     group.control = "control",
                     internal = "GAPDH",
                     value = "Ct",
                     plot.type = c(1,2)[1],
                     palette = NULL,
                     size=20,
                     label.position = c(4,5),
                     label.type = c("p.signif","p.format")[1],
                     method = "t.test",
                     x.title = "Genes",
                     y.title = "The Relative Expression of Genes",
                     legend.position = "top"){
  ## package
  nd <- c("plyr","ggplot2","ggpubr")
  Plus.library(nd)

  ###======================质量控制====================###
  data <- as.data.frame(data)
  data2 <- data
  colnames(data2)[Fastmatch(c(group,marker,value),colnames(data2))] <- c("groups","markers","Ct")
  data2$Ct <- - data2$Ct
  p.con <- ggplot(data2,aes(x = groups,y = Ct,fill = markers)) +
    geom_flat_violin(position=position_nudge(x=0.08)) +
    #geom_jitter(aes(color=markers), width=.15) +
    geom_dotplot(binaxis="y", stackdir="down", dotsize=.35,position=position_nudge(x=-0.05),binwidth = 0.35) +
    geom_boxplot(width=.1) +
    coord_flip() +
    theme_bw() +
    ggtitle("qPCR Quality Control") +
    labs(x = "Groups",y = "Reverse Ct Value") +
    theme(
      plot.title = element_text(face = "bold",size = 18*size/20,hjust = 0.5),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold",size = 15*size/20),
      axis.text = element_text(face = "bold",size = 12*size/20),
      legend.title = element_text(face = "bold",size = 12*size/20),
      legend.text =element_text(face = "bold",size = 12*size/20),
      legend.position = "bottom"
    )
 win.graph(10,10);print(p.con)

  ###======================计算有效数据====================###
  {
  ## 计算某个marker中， 某group 中，某sample中，某biorepeat中，parallel repeat的平均值
  colnames(data)[match(value,colnames(data))] <- "Ct"
  data.exp <- ddply(.data = data,
                    .variables = c(marker,group,sample,bioRepeat),
                     summarise,
                     Ct.mean = mean(Ct))
  data.internal <- data.exp[as.character(data.exp[,marker]) %in% internal,]
  data.exp2 <-data.exp[!as.character(data.exp[,marker]) %in% internal,]

  ## 对于某个sample，计算各marker与内参之间的Ct差值
  get1 <- function(x){
    # x <- data.exp2[1,]
    x <- as.character(as.matrix(x))
    sample.i <- as.character(x[3]) #sample
    Ct.mean.marker.i <- as.numeric(as.character(x[5]))#Ct.mean
    data.internal.i <- data.internal[as.character(data.internal[,sample]) %in% sample.i,] #选出某sample的内参数据
    Ct.mean.internal.i <- mean(as.numeric(as.character(data.internal.i[,"Ct.mean"]))) #某个生物学重复，其内参求均值
    dif <- Ct.mean.marker.i - Ct.mean.internal.i
    return(dif)
  }
  data.exp3 <- ddply(.data = data.exp2,
                     .variables = c(marker,group,sample,bioRepeat),
                     .fun = get1)
  colnames(data.exp3)[ncol(data.exp3)] <- "dif" #各marker与内参之间的Ct差值

  ## 计算其它group与group control的差值
  # 先计算control组中，dif的平均值
  data.exp3.control <- data.exp3[data.exp3[,group] %in% group.control,]
  data.exp3.control.Ctmean <- ddply(
    .data = data.exp3.control,
    .variables = marker,
    summarise,
    dif.mean = mean(dif))
  # 再计算2的-△t次方
  get2 <- function(x){
    # x <- data.exp3[as.character(data.exp3$markers) == "RP3",]
    ## control组的平均dif值
    marker.x <- as.character(unique(x[,marker]))
    dif.x <- data.exp3.control.Ctmean[match(marker.x,data.exp3.control.Ctmean[,marker]),"dif.mean"]
    # 计算2的-△t次方
    x$fc <- 2^(-(x[,"dif"] - dif.x))
    return(x)
  }
  newData <- ddply(.data = data.exp3,
                   .variables = marker,
                   .fun = get2)
  colnames(newData)[match(group,colnames(newData))] <- "groups"
  }
  # View(newData)

  ###=======================绘制统计图=====================###

  ## 根据组数选择颜色
  if(is.null(palette)){
    palette <- mycolor[1:length(as.character(unique(newData$groups)))]
  }

  ## 绘图
  if(plot.type == 1){
    ## 此时有多个marker(除内参外)
    p <- ggbarplot(newData, x = marker,y = "fc",
                   color = "groups",
                   fill = "groups",
                   palette = palette,
                   label = F,
                   position = position_dodge(0.8),
                   add = "mean_se",
                   error.plot = "upper_errorbar") +
      labs(x = x.title,y = y.title) +
      theme(axis.title = element_text(face = "bold",size = 15*size/20),
            axis.text = element_text(face = "bold",size = 12*size/20),
            legend.title = element_text(face = "bold",size = 12*size/20),
            legend.text =element_text(face = "bold",size = 12*size/20),
            legend.position = legend.position) +
      stat_compare_means(aes(group=groups),
                         method = method,
                         label = label.type,
                         paired = F,
                         label.x =  label.position[1],
                         label.y = label.position[2])#统计学差异
  } else {
    ## 此时只有多个marker(除内参外)
    #指标只有一个时，分组应该做横坐标
    p <- ggbarplot(newData,
                   x = "groups",
                   y = "fc",
                   fill = "groups",
                   palette = palette,
                   label = F,
                   position = position_dodge(0.8),
                   add = "mean_se",
                   error.plot = "upper_errorbar") +
      labs(x = x.title,y = y.title) +
      theme(axis.title = element_text(face = "bold",size = 15*size/20),
            axis.text = element_text(face = "bold",size = 12*size/20),
            legend.title = element_text(face = "bold",size = 12*size/20),
            legend.text =element_text(face = "bold",size = 12*size/20),
            legend.position = legend.position) +
      stat_compare_means(aes(group=groups),
                         method = method,
                         label = label.type,
                         paired = F,
                         label.x =  label.position[1],
                         label.y =  label.position[2])#统计学差异
  }

  win.graph(10,10);print(p)

  ###=======================汇总统计结果=====================###
  uniquemarker <- as.character(unique(newData[,marker]))
  uniquegroup <- as.character(unique(newData[,"groups"]))
  s <- combn(length(uniquegroup),2)

  ## 总体比较
  test <- NULL
  for(i in uniquemarker){ #i=uniquemarker[1]
    test.i <- compare_means(fc ~ groups, newData[newData[,marker] == i,], method = method)
    test.i <- as.data.frame(test.i,stringAsFactor = F)
    rownames(test.i) <- i
    test <- rbind(test,test.i)
  }
  test <- test[-2]

  ## 多个两两比较
  test2 <- NULL
  for(m in uniquemarker){ # m = uniquemarker[1]
    test2.m <- NULL
    for(i in 1:ncol(s)){ # i = 1
      d <- newData[newData[,marker] == m & newData[,"groups"] %in% uniquegroup[s[,i]],]
      test.i <- compare_means(fc ~ groups,data = d, method = "t.test")
      test.i <- as.data.frame(test.i,stringAsFactor = F)
      rownames(test.i) <- i
      test.i <- cbind(marker = m,test.i)
      test2.m <- rbind(test2.m,test.i)
    }
    if(ncol(s) >= 3){
      #校正两两比较
      test2.m$p.adj <- p.adjust(test2.m$p,method = "BH")
      test2.m$p.signif <- cut(
        test2.m$p.adj,
        breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        labels = c("****", "***", "**", "*", "ns"),
        right = F
      )
    }
    test2.m <- test2.m[-match("p.format",colnames(test2.m))]
    test2 <- rbind(test2,test2.m)
  }
  test2 <- test2[-match(".y.",colnames(test2.m))]

  ###=======================输出结果=====================###
  colnames(newData)[Fastmatch(c("dif","fc"),colnames(newData))] <- c("△Ct","2^-△△Ct")
  l <- list(
    Repeat = list(
      sample =sample,
      marker = marker,
      group = group,
      value = value,
      size = size,
      label.position = label.position,
      method = method
    ),
    Data = list(
      metadata = newData,
      statistc = list(whole = test,pair = test2)
    ),
    Plot = list(
      QCplot = p.con,
      compareplot = p
    )
  )
  return(l)
}


