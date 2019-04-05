

#' @title Fast way to draw a KM plot of ggplot style
#' @description Fast way to draw a KM plot of ggplot style
#' @param data a data frame
#' @param time the colname of time
#' @param status the colname of status
#' @param marker the marker you want to explore with KM plot
#' @param color If \code{color=NULL},it would use the color from \code{mycolor} in lucky package.You can also set self-defined colors
#' @param size the size of plot
#' @param legend the position of legend
#' @param legend.title the title of legend
#' @param title the title of the KM plot
#' @param pval whether to output P value in KM plot
#' @param pval.position the position of P value in KM plot
#' @param saveplot whether to save the KM plot
#' @param name part name of saved plot
#' @param conf.int whether to output confidence interval
#' @param risk.table whether to output risk table
#' @param linetype line type
#' @param surv.median.line  whether to show median line
#' @param ncensor.plot whether to show ncensor plot
#' @seealso \code{\link[survminer]{ggsurvplot}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(survival)
#' data("lung");
#' data = lung
#' data$status <- data$status -1
#' p <- FastSurvplot(data = data,
#'                   time = "time",
#'                   status = "status",
#'                   marker = "sex",
#'                   size = 10,
#'                   legend.title = "Sex",
#'                   pval.position = c(900,1))
#' p
#' @export
FastSurvplot <- function(data,
                         time = "OS.time",
                         status = "OS.status",
                         marker = "pT",
                         color = NULL,
                         size = 16,
                         legend = "top",
                         legend.title = "T status",
                         title = NULL,
                         pval = TRUE,
                         pval.position = c(1,1),
                         conf.int = F,
                         risk.table = F ,
                         linetype = "solid",
                         surv.median.line = "none" ,
                         ncensor.plot=F,
                         saveplot = F ,
                         name = "GC"){
  #加载包
  nd <- c("survival","ggpubr","stringr","survminer","RColorBrewer","ggplot2");Plus.library(nd)

  ## color
  if(is.null(color)){
    color = mycolor[c(1,4,6,10)]
  } else {
    color = color
  }

  ##自定义列名
  data1.FastSurvplot <- data[c(time,status,marker)]## data1.FastSurvplot一定要为全局变量
  colnames(data1.FastSurvplot)[1:3] <- c("time","status","value")
  data1.FastSurvplot$time <- as.numeric(as.character(data1.FastSurvplot$time))
  data1.FastSurvplot$status <- as.numeric(as.character(data1.FastSurvplot$status))

  ## p值
  if(pval == T){
    sdf <- survdiff(Surv(time,status)~value, data = data1.FastSurvplot)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    p1 <- round.plus(p.val,digits = 4)
    label <- ifelse(is.character(p1),p1,paste0("= ",p1))
    p.p <- annotate("text",
                    x = pval.position[1],y = pval.position[2],
                    label = paste0("p ",label),
                    size = (size/16)*8)
  } else {
    p.p <- NULL
  }

  ##画图
  fit <- survfit(Surv(time,status) ~ value, data=data1.FastSurvplot)


  ## theme
  theme1 <- theme_bw() + theme(
    panel.grid =element_blank(),
    plot.title = element_text(face = "bold",size = (size/16)*25,hjust = 0.5),
    axis.title = element_text(face = "bold",size = (size/16)*25),
    axis.text = element_text(face = "bold",size = (size/16)*20),
    legend.title = element_text(face = "bold",size = (size/16)*20),
    legend.text = element_text(face = "bold",size = (size/16)*20),
    panel.border=element_rect(fill='transparent',color='transparent'),
    axis.line = element_line(colour = "black")
  )

  ## ggsurvplot
  suvival <- ggsurvplot(fit,
                        data = data1.FastSurvplot,
                        pval = F,
                        conf.int = conf.int,
                        risk.table = risk.table,
                        linetype = linetype,
                        surv.median.line = surv.median.line,
                        palette = color,
                        ncensor.plot=ncensor.plot,
                        legend = legend,
                        legend.title = legend.title,
                        legend.labs = levels(data1.FastSurvplot$value),ggtheme = theme1);
  p2 <- suvival$plot %+% labs(title = title) + p.p

  ##展示或保存图片
  if(saveplot == F){
    return(p2)
  } else {
    #保存pdf文件
    ggsave(str_c("KMcurves_",marker,"_",name,".pdf"),p2,height = 12,width = 12)
    print(str_c("完成",marker,"_KMcurves","图的绘制"))
    return(p2)
  }

  ##End

}



