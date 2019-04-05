
#' @title  time-dependent ROC curve base on timeROC package
#' @description getTimeROC is a powerful function to do time-dependent ROC curve base on timeROC::timeROC.getTimeROC help quickly and conveniently set common parameter in a ROC curve of ggplot2 style.
#' @param data a data frame
#' @param time.col time colnames
#' @param status.col status
#' @param cause the event status.For example,if 1=death,0=alive,then cause should be 1.
#' @param markers The vector of the marker values for which we want to compute the time-dependent ROC curves. Without loss of generality, the function assumes that larger values of the marker are associated with higher risks of events. If lower values of the marker are associated with higher risks of events, then reverse the association adding a minus to the marker values.
#' @param merge.markers  a list of markers groups.Default is NULL.
#' @param weighting.univariate Default is "marginal"
#' @param weighting.multivariate Default is "cox"
#' @param time.raw.type one of "Day","Month","Year"
#' @param time.target.type one of "Day","Month","Year"
#' @param time.knot target time knots.Like 3 years and 5 years survival.
#' @param line.color line colors in ROC plot.
#' @param size the size of plot
#' @param legend.position the position of legend position
#' @param half.border whether show half border of the plot
#' @param width the width of plot
#' @param height the height of plot
#' @param merge.plot whether use a merge style by time.knot
#' @param names part of output file name.
#' @return a list contain model information and series ggplot object.
#' @seealso \code{\link[timeROC]{timeROC}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' library(survival)
#' data(pbc)
#' pbc<-pbc[!is.na(pbc$trt),]
#' pbc$status<-as.numeric(pbc$status==2)
#' data = pbc
#' time.col = "time"
#' status.col ="status"
#' cause = 1
#' markers = c("bili","chol")
#' merge.markers=list(c("bili","chol"))
#' time.raw.type=c("Day","Month","Year")[1]
#' time.target.type=c("Day","Month","Year")[3]
#' time.knot=c(3,5)
#' names = "love"
#' line.color = mycolor[c(4,5,3,8,23)]
#' legend.position = "top"
#' weighting.univariate = "marginal"
#' weighting.multivariate = "cox"
#' size = 1.2
#' half.border=T
#' reference = F
#' width = 10
#' height=10
#' merge.plot = T
#' names = "love"
#'
#' ## merge style
#' models.merge <- getTimeROC(data,
#'                            time.col,status.col,cause,
#'                            markers,
#'                            merge.markers,
#'                            time.raw.type=time.raw.type,
#'                            time.target.type=time.target.type,
#'                            time.knot=time.knot,
#'                            line.color = line.color,
#'                            half.border=half.border,
#'                            reference = reference,
#'                            names = "love")
#'
#' ## single style
#' models.single <- getTimeROC(data,
#'                             time.col,status.col,cause,
#'                             markers,
#'                             merge.markers,
#'                             time.raw.type=time.raw.type,
#'                             time.target.type=time.target.type,
#'                             time.knot=time.knot,
#'                             line.color = line.color,
#'                             half.border=half.border,
#'                             reference = reference,
#'                             merge.plot = F,
#'                             names = "love")
#'
#' ## Paquid
#' library(timeROC);library(survival)
#' data(Paquid) # rm(Paquid)
#' Paquid$DSST <- -Paquid$DSST
#' Paquid$status <- ifelse(Paquid$status == 1,1,0)
#' models.merge <- getTimeROC(data = Paquid,
#'                            time.col = "time",
#'                            status.col = "status",
#'                            cause = 1,
#'                            markers = "DSST",
#'                            merge.markers = list(
#'                              merge = c("DSST","MMSE")
#'                            ),
#'                            time.raw.type="Year",
#'                            time.target.type="Year",
#'                            time.knot= c(3,5),
#'                            line.color = mycolor[c(1,3,4)],
#'                            names = "Paquid1")
#' @export
getTimeROC <- function(
  data,
  time.col,status.col,cause,
  markers,
  merge.markers=NULL,
  weighting.univariate = c("marginal","cox","aalen")[2],
  weighting.multivariate = c("marginal","cox","aalen")[2],
  time.raw.type=c("Day","Month","Year")[1],
  time.target.type=c("Day","Month","Year")[3],
  time.knot=c(3,5),
  line.color = mycolor,
  size = 1.2,half.border=F,
  legend.position = c(0.8,0.2),
  show.auc = T,auc.digits = 2,
  reference = T,
  width = 10,height=10,
  merge.plot = T,
  names = "love"
){
  ## 加载必要包
  need <- c("timeROC","pROC","survival","ggplot2","stringr")
  Plus.library(need)

  ## 选择data
  select = c(time.col,status.col,unique(c(markers,unlist(merge.markers))))
  data1 <- subset(data,select = select)
  logic1 <- apply(data1,1,is.one.na);#table(logic1)
  data1 <- data1[!logic1,]

  ## 去因子化
  for(i in 1:ncol(data1)){
    if(!is.null(levels(data1[,i]))){
      #因子型向量
      print(paste0(colnames(data1)[i]," is a factor and converted to numeric."))
      data1[,i] <- as.numeric(as.character(data1[,i]))
    }
  }
  #str(data1)

  ## 统一time的单位:
  #1 year = (12+1/6) month = 365 day;1 month = 30 day
  convert.time <- function(time,from.ts,to.ts){
    if(from.ts == "Day"){
      day <- time
      month <- time/30
      year <- time/365
    } else {
      if(from.ts == "Month"){
        day <- time*30
        month <- time
        year <- time/(12+1/6)
      } else {
        if(from.ts == "Year"){
          day <- time*365
          month <- time*(12+1/6)
          year <- time
        } else {
          stop("Error!not a right type of time!")
        }
      }
    }
    #输出结果
    time1 <- data.frame(
      Day=day,
      Month=month,
      Year=year
    )
    time2 <- time1[,to.ts]
    return(time2)
  }
  data1[,time.col] <- convert.time(data1[,time.col],
                                   from.ts = time.raw.type,
                                   to.ts = time.target.type)
  #report
  print(paste0("The max of time is ",round(max(data1[,time.col]),2)," ",time.target.type,"s,and the min is ",round(min(data1[,time.col]),2)," ",time.target.type,"s.The time.knot parameter must be in this range."))

  ## 单变量
  single.model <- NULL #i=3
  for(i in 1:length(markers)){
    marker.i <- markers[i]
    ROC.i <-timeROC(T=data1[,time.col],
                    delta=data1[,status.col],
                    marker=data1[,marker.i],
                    cause=cause,
                    weighting=weighting.univariate,
                    times=time.knot,
                    ROC=T,
                    iid=TRUE)
    single.model <- c(single.model,list(ROC.i))
    names(single.model)[i] <- marker.i
  }

  ## 多变量
  if(is.null(merge.markers)){
    multi.model <- NULL
  } else {
    multi.model <- NULL # i=1
    m.names <- names(merge.markers)
    for(i in 1:length(merge.markers)){
      marker.i <- merge.markers[[i]]
      marker.i1 <- marker.i[1]
      marker.i2 <- marker.i[2:length(marker.i)]
      ROC.i <-timeROC(T=data1[,time.col],
                      delta=data1[,status.col],
                      marker=data1[,marker.i1],
                      other_markers = as.matrix(data1[,marker.i2]),
                      cause=cause,
                      weighting=weighting.multivariate,
                      times=time.knot,ROC=T)
      multi.model <- c(multi.model,list(ROC.i))
      if(is.null(m.names[i])){
        #未指定merge.marker的名字
        names(multi.model)[i] <- paste0(marker.i,collapse = "&")
      } else {
        names(multi.model)[i] <- m.names[i]
      }
    }
  }

  ## 合并单/多变量
  model <- c(single.model,multi.model)
  level.all <- names(model)


  ## 根据某个时间和某个model进行绘图
  #model=model  #times = time.knot[1]  #model.names = names(model)
  getroc1 <- function(model,model.names,times,
                      color=mycolor,size = 1.2,
                      half.border=F,
                      legend.position = "top",
                      show.auc = T,auc.digits =2,
                      reference = T){
    ## 提取数据
    nm <- NULL;FP <- NULL;TP <- NULL;tm<- NULL;
    for(z in 1:length(model)){
      model.z <- model[[z]]
      model.names.z <- model.names[z]
      FP.z <- NULL;TP.z <- NULL;tm.z<- NULL;
      for(i in 1:length(times)){
        p <- match(times[i],time.knot)
        TP.i <- model.z[["TP"]][,p]
        TP.z <- c(TP.z,TP.i)
        FP.i <- model.z[["FP"]][,p]
        FP.z <- c(FP.z,FP.i)
        time.z <- as.character(rep(times[i],length(TP.i)))
        tm.z <- c(tm.z,time.z)
      }
      nm.z <- rep(model.names.z,length(FP.z))
      ## 输出
      nm <- c(nm,nm.z);FP <- c(FP,FP.z);
      TP <- c(TP,TP.z);tm <- c(tm,tm.z);
    }

    ## 假阳横，真阳纵
    df1 <- data.frame(FP=FP,TP=TP,time=tm,names=nm)
    level = level.all[level.all %in% unique(df1$names)]
    df1$names <- factor(df1$names,levels = level)
    ncol <- length(unique(df1$names))

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
      for(i in 1:length(model)){ # i=1
        auc.i <- model[[i]]
        auc.i <- round(auc.i$AUC,auc.digits)
        time.names <- paste0("t=",times)
        auc.i <- auc.i[names(auc.i) %in% time.names]
        labels.i <- paste0(Fastextra(names(auc.i),"=",2)," ",time.target.type," = ",as.numeric(auc.i))
        auc <- c(auc,labels.i)
      }
      auc.labels <- str_c(names(model),"(AUC:",auc,")")
      gg.auc <- scale_color_manual(
        values = color[1:ncol],
        aesthetics = "colour",
        breaks = names(model),
        labels = auc.labels)
    } else {
      gg.auc <- scale_color_manual(
        values = color[1:ncol],
        aesthetics = "colour")
    }


    ## 图 #aes(linetype=names),
    colnames(df1)[match("names",colnames(df1))] <- "Groups"
    p <- ggplot(df1, aes(x=FP,y=TP,color = Groups)) +
      geom_line(size = size) +
      r + gg.auc +
      labs(x = "1 - Specificity",y = "Sensitivity") +
      theme_bw() +
      theme(panel.grid =element_blank(),
            axis.title = element_text(face = "bold",size = 15),
            axis.text = element_text(face = "bold",size = 12),
            legend.title = element_text(face = "bold",size = 12),
            legend.text =element_text(face = "bold",size = 12),
            legend.position = legend.position,
            panel.border=panel.border,
            axis.line = axis.line);
    return(p)
  }

  ## 画图 i=1
  if(merge.plot == T){
    #按时间来融合图片
    print("Merge plot by times.")
    p1 <- list()
    pdf(paste0(names,"_merge_Time-dependent ROC curve.pdf"),width = width,height = height)
    for(i in 1:length(time.knot)){
      time.i <- time.knot[[i]]
      p <- getroc1(model,names(model),time.i,
                   color=line.color,size = size,
                   half.border=half.border,
                   legend.position = legend.position,
                   show.auc = show.auc,auc.digits =auc.digits,
                   reference = reference)
      p <- p +
        ggtitle(paste0("ROC at ",time.i," ",time.target.type)) +
        theme(plot.title = element_text(hjust = 0.5,
                                        face = "bold",
                                        size = 15))
      print(p)
      p1 <- c(p1,list(p))
      names(p1)[i] <- paste0(time.i,"-",time.target.type)
    }
    dev.off()
  } else {
    print("plot every model in every time.")
    p1 <- list()
    pdf(paste0(names,"_single_Time-dependent ROC curve.pdf"),width = width,height = height)
    for(i in 1:length(time.knot)){
      time.i = time.knot[i]
      p.i <- list()
      for(z in 1:length(model)){ #z=1;i=1
        model.z <- model[z]
        plot.main = paste0(names(model)[z],"_ROC at ",time.knot[i]," ",time.target.type,",AUC=",round(model.z[[1]][["AUC"]][i],2))
        p <- getroc1(model = model.z,
                     model.names = names(model)[z],
                     times = time.i,
                     color=line.color,size = size,
                     half.border=half.border,
                     legend.position = legend.position,
                     show.auc = show.auc,auc.digits =auc.digits,
                     reference = reference) +
          ggtitle(plot.main) +
          theme(plot.title = element_text(hjust = 0.5,
                                          face = "bold",
                                          size = 15))
        print(p)
        p.i <- c(p.i,list(p))
        names(p.i)[z] <- names(model)[z]
      }
      p1 <- c(p1,list(p.i))
      names(p1)[i] <- paste0(time.i,"-",time.target.type)
    }
    dev.off()
  }
  print("ROC plots have been saved in the work space.")

  ## 输出数据
  result <- list(
    model=model,
    plot = p1
  )
  return(result)
}







