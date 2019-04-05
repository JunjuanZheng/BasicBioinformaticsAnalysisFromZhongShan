
###compareROC can help produce a time dependent KM curve,such as 1-year,3-year or 5-year KM curves.
#' @title Produce time dependent KM-ROC curves
#' @description  compareROC can help produce a time dependent KM curve,such as 1-year,3-year or 5-year KM curves.
#' @param data a data frame with survival data like time and status
#' @param i.values c("PNI1","NIHrisk").colnames of markers
#' @param merge.i.values list(merge1 = i.values).a list of merge markers selection
#' @param timecutoff c(1-365,3-365,5-365).cut-off for time series
#' @param roc.type c("ggplot","pROC")[1].ggplot recommanded
#' @param color brewer.pal(12, "Set3")[c(1,3,4)].curve colors
#' @param output.name part of PDF output name
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data(fat)
#' # markers=c("BMI","Waist","WHR","belly fat thickness") #marker colnames
#' # status =  "outcomes" #survival status
#' # merge.markers = list(c("BMI","Waist"))# the markers list that you want to merge as a co-prognostic factor
#' # return.plot = T # if T，return plot;if F，return data
#' # roc.type = c("ggplot","pROC")[1] # ggplot style is recommanded
#' # title="this is a title" # plot title
#' # color = brewer.pal(12, "Set3")[1:5] #curve colors.The lenght of colors must >= the number of curves.
#' # output.name = "ROC of something" # part of PDF name
#' compareROC2(data=fat,
#'             markers,
#'             status,
#'             merge.markers=NULL,
#'             return.plot=T,
#'             roc.type,
#'             title,
#'             color,
#'             output.name="test1")
#' @export
compareROC <- function(data,
                       i.values,
                       merge.i.values,
                       timecutoff,
                       roc.type,
                       color,
                       output.name,
                       return.plot){
  #加载包
  library(pROC);library(ggplot2);library(RColorBrewer);library(stringr)

  #去除空值
  data <- data[!is.na(data$DFS),]

  #将time/status转变为二分类
  time.convert <- function(time,status,days.cutoff){
    ifelse(time > days.cutoff,0,ifelse(status == 1,1,0))
  }

  D <- NULL;
  for (year in timecutoff) {
    data$DFS.year <- time.convert(time = data$DFS,data$status,days.cutoff = year)

    #计算单变量
    single.models <- list()
    for (j in 1:length(i.values)) {
      p.j <- match(i.values[j],colnames(data))
      roc.j <- roc(data$DFS.year,data[,p.j])
      single.models <- c(single.models,list(roc.j))
      names(single.models)[j] <- i.values[j]
    }

    #计算多变量
    multi.models <- list()
    for(j in 1:length(merge.i.values)){
      data1 <- data[,c("DFS.year",merge.i.values[[j]])]
      model.j <- glm(DFS.year ~ ., data=data1, family='binomial')
      pre.j <- predict(model.j)
      roc.j <- roc(data1$DFS.year,pre.j)
      multi.models <- c(multi.models,list(roc.j))
      names(multi.models)[j] <- paste0(merge.i.values[[j]],collapse = "_")
    }

    #合并单/多变量的model
    roc.models <- c(single.models,multi.models)

    #绘制2个ROC图
    if(return.plot == T){
      #绘制图并保存
      if(roc.type == "ggplot"){
        #ggplot风格
        g <- ggroc(roc.models) +
          theme_bw() + #白色背景
          geom_line(size = 1.2) + #线粗细
          scale_colour_manual(values = color) + #线颜色
          labs(title=paste(year/365,"Years DFS"),x = "Specificity",y = "Sensitivity") + #标题
          theme(plot.title = element_text(hjust = 0.5),#总标题居中
                title = element_text(size = 15,face = "bold"),#总标题字体/尺寸
                axis.title = element_text(size = 15,face = "bold"),#坐标标题字体/尺寸
                axis.text = element_text(size = 12,face = "bold"),#坐标标尺字体/尺寸
                legend.title=element_blank(),#去除legend.title
                legend.text =element_text(size = 12,face = "bold"),#legend字体/尺寸
                panel.grid.major = element_blank(),#去除网格
                panel.grid.minor = element_blank());
        ggsave(str_c(year/365," Year_",output.name,"_ggplot.style.pdf"),g,width = 10,height = 8)
        print(str_c("完成",str_c(year/365," Year_",output.name,"_ggplot.style.pdf"),"的绘制。"))
      } else {
        if(roc.type == "pROC"){
          #pROC的风格~研究不深哈，就随便画一个了~
          pdf(str_c(year/365," Year_",output.name,"_pROC.style.pdf"),width = 8,height = 8)
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
          print(str_c("完成",str_c(year/365," Year_",output.name,"_pROC.style.pdf"),"的绘制。"))
        }
      }
    } else {
      #保存两两比较的数据
      n = 1:length(roc.models);n.mt <- combn(n,2)
      if(length(n)==1){print("仅单个变量，无法比较")} else {
        #可进行两两比较
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
        d <- list(compare = d,roc.models = roc.models)
      }
      D <- c(D,list(d))
    }
  }

  #result
  if(return.plot ==F){
    names(D) <- paste("Year",timecutoff/365,sep = "_")
    return(D)
  } else {print("End")}
}
#data # a data frame with survival data like time and status
#i.values = c("PNI1","NIHrisk") # colnames of markers
#merge.i.values = list(merge1 = i.values) # a list of merge markers selection.
#timecutoff = c(1*365,3*365,5*365) # cut-off for time series
#roc.type = c("ggplot","pROC")[1] # ggplot recommanded
#color = brewer.pal(12, "Set3")[c(1,3,4)] # curve colors
#output.name = "test1" # part of PDF output name


