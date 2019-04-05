


# data # a data frame
# select # selected marker colnames
# legend.position = "right"# legend position
# half.border = T #if show half of black border of the plot
# width = 10 # new showing window width
# height = 10 # new showing window height
# save.file = T # whether save PDF file
# names = "love" # part of PDF file name
#' @export
Fastpoint1 <- function(data,
                       select,
                       color="#FB8072",
                       select.color=1:2,
                       legend.position = "right",
                       half.border = T,
                       alpha = 0.5,
                       width = 10,
                       height = 10,
                       save.file = T,
                       names = "love"){
  ## 加载包
  library(ggplot2);library(stringr)
  Fastpoint1.logic1 <- color %in% colnames(data)
  if(Fastpoint1.logic1){
    #设定颜色
    library(RColorBrewer)
    c <- brewer.pal(12, "Set3")[c(4,5,1,3,2,6:12)]
    #数据列名,进行color或fill运算
    s <- c(select,color)
    position <- Fastmatch(s,colnames(data))
    colnames(data)[position] <- c(paste("v",1:length(select),sep = ""),"color")
    g <- ggplot(data,aes(x=v1,y=v2,colour = color)) + geom_point(alpha = alpha,size = 4) + scale_colour_manual(values = c[select.color],aesthetics = "colour")
  } else {
    #非数据列名，因此是颜色向量
    data1 <- subset(data,select = select)
    colnames(data1) <- paste("v",1:length(select),sep = "")
    g <- ggplot(data1,aes(x=v1,y=v2)) + geom_point(alpha = alpha,size = 4,colour = color)
  }
  # data1一定要作为一个中间变量，否则程序运行出错。
  #data1 <- data1[1:3,]
  #save(data1,file = "data1.rda")

  ## 背景边框:是否用半保留边框
  if(half.border == T){
    panel.border=element_rect(fill='transparent',color='transparent')
    axis.line = element_line(colour = "black")
  } else {
    panel.border=element_rect()
    axis.line = element_line()
  }

  ## 提取数据
  p <-  g +
    xlab(select[1]) + ylab(select[2]) +
    ggtitle(paste0("Comparision:",paste(select[1:2],collapse = "-"))) +
    theme_bw() +
    theme(panel.grid =element_blank(),
          plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
          axis.title = element_text(face = "bold",size = 15),
          axis.text = element_text(face = "bold",size = 12),
          legend.title = element_text(face = "bold",size = 12),
          legend.text =element_text(face = "bold",size = 12),
          legend.position = legend.position,
          panel.border=panel.border,
          axis.line = axis.line);
  if(save.file ==F){
    return(p)
  } else {
    win.graph(width = width,height = height);print(p)
    ggsave(paste0(names,"_","PointPlot.pdf"),p,width = width,height = height)
    print(paste0(names,"_","PointPlot.pdf已经成功保存!"))
    return(p)
  }
}


