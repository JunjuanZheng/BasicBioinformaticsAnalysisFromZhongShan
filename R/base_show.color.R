
# Usage preview color vector
# yourcolor # color vector like "#8DD3C7" or other type that can be identified by R.
#' @export
show.color <- function(yourcolor){
  ## 数据框构建
  yourcolor <- as.character(yourcolor)
  cd <- data.frame(
    y=2,
    color=yourcolor
  )
  cd$color <- factor(cd$color,levels = yourcolor)
  ## 柱状图
  library(ggpubr)
  p <- ggbarplot(cd,
                 x = "color", y = "y",
                 fill = "color",#箱体填充颜色
                 palette = yourcolor)+
    labs(x="Color Types",y="Color show") +
    ggtitle("yourcolor")+
    rotate_x_text(angle = 45)+
    annotate("text", x = 1:length(yourcolor), y = 2.05, label = as.character(1:length(yourcolor)))+
    theme(legend.position = "none",
          plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
          axis.title = element_text(face = "bold",size = 15),
          axis.text = element_text(face = "bold",size = 10))

  ## 加上文字注释
  print(p)
  return(yourcolor)
}


