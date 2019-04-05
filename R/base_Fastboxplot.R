
#' @title fast way to draw a boxplot based on ggplot2
#' @description fast way to draw a boxplot
#' @param data a data frame
#' @param x,y,fill,facet parameters of ggpubr::ggboxplot
#' @param palette fill color of box
#' @param x.lab the title of x axis
#' @param y.lab the title of y axis
#' @param legend.position legend.position
#' @param x.test.angle  the angle of text on the x axis
#' @param add.jitter  whether add points with jitter style
#' @return a boxplot
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data(mtcars)
#' data=mtcars
#' x="vs"
#' y="disp"
#' fill="vs"
#' facet=NULL
#' facet="am"
#' x.lab="vs of mtcars"
#' y.lab="disp of mtcars"
#'
#' ## quick start of Fastboxplot
#' Fastboxplot(data,
#'             x,y,fill,facet,
#'             palette=NULL,
#'             x.lab,y.lab,
#'             legend.position="none",
#'             x.test.angle=0,
#'             add.jitter=T)
#'
#' #change fill color
#' show.color(mycolor)
#' Fastboxplot(data,
#'             x,y,fill,facet,
#'             palette=c("#8DD3C7","#FB8072"),
#'             x.lab,y.lab,
#'             legend.position="none",
#'             x.test.angle=0,
#'             add.jitter=T)
#' @export
Fastboxplot <- function(data,
                        x,y,fill,facet,
                        palette=NULL,
                        x.lab,y.lab,
                        legend.position="none",
                        x.test.angle = 0,
                        add.jitter=T){

  ## 包
  library(ggpubr);library(ggplot2)

  ## 构建颜色向量
  library(lucky)
  if(is.null(palette)){
    palette = rep(mycolor,100)
    palette1 = palette[1:length(unique(data[,x]))]
  } else {
    palette1=palette
  }

  ## 是否添加jitter
  if(add.jitter==T){
    add="jitter"
  }else {
    add="none"
  }

  ## 分面
  if(is.null(facet)){
    fa <- NULL
    th <- theme(
      axis.text= element_text(size = 12,colour = "black",face = "bold"),
      axis.title = element_text(size = 15,colour = "black",face = "bold"),
      legend.position=legend.position)
  } else {
    f <- paste0(". ~ ",facet)
    f <- as.formula(f)
    fa <- facet_grid(f)
    th <- theme(
      axis.text = element_text(size = 12,colour = "black",face = "bold"),
      axis.title= element_text(size = 15,colour = "black",face = "bold"),
      legend.position=legend.position,
      strip.background = element_rect(fill="white"))
  }

  ###绘图
  p <- ggboxplot(data,
                 x = x, y = y,fill = fill,#箱体填充颜色
                 color = "black",
                 palette = palette1,
                 add = add) + #加上jitter点
    labs(x = x.lab,y = y.lab) + fa + th +
    rotate_x_text(angle = x.test.angle) #x轴的text以逆时针旋转45度角

  ### return
  return(p)
}


