

#' @title Fast Way to Get a Bar plot based on ggpubr
#' @description Fastbar help get a barplot of ggplot2 style.
#' @keywords FastBar
#' @param data a data frame
#' @param x x is usually the colnames of symbol/names
#' @param y y is usually the colnames of frequency
#' @param order whether order from upper to lower
#' @param fill Default is white.Set fill=x is a useful way to enhance visualize barplot and recommanded
#' @param color Default is NULL,and use x as color scale
#' @param palette Default is NULL,and it use a lucky::mycolor color strategy.Or you can set self-defined colors.Note that the number of color you want to set should be the same as the unique number of x
#' @param label Default is TRUE.If variables are too much,you can set FALSE to enhance visualization
#' @param size a size control for plot title,axis titles and labels
#' @param position see ggpubr::ggbarplot
#' @param title the title of the plot
#' @param legend.position the position of legend.you can select "top","left","right" and "none"
#' @return a ggplot list
#' @seealso ggpubr::ggbarplot
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and NOT RUN
#' p  <- FastBar(data = df,
#'               x = "SYMBOL",y = "Freq",
#'               fill = "SYMBOL",
#'               title = "Frequecy",
#'               legend.position = "none",
#'               size = 20)
#' @export
testFastBar <- function(data,x,y,
                        order=NULL,
                        fill = "white",
                        color = NULL,
                        palette =NULL,
                        label = TRUE,
                        size = 15,
                        position = position_stack(),
                        title = "Frequecy",
                        legend.position = "right"){

  ## 加载必要的包
  nd <- c("ggplot2","ggpubr")
  Plus.library(nd)

  ## 加载默认值
  if(is.null(color)){
    color <-  x
  }
  if(is.null(palette)){
    x.type <- length(unique(as.character(data[,x])))
    pa <- c(1,3,4,5,6)
    pa <- mycolor[c(pa,setdiff(1:70,pa))]
    pa <- rep(pa,100)
    palette <- pa[1:x.type]
  }

  ## 绘制图
  p <- ggbarplot(data,x = x,y = y,
                 order = order,
                 fill = fill, color = color,
                 palette = palette,
                 label = label,
                 position = position,
                 lab.size = (size/15)*4) +
    ggtitle(title) +
    rotate_x_text(angle = 45) +
    theme(panel.grid =element_blank(),
          plot.title = element_text(face = "bold",size = (size/15)*25,hjust = 0.5),
          axis.title = element_text(face = "bold",size = (size/15)*22),
          axis.text = element_text(face = "bold",size = (size/15)*12),
          legend.title = element_text(face = "bold",size = (size/15)*12),
          legend.text =element_text(face = "bold",size = (size/15)*10),
          legend.position = legend.position)

  ## 环境记录
  Repeat <- list(
    order=order,
    fill = fill,
    color = color,
    palette = palette,
    label = label,
    size = size,
    legend.position = legend.position
  )

  ## 输出结果
  a <- list(
    Repeat = Repeat,
    Data = NA,
    Plot = p
  )
  class(a) <- "LuckyObject"

  return(a)
}

## S3 method
print.testFastBar <- function(x){
  return(x$Plot)
}


## S4 method

# S4 method
setClass("LuckyObject", representation = list(Repeat = "list", data = "list", Plot = "list"))

# S5 method
setMethod("show", "LuckyObject", function(object){
 return(object@Plot)
})














