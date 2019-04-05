


####=========================Fastcorrplot()====================####
## Usage:help fast to draw a corrplot with multiple styles available to select.Fastcorrplot() is based on corrplot::corrplot.mixed().

# expr#gene expression matrix or data frame with patient cols and gene rows.
# select#genes your want to explore the correlations.Note that the select genes must be part of rownames of expr.
# col=NULL#if Null,use grDevices::colorRampPalette(c("blue","white","red"))(100).Also you can defined it alternatively.
# savefile=T# whethe saving plot as PDF file.
# names="test1"# part of PDF file names.
#' @export
Fastcorrplot <- function(data,
                         transposition = T,#是否转置矩阵
                         select,
                         order = T,
                         lower.col = NULL,
                         upper.col =NULL,
                         upper = NULL,
                         tl.pos = NULL,
                         tl.col=NULL,
                         tl.srt=NULL,
                         diag = NULL,
                         width=10,height=10,
                         savefile=T,
                         names="test1"){

  ## 加载相应的包
  library(corrplot)

  ## 矩阵设置
  expr1 <- as.matrix(data)
  if(transposition==T){
    # 矩阵要转置
    expr2 <- t(expr1)
  } else {
    expr2 <- expr1
  }
  expr2 <- expr2[,select]

  ## 计算矩阵相关系数并决定是否自动排序
  mycor <- cor(expr2)
  if(order == F){
    mycor2 <- mycor
  } else {
    ord <- corrMatOrder(mycor, order = "AOE")
    mycor2 <- mycor[ord,ord]
  }

  ## 加载默认或者自定义信息
  library(grDevices)
  col <- colorRampPalette(c("blue","white","red"))(100)
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

  ## corrplot
  if(savefile==T){
    win.graph(width = width,height = width);corrplot.mixed(mycor2,
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    file.names <- paste0(names,"_Corrplot.pdf",collapse = "")
    pdf(file.names,width = width,height = height)
    corrplot.mixed(mycor2,
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    dev.off()
    print(paste0(file.names,"已经保存至当前工作目录",collapse = ""))
    return(mycor2)
  } else {
    win.graph(width = width,height = width)
    corrplot.mixed(mycor2,
                   lower.col = output[[1]],
                   upper.col =output[[2]],
                   upper = output[[3]],
                   tl.pos = output[[4]],
                   tl.col=output[[5]],
                   tl.srt=output[[6]],
                   diag = output[[7]])
    return(mycor2)
  }
#End
}





