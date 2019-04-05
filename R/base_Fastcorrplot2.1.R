

###=====Fastcorrplot2.1
## Usage:give one control markers,then test the correlation between every control marker and target.markers.

## parameters
# data #a gene expression matrix or a data frame with patient cols and gene rows.Or a matrix and a dataframe with similar construction.
# transposition = T # whethe make data transposition.
# control.markers # a gene used as internal reference to explore correlations.
# target.markers=NULL #other genes differing from control markers.It must be part of data rownames(transposition = T) or colnames.If NULL,it is the other set beyond control markers.
# method # see psych::corr.test
# p.cut.off # a cut.off of significance.Default = 0.05.
# savefile=T #Whether to save a PDF plot.
# names # part of PDF file name.
# other parameters# see corrplot::corrplot().

Fastcorrplot2.1 <- function(data,
                          transposition = T,#是否转置矩阵
                          control.markers,
                          target.markers=NULL,
                          method="pearson",p.cut.off=0.05,
                          savefile=T,#corrplot()相关参数
                          names="test1",
                          lower.col = NULL,#corrplot()相关参数
                          upper.col =NULL,#corrplot()相关参数
                          upper = NULL,#corrplot()相关参数
                          tl.pos = NULL,#corrplot()相关参数
                          tl.col=NULL,#corrplot()相关参数
                          tl.srt=NULL,#corrplot()相关参数
                          diag = NULL){
  ## 加载相应的包
  library(corrplot);library(psych)

  ## 加载默认设置
  library(grDevices)
  col <- colorRampPalette(c("blue","white","red"))(100)
  l1 <- colnames(data) %in% control.markers
  default <- c(rep(list(NULL),8)) #
  input <- list(target.markers=target.markers,
                lower.col = lower.col,
                upper.col = upper.col,
                upper = upper,
                tl.pos = tl.pos,
                tl.col=tl.col,
                tl.srt=tl.srt,
                diag = diag)
  do <- list(target.markers=colnames(data)[!l1],
             lower.col = col,
             upper.col =col,
             upper = "pie",
             tl.pos = "lt",
             tl.col="black",
             tl.srt=45,
             diag = "l")
  output <- set.default(input,default,do);
  #names(output)
  #[1] "target.markers" "lower.col"      "upper.col"
  #[4] "upper"          "tl.pos"         "tl.col"
  #[7] "tl.srt"         "diag"

  ### data transposition
  expr1 <- as.matrix(data)
  if(transposition==T){
    # 矩阵要转置
    expr2 <- t(expr1)
  } else {
    expr2 <- expr1
  }
  expr2 <- expr2[,c(control.markers,output[[1]])]

  ### 计算control marker与target marker之间的相关性，进行差异性检验
  library(psych)
  cor.t1 <- corr.test(expr2,
                      adjust = "none",
                      use = "complete",
                      method = method);cor.t1
  p.matrix <- cor.t1$p
  p.matrix1 <- ifelse(p.matrix < p.cut.off,T,F)

  ## 获得对应的长型数据
  getlong <- function(matrix){
    p.list <- list()
    for(i in 1:ncol(matrix)){
      l.i <- matrix[,i]
      names(l.i) <- rownames(matrix)
      p.list <- c(p.list,list(l.i))
    }
    names(p.list) <- colnames(matrix)
    ## 生成长型data.frame
    p.long <- stack(p.list);
    p.long$pair2 <- rep(rownames(matrix),ncol(matrix))
    colnames(p.long)[1:2] <- c("logic","pair1")
    return(p.long)
  }
  p.long1 <- getlong(p.matrix)
  p.long2 <- getlong(p.matrix1)
  p.long <- merge(p.long1,p.long2,by=c("pair1","pair2"))
  colnames(p.long)[3:4] <- c("P.val","logic")
  p.long <- p.long[p.long[,"logic"] %in% T,]
  p.long <- p.long[p.long$pair1 != p.long$pair2,]#将自我比较的行去除

  ## get information about control.markers
  p.long1 <- p.long[p.long$pair1 == control.markers,] #positive correlations between controls and targets
  select = c(unique(as.character(p.long1$pair1)),
             unique(as.character(p.long1$pair2)))

  ## Use Fastcorrplot() to draw corrplot.
  cor.matrix <- Fastcorrplot(data=expr2,
                             transposition=F,
                             order = F,
                             select=select,
                             lower.col = output[[2]],
                             upper.col =output[[3]],
                             upper = output[[4]],
                             tl.pos = output[[5]],
                             tl.col=output[[6]],
                             tl.srt=output[[7]],
                             diag = output[[8]],
                             savefile=T,
                             names = paste0(names,"_",control.markers))
  return(p.long)
}

