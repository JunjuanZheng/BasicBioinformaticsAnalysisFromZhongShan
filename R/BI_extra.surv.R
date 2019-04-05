

# extra.surv help get survival data of a gene based on expression and design objects.Note that multiple genes would not be supported.

# expr.matrix # expression matrix
# design# design object
# select# a gene.Multiple genes are not supported
# event.status# a colname of status
# event.time# a colname of time
# event.lower# the lower limit of the time
# mode# Default is "median".You can also set a number of percentum or a decimal less than 1
#' @export
extra.surv <- function(expr.matrix,
                       design,
                       select,
                       event.status,
                       event.time,
                       event.lower,
                       mode=c("median","maxstat")[1]){
  ## 加载必要的包
  nd <- c("survival","maxstat")
  Plus.library(nd)

  ## 对齐
  expr.matrix <- expr.matrix[,rownames(design)]

  ##去除event和cluster对应元素的NA值。
  test1 <- design[,c(event.status,event.time)]
  logic1 <- apply(test1,1,is.one.na)
  design1 <- design[!logic1,]
  design2 <- subset(design1,select=c(event.time,event.status))
  expr.matrix1 <- expr.matrix[,rownames(design1)]

  ## 基因表达量的二分类函数
  fun1 <- function(expr.matrix1,design1,select,mode){
    gene1.expr <- expr.matrix1[select,]
    two.tie <- function(vt,mode){
      ## mode setting
      if(mode != "maxstat"){
        if(mode %in% c("median","mean")){
          vt1 <- ifelse(mode == "median",median(vt),mean(vt))
        } else {
          #mode is a portion number
          vt <- sort(vt,decreasing = F)
          vt.count <- length(vt)
          vt.tie.count <- floor(vt.count*mode)
          vt1 <- vt[vt.tie.count]
        }
      } else {
        ##maxstat运算
        print("maxstat auto cut-off...")
        time <- as.numeric(as.character(design1[,event.time]))
        status <- as.numeric(as.character(design1[,event.status]))
        expr <- gene1.expr
        df1 <- data.frame(time = time,status=status,expr = expr)
        mtHL <- maxstat.test(Surv(time, status) ~ expr, data=df1, smethod="LogRank", pmethod="none")
        vt1 <- as.numeric(mtHL[["estimate"]])
      }
      return(vt1)
      #End
    }
    cut.off <- two.tie(gene1.expr,mode = mode)
    sur.vt <- ifelse(gene1.expr > cut.off,1,0)
    l <- list(sur.vt = sur.vt,cut.off = cut.off)
    return(l)
  }

  ##构建生存矩阵
  list.fun1 <- fun1(expr.matrix1,design1,select,mode=mode)
  surv.matrix1 <- as.data.frame(list.fun1$sur.vt)
  surv.matrix1 <- as.data.frame(surv.matrix1[rownames(design1),])
  colnames(surv.matrix1) <- select
  surv.matrix1 <- cbind(design2,surv.matrix1)
  lg1 <- as.numeric(as.character(surv.matrix1[,event.time])) > event.lower
  surv.matrix1 <- surv.matrix1[lg1,]

  ## 算p值
  ##判断有无差异
  colnames(surv.matrix1) <- c("time","status","gene")#某个基因的某个事件及混杂因素的数据
  sdf <- survdiff(Surv(time, status) ~ gene, data = surv.matrix1)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

  ## 输出结果
  l1 <- list(
    P.Value = p.val,
    CutOff = list.fun1$cut.off,
    Metadata = surv.matrix1
  )
  return(l1)
}











