

## cox.optimized4 help get survival significant markers based on KM-Curve data

#expr.matrix # expression matrix
#design # design object
#time.col # colname of time value
#status.col # colname of status value
#strategy # one of "median" and "maxstat"
#smethod # Default is "LogRank"
#pmethod# Default is "HL",
#cutoff.p#p value cut-off
#digits # digits of result
#parallel#whether to use parallel strategy when strategy is median.Note that the parallel would always be applied with maxstat method
#' @export
cox.optimized4 <- function(expr.matrix,
                           design,
                           time.col,
                           status.col,
                           strategy  = c("median","maxstat")[1],
                           smethod="LogRank",
                           pmethod="HL",
                           cutoff.p = 0.05,
                           digits = 5,
                           parallel = F){
  ## 加载必要的包
  nd <- c("plyr","survival","maxstat")
  Plus.library(nd)

  ## 对齐
  expr.matrix <- expr.matrix[,rownames(design)]

  ## 过滤掉表达为0的基因
  lg1 <- rowMeans(expr.matrix) > 0;
  expr.matrix <- expr.matrix[lg1,]


  ## time & status
  time <- as.numeric(as.character(design[,time.col]))
  status <- as.numeric(as.character(design[,status.col]))

  ## median
  if(strategy == "median"){
    #中位数作为cut-off值
    getsig <- function(expr.matrix.i){
      expr.matrix.i <- as.numeric(as.character(expr.matrix.i))
      m.i <- median(expr.matrix.i)
      expr <- ifelse(expr.matrix.i > m.i,1,0)
      df1 <- cbind(status,time,expr)
      df1 <- as.data.frame(df1)
      #计算p值
      sdf <- survdiff(Surv(time, status) ~ expr, data = df1)
      p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      p.val <- round(p.val,digits)
      df1 <- data.frame(cutoff = m.i,P.val = p.val)
      return(df1)
    }
    print("Binary P value of every gene is on calculation,please wait...")
    if(parallel == F){
      df1 <- adply(expr.matrix,1,getsig,.id = "genes")
    } else {
      print("Use parallel algrithm...")
      nd <- c("parallel","foreach","doParallel");Plus.library(nd)
      ncore <- detectCores(logical = F)
      cl <- makeCluster(mc <- getOption("cl.cores",ncore))
      registerDoParallel(cl)
      clusterExport(cl=cl, c("getsig", "expr.matrix", "status", "time","maxstat.test","Surv","survdiff","digits"), envir=environment())
      df1 <- adply(expr.matrix,1,getsig,.parallel = T,.id = "genes")
      stopCluster(cl)
    }
    colnames(df1)[1] <- "genes"
    print("Calculation is done!")
    lg2 <- df1$P.val < cutoff.p
    df3 <- subset(df1,subset = lg2)
    ## 排序
    df3 <- df3[order(df3$P.val,decreasing = F),]
  } else {
    ## maxstat运算自动优化
    getsig2 <- function(expr.matrix.i){
      expr <- as.numeric(as.character(expr.matrix.i))
      df1 <- cbind(status,time,expr)
      df1 <- as.data.frame(df1)
      #取cutoff值
      mtHL <-maxstat.test(Surv(time, status) ~ expr, data=df1, smethod=smethod, pmethod=pmethod)
      cutoff1 <- as.numeric(mtHL$estimate)
      cutoff1 <- round(cutoff1,digits = digits)
      ##0或1
      df1$expr <- ifelse(df1$expr > cutoff1,1,0)
      #计算p值
      sdf <- survdiff(Surv(time, status) ~ expr, data = df1)
      p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      df2 <- data.frame(cutoff = cutoff1,P.val = p.val)
      return(df2)
    }
    names <- rownames(expr.matrix)
    names <- data.frame(genes = names)
    print("maxstat:binary P value of every gene is on calculation,please wait...")
    nd <- c("parallel","foreach","doParallel");Plus.library(nd)
    ncore <- detectCores(logical = F)
    cl <- makeCluster(mc <- getOption("cl.cores",ncore))
    registerDoParallel(cl)
    clusterExport(cl=cl, c("getsig2", "expr.matrix", "status", "time","maxstat.test","Surv","survdiff","smethod","pmethod","digits"), envir=environment())
    df1 <- adply(expr.matrix,1,getsig2,.parallel = T,.id = "genes")
    stopCluster(cl)
    colnames(df1)[1] <- "genes"
    print("Calculation is done!")
    lg2 <- df1$P.val < cutoff.p
    df3 <- subset(df1,subset = lg2)

    ## 排序
    df3 <- df1[order(df1$P.val,decreasing = F),]
  }

  ## 输出结果
  rownames(df3) <- df3$genes
  return(df3)

}











