
## cox.optimized3 help select prognosis-predictive genes base on expression matrix and related design object

# expr.matrix# expression matrix
# design# design object
# p.cutoff # the cut-off of P value.Default is 0.0001
# FDR.cutoff# the cut-off of BH adjusted P value.Default is 0.001.If NULL,the FDR filter would not be applied
# event.status#character, the colnames of event status
# event.time#]character.the colnames of event time
# event.lower# numeric.the lower border of time
#' @export
cox.optimized3 <- function(expr.matrix,
                           design,
                           sd.cutoff = NULL,
                           p.cutoff = 0.0001,
                           FDR.cutoff = NULL,
                           event.status=c("TTR.status","DFS.status","OS.status"),
                           event.time=c("TTR.time","DFS.time","OS.time"),
                           event.lower =c(89,89,89),
                           show.music = T){

  ## 加载必要的包
  need <- c("survival","plyr","foreach","doParallel","tuneR")
  Plus.library(need)

  ## 对齐
  expr.matrix <- as.matrix(expr.matrix)
  expr.matrix <- expr.matrix[,rownames(design)]

  ## sd filter
  if(!is.null(sd.cutoff)){
    print(paste0("Delete genes with standard deviation less than ",sd.cutoff,"..."))
    sd <- apply(expr.matrix,1,sd)
    expr.matrix1 <- expr.matrix[sd >= sd.cutoff,]
    design1 <- design
    sd <- sd[sd >= sd.cutoff]
  } else {
    expr.matrix1 <- expr.matrix
    design1 <- design
    sd <- 0
  }

  ## base.environment:cox.onegene/myfun/cox.gene
  ## 批量计算函数
  cox.gene <- function(expr.matrix1,
                       design1,
                       event.status.i,
                       event.time.i,
                       event.lower.i,
                       names){
    ## 加载必要的包
    need <- c("plyr","foreach","doParallel")
    Plus.library(need)

    ## 对于某个基因，进行单因素cox风险回归计算
    cox.onegene <- function(gene,
                            expr.matrix1,
                            design1,
                            event.status.i,
                            event.time.i,
                            event.lower.i){
      need <- c("survival")
      Plus.library(need)

      ##提取生存数据和生存时间
      time <- as.numeric(as.character(design1[,event.time.i]))
      status <- as.numeric(as.character(design1[,event.status.i]))

      ## 提取基因信息
      expr <- as.numeric(as.character(expr.matrix1[gene,]))

      ## 数据框
      df1 <- cbind(time,status,expr)
      df1 <- as.data.frame(df1)
      logic <- df1$time > event.lower.i
      df1 <- subset(df1,subset = logic)

      ## cox回归
      fit1 <- coxph(Surv(time,status)~., data = df1)
      result1 <- summary(fit1)[["coefficients"]]
      rownames(result1) <- gene

      ## 输出结果
      result1 <- as.data.frame(result1)
      colnames(result1) <- c("coef","exp.coef","se.coef","z","P")
      return(result1)
    }

    ## 组装函数
    myfun <- function(x){cox.onegene(
      gene = x[,1],
      expr.matrix1=expr.matrix1,
      design1=design1,
      event.status.i = event.status.i,
      event.time.i = event.time.i,
      event.lower.i = event.lower.i
    )}

    ## 批量计算多个基因的cox模型相关参数
    genes <- as.data.frame(rownames(expr.matrix1));
    colnames(genes) <- "genes"
    print(paste0(names,":The parameters of cox model for every genes are on calculation..."))
    ncore <- detectCores(logical = F)
    cl <- makeCluster(mc <- getOption("cl.cores",ncore))
    registerDoParallel(cl)
    clusterExport(cl,c("expr.matrix1","design1","event.status.i","event.time.i","event.lower.i","cox.onegene","coxph","Surv","Plus.library","myfun"),envir = environment())
    df1 <- ddply(.data = genes,"genes",myfun,.parallel = T)
    stopCluster(cl)
    print(paste0(names,":calucalation completed!"))
    return(df1)
  }

  ## 计算某种time和status
  L1 <- list()
  for(i in 1:length(event.status)){ # i=3
    event.status.i <- event.status[[i]];
    event.time.i <- event.time[[i]];
    event.lower.i <- event.lower[[i]];
    a = Fastextra(event.status[i],"[.]",1)
    df2 <- cox.gene(expr.matrix1,
                    design1,
                    event.status.i = event.status[i],
                    event.time.i = event.time[i],
                    event.lower.i = event.lower[i],
                    names = a)

    ## 整理数据
    df2.1 <- cbind(sd,df2)
    a <- c("genes","sd")
    df2.1 <- subset(df2.1,select = c(a,setdiff(colnames(df2.1),a)))
    logic2 <- as.numeric(as.character(df2.1$P)) <= p.cutoff
    df3 <- subset(df2.1,subset = logic2)
    df3$FDR <- p.adjust(as.numeric(as.character(df3$P)),method = "BH")
    if(is.null(FDR.cutoff)){
      df4 <- df3
    } else {
      logic3 <- df3$FDR <= FDR.cutoff
      df4 <- subset(df3,subset = logic3)
    }

    ## 输出数据
    result <- list(
      cox.genes = as.character(df4$genes),
      cox.result = df4
    )
    L1 <- c(L1,list(result));
    names(L1)[i] <- Fastextra(event.status[i],"[.]",1)
  }

  ## 输出结果
  if(show.music == T){play(music)}
  return(L1)

}





