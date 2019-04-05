
# cox.threshold3 help get a comprimized cut-off between mutiple design object


# expr.matrix# expression object
# design.list# a list of design object.The names of them must be provided
# model # model from List.ModelFromOptimized2
# event.status#colname of status
# event.time#colname of time
# event.lower # lower limit of time
# interval.n  # the cut-off selected parameter.Default is 1000
# cutoff.pval # cut-off of p value
#' @export
cox.threshold3 <- function(expr.matrix,
                           design.list,
                           model,
                           event.status=status,
                           event.time=time,
                           event.lower = event.lower,
                           interval.n  = 1000,
                           cutoff.pval = 0.05){
  ## 加载必要的包
  nd <- c("survival","broom","plyr")
  Plus.library(nd)

  ## 获得risk.score矩阵和cutoff节点
  print("Getting score matrixes and cutoff limits...")
  L1 <- NULL
  for(i in 1:length(design.list)){ #i=1
    design.i <- design.list[[i]]

    ## 对齐
    expr.matrix.i <- expr.matrix[,rownames(design.i)]

    ##获得患者的risk score
    model.i <- model #仅1个model
    genes <- as.character(model.i$ENSEMBL[is.na(model.i$SYMBOL)==F])
    lg1 <- rownames(expr.matrix.i) %in% genes
    expr.matrix.i1 <- subset(expr.matrix.i,subset = lg1)
    if(length(genes) >= 2){expr.matrix.i1 <- expr.matrix.i1[genes,]} #基因一定要对齐！
    #expr.matrix.i1 <- expr.matrix.i[genes,]
    cor <- as.numeric(as.character(model.i$Cor[is.na(model.i$SYMBOL)==F]))
    expr.matrix.i1 <- cbind(cor,expr.matrix.i1)
    expr.matrix.i2 <- apply(expr.matrix.i1,1,function(x) x[1]*x[2:length(x)])
    expr.matrix.i3 <- apply(expr.matrix.i2,1,sum)
    expr.matrix.i3 <- as.data.frame(expr.matrix.i3);
    colnames(expr.matrix.i3)[1] <-"risk.score"

    ##组合数据
    prognosis <- c(event.time,event.status)
    design.i1 <- design.i[,prognosis]
    design.i1 <- cbind(design.i1,expr.matrix.i3)
    colnames(design.i1) <- c("time","status","risk.score")
    design.i1 <- design.i1[!is.na(design.i1[,1]),] #去NA值
    lg1 <- design.i1[,1] > event.lower
    design.i1 <- subset(design.i1,subset = lg1)

    ##cut-off节点
    max.i <- max(design.i1$risk.score)
    min.i <- min(design.i1$risk.score)

    ##
    l <- list(max=max.i,min = min.i,design.i = design.i1)
    L1 <- c(L1,list(l))
    names(L1)[i] <- names(design.list)[i]

  }

  ## 获得cutoff向量
  min1 <- NULL;max1 <- NULL
  for(i in 1:length(L1)){
    min1 <- c(min1,L1[[i]][["min"]])
    max1 <- c(max1,L1[[i]][["max"]])
  }
  min2 <- max(min1);max2 <- min(max1)
  interval <- max2 - min2
  int1 <- interval/interval.n
  cutoff <- c(min2,rep(0,interval.n))#length(cutoff)
  for(j in 1:interval.n){cutoff[j+1] <- cutoff[j]+int1}
  cutoff <- cutoff[c(-1,-(interval.n+1))]

  ##计算某个cutoff值的p值
  P <- NULL
  print("Getting p statistics for every design object...")
  for(i in 1:length(L1)){ # i=1
    design.i2 <- L1[[i]][["design.i"]]
    get1 <- function(cutoff.i){
      rs <- ifelse(design.i2$risk.score >= cutoff.i,1,0);
      design.i3 <- cbind(design.i2,rs)
      sdf <- survdiff(Surv(time, status) ~ rs, data = design.i3)
      p.test <- glance(sdf)
      return(p.test)
    }
    cutoff.a <- as.matrix(cutoff)
    p.test <- adply(cutoff.a,1,get1,.id = "cutoff")
    p.test$cutoff <- cutoff
    p.test <- subset(p.test,select = c("cutoff","statistic","p.value"))
    p.test <- p.test[order(p.test$p.value,decreasing = F),]
    print(paste0(names(design.list)[i],": p statistics has been gotten!"))
    p.test$group <- names(design.list)[i]
    P <- rbind(P,p.test)
  }

  ## 获取系列cut-off值
  get2 <- function(P1){
    p.val <- P1[,"p.value"]
    lg1 <- all(p.val < cutoff.pval)
    d <- data.frame(logic = lg1,mean.p= mean(p.val))
    return(d)
  }
  df1 <- ddply(P,"cutoff",.fun = get2)

  ## 优化cutoff值的结果整理
  cutoff.ok <- df1[df1$logic == T,]
  cutoff.ok <- cutoff.ok[order(cutoff.ok$mean.p,decreasing = F),]

  ## 输出结果
  print(paste0("Get ",length(cutoff.ok$cutoff)," optimized cut-off values for design list."))
  l <- list(
    cutoff = cutoff.ok$cutoff,
    data.cutoff = cutoff.ok,
    meta.data = P
  )
  print("All done!")
  return(l)

}







