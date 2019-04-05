

## FastCoxlasso help do a plot-based penalized Cox regression model via lasso algrithm.

# expr.matrix# expression matrix
# design# design object
# select#selected markers or genes
# event.status# colnames of event status
# event.time# colnames of event time
# event.lower#the lower limit of event time
# lambda.n# the counts of genes with parameters you want
# ln.lambda#  ln(lambda)
# label# whether to show the genes label on the plot.Default is T
# digits#the digits of model coefficient
# names#part of saved files name

#' @export
FastCoxlasso <- function(
  expr.matrix,
  design,
  select,
  event.status=c("TTR.status","DFS.status","OS.status")[2],
  event.time=c("TTR.time","DFS.time","OS.time")[2],
  event.lower =c(89,89,89),
  lambda.n = 3,ln.lambda = 0,label=TRUE,
  digits = 5,
  names = "love"
){
  ##加载包
  need <- c("MASS","survival","glmnet","stringr")
  Plus.library(need)

  ### 矩阵对齐
  expr.matrix1 <- expr.matrix[select,rownames(design)]
  expr.matrix1 <- t(expr.matrix1)

  ###形成cox模型用的数据
  cox.data <- cbind(design,expr.matrix1)

  ###提取某种预后的信息，step逐步筛选，得到最佳模型。
  l <- list();#i=1
  for(i in 1:length(event.status)){
    n.i <- Fastextra(event.time[i],"[.]",1)
    prognosis <- c(event.time[i],event.status[i])
    ## 去除非空值
    cox.data1 <- cox.data[!is.na(cox.data[,prognosis[1]]),]
    cox.data1 <- cox.data1[!is.na(cox.data1[,prognosis[2]]),]
    ## 数值型
    time <- as.numeric(as.character(cox.data1[,prognosis[1]]))
    status <- as.numeric(as.character(cox.data1[,prognosis[2]]))
    cox.data1 <- subset(cox.data1,select=select)
    cox.data1 <- cbind(time,status,cox.data1)
    ## 过滤数据
    cox.data1 <- cox.data1[time > event.lower[i],]
    ## 建立formula
    y <- subset(cox.data1,select = c("time","status"))
    x <- subset(cox.data1,select = setdiff(colnames(cox.data1),c("time","status")))
    y <- as.matrix(y)
    x <- as.matrix(x)

    ###通过lasso系列获得优化的模型
    cv.fit <- glmnet(x,y,family="cox")

    ### 图片展示
    win.graph(width = 10,height = 10)
    par(cex=1.2) #全局设置
    plot(cv.fit,
         xvar = "lambda",
         label=label,lwd=2,
         xlab = paste0(n.i,"_log lambda"))
    abline(v=ln.lambda,lty=1,lwd=2,col=mycolor[1])
    par(las = 0) #默认设置

    ### 图片保存
    pdf(paste0(names,"_",n.i,"_Penalized Cox regression model.pdf"),width = 8,height = 8)
    par(cex=1.2) #全局设置
    plot(cv.fit,
         xvar = "lambda",
         label=label,lwd=2,
         xlab = paste0(n.i,"_log lambda"))
    abline(v=ln.lambda,lty=1,lwd=2,col=mycolor[1])
    par(las = 0) #默认设置
    plot(cv.fit,
         xvar = "norm",
         label=label,lwd=2,
         xlab = paste0(n.i,"_Norm"))
    abline(v=lambda.n,lty=1,lwd=2,col=mycolor[1])
    dev.off()

    ###提取模型基因和相关系数
    Coefficients <- coef(cv.fit, s = exp(ln.lambda))#coef函数来提取回归模型的系数
    m1 <- as.matrix(Coefficients)
    Active.Index <- which(m1 != 0)#输出基因所在的位置
    Active.Coefficients <- Coefficients[Active.Index,]
    Active.Coefficients <- Active.Coefficients[order(abs(Active.Coefficients),decreasing = T)]
    Active.Coefficients <- Active.Coefficients[1:lambda.n]
    Active.Coefficients <- sort(Active.Coefficients,decreasing = T)
    ensembl <- names(Active.Coefficients);ensembl
    type <- convert(ensembl,totype = "gene_type");type
    genesymbol <- convert(ensembl);genesymbol
    model <- data.frame(
      ENSEMBL = ensembl,
      SYMBOL = genesymbol,
      GENE_TYPE = type,
      Cor = round(Active.Coefficients,digits = digits)
    )

    ## 结果整理1
    all1 <- paste0(paste(model$Cor,model$ENSEMBL,sep = " × expression of "),collapse = " + ")
    all2 <- paste0(paste(model$Cor,model$SYMBOL,sep = " × expression of "),collapse = " + ")
    all <- data.frame(
      ENSEMBL = c("summary_ENSEMBL","summary_SYMBOL"),
      SYMBOL = c(NA,NA),
      GENE_TYPE = c(NA,NA),
      Cor = c(all1,all2)
    )
    model2 <- rbind(model,all)
    rownames(model2) <- 1:nrow(model2)
    model2 <- subset(model2,select = c("ENSEMBL","SYMBOL","GENE_TYPE","Cor"))

    #汇总并输出
    l.i <- list(
      modeldata = model2,
      coxdata=cox.data1)
    l <- c(l,list(l.i));
    names(l)[i] <- n.i
  }

  ### 输出结果
  return(l)
}





