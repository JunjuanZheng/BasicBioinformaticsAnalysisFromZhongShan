


###==============cox.optimized()=================###
##作用：对给定的基因进行模型的优化。
##原理：向后逐步回归（backward stepwise）从模型包含所有预测变量开始，一次删除一个变量直到会降低模型质量为止。MASS包中的steAIC（）函数可实现逐步回归模型，依据的是精确AIC准则。
##coxscreen()最终将输出包含多种优化信息的list类型。

#expr.matrix# expression matrix
#design#design object
#select# selected markers or genes
#event.status#colnames of event status
#event.time#colnames of event time
#event.lower =c(89,89,89)#the lower limit of event time
#direction# parameter of MASS::stepAIC.One of "both", "backward", "forward".Default is "both"
#method# whether use PCA strategy to avoid multicollinearity.Default is pca
#digits# the digits of Cor
#' @export
cox.optimized <- function(
  expr.matrix,
  design,
  select,
  event.status=c("TTR.status","DFS.status","OS.status"),
  event.time=c("TTR.time","DFS.time","OS.time"),
  event.lower =c(89,89,89),
  direction="both",
  method = "pca",
  digits = 5
){
  ##加载包
  need <- c("MASS","survival")
  Plus.library(need)

  ###提取基因信息
  expr.matrix1 <- expr.matrix[select,rownames(design)]
  expr.matrix1 <- t(expr.matrix1)

  ###形成cox模型用的数据
  cox.data <- cbind(design,expr.matrix1)

  ###提取某种预后的信息，step逐步筛选，得到最佳模型。
  l <- list()
  for(i in 1:length(event.status)){ #i=3
    event.names <- Fastextra(event.status[i],"[.]",1)
    prognosis <- c(event.time[i],event.status[i])
    ## 去除非空值
    cox.data1 <- cox.data[!is.na(cox.data[,prognosis[1]]),]
    cox.data1 <- cox.data1[!is.na(cox.data1[,prognosis[2]]),]
    ## 数值型
    time <- as.numeric(as.character(cox.data1[,prognosis[1]]))
    status <- as.numeric(as.character(cox.data1[,prognosis[2]]))
    cox.data1 <- subset(cox.data1,select=select)
    cox.data2 <- cbind(time,status,cox.data1)
    ##过滤数据
    cox.data2 <- cox.data2[time > event.lower[i],]

    ### 获得公式
    if(method != "pca"){
      ## 建立formula
      fit <- coxph(Surv(time,status) ~.,data = cox.data2)
      step1 <- stepAIC(fit,direction=direction)
      fit2 <- coxph(step1[["formula"]],data = cox.data2)

      ## 结果整理1
      coef <- coef(fit2)
      coef <- round(coef,digits = digits)
      coef <- as.data.frame(coef)
      coef <- cbind(rownames(coef),coef)
      colnames(coef)[1] <- "ENSEMBL"
      coef$SYMBOL <- convert(coef$ENSEMBL)
      coef$GENE_TYPE <- convert(coef$ENSEMBL,totype = "gene_type")

      ## 结果整理2
      coef2 <- coef
      all1 <- paste0(paste(coef2$coef,coef2$ENSEMBL,sep = " × expression of "),collapse = " + ")
      coef2$SYMBOL <- convert(coef2$ENSEMBL)
      all2 <- paste0(paste(coef2$coef,coef2$SYMBOL,sep = " × expression of "),collapse = " + ")
      all <- data.frame(
        ENSEMBL = c("summary_ENSEMBL","summary_SYMBOL"),
        SYMBOL = c(NA,NA),
        GENE_TYPE = c(NA,NA),
        coef = c(all1,all2)
      )
      coef3 <- rbind(coef2,all)
      rownames(coef3) <- 1:nrow(coef3)
      coef3 <- subset(coef3,select = c("ENSEMBL","SYMBOL","GENE_TYPE","coef"))
      pca <- NULL
      colnames(coef3)[ncol(coef3)] <- "Cor"
    } else {
      print("Use PCA optimized...")
      pca <- princomp(~., data=cox.data1,cor=T);
      print(paste0(event.names,":PCA optimized summary..."))
      print(summary(pca, loadings=TRUE))

      ## 选择主成分
      pre <- predict(pca)
      cox.data3 <- cbind(time,status,pre)
      cox.data3 <- as.data.frame(cox.data3)

      ## 方程优化
      fit <- coxph(Surv(time,status) ~.,data = cox.data3)
      step1 <- stepAIC(fit,direction=direction)
      fit2 <- coxph(step1$formula,data = cox.data3)
      comp <- as.character(fit2$formula)
      comp <- comp[3:length(comp)]

      ### 进行变换，得到线性回归方程中的β值。
      A <- loadings(pca) #特征矩阵
      beta <- coef(fit2) #提取回归系数
      x.bar <- pca$center; x.sd <- pca$scale#中心化和标准化转换；
      ## β值
      coef <- 0
      for(b in 1:length(beta)){
        coef.i <- beta[b]*A[,b]
        coef <- coef + coef.i
      }
      coef <- coef/x.sd
      #beta0 <- 0 - sum(x.bar*coef)

      ## 结果整理1
      coef <- round(coef,digits = digits)
      coef1 <- as.data.frame(coef)
      #beta0 <- data.frame(coef=round(beta0,digits = digits),
      #                    row.names = "(interpret)")
      #coef2 <- rbind(beta0,coef1)
      coef2 <- cbind(rownames(coef1),coef1)
      colnames(coef2)[1] <- "ENSEMBL"
      coef2$GENE_TYPE <- convert(coef2$ENSEMBL,totype = "gene_type")

      ## 结果整理2
      all1 <- paste0(paste(coef2$coef,coef2$ENSEMBL,sep = " × expression of "),collapse = " + ")
      coef2$SYMBOL <- convert(coef2$ENSEMBL)
      all2 <- paste0(paste(coef2$coef,coef2$SYMBOL,sep = " × expression of "),collapse = " + ")
      all <- data.frame(
        ENSEMBL = c("summary_ENSEMBL","summary_SYMBOL"),
        SYMBOL = c(NA,NA),
        GENE_TYPE = c(NA,NA),
        coef = c(all1,all2)
      )
      coef3 <- rbind(coef2,all)
      rownames(coef3) <- 1:nrow(coef3)
      coef3 <- subset(coef3,select = c("ENSEMBL","SYMBOL","GENE_TYPE","coef"))
      colnames(coef3)[ncol(coef3)] <- "Cor"
    }

    ## 汇总并输出
    optimized.genes = as.character(coef3$ENSEMBL[!is.na(coef3$SYMBOL)])
    l.i <- list(optimized.genes = optimized.genes,
                model = coef3,
                optimized.cox = fit2,
                stepAIC.details = step1,
                PCA = pca)
    l <- c(l,list(l.i));
    names(l)[i] <- Fastextra(prognosis[1],"[.]",1)
  }

  ###输出结果
  return(l)

}













