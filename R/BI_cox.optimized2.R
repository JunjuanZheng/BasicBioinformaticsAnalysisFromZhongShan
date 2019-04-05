
###==============cox.optimized2()=================###
##作用：对给定的基因进行模型的优化。
##原理：利用lasso回归进行筛选预后相关模型。输出结果包含modeldata和coxdata,其中modeldata可以用lasso系列的uniqueModel/oneModel函数进行后续分析，coxdata则包含生存资料及基因表达矩阵的信息。结果汇总于List文件中。

#' @title A pipline of lasso cox modeling
#' @description  cox.optimized2 use glmnet::cv.glmnet to do optimized selection of lasso cox models based on random seeds
#' @param expr.matrix gene expression with sample cols and genes/markers rows.
#' @param design design object with characters cols and sample rows.
#' @param select intersting genes/markers
#' @param event.status a vector of event status names
#' @param event.time a vector of event time names
#' @param event.lower the cutoff(>) ofevent.time
#' @param k number of folds.Default is 10.
#' @param seed a number for randomization
#' @param seed.range the range of randomization
#' @param R the round of randomization
#' @param optimize.method Default is "min".You can also use "1se"
#' @param show.music whether to show music at the end of the process
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom glmnet cv.glmnet predict.cv.glmnet
#' @seealso \code{\link[glmnet]{glmnet}};\code{\link[glmnet]{cv.glmnet}};
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#'
#' # Get lots of optimized model
#' model <- cox.optimized2(expr.matrix = x,
#'                         design = y,
#'                         select = cluster,
#'                         event.status="DFSoutcome",
#'                         event.time="DFS",
#'                         event.lower =0,
#'                         k=10,
#'                         seed = 2018,
#'                         seed.range = 1:2000,R=100,
#'                         optimize.method = "min",
#'                         show.music = F)
#'
#'# select models of given status and time
#' model1 <- model$DFS$modeldata;View(model1)
#'
#'# models statistics
#' model2 <- uniqueModel(model1) ;View(model2)
#'
#'# Visualize the best model
#' res <- exampleLassoCox(expr.matrix = x,
#'                       design = y,
#'                       select = cluster,
#'                       event.status="DFSoutcome",
#'                       event.time="DFS",
#'                       event.lower =0,
#'                       k=10,
#'                       seed = 1601,#某个seed.来自model2
#'                       optimize.method = "min",
#'                       verbose = T,
#'                       save.file = T,
#'                       names = "总GIST")
#'
#'# select best model
#' mod <- oneModel(model2[1,],dig=5)
#'
#'# KM1
#'th1 <- cox.threshold(expr.matrix = tpm,
#'                     design = design.train,
#'                     model = list(mod),
#'                     event.status="DFS.status",
#'                     event.time="DFS.time",
#'                     event.lower = 89,
#'                     smethod="LogRank",
#'                     pmethod="HL",
#'                     dig = 5,
#'                     cancertype.cv="STAD",
#'                     file.name = project)
#'
#'cf <- th1$DFS$cunoff
#'
#'# KM2: Giving a specified cut off
#'th2 <- cox.threshold2(expr.matrix = tpm,
#'                      design = design.train,
#'                      model = rep(list(mod),3),
#'                      cut.off = c(cf,cf,cf),
#'                      event.status=c("TTR.status","DFS.status","OS.status"),
#'                      event.time=c("TTR.time","DFS.time","OS.time"),
#'                      event.lower =c(89,89,89),
#'                      dig = 5,
#'                      cancertype.cv="STAD",
#'                      file.name = project)
#' @export
cox.optimized2 <- function(expr.matrix,
                           design,
                           select,
                           event.status=c("TTR.status","DFS.status","OS.status"),
                           event.time=c("TTR.time","DFS.time","OS.time"),
                           event.lower =c(89,89,89),
                           k=10,
                           seed = 2018,seed.range = 1:2000,R=100,
                           optimize.method = "min",
                           show.music = T){
  ##加载包
  #need <- c("MASS","survival","glmnet","stringr","tuneR","survminer","parallel","doParallel","boot","pasilla","BiocParallel")
  #Plus.library(need)

  ##并行运算
  n_Cores <- detectCores(logical = F)#检测你的电脑的CPU核数
  cluster_Set <- makeCluster(n_Cores)#进行集群
  registerDoParallel(cluster_Set)

  ### 矩阵对齐
  expr.matrix1 <- expr.matrix[select,rownames(design)]
  expr.matrix1 <- t(expr.matrix1)

  ###形成cox模型用的数据
  cox.data <- cbind(design,expr.matrix1)

  ###随机种子
  set.seed(seed);seed1 <- sample(seed.range,R,replace = F)

  ###提取某种预后的信息，step逐步筛选，得到最佳模型。
  l <- list();
  for(i in 1:length(event.status)){ # i=1
    prognosis <- c(event.time[i],event.status[i])
    ## 去除空值
    cox.data1 <- cox.data[!is.na(cox.data[,prognosis[1]]),]
    cox.data1 <- cox.data1[!is.na(cox.data1[,prognosis[2]]),]
    ## 数值型
    time <- as.numeric(as.character(cox.data1[,prognosis[1]]))
    status <- as.numeric(as.character(cox.data1[,prognosis[2]]))
    cox.data1 <- subset(cox.data1,select=select)
    cox.data1 <- cbind(time,status,cox.data1)
    # View(cox.data1[,c(1,2)])
    ## 过滤数据
    cox.data1 <- cox.data1[time > event.lower[i],]
    ## 建立formula
    y <- subset(cox.data1,select = c("time","status"))
    x <- subset(cox.data1,select = setdiff(colnames(cox.data1),c("time","status")))
    y <- as.matrix(y)
    x <- as.matrix(x)


    ###通过lasso系列获得优化的模型
    f2 <- NULL;
    for(j in seed1[1:length(seed1)]){# j=673
      set.seed(j);foldid <- sample(rep(seq(k),length=nrow(x)))#用以进行交叉验证的k向量。

      cv.fit <- cv.glmnet(x,y,family="cox",foldid=foldid,parallel=TRUE);
      # cv.fit <- cv.glmnet(x,y,family="cox",foldid=foldid)
      plot(cv.fit)
      ##提取最佳模型的基因
      if(optimize.method == "1se"){s = cv.fit$lambda.1se} else {if(optimize.method == "min"){s = cv.fit$lambda.min} else print("error:not such optimiaze methods!")}#lambda.min是指交叉验证中使得均方误差最小的那个值λ，lambda.1se为离最小均方误差一倍标准差的λ值。
      Coefficients <- coef(cv.fit, s = s)#coef函数来提取回归模型的系数
      m1 <- as.matrix(Coefficients)
      Active.Index <- which(m1 != 0)#输出基因所在的位置
      Active.Coefficients <- m1[Active.Index,]
      names(Active.Coefficients) <- rownames(m1)[Active.Index]
      ensembl <- names(Active.Coefficients);ensembl
      type <- convert(ensembl,fromtype = "ENSEMBL",totype = "gene_type");type
      genesymbol <- convert(ensembl);genesymbol
      p <- match(j,seed1)
      fi <- data.frame(Seeds = j,
                       Model_type = paste(type,collapse = "_"),
                       Model_SYMBOL = paste(genesymbol,collapse = "_"),
                       Model_ENSEMBL=paste(ensembl,collapse = "_"),
                       Cor=paste(Active.Coefficients,collapse = "_"),
                       stringsAsFactors = F
      )
      f2 <- rbind(f2,fi)
      print(str_c("完成第",as.character(p),"次筛选","种子为",j,",Model:",paste(genesymbol,collapse = "_")))
    }

    #汇总并输出
    l.i <- list(modeldata = f2,coxdata=cox.data1)
    l <- c(l,list(l.i));
    names(l)[i] <- unlist(strsplit(prognosis[1],"[.]"))[1]
  }

  ###输出结果
  stopCluster(cluster_Set)
  if(show.music ==T){mymusic(2);return(l)} else {return(l)}
}


###==============cox.threshold()=================###
##作用：对给定的oneModel()获得最优cut-off值。此函数输出包含cutoff及survdata，其中survdata为每个患者的riskscore及生存数据。结果包装在list里。
##原理：利用maxstat包的方法。
# expr.matrix
# design
# event.status=c("TTR.status","DFS.status","OS.status")
# event.time=c("TTR.time","DFS.times","OS.time")
# event.lower =c(89,89,89)
# smethod="LogRank",为maxstat::maxstat.test函数的参数
# pmethod="HL",为maxstat::maxstat.test函数的参数
# model为oneModel()函数的输出结果。
#' @export
cox.threshold <- function(expr.matrix,
                          design,
                          model,
                          event.status=c("TTR.status","DFS.status","OS.status"),
                          event.time=c("TTR.time","DFS.times","OS.time"),
                          event.lower =c(89,89,89),
                          smethod="LogRank",
                          pmethod="HL",
                          dig = 5,
                          cancertype.cv="STAD",
                          file.name = "1"){
  ###加载包
  library("maxstat");library("survival");library(ggpubr);library(stringr);library(survminer)

  ### 对齐
  expr.matrix <- expr.matrix[,rownames(design)]

  ###对于某种event，提取相应的生存数据。
  l <- list()
  for(i in 1:length(event.status)){ # i=1
    ##获得患者的risk score
    model.i <- model[[i]]
    genes <- as.character(model.i$ENSEMBL[is.na(model.i$SYMBOL)==F])
    lg1 <- rownames(expr.matrix) %in% genes
    expr.matrix1 <- subset(expr.matrix,subset = lg1)
    if(length(genes) >= 2){expr.matrix1 <- expr.matrix1[genes,]} #基因一定要对齐！
    #expr.matrix1 <- expr.matrix[genes,]
    cor <- as.numeric(as.character(model.i$Cor[is.na(model.i$SYMBOL)==F]))
    expr.matrix1 <- cbind(cor,expr.matrix1)
    expr.matrix2 <- apply(expr.matrix1,1,function(x) x[1]*x[2:length(x)])
    expr.matrix3 <- apply(expr.matrix2,1,function(x)sum(x))
    expr.matrix3 <- as.data.frame(expr.matrix3);
    colnames(expr.matrix3)[1] <-"risk.score"

    ##组合数据
    design.i <- design[rownames(expr.matrix3),]
    prognosis <- c(event.time[i],event.status[i])
    design.i1 <- design.i[,prognosis]
    design.i1 <- cbind(design.i1,expr.matrix3)
    colnames(design.i1) <- c("time","status","risk.score")
    design.i1 <- design.i1[is.na(design.i1[,1])==F,]#去除低于event.lower的数据。
    design.i1 <- design.i1[as.numeric(design.i1[,1])>event.lower[i],]#保存数据1

    ##maxstat运算
    mtHL <- maxstat.test(Surv(time, status) ~ ., data=design.i1, smethod=smethod, pmethod=pmethod)
    n1 <- unlist(strsplit(prognosis[1],"[.]"))[1]
    threshold.i <- mtHL[["estimate"]];
    threshold.i <- round(threshold.i,dig=dig)#保存数据2
    score.status <- ifelse(design.i1$risk.score >= threshold.i,"high.risk","low.risk")
    x1 <- data.frame(score.status = score.status)
    design.i1 <- cbind(design.i1,x1)

    ##画图
    ##cut-off图
    pdf(str_c(file.name,"_threshold of risk scores_",n1,"_",cancertype.cv,".pdf"),width = 8,height = 8)
    plot(mtHL);
    #text(x=threshold.i,y=1.03*max(mtHL[["stats"]]),labels = paste("cut off =",as.character(threshold.i),collapse = ""),pos = 1,cex=1)
    text(x=threshold.i,y=max(mtHL[["stats"]]),labels = as.character(threshold.i),pos = 2,cex=1)
    dev.off()
    print(str_c("完成",n1,"_threshold of risk scores","图的绘制"))

    #频数分布直方图
    p1 <- gghistogram(data=design.i1,
                      x = "risk.score",
                      rug = TRUE,
                      title = "Distribution", xlab = "Score", ylab = "Number of Patients",
                      fill = "score.status", palette = c("#00AFBB", "#E7B800")) +
      theme(axis.title= element_text(size = 14,face = "bold"),
            axis.text = element_text(size = 12,face = "bold"),
            legend.title=element_blank(),
            #legend.position = "right",
            legend.text = element_text(size = 12,face = 'bold'),
            plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
    ggsave(str_c(file.name,"_distribution of risk scores_",n1,"_",cancertype.cv,".pdf"),p1,height = 8,width = 8)
    print(str_c("完成",n1,"_distribution of risk scores","图的绘制"))

    #预后K-M图
    fit <- survfit(Surv(time,status)~score.status, data=design.i1)
    p.value <- survdiff(Surv(time,status)~score.status,data=design.i1)
    library(broom)
    p1 <- glance(p.value)$p.value
    label <- paste0("p = ",as.character(signif(p1,3)))
    suvival <- ggsurvplot(fit,
                          pval = F, conf.int = F,
                          risk.table = F,
                          linetype = "strata",
                          surv.median.line = "hv",
                          palette = c("#E7B800", "#2E9FDF"),
                          ncensor.plot=F);
    p2 <- suvival$plot +
      labs(title = n1) +
      annotate("text",x=max(design.i1$time),y=1,label =label)+
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(str_c(file.name,"_KMcurves_",n1,"_",cancertype.cv,".pdf"),p2,height = 8,width = 8)
    print(str_c("完成",n1,"_KMcurves","图的绘制"))

    ##输出数据。
    l.i <- list(cunoff = threshold.i,survdata = design.i1)
    l <- c(l,list(l.i))
    names(l)[i] <- n1
  }
  return(l)
}


##作用：对给定的模型基因和cut.off值，进行分布分析和生存分析。
#' @export
cox.threshold2 <- function(expr.matrix,
                           design,
                           model,
                           cut.off,
                           event.status=c("TTR.status","DFS.status","OS.status"),
                           event.time=c("TTR.time","DFS.times","OS.time"),
                           event.lower =c(89,89,89),
                           dig = 5,
                           cancertype.cv="STAD",
                           file.name = "test1"){
  ###加载包
  library("survival");library(ggpubr);library(stringr);library(survminer)
  ###对齐
  expr.matrix <- expr.matrix[,rownames(design)]

  ###对于某种event，提取相应的生存数据。
  l <- list()
  for(i in 1:length(event.status)){ # i=2
    ##获得患者的risk score
    model.i <- model[[i]]
    genes <- as.character(model.i$ENSEMBL[is.na(model.i$SYMBOL)==F])
    lg1 <- rownames(expr.matrix) %in% genes
    expr.matrix1 <- subset(expr.matrix,subset = lg1)
    if(length(genes) >= 2){expr.matrix1 <- expr.matrix1[genes,]} #基因一定要对齐！
    #expr.matrix1 <- expr.matrix[genes,]
    cor <- as.numeric(as.character(model.i$Cor[is.na(model.i$SYMBOL)==F]))
    expr.matrix1 <- cbind(cor,expr.matrix1)
    expr.matrix2 <- apply(expr.matrix1,1,function(x) x[1]*x[2:length(x)])
    expr.matrix3 <- apply(expr.matrix2,1,function(x)sum(x))
    expr.matrix3 <- as.data.frame(expr.matrix3);
    colnames(expr.matrix3)[1] <-"risk.score"

    ##组合数据
    design.i <- design[rownames(expr.matrix3),]
    prognosis <- c(event.time[i],event.status[i])
    design.i1 <- design.i[,prognosis]
    design.i1 <- cbind(design.i1,expr.matrix3)
    colnames(design.i1) <- c("time","status","risk.score")
    design.i1 <- design.i1[is.na(design.i1[,1])==F,]#去除低于event.lower的数据。
    design.i1 <- design.i1[as.numeric(design.i1[,1])>event.lower[i],]#保存数据1

    ##根据cut-off值进行画图。
    n1 <- Fastextra(event.status[i],"[.]",1)
    threshold.i <- cut.off[i]
    score.status <- ifelse(design.i1$risk.score >= threshold.i,"high.risk","low.risk")
    x1 <- data.frame(score.status = score.status)
    design.i1 <- cbind(design.i1,x1)

    ##画图
    #频数分布直方图
    p1 <- gghistogram(data=design.i1,
                      x = "risk.score",
                      rug = TRUE,
                      title = "Distribution", xlab = "Score", ylab = "Number of Patients",
                      fill = "score.status", palette = c("#00AFBB", "#E7B800")) +
      theme(axis.title= element_text(size = 14,face = "bold"),
            axis.text = element_text(size = 12,face = "bold"),
            legend.title=element_blank(),
            #legend.position = "right",
            legend.text = element_text(size = 12,face = 'bold'),
            plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
    ggsave(str_c(file.name,"_distribution of risk scores_",n1,"_",cancertype.cv,".pdf"),p1,height = 8,width = 8)
    print(str_c("完成",n1,"_distribution of risk scores","图的绘制"))

    #预后K-M图
    fit <- survfit(Surv(time,status)~score.status, data=design.i1)
    p.value <- survdiff(Surv(time,status)~score.status,data=design.i1)
    library(broom)
    p1 <- glance(p.value)$p.value
    label <- paste0("p = ",as.character(signif(p1,3)))
    suvival <- ggsurvplot(fit,
                          pval = F, conf.int = F,
                          risk.table = F,
                          linetype = "strata",
                          surv.median.line = "hv",
                          palette = c("#E7B800", "#2E9FDF"),
                          ncensor.plot=F);
    p2 <- suvival$plot +
      labs(title = n1) +
      annotate("text",x=max(design.i1$time),y=1,label =label)+
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(str_c(file.name,"_KMcurves_",n1,"_",cancertype.cv,".pdf"),p2,height = 8,width = 8)
    print(str_c("完成",n1,"_KMcurves","图的绘制"))

    ##输出数据。
    l.i <- list(cunoff = threshold.i,survdata = design.i1)
    l <- c(l,list(l.i))
    names(l)[i] <- n1
  }
  return(l)
}




