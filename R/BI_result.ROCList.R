


###===================result.ROCList(lasso.by.roc)
##对于bootstrap.lasso的结果进行汇总分析，计算基因与roc面积的关系。
# train.expr.matrix、valid.expr.matrix为基因表达矩阵
# design.train,design.valid为对应的design
# ROCList 上一步的结果
# contrast = "N.status",
# contrast.control指的是对照组的名字
# return.freq =F输出roc筛选图；return.freq =T输出多种基础数据。
# palette为数字向量，代表brewer.pal(12,"Set3")颜色的位置。

## show the progonis value of multiple model genes from the result of bootstrap.lasso function
# expr.matrix.list # the list of comparing expression matrix
# design.list# the list of design object of comparing expression matrix
# ROCList # the result of lucky::bootstrap.lasso function
# contrast # the colnames of contrast
# contrast.control# the names of control object in the contrast.
# return.freq#If F,return roc related plot;If T,return data.
# palette # a integer vector to select color from lucky::mycolor
# save.file # whether save a PDF plot
# width,height# the width and height of saved plot
# names # part of file name

#' @export
result.ROCList <- function(expr.matrix.list,
                           design.list,
                           ROCList,
                           contrast="N.status",
                           contrast.control="N0",
                           palette=c(1,4),
                           save.file = T,
                           width = 12,height=10,
                           names = "love"){
  ###加载包
  library(pROC);

  ###获得模型中基因的频数信息v1
  print("get frequence information from ROCList...")
  v <- NULL
  ROCList <- ROCList$modelgenes
  for(i in 1:length(ROCList)){
    vi <- ROCList[[i]][["genes"]]
    v <- c(v,vi)
  }
  p.intercept <- grep("(Intercept)",v)
  v <- v[-p.intercept]
  v1 <- as.data.frame(table(v))
  colnames(v1)[1] <- "Modelgenes"
  v1$SYMBOL <- convert(v1$Modelgenes)
  v1 <- subset(v1,select=c("Modelgenes","SYMBOL","Freq"))
  v1 <- v1[order(v1$Freq,decreasing = T),]
  rownames(v1) <- 1:nrow(v1)

  ###构建表达矩阵和design,获得valid和train的不同基因的aoc值。
  print("get aoc based on expression matrix and design object,please wait...")
  roc.df <- NULL;#i=1;j=1;
  for(j in 1:length(expr.matrix.list)){
    ## 矩阵对齐
    design.j <- design.list[[j]]
    mt.j <- expr.matrix.list[[j]][,rownames(design.j)]

    #某个集里的信息
    roc.df.j <- NULL;
    for(i in 1:nrow(v1)){
      ### 单基因模型
      # 表达矩阵
      tgene <- as.character(v1$Modelgenes[i])
      x.mt <- as.matrix(mt.j[tgene,])
      colnames(x.mt)[1] <- tgene
      x2 <- x.mt
      #design
      design.i <- design.j[rownames(x2),]
      p.contrast <- match(contrast,colnames(design.i))
      y2 <- ifelse(design.i[,p.contrast] == contrast.control,0,1)
      g1 <- cbind(y2,x2)
      g1 <- as.data.frame(g1)
      #roc
      library(pROC)
      f <- paste0("y2 ~ ",tgene);f <- as.formula(f)
      roc.i <- roc(f,data=g1)
      auc.i <- as.numeric(roc.i$auc)

      ### 多基因模型
      testgenes <- as.character(v1$Modelgenes[1:i])
      if(i==1){
        x.mt <- as.matrix(mt.j[testgenes,])
        colnames(x.mt)[1] <- as.character(testgenes)
        x2 <- x.mt
      } else {x.mt <- mt.j[testgenes,];x2 <- t(x.mt)}
      #design
      design.i <- design.j[rownames(x2),]
      p.contrast <- match(contrast,colnames(design.i))
      y2 <- ifelse(design.i[,p.contrast] == contrast.control,0,1)
      g1 <- cbind(y2,x2)
      g1 <- as.data.frame(g1)
      #logistic回归
      modeROCList <- glm(y2 ~ ., data=g1, family='binomial')#可能是有错的
      pre <- predict(modeROCList)

      modelroc <- roc(g1$y2,pre)
      auc.i2 <- as.numeric(modelroc[["auc"]])
      roc.df.i <- data.frame(genecounts=i,
                             auc=auc.i2,sigle.auc = auc.i,
                             group=names(expr.matrix.list)[j])
      roc.df.j <- rbind(roc.df.j,roc.df.i)
    }
    roc.df.j$genes <- as.character(v1$Modelgenes)
    colnames(roc.df.j)
    roc.df.j <- subset(roc.df.j,select = c("genes","genecounts","sigle.auc","auc","group"))
    roc.df <- rbind(roc.df,roc.df.j)
  }

  ###  auc的可视化
  print("plot auc via ggpubr::ggline...")
  library(ggpubr);library(RColorBrewer)
  plot <- ggline(roc.df, x = "genecounts", y = "auc",
              linetype = "group", shape = "group",
              color = "group",palette = mycolor[palette],
              point.size = 3,
              size = 1)
  print(plot)
  if(save.file == T){
    ggsave(paste0(names,"_ROCfileter.plot.pdf"),plot,width = width,height=height)
    print("AUC plot had been saved at present work space!")
  }

  #找到全部auc达到最高auc的第一个counts。
  count1 <- table(roc.df$genecount[roc.df$auc == max(roc.df$auc)])
  count2 <- NULL
  for(i in 1:length(expr.matrix.list)){
    p <- min(grep(i,count1))
    count2.i <- as.numeric(names(count1[p]))
    count2.i <- data.frame(co = i,oc = count2.i)
    count2 <- rbind(count2,count2.i)
  }
  #l <- list(optimized.count = count2,freq=v1,rocdata = roc.df)

  # 输出结果
  l <- list(optimized.count = count2,
            freq = v1,
            rocdata = roc.df,
            plot = plot)
  print("All done!")
  return(l)

}






