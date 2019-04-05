



###=============uniqueModel==============###
#概览模型。modeldata:f2类
# model1 <- model$DFS$modeldata;View(model1)
# model2 <- uniqueModel(model1) ;View(model2)
#' @export
uniqueModel <- function(modeldata){
  # str(modeldata)
  modelfreq <- as.data.frame(table(modeldata$Model_ENSEMBL,modeldata$Cor))
  modelfreq <- modelfreq[order(modelfreq$Freq,decreasing = T),]
  modelfreq <- modelfreq[modelfreq$Freq != 0,] #

  #种子记录
  seeds <- NULL
  for(i in 1:nrow(modelfreq)){ # i=1
    id.i <- as.character(modelfreq[i,1])
    cor.i <- as.character(modelfreq[i,2])
    a1 <- as.character(modeldata$Model_ENSEMBL) %in%  id.i
    a2 <- as.character(modeldata$Cor) %in% cor.i
    a <- apply(cbind(a1,a2),1,is.all.true)
    seed.i <- paste0(modeldata$Seeds[a],collapse = "_")
    seeds <- c(seeds,seed.i)
  }
  modelfreq$seeds <- seeds

  #Model基因数
  genescounts <- function(vt,split="_"){
    gc <- NULL
    for (i in vt) {
      i.c <- length(unlist(strsplit(i,split)))
      gc <- c(gc,i.c)
    }
    gc
  }
  modelfreq$genecounts <- genescounts(modelfreq$Var1)

  #Model相关系数
  #p <- NULL
  #for (i in as.character(modelfreq$Var1)) {
  #  p.i <- match(i,modeldata$Model_ENSEMBL)
  #  p <- c(p,p.i)
  #};p
  #modelfreq$Cor <- modeldata$Cor[p]

  #改名
  colnames(modelfreq)[1:2] <- c("Models_ENSEMBL","Cor")
  return(modelfreq)
}#汇总各模型的数据


###=============oneModel==============###
#观察某个模型的具体信息。输入为uniqueModel结果的某列数据。
#' @export
oneModel <- function(modeldata,dig=5){
  model1 <- NULL;
  model1$ENSEMBL <- unlist(strsplit(as.character(modeldata[1,1]),"_"))
  if(length(model1$ENSEMBL) == 0){model1$ENSEMBL<- NA}
  model1$SYMBOL <- convert(model1$ENSEMBL)
  model1$GENE_TYPE <- convert(model1$ENSEMBL,totype = "gene_type")
  if(length(model1$SYMBOL) == 0){model1$SYMBOL<- NA;model1$GENE_TYPE <- NA}
  model1$Cor <- as.numeric(unlist(strsplit(as.character(modeldata[1,"Cor"]),"_")));model1$Cor <- round(model1$Cor,dig)
  if(length(model1$Cor) == 0){model1$Cor<- NA}
  model1 <- as.data.frame(model1)
  model1 <- model1[order(abs(model1$Cor),decreasing = T),]
  all1 <- paste0(paste(model1$Cor,model1$ENSEMBL,sep = " × expression of "),collapse = " + ")
  all2 <- paste0(paste(model1$Cor,model1$SYMBOL,sep = " × expression of "),collapse = " + ")
  all <- data.frame(ENSEMBL = c("summary_ENSEMBL","summary_SYMBOL"),
                    SYMBOL = c(NA,NA),
                    GENE_TYPE = c(NA,NA),
                    Cor = c(all1,all2))
  model1 <- rbind(model1,all)
  model1 <- subset(model1,select = c("ENSEMBL","SYMBOL","GENE_TYPE","Cor"))
  return(model1)
}


## get results from lasso model function

# lassoModel # the result from cox.optimized2 and bootstrap.lasso series
# position # the position of one model you want
# names # the name of one model like DFS
# dig # digits of correlation of model genes
#' @export
result.lassoModel <- function(lassoModel,
                              position = c(2,1,1),
                              names = c("TTR","DFS","OS"),
                              dig = 5){
  ### 判断model属于bootstrap.lasso还是cox.optimized2
  lg1 <- "coxdata" %in% names(lassoModel[[1]])
  if(lg1){
    #说明是cox.optimized的结果
    print("Results from cox.optimized2...")
    model0 <- list()
    for(i in 1:length(lassoModel)){
      l.i <- lassoModel[[i]][["modeldata"]]
      model0 <- c(model0,list(l.i))
    }
  } else {
    #说明是bootstrap.lass的结果
    print("Results from bootstrap.lasso...")
    model0 <- c(lasso = list(lassoModel[["modeldata"]]))
  }

  model1 <- NULL;
  for(i in 1:length(model0)){
    model1.i <- model0[[i]]
    model1.i <- uniqueModel(model1.i)
    model1.i <- oneModel(model1.i[position[i],],dig = dig)
    model1 <- c(model1,list(model1.i))
  }
  names(model1) <- names
  print("All done!")
  return(model1)

}


## get model matrix via expression matrix,design object and lasso.model
# expr.matrix#a log scale expression matrix
# design# a design object
# model# result from result.lassoModel
# position # the position of one model you want
#' @export
getModelMatrix <- function(expr.matrix,
                           design,
                           model,
                           cut.off = NULL,
                           position = 1){
  ## 获取model genes
  model.i <- model[[position]]
  if("(Intercept)"  %in%  model.i$ENSEMBL){
    p1 <- match("(Intercept)",model.i$ENSEMBL)
    p2 <- setdiff(1:nrow(model.i),p1)
    p <- c(p1,p2)
    model.i <- model.i[p,]
  }
  genes <- as.character(model.i$ENSEMBL[!is.na(as.character(model.i$SYMBOL))])

  ## model genes matrix
  expr.matrix <- expr.matrix[,rownames(design)]
  lg1 <- rownames(expr.matrix) %in% genes
  expr.matrix1 <- subset(expr.matrix,subset = lg1)
  if(length(genes) >= 2){expr.matrix1 <- expr.matrix1[genes,]} #基因一定要对齐！
  if("(Intercept)"  %in%  model.i$ENSEMBL){
    intercept <- rep(1,ncol(expr.matrix1))
    x2 <- rbind(intercept, expr.matrix1)
    rownames(x2)[1] <- "(Intercept)"
  } else {
    x2 <-  expr.matrix1
  }

  ## model scores
  cor <- as.numeric(as.character(model.i$Cor[1:(nrow(model.i)-2)]))
  x3 <- cbind(cor,x2)
  x4 <- apply(x3,1,function(x) x[1]*x[2:length(x)])
  x5 <- apply(x4,1,function(x)sum(x))
  x5 <- as.data.frame(x5);colnames(x5)[1] <-"score"

  ## score status
  y <- t(expr.matrix1)
  if(is.null(cut.off)){
    result <- cbind(design,y,x5)
  } else {
    result <- cbind(design,y,x5)
    result$score.status <- ifelse(as.numeric(as.character(x5$score)) >= cut.off,"high.risk","low.risk")
  }

  ## 输出结果
  result1 <- list(
    modelgenes = genes,
    metadata = result
  )

  return(result1)
}




