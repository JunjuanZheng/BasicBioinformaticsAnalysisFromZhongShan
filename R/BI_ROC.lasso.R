


## Usage :ROC.lasso bulid a binomial lasso model via exprs and design objects.

# expr.matrix#a log scale expression matrix
# design# a design object
# select# select genes names
# contrast#the colnames of the contrast
# contrast.control#the control of the contrast
# k#nfolds
# R#the number of bootstrap
# seed #random seed
# seed.range #the range of random seeds
# optimize.method#one of "min" and "1se".Default is "min"
# whether show music at the end of the programe.Default is FALSE
#' @export
ROC.lasso <- function(expr.matrix,
                      design,
                      select,
                      contrast = "N.status",
                      contrast.control = "N0",
                      k=10,
                      R=100,
                      seed = 2018,seed.range = 2500:5000,
                      optimize.method = "min",
                      show.music = F){
  ###加载包
  need <- c("survival","survminer","glmnet","parallel","doParallel","stringr","boot","pasilla","BiocParallel","tuneR")
  Plus.library(need)

  ## 并行运算
  n_Cores <- detectCores(logical = F)#检测你的电脑的CPU核数
  cluster_Set <- makeCluster(n_Cores)#进行集群
  registerDoParallel(cluster_Set)

  ## 矩阵对齐
  mt <- expr.matrix[select,rownames(design)]

  ## 随机种子
  set.seed(seed);seeds <- sample(seed.range,R,replace = F)

  ## lasso矩阵的构建
  x2 <- t(mt)
  p.contrast <- match(contrast,colnames(design))
  y2 <- ifelse(design[,p.contrast] == contrast.control,0,1)

  ## lasso回归
  f2 <- NULL;modelgenes.list <- NULL
  for (j in 1:R) {
    #lasso回归
    set.seed(seeds[j]);cv.fit <- cv.glmnet(x2,y2,family="binomial",nfolds=k,parallel=TRUE);
    plot(cv.fit)
    ##提取最佳模型的基因
    if(optimize.method == "1se"){s = cv.fit$lambda.1se} else {if(optimize.method == "min"){s = cv.fit$lambda.min} else print("error:not such optimiaze methods!")}#lambda.min是指交叉验证中使得均方误差最小的那个值λ，lambda.1se为离最小均方误差一倍标准差的λ值。
    Coefficients <- coef(cv.fit, s = s)#coef函数来提取回归模型的系数
    m1 <- as.matrix(Coefficients)
    Active.Index <- which(m1 != 0)#输出基因所在的位置
    Active.Coefficients <- Coefficients[Active.Index,]
    ensembl <- names(Active.Coefficients);ensembl
    genesymbol <- convert(ensembl,totype = "SYMBOL")
    type <- convert(ensembl,totype = "gene_type")
    # model.list
    modelgenes.list.j <-list(seed = seed[j],
                             random.row = NULL,
                             genes = ensembl)
    modelgenes.list <- c(modelgenes.list,list(modelgenes.list.j))
    names(modelgenes.list)[j] <- paste0("seed",seeds[j])
    print(paste0("完成第",j,"次bootstrap,Model:",paste(genesymbol,collapse = "_")))

    ## 模型输出
    fi <- data.frame(Seeds = seeds[j],
                     Model_type = paste(type,collapse = "_"),
                     Model_SYMBOL = paste(genesymbol,collapse = "_"),
                     Model_ENSEMBL=paste(ensembl,collapse = "_"),
                     Cor=paste(Active.Coefficients,collapse = "_")
    )
    f2 <- rbind(f2,fi)
  }

  #结束并行
  stopCluster(cluster_Set)

  #输出结果
  if(show.music == T){play(music)}
  l <- list(
    modelgenes = modelgenes.list,
    modeldata = f2
  )
  return(l)
}








