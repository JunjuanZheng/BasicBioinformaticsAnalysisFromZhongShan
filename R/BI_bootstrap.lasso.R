


## Usage :bulid a bootstrap lasso model via exprs and design objects.
# dds # DESeq2 dds object.If it was NULL,the expr.matrix should not be null.
# transformation #one of "vst" and "normTransform".
# expr.matrix#a log scale expression matrix.
# design# a design object.
# select# select genes names.
# contrast#the colnames of the contrast
# contrast.control#the control of the contrast
# k=10 #nfolds
# R=2000 #bootstrap
# seed #random seed
# seed.range #the range of random seeds
# each.size #the size of each group in one bootstrap
# optimize.method#one of "min" and "1se".Default is "1se".
#' @export
bootstrap.lasso <- function(dds,
                            transformation = "normTransform",
                            expr.matrix=NULL,
                            design,
                            select,
                            contrast = "N.status",
                            contrast.control = "N0",
                            k=10,
                            R=2000,
                            seed = 2018,seed.range = 1:15000,
                            each.size = 15,
                            optimize.method = "1se",
                            show.music = T){
  ###加载包
  need <- c("survival","survminer","glmnet","parallel","doParallel","stringr","boot","pasilla","BiocParallel","tuneR")
  Plus.library(need)

  ## 并行运算
  n_Cores <- detectCores(logical = F)#检测你的电脑的CPU核数
  cluster_Set <- makeCluster(n_Cores)#进行集群
  registerDoParallel(cluster_Set)

  ###构建淋巴结状态~基因表达矩阵
  if(is.null(expr.matrix) ==T){
    #依赖dds提供基因表达矩阵
    Plus.library("DESeq2")
    if(transformation == "vst"){mt <- assay(vst(dds))} else {if(transformation == "normTransform"){mt <- assay(normTransform(dds))} else print("error!Not such type of transformation!")}
  } else {mt <- expr.matrix}

  ## 矩阵对齐
  mt <- mt[select,rownames(design)]

  ###构建记录某次bootstrap的基因ensembl
  l1 <- c(list(seed=c(1,2,3)),list(random.row=c(1,2,4)),list(genes=c(1,2,4)));
  modelgenes.list <- rep(list(l1),R)
  x.mt <- t(mt);x.mt <- x.mt[rownames(design),] #对齐

  #design的contrast的位置向量
  control.p <- grep(contrast.control,design[,contrast])
  treat.p <- setdiff(1:nrow(design),control.p)
  set.seed(seed);sam <- sample(seed.range,R+1,replace = FALSE);#生成R+1个不重复随机数。

  ##bootstrap
  f2 <- NULL
  for (j in 1:R) {
    #从两组各选出特定数量的病例。
    modelgenes.list[[j]][[1]] <- sam[j]
    set.seed(sam[j]);i <- c(sample(control.p,each.size,replace = F),sample(treat.p,each.size,replace = F));

    #lasso矩阵的构建
    modelgenes.list[[j]][[2]] <- i
    x2 <- x.mt[i,]
    design.i <- design[rownames(x2),]
    p.contrast <- match(contrast,colnames(design.i))
    y2 <- ifelse(design.i[,p.contrast] == contrast.control,0,1)

    #lasso回归
    set.seed(sam[R+1]);cv.fit <- cv.glmnet(x2,y2,family="binomial",nfolds=k,parallel=TRUE);
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
    modelgenes.list[[j]][[3]] <- ensembl
    print(paste0("完成第",j,"次bootstrap,Model:",paste(genesymbol,collapse = "_")))

    ## 模型输出
    fi <- data.frame(Seeds = sam[j],
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







