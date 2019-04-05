


#' @export
ROC.optimized <- function(expr.matrix,
                         design,
                         contrast.col,
                         contrast.control,
                         parallel = F){

  ## 加载必要的包
  nd <- c("pROC","plyr")
  Plus.library(nd)

  ## auc计算
  if(parallel == F){
    ## 对齐
    expr.matrix <- expr.matrix[,rownames(design)]

    ## 提取相关信息
    status <- as.character(design[,contrast.col])
    status <- ifelse(status == contrast.control,0,1)

    ## 对某个基因进行数据提取
    #gene.exprs  = expr.matrix[1,]
    extra.gene <- function(gene.exprs){
      gene.exprs <- as.numeric(as.character(gene.exprs))
      gene.exprs <- as.data.frame(gene.exprs)
      a <- cbind(status,gene.exprs)
      roc.a <- roc(status~gene.exprs,data = a)
      auc.a <- as.numeric(roc.a$auc)
      return(auc.a)
    }
    print("The AUC of every gene is on calculation,please wait...")
    ## 对多个基因进行数据提取
    auc <- apply(expr.matrix,1,extra.gene)
    auc <- as.data.frame(auc)
    auc <- cbind(auc,test=0)
    auc <- auc[order(auc$auc,decreasing = T),]
    auc <- subset(auc,select = "auc")
  } else {
    ## 对齐
    RO.expr.matrix <- expr.matrix[,rownames(design)]

    ## environment:RO.extra.gene
    RO.extra.gene <- function(gene.exprs){
      gene.exprs <- as.numeric(as.character(gene.exprs))
      gene.exprs <- as.data.frame(gene.exprs)
      a2 <- cbind(RO.status,gene.exprs)
      roc.a <- roc(RO.status~gene.exprs,data = a2)
      auc.a <- as.numeric(roc.a$auc)
      return(auc.a)
    }

    ## 提取相关信息
    RO.status <- as.character(design[,contrast.col])
    RO.status <- ifelse(RO.status == contrast.control,0,1)

    ## parApply并行
    print("Use parallel calculation...")
    nd <- c("parallel");Plus.library(nd)
    ncore <- detectCores(logical = F)
    cl <- makeCluster(getOption("cl.cores",ncore))
    clusterExport(cl,c("RO.extra.gene","RO.expr.matrix","RO.status","roc"),envir = environment())
    auc <- parApply(cl = cl,RO.expr.matrix,1,RO.extra.gene)
    stopCluster(cl)

    ## 结果整理
    auc <- as.data.frame(auc)
    auc <- cbind(auc,test=0)
    auc <- auc[order(auc$auc,decreasing = T),]
    auc <- subset(auc,select = "auc")
  }

 ## 返回列表
  print("All done!")
  return(auc)

}
