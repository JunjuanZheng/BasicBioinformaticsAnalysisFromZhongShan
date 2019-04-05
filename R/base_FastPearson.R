

## FastPearson help get cor from expression matrix

# data# a expression matrix
# transposition #Whether transpose the matrix.IF you input a gene expression matrix with rownames gene,you should set transposition=T
# control.markers#If NULL,use all pairs strategy,else use one pair strategy.Only one control marker is supported
# target.markers#when control.marker is not a NULL value,it specified the target variates you want.Default is NULL,which means you wan correlation between control.markers to other variates
# method #One of pearson,kendall and spearman.Default is pearson
# order#whether order the cormatrix with corrplot::corrMatOrder
# parallel# whether use parallel strategy
#' @export
FastPearson <- function(data,
                        transposition = F,
                        control.markers=NULL,
                        target.markers=NULL,
                        method = "pearson",
                        order = F,
                        parallel = F){
  ### 加载必要包
  nd <- c("parallel")
  Plus.library(nd)

  ### data transposition
  expr1 <- as.matrix(data)
  if(transposition==T){
    # 矩阵要转置
    expr2 <- t(expr1)
  } else {
    expr2 <- expr1
  }

  ### 计算某两个基因的相关系数
  get1 <- function(genes){ #genes=n2[1,]
    df1 <- expr2[,genes]
    p1 <- cor(df1,method = method)
    p2 <- p1[genes[1],genes[2]]
    return(p2)
  }

  ### 计算相关系数
  if(is.null(control.markers)){
    print("calculate all pairs...")
    ## 计算all pair
    n1 <- colnames(expr2)
    n2 <- combn(n1,2);n2 <- t(n2)
    colnames(n2) <- c("source","target")

    ## 并行运算计算长型数据的相关系数
    if(parallel == F){
      n3 <- apply(n2,1,get1)
    } else {
      print("Use parallel strategy...")
      ncore <- detectCores()
      cl <- makeCluster(mc <- getOption("cl.cores",ncore));
      clusterExport(cl,c("get1","n2","expr2","method"),envir = environment())
      n3 <- parApply(cl = cl,n2,1,get1)
      stopCluster(cl)
    }
    n2 <- as.data.frame(n2)
    n4 <- cbind(n2,Pearson = n3)
    n4 <- n4[order(n4$Pearson,decreasing = T),]

    ## 获取相关矩阵
    print("Get cormatrix data...")
    options(warn = -1)
    n5 <- Make.cormatrix(data = n4,
                         source.col="source",
                         target.col="target",
                         value="Pearson",
                         report = F)
    options(warn = 0)

    ## 是否排序
    if(order == T){
      print("Order the cormatirx...")
      Plus.library("corrplot")
      ord <- corrMatOrder(n5,order = "AOE")
      n5 <- n5[ord,ord]
    }

    ## 输出结果
    l <- list(
      longdata = n4,
      cormatrix = n5
    )

  } else {
    ## one pair
    if(is.null(target.markers)){
      print(paste("Calculate correlation between",control.markers,"to other variates..."))
      n1 <- cbind(source = control.markers,target = setdiff(colnames(expr2),control.markers))
    } else {
      print(paste("Calculate correlation between",control.markers,"to",paste0(target.markers,collapse = " | ")))
      n1 <- cbind(source = control.markers,target = target.markers)
    }

    ## 并行运算计算长型数据的相关系数
    if(parallel == F){
      n2 <- apply(n1,1,get1)
    } else {
      print("Use parallel strategy...")
      ncore <- detectCores()
      cl <- makeCluster(mc <- getOption("cl.cores",ncore));
      clusterExport(cl,c("get1","n1","expr2","method"),envir = environment())
      n2 <- parApply(cl = cl,n1,1,get1)
      stopCluster(cl)
    }
    n1 <- as.data.frame(n1)
    n3 <- cbind(n1,Pearson = n2)
    n3 <- n3[order(n3$Pearson,decreasing = T),]

    ## 输出结果
    l <- list(
      longdata = n3,
      cormatrix = "Not Available"
    )

  }

  ### 输出结果
  print("All done!")
  return(l)

}







