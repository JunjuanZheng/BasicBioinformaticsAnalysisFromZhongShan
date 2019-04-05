
# getexpr help get long data frame from a gene expression matrix.

# exprs.matirx # genes expression matrix or a matrix with individual colnames and marker rownames.
# design#design object
# select#selected markers or genes
# contrast.col#the colnames of the contrast
# contrast.list #the elements of a contrast like Np and N0. If NULL,use all the contrast element.
#' @export
getexpr <- function(exprs.matirx,
                    design,
                    select,
                    contrast.col,
                    contrast.list=NULL){

  ## 提取数据
  mt1 <- exprs.matirx[select,]
  ind.names <- colnames(mt1) #记录个体名称
  colnames(mt1) <- as.character(design[,contrast.col])
  mt1 <- as.matrix(mt1)

  ## 按contrast.list缩小范围
  if(is.null(contrast.list)){
    contrast.list = as.character(unique(design$condition))
  } else {
    contrast.list = contrast.list
  }
  logic1 <- colnames(mt1) %in% contrast.list
  mt2 <- mt1[,logic1]
  ind.names <- ind.names[logic1]

  ## 根据select将宽型矩阵变成长型数据
  l1 <- NULL
  for(i in 1:ncol(mt2)){
    l1.i <- mt2[,i]
    l1 <- c(l1,list(l1.i))
    names(l1)[i] <- colnames(mt2)[i]
  }
  df1 <- stack(l1)
  df1$genes <- rep(select,ncol(mt2))

  ## 加入个体名
  ind1 <- NULL
  for(i in 1:length(ind.names)){
    ind.i <- rep(ind.names[i],length(select))
    ind1 <- c(ind1,ind.i)
  }
  df1$ind1 <- ind1
  colnames(df1) <- c("values","contrast","genes","indivitual")

  ## 输出结果
  return(df1)

}










