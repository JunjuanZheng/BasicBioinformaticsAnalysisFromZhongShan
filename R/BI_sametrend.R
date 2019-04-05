
## Usage:find out genes that can be plot appropriately in boxplot with student test.

# expr.matrix.list# a list of different gene expression matrix.
# res.list # a list of different results which contain logFC columns.
# design# the design of gene expression matrix
# contrast# the colname of design representing contrast.
# sigdif.co.id #  the range of test genes
# logFC.col# the colnames of results
# count.test# whether to make test like t test and so on.Always T.
# log.convert#whether to make log scale in gene expression
# method# method to do test
# p.cutoff# cut-off of P value

# expr.matrix.list:基因表达矩阵的list
# res.list:基因表达差异的results(dds)
# sigdif.co.id:交叉基因或待测基因
# count.test：F表示不进行counts数的t.test/wilcox.test
# log.convert：T表示对expr.matrix.list进行log2变换。
# method：counts数检验的统计学方法,
# p.cutoff：counts数检验的p值最大值0.05
# design:包含患者的临床信息。某一列的列名为contrast
# contrast：分组变量（如N0与Np）。仅支持1对contrast的对比。
#' @export
sametrend <- function(expr.matrix.list,
                      res.list,
                      logFC.col = "logFC",
                      design,
                      contrast="N.status",
                      sigdif.co.id,
                      count.test = T,
                      log.convert = T,
                      method="wilcox.test",
                      p.cutoff = 0.05){
  p.contrast <- match(contrast,colnames(design))
  library(ggpubr);library(plyr)
  #生成数据记录list
  x <- vector();l1 <- rep(list(x),length(res.list));names(l1) <- names(res.list);l2 <- l1#l1记录logFC,l2记录logFC的符号
  #富集趋势相同的genes
  for (i in 1:length(res.list)) {
    l1[[i]]<- res.list[[i]][sigdif.co.id,logFC.col]
    l2[[i]] <- l1[[i]]/abs(l1[[i]])
  }
  sum <- 0;
  for (i in 1:length(res.list)) {
    sum.i <- l2[[i]]
    sum <- sum+sum.i
  }
  p <- which(sum !=0);p
  #输出趋势相同的基因或进行下一步操作
  select0 = sigdif.co.id[p];
  if(count.test ==F){select <- list(select0)} else {
    ##需要进行log数据的t.test或者是wilcox.test检验等。
    #是否要进行矩阵的log2变换
    if(log.convert == F){l3 <- expr.matrix.list} else {
      as.log2 <- function(expr.matrix.list){
        for(i in 1:length(expr.matrix.list)){
          expr.matrix.list[[i]] <- log2(expr.matrix.list[[i]]+1)
        }
        return(expr.matrix.list)
      }
      l3 <- as.log2(expr.matrix.list)}
    #对于l3中的某res中的某基因，进行contrast的对比并记录p值
    as.df <- function(l3,gene.i){
      df <- NULL;
      for (j in 1:length(l3)) {
        ID=colnames(l3[[j]])
        df.i <- data.frame(ID=ID,value = as.numeric(l3[[j]][gene.i,]),gene = rep(gene.i,length(l3[[j]][gene.i,])),group = design[ID,contrast])
        x <- compare_means(value~group, data=df.i, method = method)
        df.i2 <- data.frame(gene = gene.i,p.value=as.numeric(x[1,4]),method = method,res = names(l3)[j])
        df <- rbind(df,df.i2)
      }
      return(df)
    }
    df <- NULL;
    for(i in 1:length(select0)){
      gene.i <- select0[i]
      df.i <- as.df(l3,gene.i)
      df <- rbind(df,df.i)
    }
    #将df中的每个res对应的基因的p值找到，判断是否都低于p.cutoff值。如果是，则输出此基因
    select.gene <- function(df,p.cutoff){
      #将df按每个res进行分组
      aslist <- function(df){
        l4 <- list()
        for (i in 1:length(unique(df$res))){
          res.i <- as.character(unique(df$res)[i])
          df.i <- df[df$res == res.i,]
          l4 <- c(l4,list(df.i))
        }
        return(l4)
      }
      l4 <- aslist(df)

      l5 <- list();
      for (i in 1:length(l4)) {
        l5.i <- l4[[i]]$p.value
        l5 <- c(l5,list(l5.i))
      }
      df.2 <- as.data.frame(l5)
      is <- apply(df.2,1,function(x)ifelse(F %in% (as.numeric(x)<p.cutoff),F,T))
      select <- as.character(l4[[1]]$gene[is])
      return(select)
    }
    select1 <- select.gene(df,p.cutoff)
    select <- c(list(select1),list(df))
    names(select) <- c("selectgene","pvalue.data")
  }
  return(select)
}



