

####=====================getAnnotation2()=====================####
## getAnotation2()在getAnotation1()的基础上，利用org.Hs.eg.db的注释，获得类似于ENTREZID/UNIGENE等的注释信息。一般在不同的数据取交集基因时，geo.type=T有助于将UNIGENE当作SYMBOL来处理。
#' Get further annotation after getAnnotation()
#' @description getAnotation2()在getAnotation1()的基础上，利用org.Hs.eg.db的注释，获得类似于ENTREZID/UNIGENE等的注释信息。一般在不同的数据取交集基因时，geo.type=T有助于将UNIGENE当作SYMBOL来处理。
#' @param data a result from getAnnotation()
#' @param geo.type Default=F.If T,it considers UNIGENE as SYMBOL,which is not recommanded
#' @return a data frame
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @seealso \code{\link{getAnnotation}}
#' @examples
#' df2.1 <- getAnnotation2(gtf.annotation,geo.type = T)
#' df2.2 <- getAnnotation2(gtf.annotation,geo.type = F)
#' @export
getAnnotation2 <- function(data,geo.type=F){
  ##加载包
  nd <- c("annotate","org.Hs.eg.db","dplyr","stringr","AnnotationDbi")
  Plus.library(nd)

  ##网络信息
  test.ano <- data
  ensids <- as.character(test.ano$ENSEMBL)
  col = c("SYMBOL","ENTREZID","UNIGENE","GENENAME")
  print(str_c("从org.Hs.eg.db获得",paste0(col,collapse = "/"),"的对应信息。"))
  df <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = ensids,
                              keytype = "ENSEMBL",
                              columns = col);
  df1 <- dplyr::left_join(df,test.ano,by = "ENSEMBL")

  ##此时ENSEMBL是重复的。进行智能合并。
  #merge.anno执行两个步骤：对于2元素的向量：相同，输出其中一个。不同：有空值，输出非空值；均为非空值，则输出“error”
  merge.anno <- function(vt,gold.key = 2){
    if(length(unique(vt))==1){
      #元素全部相同
      x <- vt[1]
    } else {
      #元素不全相同
      if(NA %in% vt){
        #有空值，将非空值输出
        x <- vt[!is.na(vt)]
      } else {
        #不同的非空值
        x <- vt[gold.key]
      }
    }
    #End
  }
  df1$SYMBOL <- apply(df1[,grep("SYMBOL",colnames(df1))],1,function(x)merge.anno(x,gold.key = 2))
  col2 = c("ENSEMBL","ENSEMBL.Versions","SYMBOL","GENENAME","ENTREZID","UNIGENE","gene_type")
  df2 <- subset(df1,select = col2)

  ## 按需输出
  if(geo.type==T){
    print("整合UNIGENE和SYMBOL。")
    ##在不少geo芯片中，UNIGENE和SYMBOL是混用的。因此，出于实用，我们可以将UNIGENE和SYMBOL统一为SYMBOL信息。
    cols.1 <- colnames(df2) %in% "SYMBOL"
    cols.1 <- colnames(df2)[!cols.1]
    df2.1 <- subset(df2,select = cols.1)
    cols.2 <- colnames(df2) %in% "UNIGENE"
    cols.2 <- colnames(df2)[!cols.2]
    df2.2 <- subset(df2,select = cols.2)
    colnames(df2.1)[match("UNIGENE",colnames(df2.1))] <- "SYMBOL"
    db.annotation2 <- rbind(df2.1,df2.2)
    return(db.annotation2)
  } else {
    print("分离UNIGENE和SYMBOL。")
    #将UNIGENE和SYMBOL分开。
    db.annotation2 <- df2
    return(db.annotation2)
  }
  #End
}



