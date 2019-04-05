

#' @title Data Clean for starBase
#' @description Data Clean for starBase
#' @keywords straBase1
#' @param path the absolute path of seriese file(.xls) download from starbase
#' @param mode Only available for "RNA-RNA".Other is on schedule
#' @return a data frame
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
starBase1 <- function(path,mode = "RNA-RNA"){
  ## 加载必要的包
  nd <- c("plyr")
  Plus.library(nd)

  if(mode == "RNA-RNA"){
    ## RNA-RNA数据融合
    print("staBase_RNA-RNA interaction: Merge data...")
    file_1 <- list.files(path,pattern = "starBase",full.names = T)
    get1 <- function(file_1_i){
      a <- read.table(file_1_i,sep = "\t",header = T,check.names = F)
      return(a)
    } #file_1_i <- file_1[1]
    df1 <- plyr::adply(as.matrix(file_1),1,get1)
    df1 <- df1[,-1]

    ## 去除重复冗余
    get2 <- function(vt){
      vt <- factor(vt)
      vt <- levels(vt)
      vt1 <- paste(vt,collapse = "_")
      return(vt1)
    }
    #vt <- c("a","b");vt <- c("b","a");  get2(vt)
    lg <- apply(cbind(as.character(df1$geneID),
                      as.character(df1$pairGeneID)),
                1,get2)
    test1 <- table(duplicated(lg))
    if(T %in% names(test1)){
      print(paste0("staBase_RNA-RNA interaction: There are ",test1[names(test1) %in% T]," repeat records and would be remove..."))
      df2 <- df1[!duplicated(lg),]
    } else {
      print("staBase_RNA-RNA interaction: No repeat detected...")
      df2 <- df1
    }

    ## 输出结果
    print("All done!")
    return(df2)
  } else {
    print("Please input a right mode.")
  }

}









