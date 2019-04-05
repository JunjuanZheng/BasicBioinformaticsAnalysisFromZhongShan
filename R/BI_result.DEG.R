

## result.DEG help get a summary result from multiple DEGList
#DEGList # a list of result from DESeq1 or edgeR1
#logFC.list# a list of logFC colnames
#select#selected markers or genes
#filter#whether to do a same trend filter.If you want all results,set FALSE to filter
#' @export
result.DEG <- function(DEGList,
                       logFC.list,
                       select,
                       filter=T){
  ## 加载必要的包
  nd <- c("DESeq2","edgeR","limma")
  Plus.library(nd)

  ## 提取difsig的数据
  dif <- NULL
  for(i in 1:length(DEGList)){ # i=1
    D1 <- DEGList[[i]]
    n.i <- names(DEGList)[i]
    diffSig1 <- D1$diffSig[select,]

    ## 找到logFC的列名
    logFC.col <- logFC.list[[n.i]]
    difSig2 <- subset(diffSig1,select = logFC.col)
    colnames(difSig2)[1] <- paste0(n.i,".logFC")
    difSig3 <- t(difSig2)
    dif <- rbind(dif,difSig3)

  }

  ## 输出结果
  dif2 <- t(dif)
  dif2 <- as.data.frame(dif2)
  if(filter == F){
    # 显示全部结果
    return(dif2)
  } else {
    ## 需要找到趋势相同的基因
    result.DEG2 <- function(resultDEG){
      ## 判断是否趋势相同
      get1 <- function(resultDEG1){
        resultDEG1 <- resultDEG1 > 0
        resultDEG1 <- as.numeric(resultDEG1)
        if(length(unique(resultDEG1))==1){
          ## 趋势相同
          return(T)
        } else {
          return(F)
        }
      }
      genes <- rownames(resultDEG)
      Plus.library("plyr")
      df1 <- adply(resultDEG,1,get1)
      rownames(df1) <- genes;colnames(df1)[3] <- "logic"

      ## 输出结果
      lg1 <- df1$logic %in% T
      df2 <- df1[lg1,]
      print(paste0("The number of all genes is ",nrow(resultDEG),",but only ",nrow(df2)," genes enjoy same trend."))
      return(df2)

    }
    rd2 <- result.DEG2(dif2)
    return(rd2)
  }

}

# resultDEG # the result of result.DEG function
#' @export
result.DEG2 <- function(resultDEG){
  ## 判断是否趋势相同
  get1 <- function(resultDEG1){
    resultDEG1 <- resultDEG1 > 0
    resultDEG1 <- as.numeric(resultDEG1)
    if(length(unique(resultDEG1))==1){
      ## 趋势相同
      return(T)
    } else {
      return(F)
    }
  }
  genes <- rownames(resultDEG)
  Plus.library("plyr")
  df1 <- adply(resultDEG,1,get1)
  rownames(df1) <- genes;colnames(df1)[3] <- "logic"

  ## 输出结果
  lg1 <- df1$logic %in% T
  df2 <- df1[lg1,]
  print(paste0("The number of all genes is ",nrow(resultDEG),",but only ",nrow(df2)," genes enjoy same trend."))
  return(df2)

}
