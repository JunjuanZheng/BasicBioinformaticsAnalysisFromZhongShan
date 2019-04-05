
####=====================get.ISSN2==========================####
#'After get.ISSN
#'@description get.ISSN2 is after get.ISSN and produce more completeISSN code information.
#'@param data a data containing journal,ISSN.code and search.names,which is always produce by luckey::get.ISSN
#'@param journal.name.col the colname of journals
#'@param ISSN.col the colname of ISSN codes
#'@param search.col the colname of search names
#'@details Because the correction of error records need get.journal.ISSN,the speed of get.ISSN2 would be slow and please give more patience.
#'@return data frame contain ISSN information
#'@author Weibin Huang<\email{654751191@@qq.com}>
#'@export
get.ISSN2 <- function(data,
                      journal.name.col="journal",
                      ISSN.col="ISSN.code",
                      search.col="search.names"){
  ## 加载必要的包
  need <- c("plyr","stringi")
  Plus.library(need)

  ## 提取异常数据
  logic1 <- data[,search.col] %in% "error"
  df.error <- data[logic1,]
  logic2 <- data[,search.col] %in% "ISSN.NULL"
  df.NULL <- data[logic2,]
  logic3 <- data[,search.col] %in% "NON.Online"
  df.NOtOnline <- data[logic3,]
  logic4 <- data[,search.col] %in% c("NON.Online","ISSN.NULL","error")
  df.normal <- data[!logic4,]

  ## 处理error类ISSN
  print("correct error records...") #i=1
  deal.error <- function(df.error){
    df.error1 <- NULL
    for(i in 1:nrow(df.error)){
      e1 <- as.character(df.error[i,journal.name.col])
      a <- Fastextra(e1,"-")
      df.a <- get.journal.ISSN(a,show.music = F)
      logic1 <- df.a[,search.col] != "error"
      if(all(logic1) == F){
        df.error.i <- df.error[i,]
      } else {
        p.i <- match(T,logic1)
        df.error.i <- df.a[p.i,]
        df.error.i <- as.matrix(df.error.i)
        df.error.i[1,journal.name.col] <- as.character(df.error[i,journal.name.col])
        df.error.i <- as.data.frame(df.error.i)
      }
      df.error1 <- rbind(df.error1,df.error.i)
    }
    return(df.error1)
  }
  df.error1 <- deal.error(df.error)
  df.error1 <- as.matrix(df.error1)
  df.error1[grep("NON.Online",df.error1[,search.col]),ISSN.col] <- Fastextra(df.error1[grep("NON.Online",df.error1[,search.col]),ISSN.col],"_",1)
  df.error1[grep("ISSN.NULL",df.error1[,search.col]),ISSN.col] <-  df.error1[grep("ISSN.NULL",df.error1[,search.col]),journal.name.col]
  df.error1 <- as.data.frame(df.error1)

  ## 处理ISSN.NULL类
  print("correct ISSN.NULL records...")
  df.NULL[,ISSN.col] <- df.NULL[,journal.name.col]

  ## 处理NON.Online类
  print("correct NON.Online records...")
  df.NOtOnline[,ISSN.col] <- Fastextra(df.NOtOnline[,ISSN.col],"_",1)

  ## 输出结果
  data2 <- rbind(df.normal,df.NULL,df.NOtOnline,df.error1)
  p <- Fastmatch(as.character(data[,journal.name.col]),as.character(data2[,journal.name.col]))
  data2 <- data2[p,]
  rownames(data2) <- 1:nrow(data2)
  print("All done!")
  return(data2)
}
