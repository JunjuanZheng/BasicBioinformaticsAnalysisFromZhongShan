

####=======================SearchIF系列===========================####
###SearchIF():把某个范围的IF进行提取，并输出Pubmed检索类型的文件“txt”和概览“csv”
# data = IF.data#包含影响因子、杂志全名的数据
# col.journal.names = "Full Journal Title"#名称列
# col.IF = "Journal Impact Factor"#影响因子列
# lower.limit = -Inf#下限
# upper.limit = Inf#上限

# data = IF.data# a data frame containing IF and journal names information
# col.journal.names = "Full Journal Title"# the colnames of journal names
# col.IF = "Journal Impact Factor"#the colnames of IF
# lower.limit = -Inf#the lower CI of IF
# upper.limit = Inf#the upper CI of IF
#' @export
SearchIF <- function(data,
                     col.journal.names= "Full Journal Title",
                     col.IF= "Journal Impact Factor",
                     lower.limit=-Inf,
                     upper.limit=Inf){
  library(stringr);library(stringi)
  ##获取相应的列
  data1 <- data[,c(col.journal.names,col.IF)]
  data1 <- as.data.frame(data1)
  IF1 <- as.numeric(data1[,col.IF])
  data1 <- data1[which(!is.na(IF1)),]

  ##获得某个范围的文献；
  IF2 <- as.numeric(data1[,col.IF])
  a <- which(IF2 < upper.limit)
  b <- which(IF2 >= lower.limit)
  pos <- intersect(a,b)
  rows.IF <- data1[pos,]
  write.csv(rows.IF,str_c("Information.IF_",lower.limit,"_",upper.limit,".csv"))
  print(str_c("生成：","Information.IF_",lower.limit,"_",upper.limit,".csv"))

  #构建数据框，给杂志名配引号，加上[Journal]后缀
  titles <- as.character(rows.IF[,col.journal.names])
  df <- NULL
  for (ti in titles) {
    library(stringr)
    df.i <- data.frame(a = I(str_c("\"",ti,"\"")),b = "[Journal]")
    df <- rbind(df,df.i)
  }
  x <- paste0(df$a,df$b,collapse = " OR ")
  write.table(x,str_c("search.IF_",lower.limit,"_",upper.limit,".txt"),sep = "",row.names =F, col.names =F, quote = F)
  print(str_c("生成：","search.IF_",lower.limit,"_",upper.limit,".csv"))

  #结束
  return(rows.IF)

}


###=====SearchIF2()
##检索某个字段的杂志，并输出Pubmed检索类型的文件“txt”和概览“csv”
# data = IF.data#包含影响因子、杂志全名的数据
# col.journal.names = "Full Journal Title"#名称列
# col.IF = "Journal Impact Factor"#影响因子列
# string = "Nature" #字符。不用区分大小写

# data = IF.data#a data frame containing IF and journal names information
# col.journal.names#the colnames of journal names
# col.IF#the colnames of IF
# string# part of journal names.Always a string like nature and cell.
# lower.limit #the lower CI of IF in selected journals.


#' SearchIF2:get journal list with specified strings of journal names
#' @description get journal list with specified strings of journal names
#' @param data IF.data.a data frame containing IF and journal names information
#' @param col.journal.names the colnames of journal names
#' @param col.IF the colnames of IF
#' @param string part of journal names.Always a string like "nature" or "cell"
#' @param lower.limit the lower CI of IF in selected journals
#' @return a txt file with Pubmed search fomula and a csv with journal names and IF.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
SearchIF2 <- function(data,
                      col.journal.names="Full Journal Title",
                      col.IF= "Journal Impact Factor",
                      string="Nature",
                      lower.limit){
  library(stringr);library(stringi)
  ##获取相应的列
  data1 <- data[,c(col.journal.names,col.IF)]
  data1 <- as.data.frame(data1)
  IF1 <- as.numeric(data1[,col.IF])
  data1 <- data1[which(!is.na(IF1)),]

  ##获得含有某个字段的文献；
  titles1 <- as.character(data1[,col.journal.names])
  logic1 <- stri_detect_regex(titles1,string,case_insensitive = TRUE)#忽略大小写
  rows.IF <- data1[logic1,]
  IF2 <- as.numeric(rows.IF[,col.IF])
  rows.IF <- rows.IF[which(IF2 >= lower.limit),]
  write.csv(rows.IF,str_c("Information.IF_",string,"类_","IF大于等于",lower.limit,".csv"))
  print(str_c("生成：","Information.IF_",string,"类_","IF大于等于",lower.limit,".csv"))

  #构建数据框，给杂志名配引号，加上[Journal]后缀
  titles <- as.character(rows.IF[,col.journal.names])
  df <- NULL
  for (ti in titles) {
    library(stringr)
    df.i <- data.frame(a = I(str_c("\"",ti,"\"")),b = "[Journal]")
    df <- rbind(df,df.i)
  }
  x <- paste0(df$a,df$b,collapse = " OR ")
  write.table(x,str_c("search.IF_",string,"类_","IF大于等于",lower.limit,".txt"),sep = "",row.names =F, col.names =F, quote = F)
  print(str_c("生成：","search.IF_",string,"类_","IF大于等于",lower.limit,".txt"))

  #结束
  return(rows.IF)
}


###=====getsearch():可以按给定的journal名称生成检索文件
#' @title get search formula with specified journal names or ISSN
#' @description get search formula with specified journal names or ISSN
#' @param journal.names a vector of journal names
#' @param output.names part of save file name
#' @return a txt file with Pubmed search fomula
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
getsearch <- function(journal.names,
                      searchstring="Journal",
                      output.names){
  library(stringr);library(stringi)
  titles <- as.character(journal.names)
  df <- NULL
  for (ti in titles) {
    library(stringr)
    b <- paste0("[",searchstring,"]")
    df.i <- data.frame(a = I(str_c("\"",ti,"\"")),b = b)
    df <- rbind(df,df.i)
  }
  x <- paste0(df$a,df$b,collapse = " OR ")
  write.table(x,str_c("journal.names_",output.names,".txt"),sep = "",row.names =F, col.names =F, quote = F)
  print(str_c("生成:","journal.names_",output.names,".txt"))
}



