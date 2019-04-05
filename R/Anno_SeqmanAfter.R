


####===================SeqmanAfter()=========================####
#' @title extra information from match_probe.txt from Seqman software
#' @description  In chip annotation,after Seqman process,there is a file with patterns like "match_probe",and SeqmanAfter() help to extra information from it.
#' @param path.match.probe the match_probe file created by Seqman software from Sangerbox.
#' @param information.col annotation information colname
#' @param save.file whether save results as rda
#' @param names part of saved file name.
#' @return a data frame containing chip annotation information.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' path.match.probe <- "E:/RCloud/database/DataDownload/GEOqueryDownload/GSE17154/GPL6699_match_probe.txt"
#' information.col="trans_id"
#' SeqmanAfter(path.match.probe,
#'             information.col,
#'             save.file=T,
#'             names ="GPL6699")
#' @export
SeqmanAfter <- function(path.match.probe,
                        information.col="trans_id",
                        save.file=T,
                        names = "test1"){
  ## 导入文件
  nd <- c("readr","plyr","foreach","doParallel")
  Plus.library(nd)
  mp1 <- read_table2(path.match.probe)

  ## 提取某行数据的tran_id的信息
  #fa.extra(见hwb.environment)
  merge.inf <- function(x){
    x1 <- unique(as.character(x))
    x1 <- paste0(x1,collapse = "_")
    return(x1)
  }
  fa.extra <- function(information1){
    fa1.i1 <- information1
    fa1.i2 <- Fastextra(fa1.i1,"[|]",1:8)
    df.i <- data.frame(
      ENSEMBL = Fastextra(fa1.i2[2],"[.]",1),
      ENST = Fastextra(fa1.i2[1],"[.]",1),
      OTTHUMG = Fastextra(fa1.i2[3],"[.]",1),
      OTTHUMT = Fastextra(fa1.i2[4],"[.]",1),
      chrom.loci = fa1.i2[5],
      SYMBOL = fa1.i2[6],
      NO=fa1.i2[7],
      gene.type = fa1.i2[8]
    )
    return(df.i)
  }

  ## 提取多行数据的tran_id的信息
  print("提取多行数据的tran_id的信息...")
  information <- mp1[,information.col]
  test <- 1:nrow(information)
  ncore <- detectCores(logical = F)
  cl <- makeCluster(mc <- getOption("cl.cores",ncore))
  registerDoParallel(cl)
  clusterExport(cl,c("information","fa.extra","Fastextra"),envir = environment())
  df.fa2 <- adply(test,1,
                  .fun = function(x)fa.extra(information[x,1]),
                  .parallel = T)
  #save(df.fa2,file = "df.fa2.rda")
  print("成功提取多行数据的tran_id的信息!")

  ## 加上其它信息
  df.fa2 <- df.fa2[,-1]
  select = setdiff(colnames(mp1),information.col)
  others <- subset(mp1,select = select)
  df.fa2 <- cbind(df.fa2,others)
  #df.fa2 <- information
  #table(duplicated(df.fa2[,"probe_seq"]))

  ## 按sequence重复的情况合并行
  print("按sequence重复的情况合并行...")
  #merge.inf(见hwb.environment)
  geneidfactor <- factor(df.fa2$probe_seq)
  clusterExport(cl,c("geneidfactor","merge.inf"),envir = environment())
  df.fa3 <- parApply(cl = cl,df.fa2,2,function(x)tapply(x,geneidfactor,merge.inf))
  stopCluster(cl)
  df.fa3[,"probe_seq"] <- levels(geneidfactor)
  rownames(df.fa3) <- 1:nrow(df.fa3)
  print("完成按sequence重复的情况合并行!")

  ## 保存文件
  if(save.file == T){
    information <- df.fa3
    save(information,file = paste0(names,"_information.rda"))
  } else {
    print("不保存结果为本地rda.")
  }

  ## 输出
  library(tuneR);play(music)
  return(information)
}



