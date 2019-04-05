

####====================get.fa.annoatation=================####

#'create a powerful annotation object for chips containing probe suqueces information
#'@description get.fa.annoatation() help get suquece associated annotation information from gencode transcripts .fa files.Before entering get.fa.annoatation(),the .fa file needs to be converted to .txt via Notpad++ software.
#'@param fa.path the path of fa file.It should be the .txt file from .fa via Notepad++ software.
#'@return a fa.annoatation object
#'@author Weibin Huang<\email{654751191@@qq.com}>
#'@examples
#'fa.path = "E:/RCloud/database/DataDownload/annotation/unzip data/gencode.v28.transcripts.txt"
#'fa.annotation <- get.fa.annoatation(fa.path)
#'save(fa.annotation,file = "fa.annotation.rda")
#'@export
get.fa.annoatation <- function(fa.path){
  ## 加载必要的包
  need.pac <- c("readr","plyr")
  Plus.library(need.pac)

  ## 初步整理
  print("通过readr::read_table()的快速通道读入数据...")
  rna.anno <- read_table(fa.path)
  print("完成fa类的txt读入!")
  a <- colnames(rna.anno);
  a <- as.data.frame(a)
  colnames(rna.anno) <- "a"
  rna.anno <- rbind(a,rna.anno)
  fa <- as.matrix(rna.anno)
  # save(fa,file = "fa.rda")
  # View(fa[1:100,])

  ## 合并散在序列为一行
  print("正在合并散在序列文件...")
  position <- grep(">",fa[,1])
  position.x <- c(position[2:length(position)],(length(position)+1))
  p.mt <- cbind(position,position.x)
  #position1=position[1];position2=position[2]
  fa.merge <- function(fa,position1,position2){
    fa.i <- fa[(position1+1):(position2-1)]
    fa.i1 <- paste0(fa.i,collapse = "")
    return(fa.i1)
  }
  biostring <- apply(p.mt,1,function(x)fa.merge(fa,x[1],x[2]))
  print(paste0("完成散在序列文件的合并!合并后总数为",length(biostring)))

  ## 将第一列的数据生成多列信息
  fa1 <- fa[position]
  fa.extra <- function(fa1,p1){
    fa1.i <- fa1[p1]
    fa1.i1 <- Fastextra(fa1.i,">",2)
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
  p.fa1 <- data.frame(position = 1:length(fa1),e = 0)
  library(plyr)
  print("ddply正在整合资料，请耐心等待...")
  #fa1[p.i]
  df.fa1 <- ddply(.data = p.fa1,"position",.fun = function(x)fa.extra(fa1,x[1,1]))
  #save(df.fa1,file = "df.fa1.rda")
  print("ddply完成资料整合资料!")

  ## 输出Annotation List
  l1 <- list(
    sequece = biostring,
    information = df.fa1
  )
  print("输出Annotation List完成!")
  return(l1)
}


