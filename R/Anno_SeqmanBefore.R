


####======================SeqmanBefore()=======================####
#' @title prepare chip probe information(.txt) before seqman pipelin
#' @description When we want to annotate a chip with probe sequence,Seqman software is a good choice.Before use Seqman,a probe information should be prepared,which SeqmanBefore() can help do.
#' @param eset an eset object containing feature data.
#' @param probe.id.col the colname of the probe id
#' @param sequence.col the colname of the sequence information
#' @param save.file whether save results as rda
#' @param names part of saved file name
#' @return  a .txt file containing probe ids and sequeces.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' load("E:/RCloud/database/DataDownload/GEOqueryDownload/GSE17154/GSE17154_rawEset.rda")
#' probe.id.col="ID"
#' sequence.col = "SEQUENCE"
#' test1 = SeqmanBefore(eset,
#'                      probe.id.col,
#'                      sequence.col,
#'                      save.file=T,
#'                      names = "test1")
#' @export
SeqmanBefore <- function(eset,
                         probe.id.col="ID",
                         sequence.col = "SEQUENCE",
                         save.file=T,
                         names = "test1"){
  ## 从eset中获取注释信息
  library(Biobase)
  fdata <- fData(eset)

  ## 提取
  probe.ids <- paste(">",as.character(fdata[,"ID"]),sep = "")
  seq1 <- as.character(fdata[,sequence.col])
  mt <- matrix(rep(0,2*length(probe.ids)),ncol = 1)
  mt[seq(1,(nrow(mt)-1),by=2),] <- probe.ids
  mt[seq(2,(nrow(mt)),by=2),] <- seq1

  ## 保存文件
  library(readr)
  mt <- as.data.frame(mt);
  if(save.file==T){
    write.table(mt,paste0(names,"_probesequence.txt"),
                sep = "\t",quote = F,
                col.names = F,row.names = F)
    return(mt)
  } else {
    return(mt)
  }
  ## End
}


