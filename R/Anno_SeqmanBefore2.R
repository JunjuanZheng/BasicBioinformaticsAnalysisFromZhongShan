


####======================SeqmanBefore2()=======================####
# prepare chip probe information(.txt) before seqman pipeline via GPL names

## Usage:When we want to annotate a chip with probe sequence,Seqman software is a good choice.Before use Seqman,a probe information should be prepared,which SeqmanBefore() can help do.

# gpl.name # GPL object name
# probe.id.col# the colname of the probe id
# sequence.col# the colname of the sequence information
# save.file# whether save results as rda
# names# part of saved file name.

SeqmanBefore2 <- function(gpl.name,
                          probe.id.col="ID",
                          sequence.col = "SEQUENCE",
                          save.file=T,
                          names = "test1",
                          path = "E:/RCloud/database/DataDownload/GEOqueryDownload/GPL Platform"){
  ## 加载必要的包
  need <- c("Biobase","GEOquery")
  Plus.library(need)

  ## 下载GPL文件
  dl <- paste0(path,"/",gpl.name)
  dir.create(dl,recursive = T)
  print(paste0("Download ",gpl.name," object from GEO database..."))
  gpl <- getGEO(gpl.name,destdir = dl)
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


