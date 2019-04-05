


####===================protein annotation========================####
#' @title Get protein annotation from gpl file
#' @description Get protein annotation from gpl file
#' @keywords  ProteinAnnotation
#' @param gtf.path the path of .gtf
#' @param gene.col the colname of gene_id(ENSG)
#' @param protein.col the colname of protein_id(ENSP)
#' @return a data frame of annotation for protein(ENSP id)
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' gtf.path='E:/RCloud/database/DataDownload/annotation/unzip data/gencode.v28.chr_patch_hapl_scaff.annotation.gtf'
#' ProteinAnnotation(gtf.path)
#' @export
ProteinAnnotation <- function(gtf.path,
                              gene.col = "gene_id",
                              protein.col = "protein_id"){

  #gft注释文件从官网下载:https://www.gencodegenes.org/releases/current.html

  ## 加载包
  nd <- c("rtracklayer","stringr","parallel")
  Plus.library(nd)

  ## 加载gpl文件
  print("loading gpl file,please wait...")
  gtf <- rtracklayer::import(gtf.path)
  gtf_df <- as.data.frame(gtf)
  print("gpl file had been loaded!")
  #View(gtf_df[1:1000,])
  rm(gtf,envir = environment());gc()

  ## 一般从官网上下载的gft文件有以下列。判断是否有这些列。
  #all.cols <- c("seqnames","start","end","width","strand","source","type","score","phase","gene_id","gene_type","gene_name","level","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid")
  select.cols <- c(gene.col,"gene_type",protein.col,"gene_name")
  gtf_df2 <- gtf_df[,select.cols]
  gtf_df2 <- gtf_df2[!is.na(gtf_df2$protein_id),]
  rm(gtf_df,envir = environment());gc()

  ## merge.inf:对于某个向量，如果一样，只输出其中一个；如果不一样，则以_连接起来
  merge.inf <- function(x){
    x1 <- unique(as.character(x))
    x1 <- paste0(x1,collapse = "_")
    return(x1)
  }

  #并行规则
  print("Parallel: merge variables based on unique protein_id...")
  ncore <- detectCores()
  if(ncore > 6){
    ncore <- ncore - 1
  }
  cl <- makeCluster(mc <- getOption("cl.cores",ncore))
  geneidfactor <- factor(gtf_df2$protein_id)
  clusterExport(cl,c("geneidfactor","merge.inf"),envir = environment())
  gtf_df3 <- parApply(cl = cl,gtf_df2,2,function(x)tapply(x,geneidfactor,merge.inf))
  stopCluster(cl)
  print("Merge suceed!")
  gtf_df3 <- as.data.frame(gtf_df3)
  #table(duplicated(gtf_df3$gene_id))

  ## id版本号控制
  colnames(gtf_df3) <- c("ENSEMBL","gene_type","protein_id","SYMBOL")
  gtf_df3$ENSEMBL.Versions <- Fastextra(gtf_df3$ENSEMBL,"[.]",2)
  gtf_df3$ENSEMBL <- Fastextra(gtf_df3$ENSEMBL,"[.]",1)
  gtf_df3$protein_id.Versions <- Fastextra(gtf_df3$protein_id,"[.]",2)
  gtf_df3$protein_id <- Fastextra(gtf_df3$protein_id,"[.]",1)

  ## 输出结果
  print("All done!")
  return(gtf_df3)
}

