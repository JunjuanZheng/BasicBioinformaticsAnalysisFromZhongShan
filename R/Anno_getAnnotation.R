

####======================getAnnotation()=====================####
## 对于已有的从GENECODE官网下载的gtf注释文件，可以直接生成可注释数据框。一般使用默认参数，仅确定gft.path即可。输出结果为包含注释信息的数据框。

#' Get annotation from a gtf annotation file
#' @description get annotation from a gtf annotation file.
#' @param gtf.path the path of .gtf
#' @param select.col the selected colnames.Default is "gene_id","gene_name" and "gene_type"
#' @param Newcolnames the new colnames.Default is "ENSEMBL","SYMBOL","gene_type"
#' @param control.id the onlyone id colname
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' getAnnotation(gtf.path,
#'               select.col=NULL,
#'               Newcolnames=NULL,
#'               control.id=NULL)
#' save(gtf.annotation,file="gtf.annotation.rda")
#' @export
getAnnotation <- function(gtf.path,
                          select.col=NULL,
                          Newcolnames=NULL,
                          control.id=NULL){

  #gft注释文件从官网下载:https://www.gencodegenes.org/releases/current.html

  ## 加载包
  nd <- c("rtracklayer","stringr")
  Plus.library(nd)

  ## 定义可能为空的值：select.col/Newcolnames/control.id
  if(is.null(select.col)){
    select.col=c("gene_id","gene_name","gene_type")
    print("select.col参数为空。加载默认值。")
  } else {
    select.col = select.col
  }
  if(is.null(Newcolnames)){
    Newcolnames=c("ENSEMBL","SYMBOL","gene_type")
    print("Newcolnames参数为空。加载默认值。")
  } else {
    Newcolnames=Newcolnames
  }
  if(is.null(control.id)){
    control.id="ENSEMBL"
    print("control.id参数为空。加载默认值。")
  } else {
    control.id=control.id
  }

  ## 加载gpl文件
  print("加载gpl文件中，请耐心等待...")
  gtf <- rtracklayer::import(gtf.path)
  gtf_df <- as.data.frame(gtf)
  print("成功加载gpl注释文件!")
  #View(gtf_df[1:1000,])

  ## 一般从官网上下载的gft文件有以下列。判断是否有这些列。
  all.cols <- c("seqnames","start","end","width","strand","source","type","score","phase","gene_id","gene_type","gene_name","level","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid")
  l1 <- all(select.col %in% all.cols)

  ## 选择数据
  if(l1){
    ##可选择的列的正确的。
    x1 <- select.col[Newcolnames %in% control.id]
    gtf2_df <- gtf_df[duplicated(gtf_df[,x1]) == F,]
    len1 <- length(gtf2_df$gene_id)#64485个
    print(str_c("一共有",len1,"个独立的",control.id,"。"))
    df.select1 <- subset(gtf2_df,select=select.col)
    colnames(df.select1) <- Newcolnames

    ## 一般ENSEMBL后方会带有版本号，去除版本号即可。
    print("去除ENSEMBL版本号...")
    i=df.select1$ENSEMBL[1]
    extra.ensembl2 <- function(vt,n){
      vt <- matrix(vt,ncol = 1)
      vt1 <- apply(vt, 1, function(x)unlist(strsplit(x,"[.]"))[n])
      return(vt1)
    }
    df.select1$ENSEMBL.Versions <- extra.ensembl2(df.select1$ENSEMBL,2)
    df.select1$ENSEMBL <- extra.ensembl2(df.select1$ENSEMBL,1)
    print("完成ENSEMBL版本号去除!")

    ## 输出结果
    return(df.select1)
    print("完成注释文件的输出。")
  } else {
    stop("参数错误。请重新设置。")
  }

}


