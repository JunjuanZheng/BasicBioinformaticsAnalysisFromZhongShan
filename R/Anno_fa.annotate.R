



####==========================fa.annotate()=======================##

#'annotate probe sequence via fasta object
#'@description annotate probe sequence via fasta object
#'@param sequence one or multiple sequeces (like "AGCGC")
#'@param fa.annotation fa.annotation object from lucky::get.fa.annoatation()
#'@param logical see parrallel::detectCores()
#'@return a vector containing ENSEMBL ids.
#'@author Weibin Huang<\email{654751191@@qq.com}>
#'@examples
#'ensembl <- fa.annotate(sequence,fa.annotation)
#'@export
fa.annotate <- function(sequence,fa.annotation,logical = F){
  ## 提取探针对应的序列所在位置
  biostring <- fa.annotation[["sequece"]]
  information <- fa.annotation[["information"]]
  print("获取探针序列在sequence中的具体位置,请耐心等待...")
  p <- Fastgrep2(sequence,biostring,parrallel = T,logical=logical)
  print("成功获取探针序列在sequence中的具体位置！")

  ## 报告
  logic1 <- table(is.na(p) %in% T)
  if(T %in% names(logic1)){
    na.count <- logic1[names(logic1) %in% T]
    print(paste0("一共有",na.count,"条探针无法正确注释。"))
  } else {
    print("所有探针均在sequence中成功注释。")
  }

  ## 对于某个p，提取相关注释信息
  fa.extra1 <- function(information,p.i){
    if(!is.all.na(p.i)){
      in.i <- information[p.i,]
      ensg.types <- unique(as.character(in.i$ENSEMBL))
      if(length(ensg.types) == 0){
        # 说明没有成功注释
        ensg.i <- NA
      } else {
        # 说明可以成功注释
        ensg.i <- paste0(ensg.types,collapse = " // ")
      }
    } else {
      #是NA值，说明没有东西
      ensg.i <- NA
    }
    return(ensg.i)
  }
  fa.extra1(information,p[[10]])

  ## 提取多个p中的注释信息
  print("正在提取多个探针的注释信息...")
  ensembl <- apply(as.matrix(1:length(p)),1,
                   function(i)fa.extra1(information,p[[i]]))
  ensembl <- unlist(ensembl)
  print("完成多个探针的注释信息的提取!")

  ##输出结果
  print("成功输出ENSEMBL类注释。")
  return(ensembl)
}





