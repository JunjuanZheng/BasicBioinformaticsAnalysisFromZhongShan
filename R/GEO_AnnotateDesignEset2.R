



####====================AnnotateDesignEset2=========================####
# Usage:AnnotateDesignEset2() is the plus version of AnnotateDesignEset().AnnotateDesignEset2 annotation chip via sequece enhanced stragegy.
#' @export
AnnotateDesignEset2 <- function(DesignEset,gpl.path=NULL){

  ## 加载包
  library(Biobase);library(stringr)

  ## 获取DesignEset中的探针信息
  eset <- DesignEset[["ESet"]]
  db <- DesignEset[["DesignList"]][["array.annotation"]][["db.anno"]]
  probeid.col <- DesignEset[["DesignList"]][["array.annotation"]][["probeid.col"]]
  sequence.col <-  DesignEset[["DesignList"]][["array.annotation"]][["sequence.col"]]

  ## 将control探针排除在外(从QC.cluster内移植过来)
  library(Biobase)
  expr <- exprs(eset)
  annotation0 <- fData(DesignEset[["ESet"]])
  control.probes <- DesignEset[["DesignList"]][["control.probes"]]
  if(is.null(control.probes)){
    #并无control.probes
    expr <- expr
  } else {
    #有对照探针
    con.position <- NULL
    for(i in 1:length(control.probes)){
      con.i <- control.probes[[i]]
      #找到control.probes的位置
      con.position.i <- NULL
      # s="control"
      for(s in con.i$control.symbol){
        con.position.s <- grep(s,annotation0[,con.i$control.col])
        con.position.i <- c(con.position.i,con.position.s)
      }
      con.position <- c(con.position,con.position.i)
    }
    con.position <- unique(con.position)
    # 去除control探针
    all.position <- 1:length(annotation0[,con.i$control.col])
    l1 <- all.position %in% con.position
    expr <- expr[!l1,]
  }
  #expr

  ## 根据探针序列进行增强注释
  #采用sequece增强注释。此时"GPL6699_information.rda"之类的注释文件应该在地址E:\RCloud\database\DataDownload\annotation\final anotation\GEO Annotation中提前准备好。"GPL6699_information.rda"来自lucky包中的SeqmanBefore和SeqmanAfter函数。
  print("根据探针序列进行增强注释...")

  if(is.null(gpl.path)){
    # 采用默认地址
    annot.file = "E:/RCloud/database/DataDownload/annotation/final anotation/GEO Annotation"
    gpls <- list.files(path = annot.file,pattern = ".information.rda",full.names = T)
  } else {
    gpls <- list.files(path = gpl.path,pattern = ".information.rda",full.names = T)
  }

  # 找到和本芯片相关的GPL注释地址并加载
  platform <- DesignEset[["ESet"]]@protocolData@data
  platform1 <- as.character(platform$platform_id[1])
  gpl.ps <- grep(platform1,gpls)
  gpl1.path <- gpls[gpl.ps]
  print(paste0("加载",platform1,"平台的本地gpl.rda类注释文件..."))
  load(gpl1.path) # rm(information)
  print("完成加载!")

  # 根据芯片的序列来进行重注释
  print("根据芯片的序列来进行重注释...")
  colnames(information)
  annotation0$ENSEMBL <- convert(annotation0[,sequence.col],
                                 fromtype = "probe_seq",
                                 totype = "ENSEMBL",
                                 db = information)
  annotation1 <- annotation0[!is.na(annotation0[,"ENSEMBL"]),]
  r1 <- length(annotation1[,"ENSEMBL"])
  print(paste0("一共成功注释",r1,"条探针..."))

  # 将ENSEMBL中含有“_”的删除。说明它们是非特异的探针
  repeat_p <- grep("_",annotation1[,"ENSEMBL"]) #length(repeat_p)
  print(paste0("发现",length(repeat_p),"条非特异性探针，予以删除..."))
  annotation1 <- annotation1[-repeat_p,]
  r2 <- length(annotation1[,"ENSEMBL"])
  print(paste0("完成探针序列重注释!最终可用探针共",r2,"条。"))

  ## co.matrix()进行重复探针的合并
  print("进入co.matrix3()队列,请耐心等待...")
  #annotation1$ID <- rownames(annotation1)
  list.exprs.matrix=list(DesignEset=expr)#DesignEset的名字存在时报告会
  list.annotation=list(annotation1)
  list.annotation.type = list(probeid.col)
  list.co_colname = list("ENSEMBL")#
  control.genes = annotation1[,"ENSEMBL"]#table(duplicated(control.genes))
  control.genes.type = "ENSEMBL"
  control.annotation = db # 总注释库。有最全的注释信息
  output.annotation.type = "ENSEMBL" #通常为最唯一的注释类型。一般为ENSEMBL
  expr1 <- co.matrix3(list.exprs.matrix,
                      list.annotation,
                      list.annotation.type,
                      list.co_colname,
                      control.genes,
                      control.genes.type,
                      control.annotation,
                      output.annotation.type)
  expr1 <-  expr1[[1]]
  print(paste0("获得初步合并表达矩阵,基因数=",nrow(expr1),"个。",collapse = ""))
  # nrow(expr1)
  # save(expr1,file = "expr1.rda")
  # load("expr1.rda");View(expr1[,1:5])
  # View(expr1[,1:5])

  ## 根据expr1生成新的注释。
  ensembl1 <- rownames(expr1)
  symbol1 <- convert(ensembl1 ,
                     fromtype="ENSEMBL",
                     totype="SYMBOL",
                     db=common.annot)
  annotation2 <- data.frame(
    row.names = ensembl1,
    SYMBOL = symbol1
  )

  #新的fdata
  feature <- as(annotation2,"AnnotatedDataFrame")

  ##构建新的eset类
  pheno <- pData(eset);pheno <- as(pheno,"AnnotatedDataFrame")
  proData <- protocolData(eset)#View(proData@data)
  library(methods)
  eset1 <- new('ExpressionSet',
               exprs=expr1,
               phenoData=pheno,
               featureData = feature,
               protocolData = proData)
  DesignList <- DesignEset[["DesignList"]]
  DesignList[["control.probes"]] <- NULL
  print("完成构建新的eset类!")

  ##构建DesignEset
  DesignEset1 <- list(ESet = eset1,DesignList = DesignList)
  print("完成构建新的DesignEset类!")
  ##输出结果
  tuneR::play(music)
  return(DesignEset1)

}
# DesignEset1 <- AnnotateDesignEset(DesignEset)


