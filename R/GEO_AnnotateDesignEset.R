





####====================GEO.Annotation=========================####
# 根据DesignEset对可注释基因进行抽提、合并。基于co.matrix2函数，十分强大。输出结果为一个新的DesignEset。此函数对电脑配置有一定要求，很消耗内存，且耗时，耐心等待。函数结束会自动清空函数相关的内存。
##修复
#2018-9-27 修复了“feature numbers differ between assayData and featureData 类别为“ExpressionSet”的对象不对: 2: featureNames differ between assayData and featureData”的bug。方法为：annotation3 <- annotation3[rownames(expr2),]。bug源芯片：GSE17154
#' @export
AnnotateDesignEset <- function(DesignEset,
                               enhanced.annotation=T,
                               gpl.path=NULL){

  ## 加载包
  nd <- c("Biobase","stringr")
  Plus.library(nd)

  ## 获取DesignEset中的探针信息
  eset <- DesignEset[["ESet"]]
  db <- DesignEset[["DesignList"]][["array.annotation"]][["db.anno"]]
  symbol.cols <- DesignEset[["DesignList"]][["array.annotation"]][["symbol.cols"]]
  anno.cols <- DesignEset[["DesignList"]][["array.annotation"]][["anno.cols"]]
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

  ## 从注释库中获得annotation文件。
  #注意:探针注释中的注意事项。有时symbol会有日期；或者是探针特异性差，代表了多个目标mRNA，这些都应该去除。不过，利用下面的方法，这些“奇怪”的symbol应该是无法成功注释的，也就自动归为无法注释的探针了。因此不需要额外地对这些探针进行筛选。

  f <- fData(eset)
  goldstandard <- ifelse("ENTREZID" %in% anno.cols,"ENTREZID",anno.cols[1])
  f2 <- Plus.convert(data = list(DesignEset = f),
                     type = list(DesignEset = symbol.cols),
                     fromtype = anno.cols,
                     totype = "ENSEMBL",
                     goldstandard = goldstandard,
                     db=db)
  f3 <- f2$DesignEset$metadata
  f3 <- as.matrix(f3)
  #rownames(f3) <- as.character(f3[,"id"])
  #View(head(f3))
  annotation1 <- f3

  #检验报告
  test1 <- table(duplicated(annotation1[,symbol.cols[1]]))#GeneSymbol
  dup.counts <- test1[names(test1) %in% T]
  na.ensembol <- sum(!is.na(annotation1[,"id"]))#可注释的个数。
  print(str_c("共有",sum(test1),"条探针，其中有",dup.counts,"条是重复探针。成功注释探针",na.ensembol,"个。"))
  #length(annotation1[,"ENSEMBL"][ is.na(annotation1[,"ENSEMBL"])])

  ## 是否根据探针序列进行增强注释
  if(enhanced.annotation==T){
    ## 采用sequece增强注释。此时"GPL6699_information.rda"之类的注释文件应该在地址E:\RCloud\database\DataDownload\annotation\final anotation\GEO Annotation中提前准备好。"GPL6699_information.rda"来自lucky包中的SeqmanBefore和SeqmanAfter函数。
    print("采用sequece增强注释。")

    if(is.null(gpl.path)){
      # 采用默认地址
      annot.file = "E:/RCloud/database/DataDownload/annotation/final anotation/GEO Annotation"
      gpls <- list.files(path = annot.file,pattern = ".information.rda",full.names = T)
    } else {
      gpls <- list.files(path = gpl.path,pattern = ".information.rda",full.names = T)
    }

    ## 找到和本芯片相关的GPL注释地址并加载
    platform <- DesignEset[["ESet"]]@protocolData@data
    platform1 <- as.character(platform$platform_id[1])
    gpl.ps <- grep(platform1,gpls)
    gpl1.path <- gpls[gpl.ps]
    print("加载本地gpl.rda注释类文件。")
    load(gpl1.path) # rm(information)

    ## 注释annotation1$ENSEMBL目前未成功注释的
    noannote1 <- annotation1[is.na(annotation1$ENSEMBL),sequence.col]
    p1 <- Fastmatch2(noannote1,as.character(information[,"probe_seq"]))
    ensembl1 <- information[p1,"ENSEMBL"] #head(ensembl1)
    #View(information[1:10,])
    print("完成未知序列的ENSEMBL位置寻找!")

    ## 将获得信息加入annotation1中
    annotation1 <- as.matrix(annotation1)
    annotation1[,"ENSEMBL"][is.na(annotation1[,"ENSEMBL"])] <- ensembl1
    test2 <- table(is.na(annotation1[,"ENSEMBL"]))
    if(T %in% names(test2)){
      na.c <- test2[names(test2) %in% T]
      print(paste0("经序列注释后，仍有",na.c,"条探针无法正常注释。"))
    } else {
      print("经序列注释后，所有探针均成功注释。")
    }

    ## 将多重注释的探针表现为NA值，说明属于无效探针。
    annotation1[,"ENSEMBL"] <- convert(annotation1[,"ENSEMBL"],
                                       fromtype = "ENSEMBL",
                                       totype = "ENSEMBL")
    annotation1 <- as.data.frame(annotation1)
  }

  ## co.matrix()进行重复探针的合并
  print("进入co.matrix3()队列,请耐心等待...")
  #annotation1$ID <- rownames(annotation1)
  list.exprs.matrix=list(DesignEset=expr)
  list.annotation=list(annotation1)
  list.annotation.type = list(probeid.col)
  list.co_colname = list("id") #对于Plus.convert而言，id列是ENSEMBL
  control.genes = as.character(annotation1[!is.na(annotation1[,"id"]),"id"]) #table(duplicated(control.genes))
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

  ## 构建新的annotatin和表达矩阵
  annotation2 <- annotation1[!is.na(annotation1[,"id"]),]
  annotation2 <-annotation2[!duplicated(annotation2[,"id"]),]
  rownames(annotation2) <- annotation2[,"id"]
  annotation2 <- as.data.frame(annotation2,
                               stringsAsFactors=F)
  # View(head(annotation2))
  annotation2 <- annotation2[symbol.cols]
  for(i in 1:length(colnames(annotation2))){
    s.i <- colnames(annotation2)[i]
    annotation2[,s.i] <- convert(rownames(annotation2),
                                 fromtype="ENSEMBL",
                                 totype=anno.cols[i],
                                 db=db)
  }
  #table(is.na(annotation2[,1]))
  #table(is.na(annotation2[,2]))
  #table(is.na(annotation2[,3]))
  #nrow(annotation2)
  l2 <- rownames(expr1) %in% rownames(annotation2) # table(l2)应全为T
  expr2 <- expr1[l2,]# View(expr2[,1:5])
  print(paste0("完成构建新的annotation和表达矩阵!基因数=",nrow(expr2),"个。",collapse = ""))
  sample.names <- colnames(expr2)

  #新的fdata
  feature <- as(annotation2,"AnnotatedDataFrame")

  ##构建新的eset类
  pheno <- pData(eset);pheno <- as(pheno,"AnnotatedDataFrame")
  proData <- protocolData(eset)#View(proData@data)
  library(methods)
  eset1 <- new('ExpressionSet',
               exprs=expr2,
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
  return(DesignEset1)

}
# DesignEset1 <- AnnotateDesignEset(DesignEset)




