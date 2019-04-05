

####======================MApipeline()========================####
## Usage:transform a raw MA object to a comparable MA which is needed in limma:lmFit().

## Description
# Raw MAs always contains probes with the same gene,or lots of control probes which are considered not differentially expressed in every arrays,and all of this need to deal with before MA go to the lmFit pipeline.I believe MApipeline() give a fast way to do this job.


# MA # a MA List
# MA.probeid # the probeid of M matrix or A matrix in MA List.Note that it was the same compose with eset feature data.
# MA.colnames # the new colnames of M matrix or A matrix in MA List.Note that it was the subset of colnames of eset exprs data.
# control.probes # a list containing control probes.See AddDesignList()
# array.annotation # a list containing array annotation information.See AddDesignList()
# eset # a project from getElist()
# control.weight # whether to gather the expression of control probes and weight them with 2.In most time we recommand control.weight=F.You can also compared two results and decide which one is the best choice.
#' @export
MApipeline <- function(MA,
                       MA.probeid,
                       MA.colnames,
                       control.probes,
                       array.annotation,
                       eset,
                       enhanced.annotation=T,
                       gpl.path=NULL,
                       control.weight = F){
  ## 加载必要的包
  library(Biobase);library(stringr)

  ##提取eset中的基本数据
  fdata <- fData(eset)
  pdata <- pData(eset)
  protocoldata <- protocolData(eset)
  expr <- exprs(eset)

  ##提取注释相关信息
  symbol.cols <- array.annotation[["symbol.cols"]]
  sequence.col <- array.annotation[["sequence.col"]]
  probeid.col <- array.annotation[["probeid.col"]]
  db <- array.annotation[["db.anno"]]
  anno.cols <- array.annotation[["anno.cols"]]

  ## 注释探针
  getAnnotation1 <- function(fdata,
                             symbol.col1,
                             anno.col1){
    library(Biobase)
    annotation <- fdata
    p1 <- match(symbol.col1,colnames(annotation))
    symbol1 <- as.character(annotation[,p1])
    nc.symbol1 <- nchar(symbol1)
    ## 找到非空行
    l1 <- nc.symbol1 %in% 0;l1 <- !l1
    annotation1 <- annotation[l1,]
    annotation1$ENSEMBL <- convert(annotation1[,symbol.col1],
                                   fromtype=anno.col1,
                                   totype="ENSEMBL",
                                   db=db)
    table(is.na(annotation1$ENSEMBL))
    return(annotation1)

  }
  annotation0 <- NULL
  for(i in 1:length(symbol.cols)){
    a.i <- getAnnotation1(fdata=fdata,
                          symbol.col1=symbol.cols[i],
                          anno.col1=anno.cols[i])
    annotation0 <- rbind(annotation0,a.i)
  }
  annotation1 <- annotation0[!duplicated(annotation0[,probeid.col]),]
  rownames(annotation1) <- as.character(annotation1[,probeid.col])
  #检验报告
  test1 <- table(duplicated(annotation1[,symbol.cols[1]]))#GeneSymbol
  na.ensembol <- sum(!is.na(annotation1$ENSEMBL))#可注释的个数。
  dup.counts <- test1[names(test1) %in% T]
  print(str_c("共有",sum(test1),"条探针，其中有",dup.counts,"条是重复探针。初步注释探针",na.ensembol,"个。"))

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
    platform <- eset@protocolData@data
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
      yes.c <- test2[names(test2) %in% F]
      print(paste0("经序列注释后，成功注释",yes.c,"条探针。但仍有",na.c,"条探针无法正常注释。"))
    } else {
      print("经序列注释后，所有探针均成功注释。")
    }

    ## 将多重注释的探针表现为NA值，说明属于无效探针。
    annotation1[,"ENSEMBL"] <- convert(annotation1[,"ENSEMBL"],
                                       fromtype = "ENSEMBL",
                                       totype = "ENSEMBL")
    annotation1 <- as.data.frame(annotation1)
  }

  ## 处理MA对象。由于M/A是log scale数据，所以可以通过mean方式合并。因此适用于co.matrix2()的相关策略。
  list.ma <- list();t <- c("M","A") #i=1
  for(i in 1:length(t)){
    ## 行名与fdata对齐
    M <- MA[[t[i]]] #View(M[,1:5])
    rownames(M) <- MA.probeid
    eset.probeid <- rownames(fdata)
    M <- M[eset.probeid,]

    ## 对M进行矩阵列名对齐
    colnames(M) <- MA.colnames
    expr <- M #进行新M的替换

    ## 将control探针排除在外(从QC.cluster内移植过来)
    annotation0 <- fdata
    if(is.null(control.probes)){
      #并无control.probes
      expr <- expr
      expr.control <- MA[["weights"]]
    } else {
      #有对照探针
      con.position <- NULL
      for(z in 1:length(control.probes)){
        con.i <- control.probes[[z]]
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
      expr1 <- expr[!l1,]
      expr.control <- expr[l1,] #control探针的属性
      expr <- expr1
    }

    ## co.matrix()进行重复探针的合并
    print("进入co.matrix3()队列,重复探针默认以mean的形式合并。请耐心等待...")
    list.exprs.matrix=list(MA.matrix=expr)
    list.annotation=list(annotation1)
    list.annotation.type = list(probeid.col)
    list.co_colname = list("ENSEMBL")
    control.genes = annotation1[!is.na(annotation1[,"ENSEMBL"]),"ENSEMBL"] #table(duplicated(control.genes))
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
    M <-  expr1[[1]]
    print(paste0("获得初步合并",t[i],"矩阵,基因数=",nrow(M),"个。",collapse = ""))#有时探针注释里有的基因，在实际检测时可能是没有的。

    ## 构建新的annotatin和表达矩阵
    annotation2 <- annotation1[!is.na(annotation1$ENSEMBL),]
    annotation2 <-annotation2[!duplicated(annotation2$ENSEMBL),]
    rownames(annotation2) <- annotation2$ENSEMBL
    annotation2 <- subset(annotation2,select = symbol.cols)
    for(z in 1:length(colnames(annotation2))){
      s.i <- colnames(annotation2)[z]
      annotation2[,s.i] <- convert(rownames(annotation2),
                                   fromtype="ENSEMBL",
                                   totype=anno.cols[z],
                                   db=db)
    }
    l2 <- rownames(M) %in% rownames(annotation2)
    M2 <- M[l2,]# View(expr2[,1:5])
    print(paste0("完成构建新genes和",t[i],"矩阵构建，基因数=",nrow(M2),"个。",collapse = ""))

    ## 输出genes对象:有时如果symbol只有一个的时候，需要用此法进行数据框重组。
    if(ncol(annotation2)==1){
      #保留logic?
      annotation3 <- annotation2[rownames(M2),]
      annotation3 <- as.data.frame(annotation3)
      colnames(annotation3) <- colnames(annotation2)
      rownames(annotation3) <- rownames(M2)
    } else {
      annotation3 <- annotation2[rownames(M2),]
      annotation3 <- as.data.frame(annotation3)
    }
    genes <- annotation3

    ##输出结果
    m <- list(M=M2,genes = genes,cw=expr.control)
    names(m)[1] <- t[i]
    list.ma <- c(list.ma,list(m))
    names(list.ma)[i] <- t[i]
  }

  ## 构建新MAList对象
  m <- list.ma[["M"]][["M"]]
  a <- list.ma[["A"]][["A"]]
  g <- list.ma[["M"]][["genes"]]
  cw <- list.ma[["M"]][["cw"]]
  if(control.weight == F){
    #说明control类探针可以直接去除。weight全指定为1。因此暂不支持确定高表达基因的低权重赋值。
    M <- m
    A <- a
    genes <- g
    weight <- M;weight[1:nrow(weight),1:ncol(weight)] <- 1
  } else {
    ## 说明指定control probes
    if(is.null(control.probes)){
      print("并无control probes。不适用control weight策略。")
      M <- m
      A <- a
      genes <- g
      weight <- M;weight[1:nrow(weight),1:ncol(weight)] <- 1
    } else {
      ## 有control probe探针。将control probe的对应数据加权为2.
      cw1 <- cw;cw1[1:nrow(cw1),1:ncol(cw1)] <- 2
      m1 <- m;m1[1:nrow(m1),1:ncol(m1)] <- 1
      weight <- rbind(m1,cw1)
      M <- rbind(m,cw)
      A <- rbind(a,cw)
      genes <- g
    }
  }

  ## 构建新MAList对象2
  MA1 <- MA
  MA1[["weights"]] <- weight
  MA1[["M"]] <- M
  MA1[["A"]] <- A
  MA1[["genes"]] <- genes

  ## 输出结果
  return(MA1)
}




