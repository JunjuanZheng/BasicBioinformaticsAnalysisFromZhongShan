

####======================MApipeline2()========================####
## Usage:transform a raw MA object to a comparable MA which is needed in limma:lmFit().

## Description
# MApipeline2 is the plus version of MApipeline().It annotate MA object via sequence enhanced annotation strategy.
#' @export
MApipeline2 <- function(MA,
                       MA.probeid,
                       MA.colnames,
                       control.probes,
                       array.annotation,
                       eset,
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

  ## 根据探针序列进行增强注释
  {
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
  platform <- eset@protocolData@data
  platform1 <- as.character(platform$platform_id[1])
  gpl.ps <- grep(platform1,gpls)
  gpl1.path <- gpls[gpl.ps]
  print(paste0("加载",platform1,"平台的本地gpl.rda类注释文件..."))
  load(gpl1.path) # rm(information)
  print("完成加载!")

  # 根据芯片的序列来进行重注释
  print("根据芯片的序列来进行重注释...")
  #colnames(information)
  annotation0 <- fdata
  annotation0$ENSEMBL <- convert(annotation0[,sequence.col],
                                 fromtype = "probe_seq",
                                 totype = "ENSEMBL",
                                 db = information)
  annotation1 <- annotation0[!is.na(annotation0[,"ENSEMBL"]),]

  # 将ENSEMBL中含有“_”的删除。说明它们是非特异的探针
  repeat_p <- grep("_",annotation1[,"ENSEMBL"]) #length(repeat_p)
  annotation1 <- annotation1[-repeat_p,]
  print("完成探针序列重注释!")
}
  #annotation1

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
    control.genes = annotation1[,"ENSEMBL"] #table(duplicated(control.genes))
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
    ensembl1 <- rownames(M)
    symbol1 <- convert(ensembl1,
                       fromtype="ENSEMBL",
                       totype="SYMBOL",
                       db=common.annot)
    annotation2 <- data.frame(
      row.names  = ensembl1,
      SYMBOL = symbol1
    )
    #table(duplicated(annotation3$ENSEMBL))

    ##输出结果
    m <- list(M=M,genes = annotation2,cw=expr.control)
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




