


###=====co.matrix2():co.matrix的改良版
##作用：以某些基因为参照，对目的矩阵选择可以注释成功的基因。例如，拟通过valid集对train集限定基因进行差异性分析，而valid和train分别属于两个不同的annotation系统，此时co.matrix2()函数将十分快速、有效地找到valid与train的交集基因，并输出统一注释的基因矩阵。co.matrix2()将ENSEMBL作为唯一性标识进行基因的识别。推荐co.matrix2()而不是co.matrix()进行芯片数据的预处理。
##输入：logscale的原始芯片矩阵；输出：包含matrix的List文件。
##参数
#list.exprs.matrix=list(x1=gene.mt,x2=gene.mt)#
#list.annotation=list(x1=annotation,x2=annotation)#
#list.annotation.type = list("probeid","probeid")#芯片矩阵的原标注类型。一般是probeid
#list.co_colname = list(x1 = c("Symbol","ENTREZID"),x2=c("Symbol","ENTREZID"))#
#control.co_colname=c("SYMBOL","ENTREZID")#参考annotation对应的列
#control.genes = rownames(dds.valid)#参考基因。保证输出的矩阵包含的基因在此范围内。
#control.genes.type = "ENSEMBL"#参考基因对应的类型。
#control.annotation = db.anotation # 总注释库。有最全的注释信息
#output.annotation.type = "ENSEMBL" #通常为最唯一的注释类型。一般为ENSEMBL
#' @export
co.matrix2 <- function(list.exprs.matrix,
                       list.annotation,
                       list.annotation.type,
                       list.co_colname,
                       control.co_colname,
                       control.genes,
                       control.genes.type,
                       control.annotation,
                       output.annotation.type){
  library(stringr)
  l <- list()
  # i=1
  for (i in 1:length(list.exprs.matrix)) {
    ###获得表达矩阵并进行合并处理。
    exprs.matrix1 <- list.exprs.matrix[[i]];dim(exprs.matrix1)

    ###合并注释文件
    mt.anno <- list.annotation[[i]]
    id <- NULL;
    for(j in 1:length(control.co_colname)){
      p1 <- match(list.co_colname[[i]][j],colnames(mt.anno))
      id.j <- convert(as.character(mt.anno[,p1]),
                      fromtype = control.co_colname[j],
                      totype = output.annotation.type,
                      db = control.annotation)
      id <- cbind(id,id.j);colnames(id)[j] <- paste("CO",j,sep = "_")
    }
    #提取向量中的第一个非空值。
    isnot.na.value <- function(vt){
      p <- which(!is.na(vt))
      isnot1 <- vt[p[1]]
      return(isnot1)
    }
    id1 <- apply(id,1,function(x)isnot.na.value(x))
    mt.anno <- cbind(mt.anno,id,id1)

    #并行运算：完成总矩阵的计算。
    ensembl1 <- convert(rownames(exprs.matrix1),
                        fromtype = list.annotation.type[[i]],
                        totype = "id1",
                        db = mt.anno)
    exprs.matrix2 <- exprs.matrix1[which(!is.na(ensembl1)),]
    dim(exprs.matrix2)
    rownames(exprs.matrix2) <- ensembl1[which(!is.na(ensembl1))]
    #nrow(exprs.matrix2)
    geneidfactor <- factor(rownames(exprs.matrix2))
    x.apply <- function(expr.train1,geneidfactor){
      geneidfactor1 <- 0
      require(parallel)
      cl <- makeCluster(mc <- getOption("cl.cores", 4))
      clusterExport(cl,"geneidfactor")
      geneidfactor1 <- geneidfactor
      i <- parApply(cl = cl,expr.train1,2,function(x)tapply(x,geneidfactor1,mean))
      return(i)
    }
    print(str_c(names(list.exprs.matrix)[i],":正在进行合并同基因的并行运算"))
    exprs.matrix3 <- x.apply(exprs.matrix2,geneidfactor)
    print(str_c(names(list.exprs.matrix)[i],":完成合并同基因的并行运算"))
    #table(duplicated(rownames(exprs.matrix3)))#无重复
    #dim(exprs.matrix3)

    ###按control.genes选择交集
    control.genes <- convert(control.genes,
                             fromtype = control.genes.type,
                             totype = output.annotation.type,
                             db = control.annotation)
    #which(is.na(control.genes))#不会有空值
    co.id <- intersect(control.genes,rownames(exprs.matrix3))#共同基因length(co.id)

    exprs.matrix4 <- exprs.matrix3[co.id,]

    #输出数据
    l <- c(l,list(exprs.matrix4))
  }
  names(l) <- names(list.exprs.matrix)
  return(l)
}



