


###=====co.matrix3():co.matrix2的改良版
##作用：由于更完备的外部注释逐渐开展，因此旧co.matrix2的注释已经是多余的，因此此处对注释进行精简设计。
##输入：logscale的原始芯片矩阵；输出：包含matrix的List文件。
##参数
#list.exprs.matrix=list(x1=gene.mt,x2=gene.mt)#
#list.annotation=list(x1=annotation,x2=annotation)#
#list.annotation.type = list("probeid","probeid")#芯片矩阵的原标注类型。一般是probeid
#list.co_colname = list(x1="ENSEMBL",x2="ENSEMBL")#
#' @export
co.matrix3 <- function(list.exprs.matrix,
                       list.annotation,
                       list.annotation.type,
                       list.co_colname,
                       control.genes=NULL,#参考基因。保证输出的矩阵包含的基因在此范围内。
                       control.genes.type = "ENSEMBL",#参考基因对应的类型。
                       control.annotation = common.annot, # 总注释库。有最全的注释信息
                       output.annotation.type = "ENSEMBL"){
  library(stringr)
  l <- list()
  for (i in 1:length(list.exprs.matrix)) {

    ## 获得表达矩阵并进行合并处理。
    exprs.matrix1 <- list.exprs.matrix[[i]];dim(exprs.matrix1)#View(exprs.matrix1[,1:5])

    ## 对行名进行注释
    #list.annotation.type[[i]]$ID <- as.character(list.annotation.type[[i]]$ID)
    ensembl1 <- convert(rownames(exprs.matrix1),
                        fromtype = list.annotation.type[[i]],
                        totype = list.co_colname[[i]],
                        db=list.annotation[[i]])
    #table(is.na(ensembl1))
    exprs.matrix2 <- exprs.matrix1[which(!is.na(ensembl1)),]#dim(exprs.matrix2)
    rownames(exprs.matrix2) <- ensembl1[which(!is.na(ensembl1))]#nrow(exprs.matrix2)
    geneidfactor <- factor(rownames(exprs.matrix2))
    x.apply <- function(expr.train1,geneidfactor){
      require(parallel)
      ncores <- detectCores(logical = F)
      cl <- makeCluster(mc <- getOption("cl.cores", ncores))
      clusterExport(cl,"geneidfactor",envir = environment())
      i <- parApply(cl = cl,expr.train1,2,function(x)tapply(x,geneidfactor,mean))
      stopCluster(cl)
      return(i)
    }
    print(str_c(names(list.exprs.matrix)[i],":正在进行合并同基因的并行运算"))
    exprs.matrix3 <- x.apply(exprs.matrix2,geneidfactor)
    print(str_c(names(list.exprs.matrix)[i],":完成合并同基因的并行运算"))
    #table(duplicated(rownames(exprs.matrix3)))#无重复
    #dim(exprs.matrix3)

    if(is.null(control.genes)){
      exprs.matrix4 <- exprs.matrix3
    } else {
      ## 按control.genes选择交集
      control.genes <- convert(control.genes,
                               fromtype = control.genes.type,
                               totype = output.annotation.type,
                               db = control.annotation)
      #which(is.na(control.genes))#不会有空值
      co.id <- intersect(control.genes,rownames(exprs.matrix3))#共同基因length(co.id)
      exprs.matrix4 <- exprs.matrix3[co.id,]
    }

    #输出数据
    l <- c(l,list(exprs.matrix4))

  }
  names(l) <- names(list.exprs.matrix)
  return(l)
}



