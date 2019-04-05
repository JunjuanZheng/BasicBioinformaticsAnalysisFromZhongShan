


###======AddDesignList():生成DesignEset列表。
## 参数
# factor = ifelse(pdata$pnstage == "N0","N0","Np") #分组向量
# levels = c("Np","N0") # factor的level
# control.probes = list(l1=list(control.col="Probe_Type",control.symbol = c("A","I")),l1=list(control.col="Source",control.symbol = "ILMN_Controls")) #对照探针组
## array.annotation = list(probeid.col="ID",symbol.cols=c("Symbol","Unigene_ID","Entrez_Gene_ID"),anno.cols=c("SYMBOL","UNIGENE","ENTREZID"),db.anno=common.annot) #探针注释相关信息
#' @export
AddDesignList <- function(eset,
                          factor,
                          levels,
                          control.probes,
                          array.annotation,
                          NA.method = c("one","all")[1],
                          normalizeBetweenArrays=F){
  ## 加载包
  library(stringr)
  library(limma)
  library(Biobase)

  ## 构建design
  factor <- factor(factor,levels = levels)
  design <- model.matrix(~0+factor);colnames(design) <- levels
  print("成功构建design对象!")

  ## 生成contrast
  contrasts <- combn(levels,2)
  cont1 <- NULL;
  for(i in 1:ncol(contrasts)){
    e <- as.character(contrasts[,i])
    c.i1 <- paste0(e,collapse = "-")
    #c.i2 <- paste0(e[c(2,1)],collapse = "-")
    cont1 <- c(cont1,c.i1)
  }
  contrast.matrix <- makeContrasts(contrasts=cont1,levels=levels)
  print("成功构建contrast对象!")

  ## 构建condition.list
  #可能会有多个对比，因此应该生成多个condition组合为list.
  condition <- NULL
  for(i in 1:ncol(contrast.matrix)){
    level.i <- unlist(strsplit(colnames(contrast.matrix)[i],"-"))
    #level.i对应的factor.i
    factor.i <- factor[factor %in% level.i]

    ##构建condition.i
    condition.i <- as.data.frame(factor.i);
    colnames(condition.i) <- "condition"
    rownames(condition.i) <- colnames(eset)[factor %in% level.i]
    condition.i$V1 <- rep(0,nrow(condition.i))
    condition.i$condition <- factor(condition.i$condition,levels=level.i)
    condition.i <- condition.i[order(condition.i$condition,decreasing = T),]
    condition.i <- subset(condition.i,select = "condition")
    condition <- c(condition,list(condition.i))
    names(condition)[i] <- colnames(contrast.matrix)[i]
  }
  print("成功构建condition对象!")

  ## 去除矩阵有NA值的行
  expr <- exprs(eset)
  print("正在判断表达矩阵是否有NA值,请耐心等候...")
  if(NA.method == "one"){
    #依赖于hwb.base包的NA系列函数
    x <- apply(expr,1,FUN=is.one.na)
  } else {
    if(NA.method == "all"){
      x <- apply(expr,1,FUN=is.all.na)
    } else {
      print("NA.method输入有误。default = one.")
      x <- apply(expr,1,FUN=is.one.na)
    }
  }
  len <- table(x);n <-  as.numeric(len[names(len) %in% T])
  if(length(n)==0){
    #说明没有空值。矩阵完好。
    print("表达矩阵并无空值。")
    select=rep(T,nrow(expr))
    expr1 <- expr[select,];
  } else {
    #说明有空值，予以删除
    if(NA.method == "all"){
      print(str_c("表达矩阵中共发现",n,"行全为NA值!"))
      select = !x
      #输出新矩阵
      print("成功删除NA行!")
      expr1 <- expr[select,];
      print("把矩阵中的NA值转化为0...")
      expr1[expr1 %in% NA] <- 0
      print("成功把矩阵中的NA值转化为0!")
    } else {
      print(str_c("表达矩阵中共发现",n,"行至少有一个NA值!"))
      select = !x
      expr1 <- expr[select,];
      print("成功删除NA行!")
    }
  }

  ## 将矩阵中的负值全部调为0。因为在双色芯片中，如果log Cy5/Cy3<0,意味着所谓的探针强度还没有背景强度强，说明探针很弱或者基本不表达，可以认为是没有表达。但为了避免0值带来计算的某些不方便，因此设为10^-5.
  logic1 <- expr1 <= 0;
  if(T %in% logic1){
    #说明有负值，全部校正为0
    print("矩阵中有负值。由于负值是没有意义的，因此全部校正为0。")
    expr1[logic1] <- 10^(-5)
  } else {
    #说明没有负值。
    print("矩阵中没有负值。")
    expr1  <- expr1
  }

  ## 构建新Eset类列表
  ##阵列间标准化
  if(normalizeBetweenArrays==T){
    print("以limma::normalizeBetweenArrays进行标准化...")
    expr2 <- normalizeBetweenArrays(expr1)
    print("完成normalizeBetweenArrays标准化!")
  } else {
    expr2 <- expr1
  }
  print("重新构建eset类...")
  #根据新矩阵选择新fdata
  feature <- fData(eset)[select,];
  feature <- as(feature ,"AnnotatedDataFrame")
  #提取其它信息
  pheno <- pData(eset);
  pheno <- as(pheno,"AnnotatedDataFrame")
  proData <- protocolData(eset)
  #输出新eset类
  eset1 <- new('ExpressionSet',
               exprs=expr2,
               phenoData=pheno,
               featureData = feature,
               protocolData = proData)
  print("已构建无NA值的新eset!")

  ##输出DesignList
  DesignList <- list(condition = condition,
                     design=design,
                     contrast.matrix = contrast.matrix,
                     control.probes = control.probes,
                     array.annotation = array.annotation)

  ##构建DesignEset
  DesignEset <- list(ESet = eset1,DesignList = DesignList)
  print("成功构建DesignEset类对象!")
  return(DesignEset)
}

# factor #分组向量
# levels=c("Np","N0") # factor的level
# control.probes=list(l1=list(control.col="Probe_Type",control.symbol = c("A","I")),l1=list(control.col="Source",control.symbol = "ILMN_Controls")) #对照探针组
# array.annotation=list(probeid.col="ID",symbol.cols=c("Symbol","Unigene_ID","Entrez_Gene_ID"),anno.cols=c("SYMBOL","UNIGENE","ENTREZID"),db.anno=common.annot) #探针注释相关信息



