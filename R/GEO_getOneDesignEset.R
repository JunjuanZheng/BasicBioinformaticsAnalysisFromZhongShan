


####===================getOneDesignEset()=========================####
## getOneDesignEset help to get a DesignEset with only one condition.

## arguement
# DesignEset # a result from AddDesignEset()
# n # which condition and related DesignEset you want to get

#contrast.i=c("N1.plus","N1")
#contrast.i=c("N1.plus","N0")
#contrast.i=c("N1","N0")
#' @export
getOneDesignEset <- function(DesignEset,contrast.i){

  ## 获得condition.i
  condition <- DesignEset[["DesignList"]][["condition"]]
  cn <- names(condition)
  contrast.i1 <- c(paste0(contrast.i[1:2],collapse = "-"),paste0(contrast.i[2:1],collapse = "-"))
  names.i <- cn[cn %in% contrast.i1]
  condition.i <- condition[[grep(names.i,cn)]]

  ## 获得expr.i
  library(Biobase)
  eset <- DesignEset[["ESet"]]
  expr.i <- exprs(eset)
  expr.i <- expr.i[,rownames(condition.i)] #colnames(expr.i)

  ## 提取其它信息
  feature <- fData(eset)
  feature <- as(feature ,"AnnotatedDataFrame")
  pheno <- pData(eset);
  pheno <- cbind(pheno,e=0)#有时有一列，因此作此写法
  pheno <- pheno[rownames(condition.i),]
  pheno <- subset(pheno,select = setdiff(colnames(pheno),"e"))
  pheno <- as(pheno,"AnnotatedDataFrame")
  proData <- protocolData(eset)
  proData <- proData@data
  proData <- proData[rownames(condition.i),]
  proData <- as(proData,"AnnotatedDataFrame")

  ## 输出新eset类
  eset.i <- new('ExpressionSet',
                exprs=expr.i,
                phenoData=pheno,
                featureData = feature,
                protocolData = proData)

  ## 修改DesignList
  DesignEset[["DesignList"]][["condition"]] <- condition.i

  ## 构建DesignEset.i
  DesignEset.i <- list(ESet = eset.i,
                       DesignList = DesignEset[["DesignList"]])
  print(paste0("成功构建condition(",names.i ,")的DesignEset类对象!",sep=""))
  return(DesignEset.i)
}

## example
# DesignEset.i <- getOneDesignEset(DesignEset,c("N1.plus","N0"))











