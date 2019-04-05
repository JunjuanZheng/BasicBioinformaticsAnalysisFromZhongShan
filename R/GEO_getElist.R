


####======================getElist()========================####
## getElist()通过给定GSE号，自动下载标准格式的soft文件，提取其中的矩阵、注释和target信息，最后形成eset类列表。此函数依赖于GEOquery，可以获得详尽的表型、注释、处理流程内容，十分强大。不过，如果是外地下载soft文件，应该改为GSEXXXXX.soft.gz的名称才能正确识别。如果想保存eset文件为rda，可以选择savefile = T。如果有个性化GSE.soft文件保存空间，可以自定义gse.space属性。
## 升级情况
# 2018-09-26:修复了有多种格式的pdata（比如肿瘤和正常标本列数、列名不致）无法合并的bug(from GSE17154)
# 2018-09-26:修复了有GSE中原作者可能存在的GSM名字错误导致eset无法正常生成的bug(from GSE17154)
#' @export
getElist <- function(gse.names,
                     gse.space=NULL,
                     savefile = F){
  ##加载必要的包
  library(pacman)
  p_load("stringr","GEOquery","DT","Biobase")

  ####生成私人空间
  if(is.null(gse.space)){
    public.workspace = "E:/RCloud/database/DataDownload/GEOqueryDownload" # 这个因人而异
    gse.space = str_c(public.workspace,"/",gse.names)
    dir.create(gse.space,recursive = T)
  } else {
    dir.create(gse.space,recursive = T)
  }

  ##下载GSE数据
  gse <- getGEO(gse.names,GSEMatrix=F,destdir = gse.space)

  ##判断平台信息,选择其中一个平台
  gsmplatforms <- lapply(GSMList(gse),function(x){Meta(x)$platform_id})
  platforms <- NULL;
  for (i in 1:length(gsmplatforms)) {
    g.i <- gsmplatforms[[i]]
    platforms <- c(platforms,g.i)
  }
  platforms.type = unique(platforms) #全部相同
  print(str_c("platform information:",platforms.type))

  ##构建Elist
  eset2 <- list()
  #i=1
  for (i in 1:length(platforms.type)) {
    #选择其中一个平台
    gsmlist <- Filter(function(gsm){Meta(gsm)$platform_id==platforms.type[i]},GSMList(gse));length(gsmlist)
    GSM.names <- names(gsmlist) #记录名

    ## featureData(GPL平台)
    print(str_c("开始第",i,"个featureData构建..."))
    gpl <- GPLList(gse)[[i]]
    gpl1 <- Table(gpl)
    l1 <- gpl1$ID %in% NA #有时候这一列有NA值，影响rownames注释。删除。
    gpl1 <- gpl1[!l1,]
    rownames(gpl1) <- gpl1$ID
    gpl1$ID <- rownames(gpl1) #有时ID是数字，这给后面的注释带来不良影响。因此统一ID类为字符型变量。
    feature <- as(gpl1,"AnnotatedDataFrame")
    print(str_c("完成第",i,"个featureData构建!"))

    ## exprs
    print(str_c("开始第",i,"个exprs构建..."))
    probesets <- gpl1$ID
    data.matrix <- do.call('cbind',lapply(gsmlist,function(x)
    {tab <- Table(x)
    mymatch <- match(probesets,tab$ID_REF)
    return(tab$VALUE[mymatch])
    }))
    data.matrix <- apply(data.matrix,2,function(x){as.numeric(as.character(x))})
    rownames(data.matrix) <- probesets
    colnames(data.matrix) <- names(gsmlist)
    #View(data.matrix[1:10,1:10])
    data.matrix <- data.matrix[,GSM.names]
    print(str_c("完成第",i,"个exprs构建!"))

    ## phenoData
    print(str_c("开始第",i,"个phenoData构建..."))
    pdata <- GSMList(gse)
    ## 制定统一列名:有时这里有列名经常写错（比如多打一个字母）。这里统一用第一个array的信息作用列名。
    cols <- list();len <- NULL;
    for(a in 1:length(pdata)){
      pdata.i <- pdata[a]
      pdata.i2 <- pdata.i[[1]]@header[["characteristics_ch1"]]
      #vt <- pdata.i2;z="node: 0"
      extra.pdata <- function(vt,n,split="[:] "){
        vt <- as.character(vt)
        a <- NULL
        for(z in vt){
          z.i <- unlist(strsplit(z,split))
          a <- c(a,z.i[n])
        }
        return(a)
      }
      col.a <- extra.pdata(pdata.i2,1)
      len <- c(len,length(col.a))
      cols <- c(cols,list(col.a))
    }
    len.types <- unique(len)
    len.ps <- match(len.types,len)
    list.colnames <- cols[len.ps];names(list.colnames) <- len.types
    ## 提取多个芯片的信息。
    df.pdata <- data.frame()#j=2
    for (j in 1:length(pdata)) {
      pdata.i <- pdata[j]
      pdata.i2 <- pdata.i[[1]]@header[["characteristics_ch1"]]
      df.j <- data.frame(pattern = extra.pdata(pdata.i2,1),
                         value = extra.pdata(pdata.i2,2))
      colnames(df.j)[2] <- names(pdata.i)
      df.j <- as.data.frame(t(df.j))
      colnames(df.j) <- list.colnames[[match(ncol(df.j),len.types)]]
      if(ncol(df.j) == 1){
        #即只有一个变量，此时数据框变化要多费周折
        coln <-  colnames(df.j)
        df.j <- df.j[2,]
        df.j <- as.character(df.j)
        df.j <- as.data.frame(df.j)
        colnames(df.j) <- coln;rownames(df.j) <- names(pdata.i)
      } else {
        df.j <- df.j[-1,]
      }
      library(plyr)
      df.pdata <- rbind.fill(df.pdata,df.j)
    }
    rownames(df.pdata) <- GSM.names
    pheno <- as(df.pdata,"AnnotatedDataFrame")
    print(str_c("完成第",i,"个phenoData构建!"))

    ## protocalData
    print(str_c("开始第",i,"个protocalData构建..."))
    library(plyr)
    proData <- NULL
    for(z in 1:length(gsmlist)){
      gsm1 <- gsmlist[[z]]
      gsm.df <- data.frame(v1=gsm1@dataTable@columns[["Description"]]);rownames(gsm.df) <- as.character(gsm1@dataTable@columns[["Column"]])
      #提取header的信息。如果某一项有2个及以上值，合并
      header <- gsm1@header;length(header)
      gsm.df2 <- NULL
      for(a in 1:length(header)){
        header.i <- header[[a]];len=length(header.i)
        header.i <- paste0(header.i,collapse = "_")
        header.i2 <- data.frame(v1=header.i,v2=len);rownames(header.i2)[1] <- names(header[a])
        gsm.df2 <- rbind(gsm.df2,header.i2)
        #print(paste0("a",a))
      }
      gsm.df2 <- gsm.df2[order(gsm.df2$v2,decreasing = T),]
      gsm.df2 <- subset(gsm.df2,select = "v1")
      #进一步处理header信息
      gsm.df4 <- rbind(gsm.df,gsm.df2)
      colnames(gsm.df4)[1] <- names(gsmlist[z])
      gsm.df5 <- as.data.frame(t(gsm.df4))
      #合并数据
      proData <- rbind.fill(proData,gsm.df5)
      #print(paste0("z",z))
    }
    rownames(proData) <- GSM.names
    proData <- as(proData,"AnnotatedDataFrame")
    print(str_c("完成第",i,"个protocalData构建!"))


    ## 生成ExpressionSet object
    eset <- new('ExpressionSet',
                exprs=data.matrix,
                phenoData=pheno,
                featureData = feature,
                protocolData = proData)
    print(str_c("完成第",i,"个eset构建!"))

    eset2 <- c(eset2,list(eset))
    names(eset2)[i] <- names(GPLList(gse))[i]
  }
  colnames(data.matrix)
  rownames(df.pdata)


  ##输出结果
  if(length(eset2) == 1){
    eset3 <- eset2[[1]]
    print(str_c(gse.names,"为单注释平台数据。"))
  } else {
    eset3 <- eset2
    print(str_c(gse.names,"为多注释平台数据。"))
  }

  ## End
  if(savefile == F){
    tuneR::play(music)
    return(eset3)
  } else {
    eset <- eset3
    save(eset,file = str_c(gse.space,"/",gse.names,"_eset.rda"))
    tuneR::play(music)
    return(eset3)
  }
}
## 例子
# gse.names = "GSE11121"
# View(exprs(eset)[1:10,]) #矩阵信息
# View(fData(eset)[1:10,]) #注释信息
# View(pData(eset)) #表型信息
# View(protocolData(eset)@data) #流程信息


