

####=======================ensemblGo==========================####

## ensemblGo annotate the ENSEMBL ID with GO terms via ensembl online

#genes # a character vector of ENSG id
#parallel # whether use parallel strategy to accelerate the process
#' @export
ensemblGo <- function(genes,
                      parallel = F,
                      save.file = F,
                      names = "love"){
  ## 加载必要的包
  nd <- c("XML","plyr")
  Plus.library(nd)
  go.types <- c("cellular_component","molecular_function","biological_process")
  go.types.ab <- c("CC","MF","BP")

  ## unix时间戳
  unix.time <- as.numeric(as.POSIXct(Sys.Date(), format="%Y-%m-%d"))
  us <- sample(0:9,3,replace = T)
  us <- paste0(us,collapse = "")
  us <- paste(us,us,sep = ".")

  ## 获得某个基因的go term
  get1 <- function(g){

    ## 爬取3种GO分类
    t1 <- NULL
    for(i in 1:length(go.types)){ # i = 2
      go.type <- go.types[i]

      ## 网址
      url <- paste0("http://asia.ensembl.org/Homo_sapiens/Component/Gene/Ontologies/",go.type,"/go?db=core;g=",g,";time=")


      ## 爬取网址构建
      url <- paste0(url,unix.time,us)

      ## 表格
      tables <- readHTMLTable(url)
      if(length(tables) == 0){
        #说明是空值，没有对应的go值
        table <- data.frame(a1="NotAvailable",
                            a2="NotAvailable",
                            a3="NotAvailable",
                            a4="NotAvailable",
                            a5="NotAvailable",
                            a6="NotAvailable")
        colnames(table) <- c("Accession","Term","Evidence","Annotation source","Mapped using","Transcript IDs")
      } else {
        table <- tables[[1]];table <- table[,-7]
      }
      table$`GO Type` <- go.types.ab[i]
      table$Gene <- g

      ## 合并
      colnames(table)
      colnames(t1)
      t1 <- rbind.fill(t1,table)
    }

    ## 数据整理
    e <- c("Gene","GO Type")
    t2 <- subset(t1,select = c(e,setdiff(colnames(t1),e)))
    a1 <- gsub("ENST","_ENST",as.character(t2$`Transcript IDs`))
    t2$`Transcript IDs` <- apply(as.matrix(a1),1,function(x)substring(x,2,nchar(x)))

    ## 输出结果
    print(paste0("GO terms of ",g," has been downloaded."))
    return(t2)

  }

  ## 获得多个基因的go term
  g1 <- as.matrix(genes)
  if(parallel == F){
    g2 <- adply(g1,1,.fun = get1)
    g3 <- g2[,-1]
  } else {
    print("ensemblGo:use parallel strategy...")
    nd <- c("foreach","doParallel");Plus.library(nd)
    ncore <- detectCores()
    cl <- makeCluster(mc <- getOption("cl.cores",ncore));registerDoParallel(cl)
    clusterExport(cl,c("get1","go.types","go.types.ab","readHTMLTable"),envir = environment())
    options(warn = -1) #关闭警告信息
    g2 <- adply(g1,1,.fun = get1,.parallel = T)
    options(warn = 0) #默认警告信息
    stopCluster(cl)
    g3 <- g2[,-1]
  }

  ## 输出结果
  if(save.file == T){
    write.csv(g3,paste0(names,"_ensemblGo analysis.csv",rownames=F))
  }
  return(g3)

}


####======================ensemblGoEnrich=====================####
## ensemblGoEnrich help do a simple Go annotation for given genes based on ensemblGo
# genes # a character/factor vector of genes
# parallel # whether use parallel strategy to accelerate the process
#' @export
ensemblGoEnrich <- function(genes,
                            parallel = F,
                            save.file = F,
                            names = "love"){

  ## 检查是否有重复值
  if(length(genes) > length(unique(genes))){
    #有重复值
    genes <- unique(as.character(genes))
    print(paste0("Genes repeat.Only ",length(genes)," genes would go into the next process."))
  } else {
    genes <- as.character(genes)
  }

  ## 加载必要的包
  nd <- c("plyr")
  Plus.library(nd)

  ## 获得所有基因的GO注释数据
  df1 <- ensemblGo(genes = genes,parallel = T)

  ## report
  test1 <- df1[df1$Accession == "NotAvailable",]
  if(nrow(test1)==0){
    print("All genes have Go annotation.")
  } else {
    t1 <- length(unique(test1$Gene))
    print(paste0("There are ",floor(t1)," genes with no Go annotation."))
  }

  ## 获得正确注释的基因数据
  df2 <- df1[df1$Accession != "NotAvailable",]
  colnames(df2)
  set.genes <- subset(df2,select = c("Gene"))
  set.other <- subset(df2,select = setdiff(colnames(df2),c("Accession","Gene")))

  ## 对于重复的GO进行合并处理
  print("merge repeated GO term in gene levels...")
  get.genes <- function(vt){
    vt1 <- paste(vt,collapse = "_")
    return(vt1)
  }
  get.other <- function(vt){
    vt1 <- unique(vt)
    vt2 <- paste(vt1,collapse = "_")
    return(vt2)
  }
  go.id <- factor(df2$Accession)
  if(parallel == F){
    x1 <- apply(set.genes,2,function(x)tapply(x,go.id,get.genes))
    x2 <- apply(set.other,2,function(x)tapply(x,go.id,get.other))
  } else {
    print("ensemblGoEnrich:use parallel strategy")
    Plus.library("parallel")
    ncore <- detectCores()
    cl <- makeCluster(mc <- getOption("cl.cores",ncore))
    clusterExport(cl,c("go.id","get.genes","get.other"),envir = environment())
    x1 <- parApply(cl=cl,set.genes,2,function(x)tapply(x,go.id,get.genes))
    x2 <- parApply(cl=cl,set.other,2,function(x)tapply(x,go.id,get.other))
    stopCluster(cl)
  }
  df3 <- cbind(x1,x2)

  ## 整理合并信息
  df3 <- as.data.frame(df3)
  df3$Accession <- rownames(df3)
  df3$GeneCount <- apply(as.matrix(df3$Gene),1,function(x)length(Fastextra(x,"_")))
  get.symbol <- function(x){
    s1 <- convert(Fastextra(x,"_"))
    s2 <- paste0(s1,collapse = "_")
    return(s2)
  }
  df3$GeneSymbol <- apply(as.matrix(df3$Gene),1,get.symbol)
  colnames(df3)
  s1 <- c("Accession","GO Type","Term","GeneCount","GeneSymbol","Evidence","Annotation source","Mapped using","Gene","Transcript IDs")
  df3 <- subset(df3,select = s1)

  ## 按基因数目排序
  df3 <- df3[order(df3$`GO Type`,df3$GeneCount,decreasing = T),]

  ## 输出数据
  l <- list(
    merge = df3,
    metadata = df2
  )
  if(save.file == T){
    write.csv(df3,paste0(names,"_ensemblGoEnrich analysis.csv"),row.names = F)
  }
  print("All done!")
  return(l)

}
