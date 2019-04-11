

#' @title Version control for TCGA TNM stage
#' @description Version control for TCGA TNM stage
#' @param data design object
#' @param version the colname of AJCC stage version.Note that the element of version should be "6th" like.
#' @param toVersion the final version you want to change.
#' @param T.status the colnames of T stage
#' @param N.status the colnames of N stage
#' @param M.status the colnames of M stage
#' @param Np.counts the colnames of positive lymph nodes
#' @param correctFile the path of reference database.Default = \code{system.file("extdata", "STAD_AJCC stage.xlsx", package = "lucky")}.It supports self-defined dataset with same format.
#' @importFrom readxl read_xlsx
#' @importFrom plyr adply
#' @details The most part of correction is the data from \code{correctFile}. It should be kept latest.\cr
#' For gastric cancer,before 5th(including 4th),the nodal involvement is classified into pN0, pN1, and pN2, based on the site of the metastasis in relation to the primary tumor.Afterwards,the number of positive lymph nodes become the details to stage N.
#' @return a data frame
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @seealso \href{https://link.springer.com/book/10.1007/978-3-642-82982-6}{UICC 4th}
#' @examples
#' ## get your design data
#' load("E:/RCloud/database/DataDownload/TCGA-TCGABiolinks/STAD_Download/2018-12-3/design_STAD.rda")
#' colnames(design)
#' data = design
#' ## update the TNM stage to 8th version
#' data2 <- CorrectTNM(data,
#'                     version = "stageVesion",
#'                     toVersion = 8,
#'                     T.status = "pT",
#'                     Np.counts = "Np.count",
#'                     M.status = "M.status")
#' @export
CorrectTNM <- function(data,
                       version = "stageVesion",
                       toVersion = 8,
                       T.status = "pT",
                       N.status = "pN",
                       M.status = "M.status",
                       Np.counts = "Np.count",
                       correctFile = NULL){
  ## Package
  nd <- c("readxl","plyr")
  Plus.library(nd)

  ## read corrected file
  if(is.null(correctFile)){
    correctFile <-  system.file("extdata", "STAD_AJCC stage.xlsx", package = "lucky")
    #correctFile <- "E:/RCloud/RFactory/lucky/inst/extdata/STAD_AJCC stage.xlsx"
  }

  ## toVersion
  toVersion <- paste0("Ver.",toVersion)

  ## 注释文件
  db.t <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 1,na = "NA"));  db.n <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 2,na = "NA"));  db.m <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 3,na = "NA"));
  db.stage <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 4,na = "NA"))
  db.stage <- db.stage[db.stage$Version %in% toVersion,]

  ## 提取数据
  data2 <- data[,c(version,T.status,N.status,Np.counts,M.status)]


  ## 对不同的version进行转换
  data3 <- adply(data2,1,function(x)finalStage(x,toVersion,db.t,db.n,db.m,db.stage)) # colnames(data3)
  data3 <- data3[c("T.stage","N.stage","M.stage","Stage" )]
  rownames(data3) <- rownames(data2)

  ## 输出结果
  data.x <- data
  data.x$T.stage <- data3$T.stage
  data.x$N.stage <- data3$N.stage
  data.x$M.stage <- data3$M.stage
  data.x$Stage <- data3$Stage
  return(data.x)

}


####=====================辅助函数=======================####

## convert的变种，专门用于CorrectTNM的注释
convert_mutiple <- function(vt,
                            fromtype="ENSEMBL",
                            totype="SYMBOL",
                            db=common.annot,
                            merge = T){
  s <- c("Tx","TX","Nx","NX","Mx","MX",NA)
  if(vt %in% s){
    #此为信息缺失的数据，直接使用而无需注释
    t <- vt
  } else {
    ## Get information
    f <- as.character(db[,fromtype])
    t <- db[f %in% vt,totype]
    t <- unique(t)

    if(length(t)==0){
      #空值。记为NA
      t <- NA
    } else {
      if(length(t)==1){
        # 一个数值为正解
        t <- t
      } else {
        # 多个数值。是否整合？
        if(merge){
          #整合
          t <- paste0(t,collapse = "_")
        } else {
          #不整合
          t <- t
        }
      }
    }


  }



  ## Output result
  return(t)

}


## 对给定的tnm分期，计算stage
getStage <- function(x,db.stage){ # x=df.v[1,]

  ## 校正多分期数据
  MtoS <- function(x){ # m=x[2]
    l <- NULL;s <- 1
    for(i in 1:length(x)){ # i=2
      m <- x[i]
      m1 <- Fastextra(as.character(as.matrix(m)),"_")
      s <- s*length(m1)
      l <- c(l,list(m1));names(l)[i] <- names(m)
    }

    l2 <- NULL
    for(i in 1:length(l)){ # i=2
      l.i <- l[[i]]
      l2.i <- rep(l.i,s/length(l.i))
      l2 <- c(l2,list(l2.i));names(l2)[i] <- names(l)[i]
    }

    df.x <- as.data.frame(l2,stringsAsFactors = F)

    return(df.x)

  }
  x2 <- MtoS(x)

  ## 校正后的单分期数据
  s <- c("Tx","TX","Nx","NX","Mx","MX",NA)
  Stage <- NULL
  for(i in 1:nrow(x2)){ # i=1
    x.i <- x2[i,]
    ## 检验是否有空值/M1等
    test <- sum(as.vector(as.matrix(x.i)) %in% s)
    if(test > 0){
      #说明有某些值为"Tx","TX","Nx","NX","Mx","MX",NA。此为无效数值
      stage <- "UK" # unknown
    } else {
      #说明没有异常值。此时观察是否有M1存在。因为M1存在时，分期一定为IV期
      test <- sum(as.vector(as.matrix(x.i)) %in% "M1")
      if(test > 0){
        #M1期
        stage <- "IV"
      } else {
        #M0期
        p <- 1:nrow(db.stage)
        # 寻找tnm分期所在位置
        x.t <- p[db.stage$T.stage %in% as.character(x.i[1])]
        x.n <- p[db.stage$N.stage %in% as.character(x.i[2])]
        x.m <- p[db.stage$M.stage %in% as.character(x.i[3])]

        # 位置交集
        a1 <- intersect(x.t,x.n)
        a2 <- intersect(a1,x.m)

        stage <- as.character(db.stage$Stage[a2[1]])
      }

    }

    ## 合并分期
    Stage <- c(Stage,stage)
  }

  ## 合并分期
  Stage2 <- paste0(unique(Stage),collapse = "_")

  # 输出结果
  return(Stage2)
}


## 对给定的版本号、T分期、N阳性淋巴结数、M分期，计算stage
finalStage <- function(x,toVersion,db.t,db.n,db.m,db.stage){ # x=data2[234,];x=data2[22,];x=data2[grep("TCGA-CD-5799",rownames(data2)),]

  x <- as.character(as.matrix(x))

  v.type.i <- x[1]
  v.type.i2 <- Fastextra(v.type.i,"th",1)
  v.type.i2 <- paste0("Ver.",v.type.i2)

  ## 提取多种分期数据
  t <- x[2]; n <- x[3];Np.c <- as.integer(x[4]);m <- x[5]

  ## 转换T分期
  t2 <- convert_mutiple(t,v.type.i2,toVersion,db.t)

  ## 转换N分期
  n2 <- convert_mutiple(Np.c,"Np.counts",toVersion,db.n)
  if(is.na(n2) & !is.na(n)){
    #Np无法提供N分期且N.status可以提供。由N.staus提供
    #grep("N1",data2$pN)
    n2 <- convert_mutiple(n,v.type.i2,toVersion,db.n)
  }

  ## 转换M分期
  m2 <- convert_mutiple(m,v.type.i2,toVersion,db.m)

  ## 汇总结果
  df.v <- data.frame(
    T.stage = t2,
    N.stage = n2,
    M.stage = m2
  )

  df.v$Stage <- apply(df.v,1,function(x) getStage(x,db.stage))

  ## 输出结果
  return(df.v)

}







