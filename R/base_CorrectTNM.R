

#' @title Version control for TCGA TNM stage
#' @description Version control for TCGA TNM stage
#' @param data design object
#' @param version the colname of AJCC stage version.Note that the element of version should be "6th" like.
#' @param toVersion the final version you want to change.
#' @param T.status the colnames of T stage
#' @param Np.counts the colnames of positive lymph nodes
#' @param M.status the colnames of M stage
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
                       Np.counts = "Np.count",
                       M.status = "M.status",
                       correctFile = NULL){
  ## Package
  nd <- c("readxl","plyr")
  Plus.library(nd)

  ## read corrected file
  if(is.null(correctFile)){
    #correctFile <- "E:/RCloud/database/DataDownload/TCGA-TCGABiolinks/STAD_Download/STAD_AJCC stage.xlsx"
    correctFile <-  system.file("extdata", "STAD_AJCC stage.xlsx", package = "lucky")
  }

  ## toVersion
  toVersion <- paste0("Ver.",toVersion)

  ## 注释文件
  db.t <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 1,na = "NA"));  db.n <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 2,na = "NA"));  db.m <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 3,na = "NA"));
  db.stage <- as.data.frame(readxl::read_xlsx(path = correctFile,sheet = 4,na = "NA"))
  db.stage <- db.stage[db.stage$Version %in% toVersion,]

  ## 提取数据
  data2 <- data[,c(version,T.status,Np.counts,M.status)]

  ## 对不同的version进行转换
  data3 <- adply(data2,1,function(x)finalStage(x,toVersion,db.t,db.n,db.m,db.stage))
  data3 <- data3[,c(5:8)]
  rownames(data3) <- rownames(data2)

  ## 输出结果
  data.x <- data
  data.x$T.stage <- data3$T.stage
  data.x$N.stage <- data3$N.stage
  data.x$M.stage <- data3$M.stage
  data.x$Stage <- data3$Stage
  return(data.x)

}


## 对给定的tnm分期，计算stage
getStage <- function(x,db.stage){ # x=df.v[1,]
  p <- 1:nrow(db.stage)
  # 寻找tnm分期所在位置
  x.t <- p[db.stage$T.stage %in% as.character(x[1])]
  x.n <- p[db.stage$N.stage %in% as.character(x[2])]
  x.m <- p[db.stage$M.stage %in% as.character(x[3])]

  # 位置交集
  a1 <- intersect(x.t,x.n)
  a2 <- intersect(a1,x.m)

  # 输出stage
  if(length(a2) == 0) {
    stage <- NA
  } else {
    stage <- as.character(db.stage$Stage[a2[1]])
  }

  # 输出结果
  return(stage)

}


## 对给定的版本号、T分期、N阳性淋巴结数、M分期，计算stage
finalStage <- function(x,toVersion,db.t,db.n,db.m,db.stage){ # x=data2[234,]

  x <- as.character(as.matrix(x))

  v.type.i <- x[1]
  v.type.i2 <- Fastextra(v.type.i,"th",1)
  v.type.i2 <- paste0("Ver.",v.type.i2)

  ## 提取多种分期数据
  t <- x[2]; Np.c <- as.integer(x[3]); m <- x[4]

  ## 转换T分期
  t2 <- convert(t,v.type.i2,toVersion,db.t)

  ## 转换N分期
  n2 <- convert(Np.c,"Np.counts",toVersion,db.n)

  ## 转换M分期
  m2 <- convert(m,v.type.i2,toVersion,db.m)

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

