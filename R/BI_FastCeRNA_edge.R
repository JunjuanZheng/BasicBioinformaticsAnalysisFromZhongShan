

#' @title Get miRNA-lncRNA-mRNA edge data after FastCeRNA
#' @description Get miRNA-lncRNA-mRNA edge data after FastCeRNA
#' @param object the result of \code{\link{FastCeRNA}}
#' @param lncRNA.col colnames of lncRNA
#' @param pcRNA.col colnames of protein-coding RNA,
#' @param miRNA.col colnames of miRNA
#' @param Cor.col colnames of Pearson corelation
#' @importFrom plyr adply ddply
#' @return a data frame of miRNA-lncRNA-mRNA edge data
#' @seealso \code{\link{FastCeRNA}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## load data
#' load("E:/iProjects/LM_WGCNA/data/LM_WGCNA_01A/rda/lm_ceRNA.rda")
#' load("E:/iProjects/LM_WGCNA/data/LM_WGCNA_01A/rda/lm_moduleGenes.rda")
#'
#' ## filter
#' lm_ceRNA_hubModules <- dplyr::filter(lm_ceRNA,Genes %in% as.character(moduleGenes),hyperPValue <= 0.05)
#'
#' ## Get edge data
#' ceRNA_edge <- FastCeRNA_edge(object = lm_ceRNA_hubModules,
#'                              lncRNA.col = "lncRNAs",
#'                              pcRNA.col = "Genes",
#'                              miRNA.col = "miRNAs",
#'                              Cor.col = "cor")
#' ## Or easier way
#' ceRNA_edge <- FastCeRNA_edge(lm_ceRNA_hubModules)
#' View(ceRNA_edge)
#' @export
FastCeRNA_edge <- function(object,
                           lncRNA.col = "lncRNAs",
                           pcRNA.col = "Genes",
                           miRNA.col = "miRNAs",
                           Cor.col = "cor"){
  ## Package
  nd <- c("plyr")
  Plus.library(nd)

  ## select data
  s <- c(lncRNA.col,pcRNA.col,miRNA.col,Cor.col)
  df1 <- subset(object,select = s)

  ## Get single miRNA information
  df2 <- getMiRData(df1,
                    miRNA.col,
                    lncRNA.col,
                    pcRNA.col,
                    Cor.col)

  ## Delete duplicated records
  df3 <- ddply(df2,
               c(lncRNA.col,pcRNA.col,miRNA.col),
               .fun = function(x)getCor(x,Cor.col))
  colnames(df3)[4] <- Cor.col

  ## Get edge data
  df4 <- adply(df3,1,getEdge)
  df4 <- df4[4:6]

  ## Output
  return(df4)

}

####======================Assistant Function===============####

### Get single miRNA information in 'df1' dataset
getMiRData <- function(df1,
                       miRNA.col,
                       lncRNA.col,
                       pcRNA.col,
                       Cor.col){
  ## Get one row data from 'df1' dataset
  get1 <- function(x){
    mir <- Fastextra(as.character(x[,miRNA.col]),",")
    rowCounts <- length(mir)
    df_mir <- data.frame(
      lncRNAs = x[,lncRNA.col],
      Genes = x[,pcRNA.col],
      miRNAs = mir,
      cor = x[,Cor.col]
    )
    return(df_mir)
  }

  ## Multiple plyr::adply
  df_GMD <- adply(df1,1,get1)
  df_GMD <- df_GMD[-1]
  colnames(df_GMD) <- c(lncRNA.col,pcRNA.col,miRNA.col,Cor.col)

  ## Output
  return(df_GMD)
}


### Get single miRNA Cor in 'df2' dataset
getCor <- function(x,Cor.col){
  x1 <- as.character(x[,Cor.col])
  #x1 <- as.numeric(as.character(x))
  return(unique(x1))
}


### Get edge information from 'df3' dataset
getEdge <- function(x){ # x <- df3[2,]
  x1 <- data.frame(miRNAs = as.character(x[1,3]),
                   tagetRNAs = as.character(as.matrix(x[1,1:2])),
                   cor = as.numeric(as.character(x[1,4])),
                   stringsAsFactors = F)
  return(x1)

}








