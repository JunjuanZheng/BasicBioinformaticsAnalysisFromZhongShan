




#' @title a simple method for ceRNA network edge/node selection
#' @description  a simple method for ceRNA network edge/node selection
#' @inheritParams FastCeRNA_edge
#' @param filter whether to do filter
#' @param pmcor the col of protein-coding and miRNA correlation
#' @param lmcor the col of lncRNA and miRNA correlation
#' @param plcor the col of protein-coding and miRNA correlation
#' @param verbose whether do verbose report
#' @param cutoff the cutoff of correlation.Default is \code{list(pl=0.5,lm=0.2,pm=0.2)}
#' @param pmcor.Pvalue the p value of pmcor
#' @importFrom plyr adply
#' @return a data frame
#' @details FastCeRNA_edge2 select ceRNA network:positive lnc-miRNA correaltion;miR negative.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' load("E:/iProjects/LM_WGCNA/data/LM_WGCNA_02C/rda/lm_ceRNA.rda")
#' df <- FastCeRNA_edge2(object)
#'
#' ## other parameters
#' colnames(lm_ceRNA)
#' object = lm_ceRNA
#' lncRNA.col = "lncRNAs"
#' pcRNA.col = "Genes"
#' miRNA.col = "miRNAs"
#' pmcor = "pmcor"
#' lmcor="lmcor"
#' plcor = "plcor"
#' filter = T
#' verbose = T
#' @export
FastCeRNA_edge2 <- function(object,
                            lncRNA.col = "lncRNAs",
                            pcRNA.col = "Genes",
                            miRNA.col = "miRNAs",
                            pmcor = "pmcor",
                            pmcor.Pvalue = "corPValue",
                            lmcor="lmcor",
                            plcor = "plcor",
                            cutoff = list(pl=0.5,
                                          lm=0.2,
                                          pm=0.2),
                            filter = T,
                            verbose = T){
  ## needed packages
  nd <- c("plyr");Plus.library(nd)

  ## select cols
  s <- c(lncRNA.col,pcRNA.col,miRNA.col,pmcor,lmcor,plcor,pmcor.Pvalue)
  df1 <- object[s]

  ## extra correlation data
  # x <- df1[1,]
  get_all_cor <- function(x){
    ##
    data_1 <- data.frame(
      lncRNA = as.character(x[1]),
      Genes = as.character(x[2]),
      miRNAs = Fastextra(as.character(x[3]),","),
      pmcor=as.numeric(Fastextra(as.character(x[4]),",")),
      lmcor=as.numeric(Fastextra(as.character(x[5]),",")),
      plcor=as.numeric(as.character(x[6])),
      p = as.numeric(as.character(x[7])),
      stringsAsFactors = F
    )
    colnames(data_1) <- paste(s,1,sep = "_")
    return(data_1)
  }
  if(verbose == T) {LuckyVerbose("Step1: Extra correlation data...")}
  df2 <- adply(df1,1,get_all_cor)
  df2 <- df2[8:14]
  colnames(df2) <- s

  ## whether do filter
  # x <- df2[1,]
  get_filter <- function(x){
    x <- as.character(x)
    x6 <- as.numeric(as.character(x[6]))
    x4 <- as.numeric(as.character(x[4]))
    x5 <- as.numeric(as.character(x[5]))
    x7 <- as.numeric(as.character(x[7]))
    s1 <- all((x6 > 0),(x4 < 0),(x5 < 0),(x7 < 0.05))
    s2 <- all((abs(x6) >= cutoff$pl),(abs(x4) >= cutoff$pm),(abs(x5) >= cutoff$lm))
    lg <- all(s1,s2);lg
    return(lg)
  }
  if(filter == T){
   if(verbose == T) {LuckyVerbose("Step2: Use ceRNA filter...")}
   lg <- apply(df2,1,get_filter)
   df3 <-  df2[lg,]
  } else {
    if(verbose == T) {LuckyVerbose("Step2: Not use ceRNA filter...")}
   df3 <- df2
  }

  ## Delete repeat data
  del <- paste(df3[,1],df3[,2],df3[,3],sep = "_") # table(!duplicated(del))
  df4 <- df3[!duplicated(del),]
  df4 <- df4[order(df4[,6],decreasing = T),]

  ## Get edges and nodes
  #x <- df4[1,]
  get_edges_nodes <- function(x){
    x1 <- data.frame(
      sourceRNA = as.character(x[1:2]),
      targetRNA =  as.character(x[3]),
      cor = as.numeric(as.character(x[5:4])),
      pl = as.numeric(as.character(x[6])),
      pl.P =as.numeric(as.character(x[7])),
      stringsAsFactors = F
    )
    return(x1)
  }
  if(verbose == T) {LuckyVerbose("Step3: Get edges and nodes...")}
  df5 <- adply(df4,1,get_edges_nodes)
  df5 <- df5[8:12]

  ## output data
  if(verbose == T) LuckyVerbose("All done!")
  return(df5)

}









