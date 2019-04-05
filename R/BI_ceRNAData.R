

#' @title Get lnRNA-miRNA/mRNA-miRNA interaction from starBase
#' @description Get lnRNA-miRNA/mRNA-miRNA interaction from starBase
#' @param object a name list(lncRNA/mRNA) of lncRNA and mRNAs.ENSEMBL id or SYMBOL id are all supported
#' @importFrom httr GET content
#' @importFrom plyr adply alply aaply
#' @importFrom rvest %>%
#' @return a name list(lncRNA/mRNA) of their interactive miRNAs
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' object <- list(mRNA = c("TP53","TP53","NOTEXIST"),
#'                lncRNA = c("MALAT1","MALAT1","NOTEXIST"))
#' ceRNA <- ceRNAData(object)
#' @export
ceRNAData <- function(object,verbose = F){

  ## Package
  #nd <- c("httr","plyr","rvest")
  #Plus.library(nd)

  ## merge data
  merge <- NULL
  for(i in 1:length(object)){ # i=2
    o.i <- object[[i]]
    name.i <- names(object)[i]
    if(name.i=="mRNA"){
      get1 <- ceRNANet_mRNA
    } else {
      if(name.i=="lncRNA"){
        get1 <- ceRNANet_lncRNA
      } else {
        message("Please use one of 'lnRNA' and 'mRNA'.")
        break
      }
    }
    message("Start:",name.i)
    merge.i <- NULL
    for(z in o.i){
      lg <- grep("ENSG",z);lg <- length(lg) !=0
      if(lg){
        z <- convert(z,
                     fromtype = "ENSEMBL",
                     totype = "SYMBOL")
      }
      merge.z <- get1(z,verbose = verbose)
      merge.i <- c(merge.i,merge.z)
    }
    #merge.i <- aaply(as.matrix(o.i),1,get1)
    #names(merge.i) <- o.i
    merge <- c(merge,list(merge.i))
    names(merge)[i] <- name.i
  }

  ## 返回结果
  message("All done!")
  return(merge)
}


####=================other assistant function===================####
  statbase_parse <- function(sB.test){
    ## 去表头
    sb1 <- Fastextra(sB.test,"\n")
    sb1 <- sb1[-c(1,2,3)]

    ## split every row and get mutiple elements
    #sb1.i <- sb1[1];number.i = 1
    sb_split <- function(x){
      sb1.i <- x[,1]
      number.i <- x[,2]
      sb1.i.split <- Fastextra(sb1.i,"\t")
      df_1 <- data.frame(sb1.i.split,stringsAsFactors = F)
      colnames(df_1)[1] <- number.i
      df_2 <- as.data.frame(t(df_1))
      return(df_2)
    }

    ## plyr:merge data
    df_sb <- data.frame(
      sb1 = sb1,
      number = 1:length(sb1)
    )
    df_sb2 <- adply(.data = df_sb,
                          .margins = 1,
                          .fun = sb_split)
    df_sb3 <- df_sb2[-c(1,2)]
    colnames(df_sb3) <- as.character(as.matrix(df_sb3[1,]))
    df_sb3 <- df_sb3[-1,]

    ## output
    return(df_sb3)

  }
  #df1 <- statbase_parse(sB.test)

  ### Get miRNA-mRNA network
  #mRNA = "TP53";clipExpNum=1;degraExpNum=0;pancancerNum=0;programNum = 0;program = "None"   #mRNA = "ABCDGER";
  ceRNANet_mRNA <- function(mRNA,
                            clipExpNum=1,
                            degraExpNum=0,
                            pancancerNum=0,
                            programNum = 0,
                            program = "None",
                            verbose = T){
    api <- paste0("http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=mRNA&miRNA=all&clipExpNum=",clipExpNum,"&degraExpNum=",degraExpNum,"&pancancerNum=",pancancerNum,"&programNum=",programNum,"&program=",program,"&target=",mRNA)
    sB.test <- GET(api) %>% content(as = "text",encoding = "UTF-8")
    miR_mRNA <- statbase_parse(sB.test)
    lg <- grep("not available",as.character(miR_mRNA[1,1]))
    lg <- length(lg) == 0
    if(lg){
      l <- list(as.character(miR_mRNA$miRNAname))
      names(l)[1] <- unique(as.character(miR_mRNA$geneID))
      if(verbose == T) cat("Done:",mRNA,"\n")
    } else {
      l <- list(NA)
      names(l)[1] <- mRNA
      if(verbose == T) cat("Not available:",mRNA,"\n")
    }
    return(l)
  }
  # miR_mRNA <- ceRNANet_mRNA(mRNA)


  ### Get miRNA-mRNA network
  #lncRNA = "MALAT1"
  ceRNANet_lncRNA <- function(lncRNA,
                              clipExpNum=1,
                              degraExpNum=0,
                              pancancerNum=0,
                              verbose = T){
    api <- paste0("http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=all&clipExpNum=",clipExpNum,"&degraExpNum=",degraExpNum,"&pancancerNum=",pancancerNum,"&target=",lncRNA)
    sB.test <- GET(api) %>% content(as = "text",encoding = "UTF-8")
    miR_lncRNA <- statbase_parse(sB.test)
    lg <- grep("not available",as.character(miR_lncRNA[1,1]))
    lg <- length(lg) == 0
    if(lg){
      l <- list(as.character(miR_lncRNA$miRNAname))
      names(l) <- unique(as.character(miR_lncRNA$geneID))
      if(verbose == T) cat("Done:",lncRNA,"\n")
    } else {
      l <- list(NA)
      names(l)[1] <- lncRNA
      if(verbose == T) cat("Not available:",lncRNA,"\n")
    }
    return(l)
  }
  #miR_lncRNA <- ceRNANet_lncRNA("MALAT1")


