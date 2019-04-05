

#' @title Get ceRNA network from starBase 3.0 by giving miRNAs symbol
#' @description Get ceRNA network by giving miRNAs symbol
#' @param assembly genome reference of a species.Default is homo spacies("hg19").See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param geneType "mRNA","lncRNA","pseudogene","sncRNA".See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param ceRNA If you set a lncRNA as ceRNA,you should set \code{geneType="lncRNA"}
#' @param family miRNA family names.See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param miRNAnum Default is 2.See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param pval Default is 0.01.See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param fdr Default is 0.01.See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param pancancerNum Default is 0.See \href{http://starbase.sysu.edu.cn/}{starbase}
#' @param verbose whether to show verbose report
#' @importFrom httr GET content
#' @importFrom plyr adply
#' @importFrom rvest %>%
#' @return a data frame with ceRNA interaction information
#' @seealso \code{starbase}\url{http://starbase.sysu.edu.cn/}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## Get specified miRNA based ceRNA network
#' ceRNA = c("ENSG00000228630","ENSG00000228630","FKDKEK","FKDKEK")
#' family = "all"
#' # Or:
#' ceRNA = c("ENSG00000228630","ENSG00000228630","FKDKEK","FKDKEK")
#' family = "ABCKD" # Here,family parameter would be corrected to "all"
#'
#' ## Get specified miRNA based ceRNA network
#' ceRNA <- "all"
#' family <-  c("miR-1252-5p","miR-200bc-3p/429")
#'
#' ## Get specified lncRNA-miRNA based ceRNA network
#' ceRNA = c("ENSG00000228630")
#' family <-  c("miR-1252-5p","miR-200bc-3p/429")
#'
#' ## Get all lncRNA based ceRNA network
#' ceRNA <- "all"
#' family <-  "all"
#' net_ceRNA <- ceRNANet(assembly = "hg19",
#'                       geneType="lncRNA",
#'                       ceRNA,
#'                       miRNAnum=2,
#'                       family,
#'                       pval=0.01,
#'                       fdr=0.01,
#'                       pancancerNum=0,
#'                       verbose = T)
#' @export
ceRNANet <- function(assembly = "hg19",
                     geneType="lncRNA",
                     ceRNA="HOTAIR",
                     miRNAnum=2,
                     family = "all",
                     pval=0.01,
                     fdr=0.01,
                     pancancerNum=0,
                     verbose = T){
  ## Test parameters
  if(!"all" %in% ceRNA){
    if(!"all" %in% family){
      LuckyVerbose("Attention! ceRNA is specified so family would be set 'all'.")
      family <-  "all"
    }
  }

  ## Get data from one ceRNA
  get1 <- function(x){
    x1 <- ceRNANet_1(assembly = assembly,
               geneType=geneType,
               ceRNA=x,
               miRNAnum=miRNAnum,
               family = family,
               pval=pval,
               fdr=fdr,
               pancancerNum=pancancerNum,
               verbose = verbose)
    return(x1)
  }

  ## ddply pipeline
  df1 <- adply(as.matrix(ceRNA),1,get1)
  df1 <- df1[-1]

  ## OutPut
  if(verbose == T) message("All done!")
  return(df1)

}


####=====================Assistant function==================####
ceRNANet_1 <- function(assembly = "hg19",
                       geneType="lncRNA",
                       ceRNA="all",
                       miRNAnum=2,
                       family = "all",
                       pval=0.01,
                       fdr=0.01,
                       pancancerNum=0,
                       verbose = T){
  ## Package
  #nd <- c("httr","rvest","plyr")
  #Plus.library(nd)

  ## Parse function
  statbase_parse2 <- function(sB.test){
    lg1 <- grep("not available",sB.test[1])
    lg1 <- length(lg1) == 0
    if(lg1){
      ###  有相应的数据
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
    } else {
      ### 未获得相应数据
      df_sb3 <- data.frame()
    }

    ## output
    return(df_sb3)
  }

  ## 情况一：寻找某些mRNA或者是lncRNA的ceRNA网络，此时family = "all"
  if("all" %in% family){
    ## ceRNA
    if(!"all" %in% ceRNA){
      lg1 <- grep("ENSG",ceRNA);lg1 <- length(lg1) !=0
      if(lg1){
        # ENSEMBL id
        ceRNA <- convert(ceRNA)
      } else {
        # SYMBOL
        ceRNA <- ceRNA
      }
      #ceRNA <- ceRNA %>% paste0(.,collapse = ",") %>% paste0('"',.,'"')
    } else {
      ceRNA <-"all"
    }

    ## Get data
    api <- paste0('http://starbase.sysu.edu.cn/api/ceRNA/?assembly=',assembly,'&geneType=',geneType,'&ceRNA=',ceRNA,'&miRNAnum=',miRNAnum,'&pval=',pval,'&fdr=',fdr,'&pancancerNum=',pancancerNum)
    if(verbose == T) message(ceRNA,": Getting text from API...")
    sB.test <- GET(api) %>% content(as = "text",encoding = "UTF-8")
    ## starBase parse
    if(verbose == T) message(ceRNA,": Parse information from starBase...")
    miR_ceRNA <- statbase_parse2(sB.test)
    if(nrow(miR_ceRNA) == 0 & verbose == T){
      LuckyVerbose(ceRNA,": Data Not Available!",levels = 2)
    }
  } else {
    LuckyVerbose("Note: miRNA-based ceRNA network is always time-consuming.Be patient.")
    ## Tidy up miRNA ids
    miRNA_right <- function(miRNA.i){
      lg <- grep("hsa-",miRNA.i);lg <- length(lg) !=0
      if(lg){
        ## 有hsa-头
        miRNA.i <- Fastextra(miRNA.i,"hsa-",2)
      }
      return(miRNA.i)
    }
    test <- sapply(family,miRNA_right)
    miRNA <- as.character(test)
    miRNA_symbol <- family <- paste(miRNA,collapse = ",")
    family <- family %>% paste0('\"',.,'\"')

    ## Get data
    api <- paste0('http://starbase.sysu.edu.cn/api/ceRNA/?assembly=',assembly,'&geneType=',geneType,'&ceRNA=',ceRNA,'&miRNAnum=',miRNAnum,'&family=',family,'&pval=',pval,'&fdr=',fdr,'&pancancerNum=',pancancerNum)
    if(verbose == T) message("Getting text from API: ",miRNA_symbol)
    sB.test <- GET(api) %>% content(as = "text",encoding = "UTF-8")
    ## starBase parse
    if(verbose == T) message("Parse information from starBase: ",miRNA_symbol)
    miR_ceRNA <- statbase_parse2(sB.test)
    if(nrow(miR_ceRNA) == 0 & verbose == T){
      LuckyVerbose(miRNA_symbol,": Data Not Available!",levels = 2)
    }
  }

  ## Output data
  return(miR_ceRNA)

}


