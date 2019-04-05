


##
#' @title Fast Molecular Signatures Analysis via clusterProfilter
#' @description Fast Molecular Signatures Analysis via clusterProfilter
#' @param genes a character of genes.The ENSEMBL id is the best.Note that no replicate should exist.
#' @param geneList a named numeric of genes.The ENSEMBL id is the best.Note that no replicate should exist.
#' @param pvalueCutoff  a list of pvalue adjust method for \code{\link[clusterProfiler]{enricher}} and \code{\link[clusterProfiler]{GSEA}}
#' @inheritParams FastGO
#' @importFrom clusterProfiler enricher GSEA
#' @importFrom enrichplot dotplot emapplot goplot cnetplot gseaplot
#' @importFrom ggplot2 ggtitle
#' @importFrom stringr str_c
#' @importFrom dplyr filter select
#' @importFrom tidyr  %>%
#' @seealso \code{\link[clusterProfiler]{enricher}}; \code{\link[clusterProfiler]{GSEA}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data(geneList, package = "DOSE")
#' genes = names(geneList)[1:100]
#' geneList = geneList
#' ms <- FastMS(genes,
#'              geneList,
#'              pAdjustMethod = "BH",
#'              pvalueCutoff = list(enrichMS=0.00001,
#'                                  gseMS = 0.00001),
#'              qvalueCutoff  = 0.2,
#'              cnet.showCategory = 5,
#'              verbose = TRUE,
#'              save.path = "Molecular Signatures",
#'              names = "love")
#' @export
FastMS <- function(genes,
                   geneList,
                   default.universe = T,
                   pAdjustMethod = "BH",
                   pvalueCutoff = list(enrichMS=0.05,
                                       gseMS = 0.05),
                   qvalueCutoff  = 0.2,
                   cnet.showCategory = 5,
                   verbose = TRUE,
                   save.path = "Molecular Signatures",
                   names = "love"){
  ## 包
  #nd <- c("clusterProfiler","enrichplot","ggplot2","stringr","msigdf","dplyr","tidyr");Plus.library(nd)

  ## Test whether the msigdf package had been installed
  x <- installed.packages();x <- as.data.frame(x)
  installed <- as.character(x$Package)
  if(!"msigdf" %in%  installed){
    Plus.library("devtools")
    LuckyVerbose("Attention!You haven't intall 'msigdf' package.Please use the next code to install it:" );LuckyVerbose("devtools::install_github(\"ToledoEM/msigdf\")");LuckyVerbose("FastMS pause.→_→")
  } else {
    ## FastMS pipeline
    #devtools::install_github("ToledoEM/msigdf")
    data("msigdf.human",package = "msigdf",envir = environment())

    ## 产生储存文件夹
    old <- options()
    dir <- paste0("./",names,"/",save.path,"/")
    options(warn = -1)
    dir.create(dir,recursive = T)
    options(warn = old$warn)


    ## geneList处理和排序
    LuckyVerbose("Step1: 'geneList' process...")
    universe <- names(geneList)
    lg1 <- grep("ENSG",universe);lg1 <- length(lg1) == 0
    if(lg1){
      ## 是ENTREZID id
      lg2 <- Fastgrep(c(letters,LETTERS),universe);lg2 <- length(lg2) == 0
      if(lg2){
        ## 不是symbol。应该转为symbol
        LuckyVerbose("geneList: ENTREZ ID to SYMBOL ID...",levels = 2)
        LuckyVerbose("genes: ENTREZ ID to SYMBOL ID...",levels = 2)
        symbol <- convert(universe,
                          fromtype = "ENTREZID",
                          totype = "SYMBOL")
        na.count <- length(symbol[is.na(symbol)])
        LuckyVerbose("There are",na.count,"genes with no SYMBOL...",levels = 2)
        geneList <- geneList[!is.na(symbol)]
        genes <- genes[genes %in% names(geneList)]
        genes  <- convert(genes,
                          fromtype = "ENTREZID",
                          totype = "SYMBOL")
        names(geneList) <- convert(names(geneList),
                                   fromtype = "ENTREZID",
                                   totype = "SYMBOL")
      }
    } else {
      ## 属于ENSEMBL id.一定有SYMBOL值
      LuckyVerbose("geneList: ENSEMBL ID to SYMBOL...",levels = 2)
      LuckyVerbose("genes: ENSEMBL ID to SYMBOL...",levels = 2)
      names(geneList) <- convert(universe)
      genes <- convert(genes)
    }

    ## repeat test
    #如果转换成entrezid后有重复者，则对genList取mean值
    repeat.test <- table(duplicated(names(geneList)))
    if(T %in% names(repeat.test)){
      r.n <- repeat.test[names(repeat.test) %in% T]
      LuckyVerbose(paste0("There are ",r.n," repeats in geneList and merge them by mean stragegy."),levels = 2)
      geneList2 <- tapply(geneList,names(geneList),mean)
      geneList <- as.numeric(geneList2)
      names(geneList) <- names(geneList2)
      geneList <- sort(geneList,decreasing = T)
    } else {
      geneList <- sort(geneList,decreasing = T)
    }

    ### MS pipeline
    LuckyVerbose("Step2:Molecular Signatures pipeline...")
    ms.type <- c("c1","c2","c3","c4","c5","c6","c7","hallmark")
    MS <- NULL;MS.GSEA <- NULL
    for(i in 1:length(ms.type)){ # i=6
      msType <- ms.type[i]
      df_MS <-  dplyr::filter(msigdf.human,category_code == msType) %>% dplyr::select(geneset, symbol) %>% as.data.frame
      ## 超几何分布
      LuckyVerbose(msType,":hypergeometric distribution test...",levels = 2)
      if(universe == F){
        msObject <- enricher(genes,
                             pvalueCutoff = pvalueCutoff$enrichMS,
                             qvalueCutoff  = qvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             TERM2GENE = df_MS)
      } else {
        msObject <- enricher(genes,
                             universe = names(geneList),
                             pvalueCutoff = pvalueCutoff$enrichMS,
                             qvalueCutoff  = qvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             TERM2GENE = df_MS)
      }
      MS <- c(MS,list(msObject))
      names(MS)[i] <- msType
      if(verbose == T) LuckyVerbose(paste0(msType,": MS object completed!"),levels = 2)

      ## save result
      result.i <- as.data.frame(msObject)
      result.i.path <- str_c(dir,names,"_",msType,"_Molecular Signatures enrichment.csv")
      write.csv(result.i,result.i.path)

      ## Plot
      if(nrow(result.i) == 0){
        if(verbose == T) LuckyVerbose(paste0(msType,": No MS enrichment"),levels = 2)
      } else {
        #enrichGo plot
        plot.names <- str_c(dir,names,"_",msType,"_enrichMS plot.pdf")
        pdf(plot.names,width = 10,height = 10)
        print(barplot(msObject, showCategory=10))
        print(dotplot(msObject, showCategory=10))
        #print(emapplot(msObject)) #有时会出错。暂不使用
        dev.off()
        if(verbose == T) LuckyVerbose(str_c(msType,": Get enrichMS plot!"),levels = 2)
      }

      ### MS GSEA
      LuckyVerbose(msType,": GSEA of Molecular Signatures...")
      gsea <- GSEA(geneList = geneList,
                   TERM2GENE = df_MS,
                   pvalueCutoff = pvalueCutoff$gseMS,
                   verbose = F)
      MS.GSEA <- c(MS.GSEA,list(gsea))
      names(MS.GSEA)[i] <- msType
      result.gsea <- as.data.frame(gsea)

      ## gsea plot
      if(nrow(result.gsea) == 0){
        if(verbose == TRUE) LuckyVerbose(paste0("gseMS-",msType,": NO GSEA enrichment."),levels = 2)
      } else {
        # gseaplot
        pdf(paste0(dir,names,"_MS_",msType,"_gseaplot_Gene Set Enrichment Analysis.pdf"),10,8)
        for(z in 1:length(gsea$ID)){ # z=1
          ## gsea view
          id <- gsea$ID[z]
          title <- gsea$Description[z]
          p <- gseaplot(gsea,
                        geneSetID=id,
                        title = title)
          print(p)
          if(verbose == TRUE) LuckyVerbose(paste0(msType,": GSEAplot: ",id," completed!"),levels = 3)
        }
        dev.off()
      }

    }

    ##Cnetplot
    pdf.pathway <- str_c(dir,names,"_MS_Get Gene-Concept Network plot.pdf");
    p_MS <- NULL
    for(i in 1:length(MS)){ # i=1
      df <- MS[[i]]
      if(nrow(df) > 0){
        p <- cnetplot(MS[[i]],
                      showCategory = cnet.showCategory,
                      categorySize="pvalue",
                      foldChange = geneList)
        p <- p + ggtitle(names(MS)[i])
        p_MS <- c(p_MS,list(p))
      }
    }
    pdf(pdf.pathway,width = 20,height = 20)
    for(i in 1:length(p_MS)) print(p_MS[[i]])
    dev.off()
    LuckyVerbose(str_c("Get Gene-Concept Network!"))

    ## Output data
    l <- list(
      Repeat = list(
        genes = genes,
        geneList = geneList,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        cnet.showCategory = cnet.showCategory,
        save.path = save.path,
        names = names
      ),
      Data = list(
        enrichMS = MS ,
        gseMS = MS.GSEA
      )
    )
    LuckyVerbose("All done!")
    return(l)
}









  }





