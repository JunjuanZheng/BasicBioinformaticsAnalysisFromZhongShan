

#' @title Fast way to do GSEA for a geneList object
#' @description Fast way to do GSEA for a geneList object
#' @param pvalueCutoff a list of pvalue adjust method for \code{\link[clusterProfiler]{gseGO}},\code{\link[clusterProfiler]{gseKEGG}} and \code{\link[clusterProfiler]{gseMKEGG}}.Note that the default value is 0.01,you can set 0.05 if you want a robust exploration.
#' @inheritParams FastGO
#' @importFrom clusterProfiler gseMKEGG gseKEGG gseGO
#' @importFrom enrichplot gseaplot
#' @importFrom ggplot2 ggtitle
#' @importFrom stringr str_c
#' @seealso \code{\link[clusterProfiler]{gseGO}}; \code{\link[clusterProfiler]{gseKEGG}}; \code{\link[clusterProfiler]{gseMKEGG}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' #'import org.Hs.eg.db
#' data(geneList, package = "DOSE")
#' gsea <- FastGSEA(geneList,
#'                  OrgDb = NULL,
#'                  keyType = NULL,
#'                  pvalueCutoff = list(gseGO=0.05,
#'                                      gseKEGG = 0.05,
#'                                      gseMKEGG = 0.05),
#'                  verbose = T,
#'                  save.path = "GSEA",
#'                  names = "love")
#' @export
FastGSEA <- function(geneList,
                     OrgDb = NULL,
                     keyType = NULL,
                     pvalueCutoff = list(gseGO=0.01,
                                         gseKEGG = 0.01,
                                         gseMKEGG = 0.01),
                     verbose = T,
                     save.path = "GSEA",
                     names = "love"){
  ## Packages
  #nd <- c("clusterProfiler","org.Hs.eg.db","stringr","pathview","DOSE","ggplot2");
  #Plus.library(nd)


  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  pathview.dir <- paste0(dir,"/pathview")
  dir.create(pathview.dir,recursive = T)
  options(warn = old$warn)

  ## geneList 排序
  if(verbose == T){
    LuckyVerbose("Step1: Order 'geneList' from upper to lower...")
  }
  geneList <- sort(geneList,decreasing = T)

  ## whether is ensembl id
  if(is.null(keyType)){
    lg1 <- length(grep("ENSG",names(geneList))) == 0
    if(lg1){
      #说明不是ensembl id
      if(verbose == T) LuckyVerbose("The names of geneList is ENTREZ ID...",levels = 2)
      keyType <-  "ENTREZID"
    } else {
      #说明是ensembl id
      if(verbose == T) LuckyVerbose("The names of geneList is ENSEMBL ID and transfered to ENTREZID...",levels = 2)
      entrez.id <- convert(names(geneList),totype = "ENTREZID")
      geneList <- geneList[!is.na(entrez.id)]
      na.count <- length(geneList[is.na(entrez.id)])
      if(verbose == T) LuckyVerbose("There are",na.count,"geneList genes with no ENTREZ ID...",levels = 2)
      names(geneList) <- convert(names(geneList),totype = "ENTREZID")
      keyType <-  "ENTREZID"
    }
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

  ## reference database
  if(is.null(OrgDb)){
    OrgDb <- "org.Hs.eg.db"
    organism <- "hsa"
    LuckyVerbose("You don't specify the 'OrgDb' parameter.Here we use 'org.Hs.eg.db' and 'hsa'(homo spacies)...",levels = 2)
  }

  ###===============GO Gene Set Enrichment Analysis========###
  if(verbose == T) {
    LuckyVerbose("Step2: GO Gene Set Enrichment Analysis...")}
  ## 差异性基因的多个Categery的enrichGO analysis
  Go.categery <- c("MF","BP","CC")

  ## GO GSEA
  go.gsea <- NULL
  for(i in 1:length(Go.categery)){ # ont = "CC"
    ont <- Go.categery[i]
    if(verbose == TRUE){
      LuckyVerbose(paste0("gseGO-",ont,": Getting GSEA object..."),levels = 2)
    }
    go.gsea.i <- gseGO(geneList = geneList,
                       ont = ont,
                       OrgDb = OrgDb,
                       keyType = keyType,
                       exponent = 1,
                       nPerm = 1000,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = pvalueCutoff$gseGO,
                       verbose = F)
    go.gsea <- c(go.gsea,list(go.gsea.i));
    names(go.gsea)[i] <- ont
    go.gsea.i2 <- as.data.frame(go.gsea.i)
    if(nrow(go.gsea.i2) == 0){
      if(verbose == TRUE){
        LuckyVerbose(paste0("gseGO-",ont,": NO GSEA enrichment."),levels = 2)
      }
    } else {
      # gseaplot
      pdf(paste0(dir,names,"_GO_",ont,"_gseaplot_Gene Set Enrichment Analysis.pdf"),10,8)
      for(i in 1:length(go.gsea.i2$ID)){ # i=1
        ## gsea view
        id <- go.gsea.i2$ID[i]
        title <- go.gsea.i2$Description[i]
        p <- gseaplot(go.gsea.i,
                      geneSetID=id,
                      title = title)
        print(p)
        if(verbose == TRUE){
          LuckyVerbose(paste0(i,".GO_GSEAplot: ",id," completed!"),levels = 3)
        }
      }
      dev.off()
    }
  }
  save(go.gsea,file = paste0(dir,names,"_go.gsea.rda"))

  ###===========KEGG Gene Set Enrichment Analysis===============###
  if(verbose == T) LuckyVerbose("Step3: KEGG Gene Set Enrichment Analysis...")
  kk2 <- gseKEGG(geneList     = geneList,
                 organism     = organism,
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = pvalueCutoff$gseKEGG,
                 verbose      = FALSE)
  #head(kk2)
  d_kk2 <- as.data.frame(kk2)
  write.csv(d_kk2,paste0(dir,names,"_KEGG_Gene Set Enrichment Analysis .csv"))
  if(nrow(d_kk2) >= 1){
    # gseaplot
    pdf(paste0(dir,names,"_KEGG_gseaplot_Gene Set Enrichment Analysis.pdf"),10,8)
    for(i in 1:length(d_kk2$ID)){ # i=1
      ## gsea view
      id <- d_kk2$ID[i]
      title <- d_kk2$Description[i]
      p <- gseaplot(kk2,
                    geneSetID=id,
                    title = title)
      print(p)
      if(verbose == T) LuckyVerbose(paste0(i,".KEGG_GSEAplot: ",id," completed!"),levels = 3)
    }
    dev.off()

    # pathview
    for(i in 1:length(d_kk2$ID)){
      id <- d_kk2$ID[i]
      # pathview
      Plus.library("pathview")
      x <- pathview(gene.data = geneList,
                    pathway.id = id,
                    species = organism,
                    limit = list(gene = max(abs(geneList)),
                                 cpd = 1),
                    kegg.dir = pathview.dir,
                    kegg.native = T)
      if(verbose==T) LuckyVerbose(paste0(i,".KEGG_pathview: ",id," completed!"),levels = 2)
    }
  }

  ###========KEGG Module Gene Set Enrichment Analysis======###
  LuckyVerbose("Step3: KEGG Module Gene Set Enrichment Analysis...")
  mkk2 <- gseMKEGG(geneList = geneList,
                   organism = organism,
                   nPerm = 1000,
                   minGSSize = 10,
                   pvalueCutoff = pvalueCutoff$gseMKEGG,
                   verbose = F)
  d_mkk2 <- as.data.frame(mkk2)
  write.csv(d_mkk2,paste0(dir,names,"_KEGG_Module Gene Set Enrichment Analysis.csv"))
  if(nrow(d_mkk2) >= 1){
    # gseaplot
    pdf(paste0(dir,names,"_KEGG_Module Gene Set Enrichment Analysis.pdf"),10,8)
    for(i in 1:length(d_mkk2$ID)){ # i=1
      ## gsea view
      id <- d_mkk2$ID[i]
      title <- d_mkk2$Description[i]
      p <- gseaplot(mkk2,
                    geneSetID=id,
                    title = title)
      print(p)
      if(verbose == T) LuckyVerbose(paste0(i,".KEGG_GSEAplot: ",id," completed!"),levels = 3)
    }
    dev.off()
  }

  ###=====================Output data========================###
  l <- list(
    Repeat = list(
      geneList = geneList,
      OrgDb = OrgDb,
      keyType = keyType,
      organism = organism,
      pvalueCutoff = pvalueCutoff,
      verbose = verbose,
      save.path = save.path,
      names = names
    ),
    Data = list(
      gseGO = go.gsea,
      gseKEGG = kk2,
      gseMKEGG = gseMKEGG
    )
  )
  LuckyVerbose("All done!")
  return(l)
}




