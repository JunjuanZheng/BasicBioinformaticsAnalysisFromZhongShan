

#' @title GO analysis via clusterProfiler package
#' @description FastGOplot give fast way to draw GO series plot like barplot/dotplot/emapplot,Cnetplot and goplot for a enrichGO object based on clusterProfiler package.
#' @param genes a character of gene id.
#' @param geneList a numeric with id names.For example,logFoldChang with ENSEMBL id names.Note that the names of \code{geneList} must be the same type of \code{genes}
#' @param default.universe whether to use default universe.If \code{default.universe = F},the background genes would be provide by \code{geneList}
#' @param classlevel Default is 2.I had tested that 2:7 was still available.If many levels is set,the process is time-consuming
#' @param OrgDb OrgDb dataset.If NULL,use "org.Hs.eg.db"
#' @param keyType the type of gene.Support automatically test
#' @param pAdjustMethod the mathod of pvalue adjustment
#' @param pvalueCutoff  a list of pvalue adjust method for \code{\link[clusterProfiler]{enrichGO}}
#' @param qvalueCutoff  cutoff of q value
#' @param cnet.showCategory the number of showed cluster at cnetplot
#' @param verbose LuckyVerbose gseplot running message or not
#' @param save.path the sub path of saved files
#' @param names the main path and part names of saved files
#' @importFrom clusterProfiler groupGO enrichGO
#' @importFrom enrichplot dotplot emapplot goplot cnetplot
#' @importFrom ggplot2 ggtitle
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' need <- c("clusterProfiler","org.Hs.eg.db");Plus.library(need)
#' data(geneList, package = "DOSE")
#' data(geneList, package='DOSE')
#' universe <- unique(as.character(common.annot$ENTREZID))
#' system.time(
#' l <- FastGO(genes,
#'             geneList,
#'             classlevel = 2:2,
#'             OrgDb  = NULL,
#'             keyType = NULL,
#'             pAdjustMethod = "BH",
#'             pvalueCutoff = list(
#'                gseGO = 0.05,
#'                enrichGO = 0.05),
#'             qvalueCutoff  = 0.05,
#'             cnet.showCategory = 5,
#'             save.path = "GO",
#'             names = "love")
#' )
#' ## enhanced plot
#' g <- enrichplot::cnetplot(goList$GOEnrichment$BP,
#'                           showCategory=10,
#'                           colorEdge=T,
#'                           node_label=T,
#'                           circular = F)
#' g1 <- g + theme(legend.text = element_text(face = "bold",size = 12))
#'
#' @export
FastGO <- function(genes,
                   geneList,
                   default.universe = F,
                   classlevel = 2:2,
                   OrgDb  = NULL,
                   keyType = NULL,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   cnet.showCategory = 5,
                   verbose = TRUE,
                   save.path = "GO",
                   names = "love"){
  ## Package
  nd <- c("clusterProfiler","org.Hs.eg.db","stringr","DOSE","ggplot2");
  Plus.library(nd)
  
  ### 差异性基因的多个Categery的enrichGO analysis
  Go.categery <- c("MF","BP","CC")
  
  ### reference database
  if(is.null(OrgDb)){
    OrgDb <- "org.Hs.eg.db"
    LuckyVerbose("You don't specify the 'OrgDb' parameter.Here we use 'org.Hs.eg.db'(homo spacies)...")
  }
  
  ### whether is ensembl id
  if(is.null(keyType)){
    lg1 <- length(grep("ENSG",names(geneList))) == 0
    if(lg1){
      #说明不是ensemblid
      LuckyVerbose("The names of geneList is ENTREZ ID...")
      keyType <-  "ENTREZID"
    } else {
      LuckyVerbose("The names of geneList is ENSEMBL ID...")
      keyType <-  "ENSEMBL"
    }
  }
  
  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)
  
  ## geneList 排序
  LuckyVerbose("Step1: order 'geneList' from upper to lower...")
  geneList <- sort(geneList,decreasing = T)
  
  ###===================GO classification=================###
  {
    LuckyVerbose("Step2: GO classification...")
    go.classification <- NULL
    #for(ont in Go.categery){
    #  a <- NULL
    #  for(level in classlevel){
    #    ggo <- groupGO(gene     = genes,
    #                   keyType  = keyType,
    #                   OrgDb    = OrgDb,
    #                   ont      = ont,
    #                   level    = 7,
    #                   readable = TRUE);
    #    go.classification <- c(go.classification,list(ggo))
    #    names(go.classification)[grep(ont,Go.categery)] <- ont
    #    a.i <- as.data.frame(ggo)
    #    if(verbose == TRUE){
    #      LuckyVerbose(paste0(ont,": GO classification_level ",level," was collected!"),levels = 2)
    #    }
    #    a.i <- a.i[a.i$Count != 0,]
    #    a <- rbind(a,a.i)
    #  }
    #   a <- a[!duplicated(a$ID),]
    #  a$categery <- ont
    #
    #}
    #save(go.classification,file = paste0(dir,names,"_go.classification.rda"))
  }
  
  ###==================GO over-representation test========###
  {
    LuckyVerbose("Step3: GO over-representation test...")
    go.enrich <- NULL #记录egoList的list。
    go.enrich.table <- NULL #记录ego表格的list。
    for(i in 1:length(Go.categery)) {  # i=1
      ont <- Go.categery[i]
      if(verbose == TRUE){
        LuckyVerbose(paste0(ont,": enrichGO List is on establishment..."),levels = 2)
      }
      if(default.universe==T){
        LuckyVerbose("You try to use a default background provided by ClusterProfilter.In factor,a self-defined background is more proper.")
        ego.sig.tv <- enrichGO(
          gene = genes,
          OrgDb  = OrgDb,
          ont   = ont,
          keyType = keyType,
          pAdjustMethod = pAdjustMethod,
          pvalueCutoff  = pvalueCutoff,
          qvalueCutoff  = qvalueCutoff,
          readable = TRUE)
      }else{
        ego.sig.tv <- enrichGO(
          gene = genes,
          universe  = names(geneList),
          OrgDb  = OrgDb,
          ont   = ont,
          keyType = keyType,
          pAdjustMethod = pAdjustMethod,
          pvalueCutoff  = pvalueCutoff,
          qvalueCutoff  = qvalueCutoff,
          readable = TRUE)
      }
      go.enrich <- c(go.enrich,list(ego.sig.tv))
      names(go.enrich)[i] <- ont
      if(verbose == T){
        LuckyVerbose(paste0(ont,": enrichGO List completed!"),levels = 2)
      }
      ego.i <- as.data.frame(ego.sig.tv)
      go.enrich.table <- c(go.enrich.table,list(ego.i))
      ## save result
      ego.sig.tv.savepath <- str_c(dir,names,"_",ont,"_enrichGO analysis.csv");
      write.csv(ego.i,ego.sig.tv.savepath)
      
      if(nrow(ego.i) == 0){
        if(verbose == T){
          LuckyVerbose(paste0(i,": No GO enrichment"),levels = 2)
        }
      } else {
        sC <- ifelse(nrow(ego.i) < 10,nrow(ego.i),10)
        
        #enrichGo plot
        plot.names <- str_c(dir,names,"_",ont,"_enrichGO plot.pdf")
        pdf(plot.names,width = 10,height = 10)
        print(barplot(go.enrich[[i]], showCategory=sC))
        print(dotplot(go.enrich[[i]], showCategory=sC))
        if(nrow(ego.i) >= 5) print(emapplot(go.enrich[[i]],showCategory=sC))
        dev.off()
        if(verbose == T){
          LuckyVerbose(str_c(ont,": Get enrichGO plot!"),levels = 2)
        }
      }
    }
    save(go.enrich,file = paste0(dir,names,"_go.enrich.rda"))
    print("储存完毕")
    
    ##Cnetplot
    # if(nrow(ego.i) >= 5){
    #  pdf.pathway <- str_c(dir,names,"_Get Gene-Concept Network plot.pdf");
    #  p1 <- cnetplot(go.enrich[[1]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList)) + ggtitle(Go.categery[1])
    #  p2 <- cnetplot(go.enrich[[2]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList))+ ggtitle(Go.categery[2])
    #  p3 <- cnetplot(go.enrich[[3]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList))+ ggtitle(Go.categery[3])
    #  pdf(pdf.pathway,width = 20,height = 20)
    #  print(p1);print(p2);print(p3)
    #  dev.off()
    #  LuckyVerbose(str_c("Get Gene-Concept Network!"))
    # }
    if(nrow(go.enrich.table[[1]]) >= 5){
      p1 <- cnetplot(go.enrich[[1]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList)) + ggtitle(Go.categery[1])
    }else{
      p1 <- NULL
    }
    if(nrow(go.enrich.table[[2]]) >= 5){
      p2 <- cnetplot(go.enrich[[2]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList))+ ggtitle(Go.categery[2])
    } else{
      p2 <- NULL
    }
    if(nrow(go.enrich.table[[3]]) >= 5){
      p3 <- cnetplot(go.enrich[[3]],showCategory = cnet.showCategory,categorySize="pvalue",foldChange = as.numeric(geneList))+ ggtitle(Go.categery[3])
    }else{
      p3 <- NULL
    }
    pdf.pathway <- str_c(dir,names,"_Get Gene-Concept Network plot.pdf");
    pdf(pdf.pathway,width = 20,height = 20)
    print(p1);print(p2);print(p3)
    dev.off()
    LuckyVerbose(str_c("Get Gene-Concept Network!"))
    
    
    ##goplot
    pdf.pathway <- str_c(dir,names,"_Goplot.pdf");
    pdf(pdf.pathway,width = 12,height = 8)
    for(i in 1:3){
      go.enrich.i <- go.enrich[[i]]
      ego.i <- as.data.frame(go.enrich.i)
      if(nrow(ego.i) > 5){
        g1 <- goplot(go.enrich.i,main=Go.categery[i])
        print(g1)
      }
    }
    dev.off()
    LuckyVerbose(paste0("GOplot have been saved!"))
  }
  
  ###==================Output the result===================###
  l <- list(
    Repeat = list(
      genes = genes,
      geneList = geneList,
      classlevel = classlevel,
      OrgDb  = OrgDb,
      keyType = keyType,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff =  pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      cnet.showCategory = cnet.showCategory,
      verbose = verbose,
      save.path = save.path,
      names = names
    ),
    GOClassification = go.classification,
    GOEnrichment = go.enrich
  )
  LuckyVerbose("All done!")
  return(l)
}





