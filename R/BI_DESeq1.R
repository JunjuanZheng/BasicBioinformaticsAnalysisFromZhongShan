
#' @title Standard process for DESeq2 RNA-Seq differencial expression analysis
#' @description DESeq1 is a customed way to use DESeq2 package for differencially expression analysis based on raw counts and design object,and create a DESeqList with dds,difSig and siggenes object.
#' @param counts raw counts matrix
#' @param design design object with sample rows and clinic feature cols
#' @param contrast.col the colnames of contrast.Like "N.status"
#' @param count.filter the rowmean filter of counts
#' @param cutoff.lFC the cutoff of log foldchange
#' @param cutoff.padj the cutoff of adjusted P value
#' @param save.file whether to save .rda object for DESeqList
#' @param names part name of saved files
#' @param report whether to do quality control report.Default is F.
#' @importFrom grDevices colorRampPalette
#' @return DESeqList
#' @seealso \code{\link["DESeq2"]{DESeqDataSetFromMatrix}};\code{\link["DESeq2"]{DESeq}};
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' ## data preparation
#' data("rna.counts")
#' data("rna.design")
#' counts = rna.counts;rm(rna.counts)
#' design = rna.design;rm(rna.design)
#'
#' ## DESeq2 package standard pipeline
#' dl1 <- DESeq1(counts,
#'               design,
#'               contrast.col= "condition",
#'               contrast.level =  c("normal","tumor"),
#'               contrast.control = "normal",
#'               count.filter=10,
#'               cutoff.lFC = 2,
#'               cutoff.padj = 0.05,
#'               save.file = F,
#'               names = "love",
#'               report = T)
#' View(dl1)
#' @export
DESeq1 <- function(counts,
                   design,
                   contrast.col,
                   contrast.level = c("Control","Treat"),
                   contrast.control = "Control",
                   count.filter=10,
                   cutoff.lFC = 1,
                   cutoff.padj = 0.1,
                   report = F,
                   save.file = T,
                   names = "love"){

  # package
  need  <-  c("DESeq2","BiocParallel","doParallel","parallel","pasilla","vsn","ggplot2");Plus.library(need)

  # adjusted counts
  counts <- counts[,rownames(design)]

  # contrast.level
  contrast.level <- c(contrast.control,setdiff(contrast.level,contrast.control))

  # rename design
  colnames(design)[match(contrast.col,colnames(design))] <- "contrast"

  ## DESeq pipeline
  f <- paste0("~ contrast");f <- as.formula(f)
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = f);
  dds <- dds[rowMeans(DESeq2::counts(dds)) >= count.filter,]
  dds$contrast <- factor(dds$contrast,levels = contrast.level)

  LuckyVerbose("DESeq pipeline is processing,please be patient...")
  ncore <- detectCores(logical = F)
  register(SnowParam(ncore))
  dds <- DESeq(dds,parallel=T)
  LuckyVerbose("DESeq pipeline had been completed!")


  ## DEA
  LuckyVerbose("Defferential Expression Analysis is processing...")
  res  <- results(dds)
  res$padj[is.na(res$padj)] <- 1
  res <- res[order(res$pvalue),] #head(res)

  ## other DEA
  # rn <- resultsNames(dds)[2]
  # resLFC <- lfcShrink(dds, coef=rn, type="apeglm")
  # resNorm <- lfcShrink(dds, coef=rn, type="normal")
  # resAsh <- lfcShrink(dds, coef=rn, type="ashr")

  ## get target genes
  list.res <- list(res=res)#,
                   #resLFC=resLFC,
                   #resNorm=resNorm,
                   #resAsh=resAsh)
  res.dds <- NULL
  for(i in 1:length(list.res)){
    res.i <- list.res[[i]]
    diffSig <- subset(res.i,abs(log2FoldChange) > cutoff.lFC & padj < cutoff.padj)
    diffSig <- as.data.frame(diffSig)
    siggenes <- rownames(diffSig)
    symbols <- convert(siggenes)
    gn <- convert(siggenes,totype = "GENENAME")
    gt <- convert(siggenes,totype = "gene_type")
    diffSig <- cbind(SYMBOL = symbols,
                     GENENAME = gn,
                     GENETYPE = gt,
                     diffSig)
    res.dds.i <- list(
      siggenes = siggenes,
      diffSig = diffSig,
      res = res.i
    )
    res.dds <- c(res.dds,list(res.dds.i))
    names(res.dds)[i] <- names(list.res)[i]
  }
  LuckyVerbose("Defferential Expression Analysis had been doned!")

  ###============================Report=======================###
  ntd <- normTransform(dds)
  gc_error <- base::tryCatch(vsd <- vst(dds, blind=T), error = function(e)e, finally = NULL)
  gc_e2 <- grep("varianceStabilizingTransformation",gc_error[[1]])
  lg_e2 <- length(gc_e2) != 1
  if(lg_e2){
    vsd <- vst(dds, blind=T) #if you wish to transform the data for downstream analysis,blind = F
  } else {
    vsd <- varianceStabilizingTransformation(dds,blind = T)
  }

  ## MA plot
  LuckyVerbose("Report1: Alternative shrinkage estimators...",levels = 2)
  {
     pdf(paste0(names,"_DESeq2 QC_MAplot.pdf"),8,8)
     xlim <- c(1,1e5); ylim <- c(-5,5)
     DESeq2::plotMA(res, xlim=xlim, ylim=ylim, main="MA:common")
    # plotMA(resLFC, xlim=xlim, ylim=ylim, main="MA:apeglm")
    # plotMA(resNorm, xlim=xlim, ylim=ylim, main="MA:normal")
    # plotMA(resAsh, xlim=xlim, ylim=ylim, main="MA:ashr")
     dev.off()
     if(report == T){
     win.graph(12,12)
     #par(mfrow=c(2,2), mar=c(4,4,2,1))
     xlim <- c(1,1e5); ylim <- c(-5,5)
     DESeq2::plotMA(res, xlim=xlim, ylim=ylim, main="MA:common")
    # plotMA(resLFC, xlim=xlim, ylim=ylim, main="MA:apeglm")
    # plotMA(resNorm, xlim=xlim, ylim=ylim, main="MA:normal")
    # plotMA(resAsh, xlim=xlim, ylim=ylim, main="MA:ashr")
     par(mfrow=c(1,1))
  }
  # print("Note: the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018) had been recommended by the DESeq2 authors!")
     }

  ## Plot row standard deviations versus row means
  LuckyVerbose("Report2: Plot row standard deviations versus row means...",levels = 2)
  {
      th.msd <- theme_bw() + theme(
       panel.grid =element_blank(),
       plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
       axis.title = element_text(face = "bold",size = 15),
       axis.text = element_text(face = "bold",size = 12),
       legend.title = element_text(face = "bold",size = 12),
       legend.text =element_text(face = "bold",size = 12),
       legend.position = "right"
      )

     msd1 <- meanSdPlot(assay(ntd),plot = F)
     msd2 <- meanSdPlot(assay(vsd),plot = F)
     gg.msd1 <- msd1$gg + scale_fill_gradient(low = mycolor[1], high = mycolor[4]) + th.msd + ggtitle("Normalized counts transformation")
     gg.msd2 <- msd2$gg + scale_fill_gradient(low = mycolor[1], high = mycolor[4]) + th.msd + ggtitle("Variance stabilizing transformation")
     pdf(paste0(names,"_DESeq2 QC_meanSdPlot.pdf"),8,8)
     print(gg.msd1);print(gg.msd2)
     dev.off()

     if(report == T){
       win.graph(15,14);print(gg.msd1)
       win.graph(15,14);print(gg.msd2)}
  }

  ## sample clustering and visualization
  LuckyVerbose("Report3:sample clustering and visualization..",levels = 2)
  {
     df1 <- as.data.frame(colData(dds))
     df <- df1[order(df1$contrast),]
     cl.name <- unique(as.character(df1$contrast))
     cl <- list(contrast = mycolor[1:length(cl.name)])
     names(cl$contrast) <- cl.name

     hp1 <- NULL
     pdf(paste0(names,"_DESeq2 QC_Heatmap.sample clustering.pdf"),8,8)
     for(i in 1:length(res.dds)){ # i =1
     select <- res.dds[[i]]$siggenes
     p.name <- names(res.dds)[i]
     he <- assay(vsd)[select,]
     he <- he[,rownames(df)]

     ## heatmap
     hp <- heatmap.dds2(
       dds=NULL,
       expr.matrix = he,log.convert = F,
       select = select,
       design = df,
       contrast.col = "contrast",
       contrast.level = NULL,
       contrast.list = cl,
       rowscale=T,
       expr.name = p.name,
       cluster_rows = T,cluster_cols=T,
       row.k = 4,
       show_column_names = F,
       show_row_names = F)
     print(hp)
     hp1 <- c(hp1,list(hp))
     }
     dev.off()
   if(report == T){
     for(i in 1:length(res.dds)){
       win.graph(15,14);print(hp1[[i]])
     }}
}

  ## Heatmap of the sample-to-sample distances
  LuckyVerbose("Report4: Heatmap of the sample-to-sample distances...",levels = 2)
  {
      sampleDists <- stats::dist(t(assay(vsd)))
      sampleDistMatrix <- as.matrix(sampleDists)
      sample.id <- colnames(sampleDistMatrix)
      df2 <- subset(design,select="contrast")
      rownames(sampleDistMatrix) <- paste(vsd$contrast, sep="-")
      colnames(sampleDistMatrix) <- sample.id
      colors <- colorRampPalette(rev(brewer.pal(9, "BuGn")))(255)
      ha <- HeatmapAnnotation(df2, col = cl)
      pdf(paste0(names,"_DESeq2 QC_Sample-to-sample Distances.pdf"),8,8)
      ## heatmap
      hp <- Heatmap(
        matrix = sampleDistMatrix,
        col = colors,
        name = "StS",
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        row_dend_side = c("left", "right"),
        cluster_columns = TRUE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        show_column_names = F,
        show_row_names = F,
        top_annotation = ha
      )
      print(hp)
     dev.off()
     if(report == T){win.graph(15,12);print(hp)}
  }

  ## Principal component plot of the samples
  LuckyVerbose("Report5: Principal component plot of the samples...",levels = 2)
  {
      pcaData <- DESeq2::plotPCA(vsd,
                                intgroup=c("contrast"),
                                 returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      pca <- ggplot(pcaData, aes(PC1, PC2 , color = contrast,shape = contrast)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        ggtitle("Principal component plot") +
       th.msd
     pdf(paste0(names,"_DESeq2 QC_Principal component plot.pdf"),8,8);print(pca);dev.off()
     if(report == T){win.graph(10,10);print(pca)}
  }


  ###========================Data Output====================###

  DESeqList <- list(
    siggenes = res.dds$res$siggenes,
    diffSig = res.dds$res$diffSig,
    dds = dds,
    result = res.dds
  )

  ## save dds object
  if(save.file == T){
    LuckyVerbose("DESeqList would be saved...")
    save(DESeqList,file = paste0(names,"_DESeqList.rda"))
    LuckyVerbose("DESeqList had been saved!")
  }

  return(DESeqList)
}

