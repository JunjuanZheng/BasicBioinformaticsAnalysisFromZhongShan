

## edgeR1 help do rawcounts differencially expreesion analysis via edgeR package.
# counts# raw counts
# design# design object
# contrast.col# contrast colname in the design object
# contrast.level# the order should be treat then control
# cutoff.lFC # the cut-off of log fold change
# cutoff.padj # the cut-off of FDR
# save.file # whether to save edgeRList
# names # part name of saved file
#' @export
edgeR1 <- function(counts,
                   design,
                   contrast.col,
                   contrast.level = c("Treat","Control"),
                   contrast.control = "Control",
                   cutoff.lFC = 1,
                   cutoff.padj = 0.1,
                   show.report = F,
                   save.file = T,
                   names = "love"){
  ## 加载必要的包
  nd <- c("limma","edgeR")
  Plus.library(nd)

  ## 矩阵对齐
  counts <- counts[,rownames(design)]

  ## contrast.level
  contrast.level <- c(setdiff(contrast.level,contrast.control),contrast.control)

  ## condition向量
  design.x <- design
  group <- factor(as.character(design[,contrast.col]),levels = contrast.level)

  ## design矩阵
  design <- model.matrix(~ 0 + group)
  colnames(design) <- contrast.level

  ## contrast.matrix
  contrasts <- combn(contrast.level,2)
  cont1 <- NULL;
  for(i in 1:ncol(contrasts)){
    e <- as.character(contrasts[,i])
    c.i1 <- paste0(e,collapse = "-")
    #c.i2 <- paste0(e[c(2,1)],collapse = "-")
    cont1 <- c(cont1,c.i1)
  }
  contrast.Matrix <- makeContrasts(contrasts=cont1,levels=contrast.level)


  ###===============Pre-Differential expression=================###
  y <- DGEList(counts, group=group)
  keep <- rowSums(cpm(y) > 0.5) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("calcNormFactors:Calculating normalized factors,please wait...")
  y <- calcNormFactors(y)
  #图1
  if(show.report == T){
    win.graph(width = 12,height = 12)
    par(mfrow = c(2,2))
    plotMD(cpm(y, log=TRUE), column=1,
           main = "plotMD:After calcNormFactors")
    abline(h=0, col="red", lty=2, lwd=2)
  }
  print("estimateDis:Estimate negative binomial dispersions by weighted likelihood empirical Bayes,please wait...")
  y <- estimateDisp(y, design, robust=TRUE)
  if(show.report == T){
    plotBCV(y,main = "plotBCV:After estimateDisp") #图2
  }
  print("glmQLFit:Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests...")

  ###=====================Get matrix for Clustering===============###
  logcpm <- cpm(y,prior.count = 2 ,log=TRUE)

  ###=======================Differential expression===================###
  fit <- glmQLFit(y, design, robust=TRUE)
  if(show.report == T){
    plotQLDisp(fit,main = "plotQLDisp:After glmQLFit")#图3
  }

  qlf <- glmQLFTest(fit, contrast=contrast.Matrix)
  if(show.report == T){
    plotMD(qlf) #图4
  }
  summary(decideTests(qlf, p.value = cutoff.padj,lfc = cutoff.lFC))
  dif <- topTags(qlf,n = Inf) # colnames(dif)
  dif <- as.data.frame(dif)
  sigdif <- subset(dif,abs(logFC) > cutoff.lFC & FDR < cutoff.padj)
  print(all(rownames(sigdif) %in% rownames(counts)))
  siggenes <- rownames(sigdif)
  symbols <- convert(siggenes)
  gt <- convert(siggenes,totype = "gene_type")
  gn <- convert(siggenes,totype = "GENENAME")
  diffSig <- cbind(SYMBOL = symbols,
                   GENETYPE = gt,
                   GENENAME = gn,
                   sigdif)
  par(mfrow = c(1,1))

  ## output data
  edgeRList <- list(
    siggenes  = siggenes,
    diffSig = diffSig,
    edgeR.fit = qlf,
    design = design.x,
    logcpm = logcpm
  )
  if(save.file == T){
    print("edgeRList would be saved...")
    save(edgeRList,file = paste0(names,"_edgeRList.rda"))
    print("edgeRList had been saved!")
  }
  return(edgeRList)

}








