


#' @export
limma1 <- function(counts,
                   design,
                   contrast.col,
                   method = c("voom","trend")[1],
                   contrast.level = c("Np","N0"),
                   cutoff.lFC = 1,
                   cutoff.padj = 0.1,
                   save.file = T,
                   names = "love"){

  ## 加载必要的包
  need <- c("limma","edgeR")
  Plus.library(need)

  ## 对齐
  counts <- counts[,rownames(design)]

  ## 生成contrast
  contrasts <- combn(contrast.level,2)
  cont1 <- NULL;
  for(i in 1:ncol(contrasts)){
    e <- as.character(contrasts[,i])
    c.i1 <- paste0(e,collapse = "-")
    #c.i2 <- paste0(e[c(2,1)],collapse = "-")
    cont1 <- c(cont1,c.i1)
  }
  contrast.Matrix <- makeContrasts(contrasts=cont1,levels=contrast.level)

  ## design matrix.
  group <- factor(design[,contrast.col],levels = contrast.level)
  design.mt <- model.matrix(~ 0 + group);
  colnames(design.mt) <- contrast.level

  ## edgeR
  genelist <- DGEList(counts=counts,group = group)
  keep <- rowSums(cpm(genelist) > 0.5) >= 2
  genelist.filted <- genelist[keep,keep.lib.sizes=FALSE]
  print("limma:normalization is processing...")
  genelist.norm <- calcNormFactors(genelist.filted, method = "TMM")
  print("limma:normalization done!")

  ## 差异性分析
  if(method == "trend"){
    print("stratgy:limma-trend...")
    ## log expression
    logCPM <- cpm(genelist.norm, log=TRUE, prior.count=1)

    ## 差异性分析
    fit <- lmFit(logCPM,design.mt)
    fit <- contrasts.fit(fit, contrast.Matrix)
    fit <- eBayes(fit,trend = T)
    results <- decideTests(fit,p.value = cutoff.padj,lfc = cutoff.lFC);summary(results) #查看up和down的信息
    dif.expr <- topTable(fit, coef=paste(contrast.level,collapse = "-"), n=Inf) #输出
    sig.expr <- subset(dif.expr,abs(logFC) >= cutoff.lFC & adj.P.Val < cutoff.padj)
    siggenes <- rownames(sig.expr)

    ## 输出结果
    limmaList <- list(
      exprs = logCPM,
      diffSig = sig.expr,
      siggenes = siggenes
    )

  } else {
    ## 使用voom策略
    print("stratgy:voom...")
    Elist <- voom(genelist.norm, design.mt, plot=TRUE)

    ## 差异性分析
    fit <- lmFit(Elist,design.mt)
    fit <- contrasts.fit(fit, contrast.Matrix)
    fit <- eBayes(fit,trend = T)
    results <- decideTests(fit,p.value = cutoff.padj,lfc = cutoff.lFC);summary(results) #查看up和down的信息
    dif.expr <- topTable(fit, coef=paste(contrast.level,collapse = "-"), n=Inf) #输出
    sig.expr <- subset(dif.expr,abs(logFC) >= cutoff.lFC & adj.P.Val < cutoff.padj)
    siggenes <- rownames(sig.expr)

    ## 输出结果
    limmaList <- list(
      exprs = Elist,
      diffSig = sig.expr,
      siggenes = siggenes
    )
  }

  ## 结果输出
  if(save.file == T){
    save(limmaList,file = paste0(names,"_limmaList.rda"))
    print("limmaList has been saved in present work space!")
  }
  print("All done!")
  return(limmaList)
}









