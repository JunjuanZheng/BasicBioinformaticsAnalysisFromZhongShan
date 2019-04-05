

#' @keywords testBatchEffect
#' @title  test and remove potencial batch effects
#' @description testBatchEffect help test potencial and specified fators that may bring batch effect base on DESeq2 pipeline
#' @param SigObject the result of \code{\link{DESeq1}}, \code{\link{edgeR1}} or a expression matrix
#' @param design design object.If \code{SigObject} is a matrix,it is needed.
#' @param candidate the candidate variables you selectd like years or days
#' @inheritParams sva::num.sv
#' @param names the name of saved file
#' @param report whether to show plot reports
#' @importFrom sva num.sv svaseq f.pvalue ComBat
#' @importFrom SummarizedExperiment colData assay
#' @importFrom DESeq2 vst plotPCA
#' @importFrom limma removeBatchEffect
#' @import ggplot2
#' @seealso \code{\link{DESeq1}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' coldata <- as.data.frame(colData(SigObject$dds))
#' colnames(coldata)
#' candidate=c("gender","diognosis.year","pT")
#' t.be <- testBatchEffect(SigObject,
#'                         candidate,
#'                         names = "test")
#' View(t.be)
#' View(t.be$Unknown$SigTest)
#' @export
testBatchEffect <- function(SigObject,
                            design=NULL,
                            candidate=c("gender","diognosis.year"),
                            n.sv = NULL,
                            names = "love",
                            report =T){
  ## 加载必要的包
  nd <- c("DESeq2","sva","ggplot2","vsn","plyr")
  #,"doParallel","parallel","pasilla","BiocParallel"
  Plus.library(nd)

  if("dds" %in% names(SigObject)){

    ## 选择数据
    LuckyVerbose("Test batch effect from DESeqList object...")
    dds <- SigObject$dds
    dat <- DESeq2::counts(dds, normalized=TRUE)
    idx <- rowMeans(dat) > 10
    dat <- dat[idx,]

    ## 选择潜在的批次效应因素
    coldata <- as.data.frame(colData(dds))

    ###=============Unknown batch effect exploration==========###
    LuckyVerbose("Unknown batch effect exploration...")
    {
      ## sva移除模型
      mod <- model.matrix(~ contrast,coldata)
      mod0 <- model.matrix(~ 1,coldata)

      ##calculating the variables
      if(is.null(n.sv)){
        n.sv <- num.sv(dat = dat,
                       mod = mod,
                       method="leek") #gives 16

        ### find a appropriate  surogate variable number
        for(i in n.sv:1){ # i=16
          a <- tryCatch({lnj.corr <- svaseq(dat,mod,mod0,n.sv = i)},
                        error = function(e)e,
                        warning = function(w)w,
                        finally = NULL)
          if("message" %in% names(a)){
            LuckyVerbose(paste0("The number of surogate variables is not appropriate so another one would be tested..."),levels = 2)
            next
          } else {
            break #跳出循环
          }
        }
      } else {
        lnj.corr <- svaseq(dat,mod,mod0,n.sv = n.sv)
      }

      ## 结果
      plot(lnj.corr$sv,pch=19,col="blue")
      cmodSv = cbind(mod, lnj.corr$sv)
      cmod0Sv = cbind(mod0, lnj.corr$sv)
      cpValuesSv = f.pvalue(dat, cmodSv, cmod0Sv)
      cqValuesSv = p.adjust(cpValuesSv, method = "BH")
      df <- cbind(P.Val = cpValuesSv,FDR = cqValuesSv)
      df <- df[order(df[,"FDR"],decreasing = F),]
      write.csv(df,"DSEeq2_Unknown batch effect exploration.csv")
      LuckyVerbose("The data of Unknown batch effect exploration had been saved...",levels = 2)
    }

    ###===========Candidate batch effect exploration===========###
    LuckyVerbose("Candidate batch effect exploration...")

    ## candidate data
    if(!is.null(candidate)){
      select <- c("contrast",candidate)
      colData(dds) <- subset(colData(dds),select = select)
    }
    coldata <- as.data.frame(colData(dds))

    ## 去除colData中的NA值
    lg1 <- !apply(coldata,1,is.one.na)
    coldata <- coldata[lg1,]
    if(F %in% names(table(lg1))){
      test <- table(lg1)[names(table(lg1)) %in% F]
      cat("There are",test,"samples are deleted for NA existing.","\n")
    }

    ## 生成矩阵
    vsd <- vst(dds,blind = T)
    for(i in 1:ncol(vsd@colData)){
      vsd@colData[,i] <- factor(vsd@colData[,i],levels = unique(as.character(vsd@colData[,i])))
    }
    vsd_1 <- vsd[,rownames(coldata)]
    vsdMat <- assay(vsd_1)

    ## removeBatchEffect
    vsdMat.c <- NULL
    for(i in 1:length(candidate)){ # i =2
      c.i <- candidate[i]
      vsdMat.ci <- removeBatchEffect(vsdMat,coldata[,c.i])
      vsdMat.c <- c(vsdMat.c,list(vsdMat.ci))
      names(vsdMat.c)[i] <- c.i
    }

    ## candidata data
    pca <- NULL
    pdf(paste0(names,"_DESeq2 QC_Batch efect_PCA plot.pdf"),10,7)
    for(i in 1:length(candidate)){ # i = 1
      ## batch数据提取
      batch <- candidate[[i]]

      ## pca图
      pcaData <- plotPCA(object = vsd_1,
                         intgroup = c(batch,"contrast"),
                         returnData=TRUE)
      b.p <- match(batch,colnames(pcaData))
      colnames(pcaData)[b.p] <- "batch"
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      pca.i <- ggplot(pcaData, aes(PC1, PC2 ,
                                   color = batch,
                                   shape = contrast)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        ggtitle(paste0("Principal component plot: ",batch)) +
        theme_bw() +
        theme(
          panel.grid =element_blank(),
          plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
          axis.title = element_text(face = "bold",size = 12),
          axis.text = element_text(face = "bold",size = 12),
          legend.title = element_text(face = "bold",size = 12),
          legend.text =element_text(face = "bold",size = 12),
          legend.position = "right"
        )
      pca <- c(pca,list(pca.i))
      names(pca)[i] <- batch
      print(pca.i)
    }
    dev.off()
    if(report== T){
      for(i in pca) {win.graph(15,12);print(i)}
    }

    ## output data
    l <- list(
      Unknown = list(
        SigTest = df,
        n.sv = ncol(lnj.corr$sv)
      ),
      Candidate = list(
        AfterCorrect = vsdMat.c,
        Plot.PCA = pca
      )
    )
    LuckyVerbose("All done!")
    return(l)

  }

  if("edgeR.fit" %in% names(SigObject)){
    ## 说明是edgeR1 result
    LuckyVerbose("Test batch effect from edgeRList object...")
    expr <- SigObject$logcpm # dim(expr)
    expr <- expr[,rownames(design)]
    expr2 <- expr
    if("design" %in% names(SigObject)) {design <- SigObject$design}
    for(i in 1:length(candidate)){ # i=2
      modcombat <-  model.matrix(~1, data = design)
      cand.i <- candidate[i]
      batch <-  as.character(design[,cand.i])
      p <- which(is.na(batch))
      if(length(p) != 0){
        batch <- batch[-p]
        expr2 <- expr2[,-p]
        modcombat <- modcombat[-p,]
      }
      expr2 <-  ComBat(dat=expr2,
                       batch=batch,
                       mod=modcombat)
      design <- design[colnames(expr2),]
    }
    #View(head(expr2))

    ## test NA sample
    na.conuts <- ncol(expr) - ncol(expr2)
    if(na.conuts > 0){
      LuckyVerbose("There are ",na.conuts," sample with NA value in specified batch variables and removed!")
    }

    ## 输出结果
    l <- list(
      Unknown = list(
        SigTest = NA,
        n.sv = NA
      ),
      Candidate = list(
        AfterCorrect = expr2,
        Plot.PCA = NA
      )
    )
    LuckyVerbose("All done!")
    return(l)
  }

  x <- !("dds" %in% names(SigObject)) & !("edgeR.fit" %in% names(SigObject))
  if(x){
    ## 说明是矩阵
    LuckyVerbose("Test batch effect from expression matrix...")
    expr <- SigObject
    expr <- expr[,rownames(design)]
    expr2 <- expr
    for(i in 1:length(candidate)){ # i=1
      cand.i <- candidate[i]
      modcombat <-  model.matrix(~1, data = design)
      batch <-  as.character(design[,cand.i])
      batch <- as.factor(batch)
      #防止空值
      p <- which(is.na(batch))
      if(length(p) != 0){
        batch <- batch[-p]
        expr2 <- expr2[,-p]
        modcombat <- as.matrix(modcombat[-p,])
        names(modcombat) <- "(interpret)"
      }
      expr2 <-  ComBat(dat=expr2,
                       batch=batch,
                       mod=modcombat)#mod=modcombat
      design <- design[colnames(expr2),]
    }
    #View(head(expr2))

    ## test NA sample
    na.conuts <- ncol(expr) - ncol(expr2)
    if(na.conuts > 0){
      LuckyVerbose("There are ",na.conuts," sample with NA value in specified batch variables and removed!")
    }

    ## 输出结果
    l <- list(
      Unknown = list(
        SigTest = NA,
        n.sv = NA
      ),
      Candidate = list(
        AfterCorrect = expr2,
        Plot.PCA = NA
      )
    )
    LuckyVerbose("All done!")
    return(l)
  }

}

