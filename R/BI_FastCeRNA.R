
#' @title easy ceRNA network via GDCRNATools::gdcCEAnalysis
#' @description easy ceRNA network via GDCRNATools::gdcCEAnalysis
#' @param gene.counts the raw counts matrix of mRNA and lncRNA
#' @param miRNA.counts the raw counts matrix of miRNA
#' @param matrix.type  the type of expression matrix.Default is "counts".If you have a normalized matrix,please set it "normalized".
#' @param deMIR differencial expression miRNAs
#' @param lnc.targets a list with lcnRNA ENSEMBL names and miRNAs symbol elements.Or "starBase"
#' @param pc.targets a list with protein-coding RNA(mRNA) ENSEMBL names and miRNAs symbol elements.
#' @inheritParams FastWGCNA
#' @importFrom limma voom
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom plyr adply
#' @seealso \code{\link[GDCRNATools]{gdcCEAnalysis}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
FastCeRNA <- function(gene.counts,
                      miRNA.counts,
                      matrix.type = "counts",
                      select = NULL,
                      design,
                      deMIR = NULL,
                      lnc.targets = "starBase",
                      pc.targets = "starBase",
                      save.path = "ceRNA",
                      names = "love"){
  ## 加载包
  nd <- c("limma","edgeR","plyr")
  Plus.library(nd)

  ## 保存路径
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)

  ###===============Get lncRNA and protein-coding=============###
  LuckyVerbose("Step1: Get lncRNA and mRNA from expression matrix...")
  ## 基因类型
  genes.type <- convert(rownames(gene.counts),
                        totype = "gene_type")

  ## 获得lncRNA的id
  lncRNA.type <- c("3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","lincRNA","macro_lncRNA","non_coding","processed_transcript","sense_intronic","sense_overlapping","TEC")
  lnc <- rownames(gene.counts)[genes.type %in% lncRNA.type]

  ## 获得protein coding的id
  pc <- rownames(gene.counts)[genes.type %in% "protein_coding"]
  if(!is.null(select)){
    lg1 <- !all(select %in% pc)
    if(lg1){
      LuckyVerbose("Some of selected pcRNA not exist!",levels = 2)
    } else {
      LuckyVerbose("Right select...",levels = 2)
    }
    pc <- intersect(select,pc)
  }

  ## 对齐
  #s <- intersect(colnames(gene.counts),colnames(miRNA.counts))
  s <- rownames(design)
  gene.counts_1 <- gene.counts[c(lnc,pc),s]
  miRNA.counts_1 <- miRNA.counts[,s]

  ###==================voom transformed matrix=================###
  if(matrix.type == "counts"){
    LuckyVerbose("Step2: Get voom normalized expression matrix...")
    # from limma1 function in lucky package
    counts.list <- list(pcLnc = gene.counts_1,
                        miRNA = miRNA.counts_1)
    voom.Expr <- NULL
    for(i in 1:length(counts.list)){ # i=2
      names.i <- names(counts.list)[i]
      LuckyVerbose("limma::voom:","deal with",names.i,"expression matrix...",levels = 2)
      counts <- counts.list[[i]]
      genelist <- DGEList(counts=counts)
      genelist.norm <- calcNormFactors(genelist, method = "TMM")
      Elist <- voom(counts = genelist.norm,
                    design = NULL,
                    plot = FALSE)
      voom.i <- Elist$E
      voom.Expr <- c(voom.Expr,list(voom.i))
      names(voom.Expr)[i] <- names.i
    }
    rna.expr = voom.Expr$pcLnc
    mir.expr = voom.Expr$miRNA
  } else {
    rna.expr = gene.counts_1
    mir.expr = miRNA.counts_1
  }

  ###==================ceRNA network building=================###
  LuckyVerbose("Step3: build ceRNA network.It's time-consuming,please wait...")
  ceOutput <- gdcCEAnalysis_1(lnc = lnc,
                              pc = pc,
                              deMIR = deMIR,
                              lnc.targets = lnc.targets,
                              pc.targets = pc.targets,
                              rna.expr = rna.expr,
                              mir.expr = mir.expr)

  ###=====================Output data=======================###
  save(ceOutput,file = paste0(dir,"/",names,"_ceOutput.rda"))
  LuckyVerbose("All done!")
  return(ceOutput)

}

####====================Assistant Functions====================####
gdcCEAnalysis_1 <- function(lnc, pc,
                            deMIR=NULL,
                            lnc.targets,
                            pc.targets,
                            rna.expr,
                            mir.expr){

  hyperOutput <- hyperTestFun_1(lnc = lnc,
                              pc = pc,
                              deMIR = deMIR,
                              lnc.targets=lnc.targets,
                              pc.targets=pc.targets,
                              mir.expr = mir.expr)
  LuckyVerbose('Hypergenometric test done !',levels = 2)

  regOutput <- multiRegTestFun_1(hyperOutput = hyperOutput,
                               rna.expr=rna.expr,
                               mir.expr=mir.expr,
                               lnc.targets=lnc.targets,
                               pc.targets=pc.targets)
  LuckyVerbose('Correlation analysis done !',levels = 2)
  LuckyVerbose('Regulation pattern analysis done!',levels = 2)

  ceOutput <- data.frame(hyperOutput, regOutput, row.names=NULL)

  return(ceOutput)
}

### hypergeometric test
#lncTargets和pcTargets的名称是数据库名，内部的数据是list(ENSG = miRNA...)的形式。可以自己构造。
hyperTestFun_1 <- function(lnc, pc, deMIR,
                           lnc.targets,
                           pc.targets,
                           mir.expr){
  old0 <- options()
  options(stringsAsFactors = F)
  ## 预测的数据格式
  if (!is.character(lnc.targets)) {
    x <- NULL
    for(i in 1:length(lnc.targets)){ #i=1
      x1 <- lnc.targets[[i]]
      x <- c(x,x1)
    }
    lnc.targets <- x
  } else {
    lnc.targets <- ceRNA_lncTargets[[lnc.targets]]
  }

  if (!is.character(pc.targets)) {
    x <- NULL
    for(i in 1:length(pc.targets)){ #i=1
      x1 <- pc.targets[[i]]
      x <- c(x,x1)
    }
    pc.targets <- x
  } else {
    pc.targets <- ceRNA_pcTargets[[pc.targets]]
  }
  mir1 <- unique(unlist(lnc.targets))
  mir2 <- unique(unlist(pc.targets))
  mirs <- union(mir1,mir2)
  popTotal <- length(mirs) #总miRNA背景

  ceLNC <- lnc[lnc %in% names(lnc.targets)]
  cePC <- pc[pc %in% names(pc.targets)]

  ## get data from every lncRNA
  hyper_lncID <- function(lncID){#lncID = ceLNC[1]
    listTotal <- length(unique(unlist(lnc.targets[[lncID]])))
    ## get data from every gene
    hyper_gene <- function(gene){ #gene = cePC[1]
      ## 某个pc与lncRNA的共同miRNA,而且此miRNA必须在mir.expr中出现
      ovlp <- intersect(unique(unlist(lnc.targets[[lncID]])),unique(unlist(pc.targets[[gene]])))
      ovlp1 <- intersect(ovlp,rownames(mir.expr)) #选择miR矩阵中存在的（原代码中没有的。这是为了后方防止regSim和sppc值的计算有误）
      if(length(ovlp)==0){
        #没有共同的miRNA.此时输出空值
        hyperOutput.i <- data.frame(
          lncRNAs = lncID,
          Genes = gene,
          Counts = 0,
          listTotal = listTotal,
          popHits = NA,
          popTotal = NA,
          foldEnrichment = NA,
          hyperPvalue = NA,
          realCounts = 0,
          miRNAs = NA,
          deMIRCounts = NA,
          deMIRs = NA
        )
      } else {
        #有共同的miRNA.
        popHits <- length(unique(unlist(pc.targets[[gene]]))) #背景：某pc的靶点miRNA数目
        Counts <- length(ovlp) #pc/lnc共同的靶点数

        ## 富集倍数和p值
        foldEnrichment <- Counts/listTotal*popTotal/popHits
        pValue <- phyper(Counts-1, popHits, popTotal-popHits,
                         listTotal, lower.tail=FALSE, log.p=FALSE)

        ## 取pc/lnc共同的、且在miR.expr中存在的miRNAs
        ovlpMIRs <- paste(ovlp1, collapse = ',')
        realCounts <- length(ovlp1)

        #取差异性基因与pc/lnc共同miR的交集
        ceMIR <- Reduce(intersect, list(ovlp1, deMIR))
        deMIRCounts <- length(ceMIR)#计算有多少个
        deMIRs <- paste(ceMIR, collapse = ',')

        hyperOutput.i <- data.frame(
          lncRNAs = lncID,
          Genes = gene,
          Counts = Counts,
          listTotal = listTotal,
          popHits = popHits,
          popTotal = popTotal,
          foldEnrichment = foldEnrichment,
          hyperPvalue = pValue,
          realCounts = realCounts,
          miRNAs = ovlpMIRs,
          deMIRCounts = deMIRCounts,
          deMIRs = deMIRs
        )
      }
      return(hyperOutput.i)
    }
    df2 <- adply(as.matrix(cePC),1,hyper_gene)
    df2 <- df2[-1]
    return(df2)
  }
  df2 <- adply(as.matrix(ceLNC),1,hyper_lncID)


  ## tidy the data
  hyperOutput <- df2[-1]
  hyperOutput <- hyperOutput[as.numeric(hyperOutput$realCounts) > 0,]
  if (is.null(deMIR)) {
    hyperOutput <- hyperOutput[,! colnames(hyperOutput) %in%
                                 c('deMIRCounts','deMIRs')]
  }

  ## Output result
  options(stringsAsFactors = old0$stringsAsFactors)
  return (hyperOutput)
}

# mir = mirs[1]
#lnc.i = "ENSG00000177410";pc.i="ENSG00000077943" ;lncDa <- unlist(rna.expr[lnc.i,]);pcDa <- unlist(rna.expr[pcDa,]);mir <-  mirs
mirCorTestFun_1 <- function(lncDa,pcDa,
                            mir,
                            mir.expr){
  if(mir %in% rownames(mir.expr)){
    ## miRNA在miR矩阵中
    mirDa <- unlist(mir.expr[mir,])
    corlm <- cor.test(lncDa, mirDa, alternative='less')
    corpm <- cor.test(pcDa, mirDa, alternative='less')
    reglm <- corlm$estimate
    regpm <- corpm$estimate
  } else {
    ## miRNA在miR矩阵中不存在，无法计算所谓的Cor
    reglm <- NA
    regpm <- NA
  }
  return (c(reglm, regpm)) ## lnc then pc
}

# lnc = "ENSG00000253552";pc = "ENSG00000171791";mirs = "hsa-miR-302a-3p,hsa-miR-302b-3p,hsa-miR-302c-3p,hsa-miR-372-3p,hsa-miR-373-3p";rna.expr = rna.expr;mir.expr = mir.expr
multiRegFun_1 <- function(lnc,
                          pc,
                          mirs,
                          rna.expr,
                          mir.expr){
  lncDa <- unlist(rna.expr[lnc,])
  pcDa <- unlist(rna.expr[pc,])

  corpl <- cor.test(pcDa, lncDa, alternative='greater')
  ppl <- corpl$p.value
  regpl <- corpl$estimate
  mirs <- unique(as.character(mirs))

  ## 计算敏感度和相似度
  mirs <- unlist(strsplit(mirs, ',', fixed=TRUE))
  mirCor <- vapply(mirs, function(mir)
    mirCorTestFun_1(lncDa, pcDa, mir, mir.expr),
    numeric(2))

  reglm <- mirCor[1,]
  regpm <- mirCor[2,]
  reglm_1 <- paste0(reglm,collapse = ",")
  regpm_1 <- paste0(regpm,collapse = ",")

  regSim <- 1-mean((abs(reglm - regpm)/(abs(reglm) + abs(regpm)))
                   ^length(mirs))
  sppc <- mean(regpl-(regpl-reglm*regpm)/(sqrt(1-reglm^2)*
                                            sqrt(1-regpm^2)))
  scores <- data.frame(pmcor = regpm_1,lmcor = reglm_1,plcor=regpl, corPValue=ppl, regSim = regSim, sppc = sppc,stringsAsFactors = F)
  return (scores)
}

multiRegTestFun_1 <- function(hyperOutput,
                              rna.expr,
                              mir.expr,
                              lnc.targets,
                              pc.targets){
  samples <- intersect(colnames(rna.expr), colnames(mir.expr))
  rna.expr <- rna.expr[,samples]

  ## 去除在系统中不存在的miRNA矩阵数据
  if (length(lnc.targets) > 1) {
    lnc.targets <- lnc.targets
  } else {
    lnc.targets <- ceRNA_lncTargets[[pc.targets]]
  }
  if (length(pc.targets) > 1) {
    pc.targets <- pc.targets
  } else {
    pc.targets <- ceRNA_pcTargets[[pc.targets]]
  }
  mir1 <- unique(unlist(lnc.targets))
  mir2 <- unique(unlist(pc.targets))
  mirs <- union(mir1,mir2)
  s1 <- intersect(mirs,rownames(mir.expr))
  mir.expr <- mir.expr[s1,samples]

  ## 提取lncRNA、pcRNA和miRNA的symbol
  lncID <- hyperOutput$lncRNAs
  pcID <- hyperOutput$Genes
  mirID <- hyperOutput$miRNAs

  ## 计算相似性和敏感度
  #reg <- vapply(
  #  seq_len(nrow(hyperOutput)),
  #  function(i)multiRegFun_1(lncID[i], pcID[i], mirID[i], rna.expr, mir.expr),
  #  numeric(6))
  get1 <- function(i){
    a <- multiRegFun_1(lncID[i],
                       pcID[i],
                       mirID[i],
                       rna.expr,
                       mir.expr)
    return(a)
  }
  reg.x <- adply(as.matrix(1:nrow(hyperOutput)),1,get1)

  ## 输出结果
  reg <- reg.x[-1]
  #colnames(reg) <- c("pmcor","lmcor",'plcor','corPValue','regSim', 'sppc')
  return (reg)
}


