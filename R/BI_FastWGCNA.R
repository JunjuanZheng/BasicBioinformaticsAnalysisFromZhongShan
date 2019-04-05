
##
#' @title Standard pipeline for WGCNA
#' @description FastWGCNA give a fast and standard pipeline to do Weighted Correlation Network Analysis(WGCNA) for specified expression matrix and design object
#' @keywords FastWGCNA
#' @param expr.matrix expression matrix
#' @param design a design object
#' @param log.convert whether do log2 scale for expression matrix
#' @param contrast.col  the colname of contrast in design object like "N.status"
#' @param contrast.control  the contrast control like "N0"
#' @param check whether to check critical objecta in save-file space:"dataExpr.rda","Net of WGCNA.rda","soft-thresholding list.rda".Set check=T help continue job based on the previous efforts
#' @param mad.portion If \code{mad.portion} <= 1, it represent the portion of data you want base on median absolute deviation(mad) filter.If \code{mad.portion} >1(for example 3000),it means the first \code{mad.portion} mad genes would be selected to use.
#' @param mad.min the lower cut-off of mad filter
#' @param verbose  a named vector of verboses in multiple scenes
#' @param maxBlockSize  set as large as possible according to the internal storage of your computer
#' @param minModuleSize the min count of module.Default is 30
#' @param cutoff.pval cut-off of the p value in significant module-phenotype filter
#' @param corType One of "pearson" and "bicor".Default is "pearson"
#' @param hub_cutoffSigGM  the cut-off of significant genes-Modules in hub genes exploration.Default is 0.2
#' @param hub_MM the cut-off of Module Memberships in hub genes exploration.Default is 0.8
#' @param hub_WeightedQ the cut-off of weighted q value in hub genes exploration.Default is 0.01
#' @param report whether to do a plot report
#' @param parallel whether use enableWGCNAThreads() to accelarate cor or bicor functions
#' @param two.step whether use two step to complete calculation.If the nGenes is large,"two.step = T" is recommaned.
#' @param save.path the space of the save file.Default is "WGCNA"
#' @param names  part of saved files name
#' @note 1.Please choose a high-performance computer for better running. /n 2.properly set "maxBlockSize" and "mad.portion" parameters to avoid the different block produced.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#'## This is a simulative process and NOT RUN
#'
#'## data preparation
#'load1(c("fpkm","design.model2"))
#'set.seed(2018);select <- sample(1:nrow(data.fpkm),2000)
#'expr.matrix = data.fpkm[select,]
#'design = design.model2
#'
#'## Quick Start
#'wgcna <- FastWGCNA(expr.matrix,
#'                   design,
#'                   contrast.col = "N.status",
#'                   contrast.control = "N0",
#'                   check = T,
#'                   parallel = F,
#'                   save.path = "WGCNA",
#'                   names = "test1");;mymusic(1)
#' @export
FastWGCNA <- function(expr.matrix,
                      design,
                      log.convert = T,
                      contrast.col = "N.status",
                      contrast.control = "N0",
                      check = F,
                      mad.portion = 0.75,
                      mad.min = 0.01,
                      verbose = c(
                        goodSamplesGenes = 3,
                        pickSoftThreshold = 5,
                        blockwiseModules = 3),
                      maxBlockSize = 50000,
                      minModuleSize = 30,
                      cutoff.pval = 0.05,
                      corType = c("pearson","bicor")[1],
                      hub_cutoffSigGM=0.2,
                      hub_MM = 0.8,
                      hub_WeightedQ = 0.01,
                      parallel = F,
                      report = T,
                      two.step = F,
                      save.path = "WGCNA",
                      names = "love"){

  ## 加载包
  nd <- c("AnnotationDbi", "impute","GO.db","preprocessCore", "stringr", "reshape2","WGCNA")
  Plus.library(nd)

  ## 系统环境备份
  old <- options()

  ## 内存管理
  # memory_size <- memory.limit()
  # if(memory_size < 32*1024){
  #   print("It seems the memory size in you computer is small.Please enhance virtual memory. ")
  #   new.m <- memory.limit(60*1024)
  #   print(paste0("Use virtual meory: ",new.m," MB."))
  # }


  ## 基础参数设置
  {
    # 字符不为因子
    options(stringsAsFactors = FALSE)

    # 官方推荐 "signed" 或 "signed hybrid"。为与原文档一致，故未修改
    type = "unsigned"

    # 相关性计算：官方推荐 biweight mid-correlation & bicor。corType: pearson or bicor。为与原文档一致，故未修改
    corFnc <- ifelse(corType=="pearson", cor, bicor)

    # 对二元变量，如样本性状信息计算相关性时，或基因表达严重依赖于疾病状态时，需设置下面参数
    maxPOutliers = ifelse(corType=="pearson",1,0.05)

    # 关联样品性状的二元变量时，设置
    robustY = ifelse(corType=="pearson",T,F)

    if(parallel==T) {enableWGCNAThreads()} #打开多线程
  }

  ## 产生储存文件夹
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = 0)

  ## 是否检查当前文件夹是否有关键文件
  if(check == T){
    checkFile <- list.files(dir)
    StepDF <- length(grep("dataExpr.rda",checkFile))==0;StepDF
    StepST <- length(grep("soft-thresholding list.rda",checkFile))==0;StepST
    StepNW <- length(grep("Net of WGCNA.rda",checkFile))==0;StepNW
    StepNwHm <-length(grep("Network heatmap plot.pdf",checkFile))==0;
    Stepscaletree <- length(grep("Check Scale free topology.pdf",checkFile))==0;
    StepCyt <- length(grep("Cytoscape Network.rda",checkFile))==0;
    if(!StepDF){load1("dataExpr.rda",path = dir,envir = environment())
      print("dataExpr.rda exists and loaded!");}
    if(!StepST){load1("soft-thresholding list.rda",path = dir,envir = environment())
      print("soft-thresholding list.rda exists and loaded!");}
    if(!StepNW){load1("Net of WGCNA.rda",path = dir,envir = environment())
      print("Net of WGCNA.rda exists and loaded!");}
  } else {
    StepDF <- StepST <- StepNW <-StepNwHm <- Stepscaletree <- StepCyt <- T
  }


  ###==========================数据筛选=======================###
  print("Step1:Data filtering...")
  {
    if(StepDF){
      ## 对齐
      dataExpr <- expr.matrix[,rownames(design)]
      if(log.convert == T){
        dataExpr <- log2(dataExpr + 1)
      }

      ## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01，筛选后会降低运算量，也会失去部分信息。也可不做筛选，使MAD大于0即可
      m.mad <- apply(dataExpr,1,mad)
      if(mad.portion <= 1){
        print(str_c("Filter: mad.portion = ",mad.portion,",mad >= ",mad.min))
        if(mad.portion < 1){
          dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad,probs=seq(0, 1, (1-mad.portion)))[2],mad.min)),]
        } else {
          dataExprVar <- dataExpr[m.mad > mad.min,]
        }
      } else {
        ## 选择mad排前几位的基因
        mad.min <- NULL
        cat(str_c("The first ",mad.portion," genes would be selected into WGCNA based on Median Absolute Deviation(MAD) order...","\n","The parameter 'mad.min' is useless and converted to NULL...","\n"))
        m.mad1 <- sort(m.mad,decreasing = T)
        dataExprVar <- dataExpr[names(m.mad1),]
        dataExprVar <- dataExprVar[1:mad.portion,]
      }

      ## 转换为样品在行，基因在列的矩阵
      dataExpr <- as.data.frame(t(dataExprVar)) # View(dataExpr[,1:5])
    }
    ## 检测缺失值
    gsg <-  goodSamplesGenes(dataExpr, verbose = verbose[names(verbose) == "goodSamplesGenes"])
    if (!gsg$allOK){
      # 有不合格的基因或者是不合格样本。提取合格基因
      if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:",paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
      if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:",paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
      # Remove the offending genes and samples from the data:
      dataExpr = dataExpr[gsg$goodSamples,gsg$goodGenes]
    }

  ## 基因数与样本数
  save(dataExpr,file = str_c(dir,names,"_dataExpr.rda"))
  rm(dataExprVar,envir = environment());gc()
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
}

  ###==========================软阈值筛选=====================###
  print("Step2:Analysis of scale free topology for soft-thresholding...")
  {
  if(StepST){
    ## 查看是否有离群样品。
    sampleTree = hclust(dist(dataExpr), method = "average")
    pdf(str_c(dir,names,"_Sample clustering to detect outliers.pdf"),15,10);plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="");dev.off()
    if(report == T){
      win.graph(15,10);plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")}

    ## 软阈值筛选:使构建的网络更符合无标度网络特征。
    powers <-  c(c(1:10), seq(from = 12, to=30, by=2))
    v2 <- verbose[names(verbose) == "pickSoftThreshold"]
    print("Analysis of scale free topology for soft-thresholding,please wait...")
    sft <-  pickSoftThreshold(dataExpr,
                              powerVector=powers,
                              networkType=type,
                              verbose = v2)
    save(sft,file = str_c(dir,names,"_soft-thresholding list.rda"))
    print("The list of soft-thresholding had been saved!")

    ## Report
    pdf(str_c(dir,names,"_Soft Threshold Evaluation.pdf"),15,10);
    {
      par(mfrow = c(1,2));cex1 = 1.5
      # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
      # 网络越符合无标度特征 (non-scale)
      plot(sft$fitIndices[,1],
           -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",
           ylab="Scale Free Topology Model Fit,signed R^2",
           type="n",
           main = paste("Scale independence"))
      text(sft$fitIndices[,1],
           -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex=cex1,col="red")
      abline(h=0.85,col="red")# 筛选标准。R-square=0.85

      # Soft threshold与平均连通性
      plot(sft$fitIndices[,1],
           sft$fitIndices[,5],
           xlab="Soft Threshold (power)",
           ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1],
           sft$fitIndices[,5],
           labels=powers,
           cex=cex1, col="red")
    }
    dev.off()
    if(report == T){
      win.graph(15,10)
      par(mfrow = c(1,2));cex1 = 1.2
      # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
      # 网络越符合无标度特征 (non-scale)
      plot(sft$fitIndices[,1],
           -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",
           ylab="Scale Free Topology Model Fit,signed R^2",
           type="n",
           main = paste("Scale independence"))
      text(sft$fitIndices[,1],
           -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex=cex1,col="red")
      abline(h=0.85,col="red")# 筛选标准。R-square=0.85

      # Soft threshold与平均连通性
      plot(sft$fitIndices[,1],
           sft$fitIndices[,5],
           xlab="Soft Threshold (power)",
           ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1],
           sft$fitIndices[,5],
           labels=powers,
           cex=cex1, col="red")}
    par(mfrow = c(1,1))

  }
  # disableWGCNAThreads()

  ## power值的确定
  power <- sft$powerEstimate
  #无向网络在power小于15或有向网络power小于30内，没有一个power值可以使无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值
  if (is.na(power)){
 power <-  ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),ifelse(type == "unsigned", 6, 12))))
 print("Power is NA,so an experiential power",power," would be used.")
  }
  print(paste0("The finally soft threshold power value is ",power))

  ## 检验软阈值的有效性
  print("Test the efficiency of soft threshold...")
  # 根据β值获得临近矩阵和拓扑矩阵
  softPower = power
  #adjacency = adjacency(dataExpr, power = softPower);# 获得临近矩阵：
  #TOM = TOMsimilarity(adjacency);# 将临近矩阵转为 Tom 矩阵
  #dissTOM = 1-TOM# 计算基因之间的相异度
  #hierTOM = hclust(as.dist(dissTOM),method="average");

  # 检验选定的β值下记忆网络是否逼近 scale free
  if(Stepscaletree){
    if(ncol(dataExpr) < 5000){
      # 基因少（<5000）的时候使用下面的代码：
      ADJ1_cor <- abs(WGCNA::cor(dataExpr,use = "p" ))^softPower
      k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
    } else {
      # 基因多的时候使用下面的代码：
      k <- softConnectivity(datE=dataExpr,power=softPower)
    }

    # k与p(k)成负相关(相关性系数0.9),说明选择的β值能够建立基因无尺度网络。
    pdf(paste0(dir,names,"_Check Scale free topology.pdf"),10,5)
    par(mfrow=c(1,2))
    hist(k)
    test <- scaleFreePlot(k,main="Check Scale free topology\n")
    dev.off()
    par(mfrow=c(1,1))
    if(report == T){
      win.graph(12,6)
      par(mfrow=c(1,2))
      hist(k)
      test <- scaleFreePlot(k,main="Check Scale free topology\n")
      par(mfrow=c(1,1))
    }
    print(paste0("Soft threshold test: scaleFreeRsquared = ",test$scaleFreeRsquared,";slope = ",test$slope))
  }

}

  ###========================网络构建=======================###
  print("Step3:Automatic network construction and module detection...")
  {
  ##一步法网络构建：One-step network construction and module detection##
  # power: 上一步计算的软阈值
  # maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个。计算资源允许的情况下最好放在一个block里面。
  # corType: pearson or bicor
  # mergeCutHeight: 合并模块的阈值，越大模块越少
  v3 <- verbose[names(verbose) == "blockwiseModules"]

  ## 虽然此步慢，但不需要用并行。为节省内存，可提前终止并行要求
  if(StepNW){
    net <-  blockwiseModules(
      dataExpr,
      power = power,
      maxBlockSize = maxBlockSize,#其实没必要设置这一参数，默认即可
      TOMType = type,
      minModuleSize = minModuleSize,
      reassignThreshold = 0,
      mergeCutHeight = 0.25,#0.25以上的树保留，以下的数去除
      numericLabels = TRUE,#返回数字为模块名字，后面可转换为颜色
      pamRespectsDendro = FALSE,
      saveTOMs=TRUE,#最耗费时间的计算，存储起来，供后续使用
      corType = corType,
      maxPOutliers=maxPOutliers,
      loadTOMs=TRUE,
      saveTOMFileBase = paste0(dir,names,"_TOMFileBase.tom"),
      verbose = v3
    )
    save(net,file = str_c(dir,names,"_Net of WGCNA.rda"))
  }

  ##层级聚类树展示各个模块
  print("Print Cluster Dendrogram...")
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  pdf(str_c(dir,names,"_Cluster Dendrogram.pdf"),15,12)
  for(i in 1:length(net$dendrograms)){ # i=1
    plotDendroAndColors(
      dendro = net$dendrograms[[i]],
      colors = moduleColors[net$blockGenes[[i]]],
      groupLabels ="Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      cex.colorLabels = 1,
      main = paste0("Block",i,":Cluster Dendrogram")
    )
  }
  dev.off()
  ## Report
  if(report == T){
    win.graph(15,12)
    for(i in 1:length(net$dendrograms)){ # i=1
      plotDendroAndColors(
        dendro = net$dendrograms[[i]],
        colors = moduleColors[net$blockGenes[[i]]],
        groupLabels ="Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05,
        cex.colorLabels = 1,
        main = paste0("Block",i,":Cluster Dendrogram")
      )
    }

  }

  ## 绘制模块之间相关性图
  # module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
  MEs_col <- MEs <- net$MEs
  colnames(MEs_col) = paste0("ME",labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  net$MEs <- MEs_col #形成新的ME矩阵(列为Modules,行为患者)

  # plot
  print("Print Eigengene adjacency heatmap...")
  pdf(str_c(dir,names,"_Eigengene adjacency heatmap.pdf"),15,12)
  if(ncol(MEs_col) >= 3){
    for(i in c(T,F)){
      plotEigengeneNetworks(multiME = MEs_col,
                            setLabels = "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2),
                            plotDendrograms = i,
                            xLabelsAngle = 90)
    }
    dev.off()
    if(report == T){
      win.graph(15,10)
      plotEigengeneNetworks(multiME = MEs_col,
                            setLabels = "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2),
                            plotDendrograms = T,
                            xLabelsAngle = 90)
    }
  } else {
    print("The MEs col counts is too small to draw a Eigengene adjacency heatmap...")
  }

}

  ###==================模块与表型数据关联====================###
  print("Step4:Module and phenotype analysis...")
  {
    level <- unique(as.character(design[,contrast.col]))
    level <- c(contrast.control,setdiff(level,contrast.control))
    group <- factor(design[,contrast.col],levels = level)
    traitData <- model.matrix(~ 0 + group)
    colnames(traitData) <-  level

    ## 计算相关系数矩阵
    if (corType=="pearson") {
      modTraitCor = cor(MEs_col,traitData, use = "p")
      modTraitP = corPvalueStudent(modTraitCor, nSamples)
    } else {
      modTraitCorP = bicorAndPvalue(MEs_col,
                                    traitData,
                                    robustY=robustY)
      modTraitCor = modTraitCorP$bicor
      modTraitP   = modTraitCorP$p
    }

    ## correlation heatmap
    textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
    dim(textMatrix) <-  dim(modTraitCor)
    #尺寸
    h <- dim(modTraitCor)[1]; w <- dim(modTraitCor)[2]
    h <- h*(10/h);w <- w*(10/h)
    win.graph(15,15);
    test <- labeledHeatmap(
      Matrix = modTraitCor,
      xLabels = colnames(traitData),
      yLabels = colnames(MEs_col),
      ySymbols = colnames(MEs_col),
      colorLabels = FALSE,
      colors = blueWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = T,
      cex.lab = 1,
      cex.text =1, zlim = c(-1,1),
      main = paste("Module-Trait relationships")
    )
    pdf(str_c(dir,names,"_Module-Trait correlation.pdf"),6,15)
    test <- labeledHeatmap(
      Matrix = modTraitCor,
      xLabels = colnames(traitData),
      yLabels = colnames(MEs_col),
      ySymbols = colnames(MEs_col),
      colorLabels = FALSE,
      colors = blueWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = T,
      cex.lab = 1,
      cex.text =1, zlim = c(-1,1),
      main = paste("Module-Trait relationships")
    )
    dev.off()

  }

  ###==================寻找关键模块和中心基因================###
  print("Step5:Build network and find hub genes...")
  ## 加载TOM文件
  tomfile <- list.files(path = dir,
                        pattern = "TOMFileBase.tom",
                        full.names = T)
  print(str_c("All is ",length(tomfile)," TOM object."))

  Hub <- NULL
  for(i in 1:length(tomfile)){ # i = 1
    ## load TOM matrix
    print(str_c("Block",i,": making TOM matrix..."))
    load(tomfile[i])
    TOM <- as.matrix(TOM)
    dissTOM = 1-TOM;
    plotTOM = dissTOM^power;diag(plotTOM) = NA
    if(StepNwHm){
      print(str_c("Block",i,": Calculating Topological overlap matrix..."))
      if(nGenes > 5000){
        print(str_c("Block",i,": nGenes is too large.Randomly select 5000 genes into plot."))
        set.seed(2018);random <- sample(1:nrow(plotTOM),5000)
        plotTOM.i = dissTOM[random, random];
        dendro.tom = hclust(as.dist(plotTOM.i), method = "average")
        color.tom = moduleColors[net[["blockGenes"]][[i]]][random];#选择颜色
        plotTOM.i = plotTOM.i^power;diag(plotTOM) = NA;
        pdf(str_c(dir,names,"_Block",i,"_Network heatmap plot.pdf"),12,12)
        print(str_c("Block",i,": drawing the Network heatmap,please wait..."))
        TOMplot(dissim = plotTOM.i,
                dendro = dendro.tom,
                Colors = color.tom,
                main = "Part genes:Network heatmap plot")
        dev.off()
        print(str_c("Block",i,": Network heatmap plot had been saved!"))
        rm(plotTOM.i,envir = environment());gc()
      } else {
        dendro.tom <- net$dendrograms[[i]]
        color.tom <- moduleColors[net[["blockGenes"]][[i]]]

        ## plot Network heatmap plot for all genes
        if(ncol(dataExpr) < 3600){
          pdf(str_c(dir,names,"_Block",i,"_Network heatmap plot.pdf"),12,12)
        } else {
          png(str_c(dir,names,"_Block",i,"_Network heatmap plot.png"),width = 1440, height = 1440)
        }
        print(str_c("Block",i,": drawing the Network heatmap,please wait..."))
        TOMplot(dissim = plotTOM,
                dendro = dendro.tom,
                Colors = color.tom,
                main = "All genes:Network heatmap plot")
        dev.off()
        print(str_c("Block",i,": Network heatmap plot had been saved!"))
      }
    }
    rm(dissTOM,envir = environment());gc()

    ## creat a network object
    probes <-  colnames(dataExpr)[net[["blockGenes"]][[i]]]
    dendro.tom <- net$dendrograms[[i]]
    color.tom <- moduleColors[net[["blockGenes"]][[i]]]
    if(StepNwHm & nGenes > 5000 & two.step == T){
      stop("Attention!This is not a real error information.It appears because you set 'tow.step=T' and the nGenes is > 5000. The following cytoscape network building would spend large RAM,so We pause the process here.Please set 'check = T' if you haven't done before and RUN YOUR SCRIPT AGAIN.")
      }
    if(StepCyt){
      dimnames(TOM) <- list(probes, probes)
      print(str_c("Block",i,": making Cytoscape network object,please wait..."))
      cyt <- exportNetworkToCytoscape(
        TOM,
        #edgeFile = str_c(save.path,names,"_edges.txt"),
        #nodeFile = str_c(save.path,names,"_nodes.txt"),
        weighted = TRUE, threshold = 0,
        altNodeNames = convert(probes),
        nodeNames = probes,
        nodeAttr = color.tom)
      save(cyt,file = str_c(dir,names,"_","Block",i,"_Cytoscape Network.rda"))
      print(str_c("Block",i,": Cytoscape object of the network had been saved!"))
    }
    rm(TOM,envir = environment());gc()

    ## 根据某block里提取dataExpr数据
    p.i <- net[["blockGenes"]][[i]]
    dataExpr.i <- dataExpr[,p.i]

    ## 计算某block里模块与基因的相关性矩阵
    if (corType=="pearson") {
      geneModuleMembership.i = as.data.frame(cor(dataExpr.i, MEs_col, use = "p"))
      MMPvalue.i = as.data.frame(corPvalueStudent(
        as.matrix(geneModuleMembership.i), nSamples))
    } else {
      geneModuleMembershipA.i = bicorAndPvalue(dataExpr.i, MEs_col, robustY=robustY)
      geneModuleMembership.i = geneModuleMembershipA.i$bicor
      MMPvalue.i   = geneModuleMembershipA.i$p
    }

    ## 计算某block里基因与性状的相关系数
    if (corType=="pearson") {
      geneTraitCor.i = as.data.frame(cor(dataExpr.i, traitData, use = "p"))
      geneTraitP.i = as.data.frame(corPvalueStudent(
        as.matrix(geneTraitCor.i), nSamples))
    } else {
      geneTraitCorA.i = bicorAndPvalue(dataExpr.i, traitData, robustY=robustY)
      geneTraitCor.i = as.data.frame(geneTraitCorA.i$bicor)
      geneTraitP.i   = as.data.frame(geneTraitCorA.i$p)
    }

    ## 计算差异性模块-性状
    sig1 <- modTraitP < cutoff.pval
    if(all(sig1 == F)){
      print("No significant module-pheno memberships")
      pairs <- NULL
    } else {
      get1 <- function(sig1){
        p1 <- grep(T,sig1)
        r1 <- ifelse(p1%%nrow(sig1) != 0,p1%%nrow(sig1),nrow(sig1))
        c1 <- ifelse(p1%%nrow(sig1) == 0,p1%/%nrow(sig1),p1%/%nrow(sig1)+1)
        df1 <- cbind(r1,c1)
        # df1.1 <- df1[1,]
        get2 <- function(df1.1){
          n1 <- rownames(sig1)[df1.1[1]]
          n2 <- colnames(sig1)[df1.1[2]]
          return(c(n1,n2))
        }
        df2 <- t(apply(df1,1,get2))

        ## 改rownames
        df2[,1] <- substring(df2[,1],3)
        return(df2)
      }
      pairs <- get1(sig1)
      pdf(str_c(dir,names,"_Module membership vs. gene significance.pdf"),10,10)
      for(z in 1:nrow(pairs)){ # i = 1
        pair.i <- pairs[z,]
        # 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
        module <-  pair.i[1]
        pheno <-  pair.i[2]
        modNames <-  substring(colnames(MEs_col),3)
        # 获取关注的列
        module_column <-  match(module, modNames)
        pheno_column <-  match(pheno,colnames(traitData))
        # 获取模块内的基因
        moduleGenes <-  moduleColors == module

        par(mfrow = c(1,1))
        # 与性状高度相关的基因，也是与性状相关的模型的关键基因
        verboseScatterplot(
          x=abs(geneModuleMembership.i[moduleGenes, module_column]),
          y=abs(geneTraitCor.i[moduleGenes, pheno_column]),
          xlab = paste("Module Membership in", module, "module"),
          ylab = paste("Gene significance for", pheno),
          main = paste("Module membership vs. gene significance\n"),
          cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
          col = module
        )

      }
      dev.off()
      print(str_c("Block",i,": Module membership vs. gene significance plot had been saved!"))
    }

    ##=====================关键Module筛选=====================##
    ## 获得关键Module
    color.i <- moduleColors[p.i]
    ModuleSignificance <- apply(geneTraitCor.i,2,function(x)tapply(x,color.i,mean))
    pheno.i <- setdiff(colnames(ModuleSignificance),contrast.control)
    ModuleSignificance <- ModuleSignificance[,pheno.i] #找到非对照的相关性平均值
    max.module.i <- which.max(ModuleSignificance[names(ModuleSignificance) != "grep"])
    max.module.i <- names(max.module.i)
    print(paste0("Block",i,": The most related module is ",max.module.i,"."))

    ## 计算模块内连接度
    Alldegrees1 <- intramodularConnectivity(
      adjMat = plotTOM,
      colors = color.i)

    ## 计算module membership。与内部连接度不同,MM 衡量了基因在全局网络中的位置。
    datKME <- signedKME(dataExpr.i,MEs_col, outputColumnName="MM.")

    ## 绘图:模块内-模块间连接度的相关性
    pdf(str_c(dir,names,"_Intramodular Connectivity vs. Module Membership.pdf"),10,10)
    which.color=max.module.i
    restrictGenes <- color.i == which.color
    verboseScatterplot(
      Alldegrees1$kWithin[restrictGenes],
      (datKME[restrictGenes, paste("MM.", which.color, sep="")])^4,
      col=which.color,
      xlab="Intramodular Connectivity",
      ylab="(Module Membership)^4",
      main = paste0(pheno.i," : ",which.color,"\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
    )
    dev.off()
    if(report == T){
      win.graph(10,10)
      verboseScatterplot(
        Alldegrees1$kWithin[restrictGenes],
        (datKME[restrictGenes, paste("MM.", which.color, sep="")])^4,
        col=which.color,
        xlab="Intramodular Connectivity",
        ylab="(Module Membership)^4",
        main = paste0(pheno.i," : ",which.color,"\n"),
        cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
      )
    }
    print(paste0("Block",i,"_Module-",max.module.i,"_Intramodular Connectivity vs. Module Membership plot had been saved!"))
    rm(list = c("Alldegrees1"),envir = environment());gc()

    ##================使用3个标准来筛选枢纽基因==============##
    #基因与指定模块显著性 > 0.2
    #最大MM值 > 0.8
    #加权q值 < 0.01

    ## 计算基因-性状的biweight midcorrelations
    #GS_spe <- bicorAndPvalue(dataExpr.i,
    #                         traitData,
    #                         robustY=robustY)
    #GeneSignificance_spe <- abs(GS_spe$bicor[,pheno.i])
    GeneSignificance_spe <- abs(geneTraitCor.i[,pheno.i])

    ## 基于显著性和MM计算每个基因与指定性状的关联，结果包括p, q, cor, z
    p <- Fastmatch(rownames(dataExpr.i),rownames(design))
    trait <- design[p,contrast.col]
    trait <- ifelse(trait == contrast.control,0,1)
    NS1 <- networkScreening(y=trait,
                            datME=MEs_col,
                            datExpr=dataExpr,
                            oddPower=3,
                            blockSize=1000,
                            minimumSampleSize=4,
                            addMEy=TRUE,
                            removeDiag=FALSE,
                            weightESy=0.5)#?

    ## 根据基因与指定性状的直接相关性(biserial.cor)，模块身份，和加权相关性筛选基因：
    FilterGenes_spe <- ((GeneSignificance_spe > hub_cutoffSigGM) & (abs(datKME[,paste("MM.",max.module.i,sep="")]) > hub_MM) & (NS1$q.Weighted <= hub_WeightedQ))
    test <- table(FilterGenes_spe)
    if(T %in% names(test)){
      print(paste0("There are ",test[names(test) %in% T]," hub genes in Module ",max.module.i,"..."))
      ## 找到满足上面条件的基因：
      trait_hubGenes_spe <- colnames(dataExpr.i)[FilterGenes_spe]
    } else {
      print(paste0("There are no hub gene in Module ",names(max.module.i),"..."))
      trait_hubGenes_spe <- NULL
    }

    ## hub 基因热图：
    if(!is.null(trait_hubGenes_spe)){
      p.names <- paste0(dir,names,"_Module-",max.module.i,"_Heatmap of hub gene unsigned correlations.pdf")
      main <- paste0("Block",i,"_Module",max.module.i,"_unsigned correlations")
      pdf(p.names,15,15)
      plotNetworkHeatmap(dataExpr.i,
                         plotGenes = trait_hubGenes_spe,
                         networkType = "unsigned",
                         useTOM = TRUE,
                         power=power,
                         main=main)
      dev.off()
      win.graph(15,15)
      plotNetworkHeatmap(dataExpr.i,
                         plotGenes = trait_hubGenes_spe,
                         networkType = "unsigned",
                         useTOM = TRUE,
                         power=power,
                         main=main)
    }

    ## 结果汇总
    print(paste0("Block",i,": Output WGCNA result,please wait..."))
    load1(paste0("Block",1,"_Cytoscape Network.rda"),dir)
    l1 <- list(
      plotTOM =plotTOM,
      cyt = cyt,
      Cor = list(
        Module_Trait = list(cor = modTraitCor,
                            p.val = modTraitP),
        Module_Gene=list(cor = geneModuleMembership.i,
                         p.val = MMPvalue.i),
        Gene_Trait = list(cor = geneTraitCor.i,
                          p.val = geneTraitP.i)
      ),
      SigPairs = pairs,
      MaxModule = max.module.i,
      HubGenes = trait_hubGenes_spe
    )
    Hub <- c(Hub,list(l1))
    names(Hub)[i] <- paste0("Block",i)
    rm(list = c("cyt","l1","plotTOM"),envir = environment());gc()
  }

  ###=======================保存结果=========================###
  wgcna <- list(
    dataExpr = dataExpr,
    sft = sft,
    net = net,
    Hub = Hub,
    Repeat = list(
      log.convert = log.convert,
      contrast.col = contrast.col,
      contrast.control = contrast.control,
      mad.portion = mad.portion,
      mad.min = mad.min,
      verbose = verbose,
      maxBlockSize = maxBlockSize,
      cutoff.pval = cutoff.pval,
      corType = corType,
      hub_cutoffSigGM=hub_cutoffSigGM,
      hub_MM = hub_MM ,
      hub_WeightedQ =  hub_WeightedQ,
      save.path = save.path,
      names = names
    )
  )
  save(wgcna,file = paste0(dir,names,"_wgcna.rda"))
  on.exit(options(old), add = TRUE)
  print("All done!")
  return(wgcna)
}



#' @title Before WGCNA
#' @description Test the number of genes selected in WGCNA
#' @keywords BeforeWGCNA
#' @param expr.matrix expression matrix
#' @param design a design object
#' @param log.convert whether do log2 scale for expression matrix
#' @param mad.portion a portion of data you want base on median absolute deviation(mad) filter
#' @param mad.min the lower cut-off of mad filter
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
BeforeWGCNA <- function(expr.matrix,
                        design,
                        log.convert=T,
                        mad.portion=0.75,
                        mad.min = 0.01){

  ## 对齐
  dataExpr <- expr.matrix[,rownames(design)]

  ## log转换
  if(log.convert == T){
    dataExpr <- log2(dataExpr + 1)
  }

  ## 求mad
  m.mad <- apply(dataExpr,1,mad)

  ## 筛选
  if(mad.portion < 1){
    dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad,probs=seq(0, 1, (1-mad.portion)))[2],mad.min)),]
  } else {
    dataExprVar <- dataExpr[m.mad > mad.min,]
  }

  ## 输出结果
  a <- paste0(nrow(dataExprVar)," Genes.")
  return(a)
}



