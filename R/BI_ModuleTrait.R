

#' @keywords ModuleTrait
#' @title Enhanced Module-Trait Relationships Plot after FastWGCNA Process
#' @description After \code{\link{FastWGCNA}}, you can input the result into \code{ModuleTrait} function,and get a more powerful Module-Trait relationships visualization.
#' @param object the result of  \code{\link{FastWGCNA}}.
#' @param variable the variables you want to show in Module-Trait relationships plot.
#' @param design a trait-design object
#' @param corType one of "pearson" and "bicor".Default is pearson
#' @param report whether do a plot report
#' @param height,width the size of saved plot
#' @param save.path the space of the save file.Default is "WGCNA"
#' @param names  part of saved files name
#' @importFrom WGCNA corPvalueStudent cor bicorAndPvalue labeledHeatmap blueWhiteRed plotModuleSignificance labels2colors
#' @details trait-design object means the design should convert to a numeric one.A good method is to use \code{\link{GetWGCNATrait}}.
#' @return a Module-Trait relationships plot
#' @seealso \code{\link{GetWGCNATrait}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and NOT RUN
#' library(lucky)
#' object = wgcna;rm(wgcna);gc()
#' design = rna.design.tumor
#' variable = c("age","his1","gender","N.status","T.status")
#' save.path = "WGCNA"
#' names = "love"
#' corType = "pearson"
#' height = 10;width = 8
#' result_MT <- ModuleTrait(object,
#'                          design,
#'                          variable,
#'                          corType = "pearson",
#'                          report=T,
#'                          height = 10,width = 8,
#'                          save.path = "WGCNA",
#'                          names = "love")
#' @export
ModuleTrait <- function(object,
                        design,
                        variable,
                        corType = c("pearson","bicor")[1],
                        report=F,
                        height = 10,width = 8,
                        save.path = "WGCNA",
                        names = "love"){
  ### Background
  # There are two approaches one can follow. We discuss each below.
  #6.e.1 Correlate the module eigengenes with the trait
  #6.e.2 Measure of module significance as average gene significance

  ## load needed package
  #nd <- c("plyr","WGCNA")
  #Plus.library(nd)

  ## 产生储存文件夹
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = 0)

  ## 系统变量
  robustY <-  ifelse(corType=="pearson",T,F)

  ## get design matrix from every variable
  traitData <- design[,variable]
  #f <- paste0("~ 0 + ",paste0(variable,collapse = " + "))
  #f <- as.formula(f)
  #traitData <- model.matrix(f,data = design_1)

  ## get ME matrix from wcgna object
  object.names <- names(object)
  lg1 <- all(c("dataExpr","sft","net","Hub") %in% object.names)
  if(lg1){
    if(report == T) {LuckyVerbose("You input a LuckyWGCNA object.")}
    dataExpr <- object$dataExpr
    MEs_col <- object$net$MEs
    nSamples <- nrow(MEs_col)
  } else {
    if(report == T) {LuckyVerbose("You input a ME object.")}
    dataExpr <- NULL
    MEs_col <- object
    nSamples <- nrow(MEs_col)
  }

  ## Get gene color
  moduleLabels <-  object$net$colors
  moduleColors <-  labels2colors(moduleLabels)
  # length(moduleColors)

  ###==================Module-Trait Plot===================###
  if(report == T) LuckyVerbose("Step1:Plot Module-Trait relationship...")
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

  ## report
  if(report==T){
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
  }

  ## print the plot
  pdf(paste0(dir,names,"_Module-Trait correlation.pdf"),width,height)
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
  if(report == T) LuckyVerbose("The Module-Trait relationships plot had been saved...",levels = 2)


  ###==============Module-Significance Plot================###
  # The advantage of this second approach is that it can be used for any gene significance measure. A gene significance measure could be defined without reference to a sample trait. For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
  # 对于不同的trait是有不同的Module-Significance Plot的

  if(!is.null(dataExpr)){
    if(report == T) LuckyVerbose("Step2: Plot module significance...")
    GS <- NULL;MS <- NULL
    pdf(paste0(dir,names,"_Gene significance across modules.pdf"),nrow(modTraitCor),8)
    for(i in 1:length(variable)){ # i=1
      var.i <- variable[i]
      y <- design[,var.i]
      ## 计算某个模板里每个基因的gene significance
      grep.p <- grep("grey",moduleColors)
      colorh1 <- moduleColors[-grep.p]
      dataExpr1 <- dataExpr[,-grep.p] #去除grey module
      GS1 <- as.numeric(cor(y,dataExpr1, use="p"))
      GS <- c(GS,list(GS1));names(GS)[i] <- var.i
      GS2 <- abs(GS1)
      MS.i <- tapply(GS2,colorh1, mean, na.rm=T)
      MS <- c(MS,list(MS.i));names(MS)[i] <- var.i
      title <- paste0("Gene significance across modules_",var.i,",")

      plotModuleSignificance(geneSignificance=GS2,
                             color = colorh1,
                             main = title)

      if(report == T) LuckyVerbose(var.i,"Module-Significance Plot had been saved...",levels = 2)
      }
    dev.off()
    } else {
    ## 不是wgcna对象，无法进行绘图
    GS <- NULL
    MS <- NULL
  }

  ###=====================输出表格化结果===================###

  textMatrix2 <- paste(signif(modTraitCor, 2), "(", signif(modTraitP, 1), ")", sep = "")
  dim(textMatrix2) <-  dim(modTraitCor)
  colnames(textMatrix2) <- colnames(modTraitCor)
  rownames(textMatrix2) <- rownames(modTraitCor)
  write.csv(textMatrix2,paste0(dir,names,"_Module-Trait relationships.csv"))
  if(report == T) LuckyVerbose("Step3: Save meta data..")

  ## 输出结果
  result <- list(
    Repeat = list(variable = variable,
                  corType = corType,
                  height = height,width = width,
                  save.path = save.path,
                  names = names),
    Data = list(ModuleTrait = list(traitData = traitData,
                                   modTraitCor = modTraitCor,
                                   modTraitP = modTraitP,
                                   textMatrix =  textMatrix),
                AverageGS = list(GeneSignificance = GS,
                                 TraitAverageGS = MS)),
    Plot = NULL
  )
  if(report == T) LuckyVerbose("All done!")
  return(result)
}



