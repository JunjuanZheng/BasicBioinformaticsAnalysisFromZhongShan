

#' @keywords ModuleHub
#' @title Enhanced hub genes exploration after FastWGCNA pipeline
#' @description The speed of \code{ModuleHub} is very fast.If you want to set different parameters(always \code{hub_WeightedQ} and \code{cutoff.pval}) for hub genes exploration, \code{ModuleHub} is a nice function to do it.
#' @inheritParams ModuleTrait
#' @inheritParams FastWGCNA
#' @importFrom WGCNA cor corPvalueStudent bicorAndPvalue labels2colors
#'   verboseScatterplot intramodularConnectivity signedKME networkScreening
#' @importFrom stringr str_c
#' @details 1.Only work on wgcna result with ONE Block.   2.\code{cutoff.pval}
#'   is a useful parameter.If you want significant module,you can set
#'   \code{cutoff.pval = 0.05}(or 0.01,It depends on your custom);If you want to
#'   see all the modules regardless of significance,just set \code{cutoff.pval =
#'   1}.  3.\code{hub_WeightedQ} is a stricter filter for hub genes.If your hub
#'   genes is too much,you can set \code{hub_WeightedQ = 0.05}(or 0.01,It
#'   depends on your custom).However,most researchers do not use
#'   \code{hub_WeightedQ} to filter their hub genes and often use only
#'   \code{hub_MM=0.8} and \code{hub_cutoffSigGM=0.2}.
#' @return LuckWGCNA object
#' @seealso \code{\link{FastWGCNA}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' library(lucky)
#' load("E:/RCloud/RFactory/lucky/love/WGCNA-test/love_wgcna.rda")
#' object = wgcna;rm(wgcna);gc()
#' design = rna.design.tumor
#' variable = c("age","his1","gender","N.status","T.status")
#' result_MH <- ModuleHub(object,
#'                        design,
#'                        variable,
#'                        save.path = "WGCNA",
#'                        names = "love")
#' @export
ModuleHub <- function(object,
                      design,
                      variable,
                      corType = "pearson",
                      cutoff.pval = 1,
                      hub_cutoffSigGM=0.2,
                      hub_MM = 0.8,
                      hub_WeightedQ = 1,
                      save.path = "WGCNA",
                      names = "love"){
  ## load needed package
  #nd <- c("WGCNA","stringr")
  #Plus.library(nd)

  ## 全局变量
  old <- options()

  ## 产生储存文件夹
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = 0)

  ## 全局变量
  robustY <- ifelse(corType=="pearson",T,F)

  ## create traitData
  traitData <- design[,variable]

  ## 计算各种参数
  tomfile <- object$net$TOMFiles
  Hub <- NULL;Data <- NULL
  for(i in 1:length(tomfile)){ # i=1

    ## 加载有用的结果
    net <- object$net
    dataExpr <- object$dataExpr
    nSamples <- nrow(dataExpr)
    MEs_col <- net$MEs
    moduleLabels <-  net$colors
    moduleColors <-  labels2colors(moduleLabels)
    plotTOM <- object$Hub$Block1$plotTOM

    ## modulegenes
    genes <- colnames(dataExpr)
    names(genes) <- moduleColors
    module.type <- unique(as.character(moduleColors))
    modulegenes <- NULL
    for(z in 1:length(module.type)){
      module.z <- module.type[z]
      genes.i <- genes[names(genes) %in% module.z]
      modulegenes <- c(modulegenes,list(genes.i))
      names(modulegenes)[z] <- module.z
    }

    ## 根据某block里提取dataExpr数据
    p.i <- net[["blockGenes"]][[i]]
    dataExpr.i <- dataExpr[,p.i]

    ## 计算Module-Trait相关系数矩阵
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

    ## 计算Module-Gene相关性矩阵
    if (corType=="pearson") {
      geneModuleMembership.i = as.data.frame(cor(dataExpr.i, MEs_col, use = "p"))
      MMPvalue.i = as.data.frame(corPvalueStudent(
        as.matrix(geneModuleMembership.i), nSamples))
    } else {
      geneModuleMembershipA.i = bicorAndPvalue(dataExpr.i, MEs_col, robustY=robustY)
      geneModuleMembership.i = geneModuleMembershipA.i$bicor
      MMPvalue.i   = geneModuleMembershipA.i$p
    }

    ## 计算Gene-Trait相关系数矩阵
    if (corType=="pearson") {
      geneTraitCor.i = as.data.frame(cor(dataExpr.i, traitData, use = "p"))
      geneTraitP.i = as.data.frame(corPvalueStudent(
        as.matrix(geneTraitCor.i), nSamples))
    } else {
      geneTraitCorA.i = bicorAndPvalue(dataExpr.i, traitData, robustY=robustY)
      geneTraitCor.i = as.data.frame(geneTraitCorA.i$bicor)
      geneTraitP.i   = as.data.frame(geneTraitCorA.i$p)
    }

    ## Module_Gene相关性的可视化
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

      ## Module_Gene
      pdf(str_c(dir,names,"_ModuleHub_MM-GS relationship.pdf"),10,10)
      for(z in 1:nrow(pairs)){ # z = 1
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
        abline(v=hub_MM,h=hub_cutoffSigGM,lty=2,col="red")
      }
      dev.off()
      cat(str_c("MM-GS relationship plot had been saved!"),"\n")
      cat("Note: In MM-GS plot,the points in top right corner are target hub genes.","\n")
      }

    ##================使用3个标准来筛选枢纽基因==============##
    #基因与指定模块显著性 > 0.2
    #最大MM值 > 0.8
    #加权q值 < 0.01

    ## 计算intramodular connectivity
    Alldegrees1 <- intramodularConnectivity(
      adjMat = plotTOM,
      colors = moduleColors[p.i])

    ## MM 衡量了基因在全局网络中的位置
    datKME <- signedKME(dataExpr.i,MEs_col, outputColumnName="MM.")

    ## Hub Genes
    hubgenes <- NULL;NS <- NULL
    for(z in 1:length(variable)){ #i=1

      pheno.i <- variable[z]

      ## 计算基因-性状的biweight midcorrelations
      GeneSignificance_spe <- abs(geneTraitCor.i[,pheno.i])

      ## 基于显著性和MM计算每个基因与指定性状的关联，结果包括p, q, cor, z
      p <- Fastmatch(rownames(dataExpr.i),rownames(design))
      trait <- design[,pheno.i]
      NS1 <- networkScreening(y=trait,
                              datME=MEs_col,
                              datExpr=dataExpr.i,
                              oddPower=3,
                              blockSize=1000,
                              minimumSampleSize=4,
                              addMEy=TRUE,
                              removeDiag=FALSE,
                              weightESy=0.5)#?
      NS1 <- NS1[colnames(dataExpr.i),]
      NS <- c(NS,list(NS1));names(NS)[z] <- pheno.i

      ## 在不同的模块中选出所谓的hub genes
      hubgenes.i <- NULL
      for(j in 1:length(unique(moduleColors))){ # j=1
        module.j <- unique(moduleColors)[j]
        FilterGenes_spe <- ((GeneSignificance_spe > hub_cutoffSigGM) & (abs(datKME[,paste("MM.",module.j,sep="")]) > hub_MM) & (NS1$q.Weighted <= hub_WeightedQ))
        test <- table(FilterGenes_spe)
        if(T %in% names(test)){
          cat(paste0("There are ",test[names(test) %in% T]," hub genes in Module ",module.j,"..."),"\n")
          ## 找到满足上面条件的基因：
          trait_hubGenes_spe <- colnames(dataExpr.i)[FilterGenes_spe]
        } else {
          cat(paste0("There are no hub gene in Module ",module.j,"..."),"\n")
          trait_hubGenes_spe <- NULL
        }
        hubgenes.i <- c(hubgenes.i,list(trait_hubGenes_spe))
        names(hubgenes.i)[j] <- module.j
      }
      hubgenes <- c(hubgenes,list(hubgenes.i))
      names(hubgenes)[z] <- pheno.i
    }

    ## 输出结果
    l1 <- list(
      pairs = pairs,
      datKME  = datKME,
      networkScreening = NS,
      modulegenes = modulegenes,
      hubgenes =  hubgenes
    )
    Hub <- c(Hub,list(l1))
    names(Hub)[i] <- paste0("Block",i)

    Data <- c(Data,list(dataExpr.i))
    names(Data)[i] <- paste0("Block",i)

  }

  ## output result
  result <- list(
    Repeat = list(
      traitData = traitData,
      corType = corType,
      cutoff.pval = cutoff.pval,
      hub_cutoffSigGM=hub_cutoffSigGM,
      hub_MM = hub_MM,
      hub_WeightedQ = hub_WeightedQ,
      save.path = save.path,
      names = names
    ),
    Data = list(Hub=Hub,Expr=Data),
    Plot = NULL)
  on.exit(options(old), add = TRUE)
  print("All done!")
  return(result)
}














