

#' @keywords ModulePreservation
#' @title a easy method for Module Preservation in FastWGCNA pipeline
#' @description If you want to test the preservation of modules created by
#'   FastWGCNA and you have also another datasets as validataion
#'   cohorts,\code{ModulePreservation} would be a nice function for this
#'   mission.
#' @param object a result from \code{\link{FastWGCNA}}
#' @param valid.matrix a list of expression matrix with genes rows and sample
#'   cols.
#' @param setLabels the names of expression matrix in \code{valid.matrix}.
#' @param nPermutations specifies the number of permutations that will be
#'   calculated in the permutation test.
#' @param randomSeed seed for the random number generator. If NULL, the seed
#'   will not be set. If non-NULL and the random generator has been initialized
#'   prior to the function call, the latter's state is saved and restored upon
#'   exit
#' @param quickCor number between 0 and 1 specifying the handling of missing
#'   data in calculation of correlation. Zero means exact but potentially slower
#'   calculations; one means potentially faster calculations, but with
#'   potentially inaccurate results if the proportion of missing data is large.
#'   See cor for more details.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#'   make the output progressively more and more verbose.
#' @param cutoff.Z Default is 10.10 is efficient enought and should not be
#'   changed.
#' @inheritParams FastWGCNA
#' @importFrom WGCNA labels2colors modulePreservation
#' @return LuckyWGCNA object
#' @details I.The expression matrix in \code{object} would be used as a training
#'   cohorts in module preservation evaluation.       II.\code{nPermutations = 200} is a usual method but it would be time-consuming.It seems that \code{nPermutations = 50} would be stable enough for module presevation exploration.
#' @seealso \code{\link{FastWGCNA}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' object = wgcna;rm(wgcna);gc()
#' valid.matrix = list(Valid1 = rna.fpkm.tumor,
#'                     Valid2 = rna.fpkm.tumor)
#' result_MP <- ModulePreservation(object,
#'                                 valid.matrix,
#'                                 nPermutations = 20,
#'                                 cutoff.Z = 10,
#'                                 report = T,
#'                                 width = 15,height = 7.5,
#'                                 save.path = "WGCNA-test",
#'                                 names = "love")
#' @export
ModulePreservation <- function(object,
                               valid.matrix = list(),
                               nPermutations = 50,
                               randomSeed = 1,
                               quickCor = 0,
                               verbose = 3,
                               cutoff.Z = 10,
                               report = T,
                               width = 15,height = 7.5,
                               save.path = "WGCNA",
                               names = "love"){

  ###======================Preparation=====================###
  ## 包
  #nd <- c("WGCNA")
  #Plus.library(nd)

  ## 产生储存文件夹
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = 0)

  ###=====================data processing====================###

  # 将train和test合成list的形式
  setLabels = c("Train", names(valid.matrix))

  ## select genes
  genes <- colnames(object$dataExpr)

  ## Valid matrix transformation
  vm_2 <- NULL
  for(i in 1:length(valid.matrix)){ # i = 1
    m.i <- valid.matrix[[i]]
    m.i <- m.i[genes,]
    m.i_2 <- as.data.frame(t(m.i))
    l.i <- list(data = m.i_2)
    vm_2 <- c(vm_2,list(l.i));names(vm_2)[i] <- names(valid.matrix)[i]
  }
  # rm(valid.matrix,envir = environment());gc()

  ## moduleColors
  moduleLabels <-  object$net$colors
  moduleColors <-  labels2colors(moduleLabels)

  ## create matrix list
  train <- list(data = object$dataExpr)
  multiExpr <-  c(list(train),vm_2);names(multiExpr)[1] <- "Train"
  multiColor <- list(Train = moduleColors)


  ###===============modulePreservation process==============###
  # Z大于10，代表strong preserved，好的module
  # 大于2小于10代表weak preserved
  # 小于2代表not preserved，不好的module
  mp <-  modulePreservation(multiExpr,
                            multiColor,
                            referenceNetworks = 1,
                            nPermutations = nPermutations,
                            randomSeed = randomSeed,
                            quickCor = quickCor,
                            verbose = verbose)
  save(mp,file = paste0(dir,names,"_module Preservation.rda"))
  cat("The module Preservation analysis had been save.","\n")

  ###================test every dataset=====================###
  ref = 1
  for(i in 2:length(multiExpr)){
    test <-  i
    names.i <- names(multiExpr)[i]
    statsObs <-  cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
    statsZ <-  cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

    ## save Z summary data
    data.z <- cbind(statsObs[,c("medianRank.pres", "medianRank.qual")],
                    signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
    write.csv(data.z,paste0(dir,names,"_",names.i,"_module Preservation result.csv"))

    # 划分训练集后有多少module，每个module的大小
    modColors <-  rownames(mp$preservation$observed[[ref]][[test]])
    moduleSizes <-  mp$preservation$Z[[ref]][[test]][, 1]

    ## 去除Z summary值太小的Modules
    module.out <- row.names(statsZ[statsZ$Zsummary.pres < cutoff.Z,])
    plotMods <-  !(modColors %in% module.out)
    text <-  modColors[plotMods]
    plotData <-  cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])

    ###==================vitualization=======================###
    mains = c("Preservation Median rank", "Preservation Zsummary")
    pdf(paste0(dir,names,"_",names.i,"_module Preservation plot.pdf"),15,7.5)
    par(mfrow = c(1,2));par(mar = c(4.5,4.5,2.5,1))
    for (p in 1:2){
      min = min(plotData[, p], na.rm = TRUE);
      max = max(plotData[, p], na.rm = TRUE);
      # Adjust ploting ranges appropriately
      if (p==2){
        if (min > -max/10) min = -max/10
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
      } else
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
      #bg 颜色 pch 圆圈的种类
      plot(moduleSizes[plotMods], plotData[plotMods, p],
           col = 1, bg = modColors[plotMods],
           pch = 21,
           main = mains[p],
           cex = 2.4,##圆圈的大小
           ylab = mains[p], xlab = "Module size", log = "x",
           ylim = ylim,
           xlim = c(10, 2000),
           cex.lab = 1.2,
           cex.axis = 1.2,
           cex.main =1.4)
      ##贴上标签
      labelPoints(moduleSizes[plotMods],
                  plotData[plotMods, p],
                  text,
                  cex = 1,
                  offs = 0.08);
      # For Zsummary, add threshold lines
      if (p==2){
        abline(h=0)
        abline(h=2, col = "blue", lty = 2)
        abline(h=cutoff.Z, col = "darkgreen", lty = 2)
      }
    }
    dev.off()

    ## Plot report
    if(report == T){
      win.graph(15,7.5)
      par(mfrow = c(1,2));par(mar = c(4.5,4.5,2.5,1))
      for (p in 1:2){
        min = min(plotData[, p], na.rm = TRUE);
        max = max(plotData[, p], na.rm = TRUE);
        # Adjust ploting ranges appropriately
        if (p==2){
          if (min > -max/10) min = -max/10
          ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
        } else
          ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
        #bg 颜色 pch 圆圈的种类
        plot(moduleSizes[plotMods], plotData[plotMods, p],
             col = 1, bg = modColors[plotMods],
             pch = 21,
             main = mains[p],
             cex = 2.4,##圆圈的大小
             ylab = mains[p], xlab = "Module size", log = "x",
             ylim = ylim,
             xlim = c(10, 2000),
             cex.lab = 1.2,
             cex.axis = 1.2,
             cex.main =1.4)
        ##贴上标签
        labelPoints(moduleSizes[plotMods],
                    plotData[plotMods, p],
                    text,
                    cex = 1,
                    offs = 0.08);
        # For Zsummary, add threshold lines
        if (p==2){
          abline(h=0)
          abline(h=2, col = "blue", lty = 2)
          abline(h=cutoff.Z, col = "darkgreen", lty = 2)
        }
      }
    }
    par(mfrow = c(1,1))

    cat("Module Preservation plot of",names.i,"had been save...","\n")
  }

  ###====================Result Output======================###
  result <- list(
    Repeat = list(
      nPermutations = nPermutations,
      randomSeed = randomSeed,
      quickCor = quickCor,
      verbose = verbose,
      cutoff.Z = cutoff.Z,
      width = width,height = height,
      save.path = save.path,
      names = names
    ),
    Data = list(mp = mp),
    Plot = NULL
  )
  cat("All done!","\n")
  return(result)
}





