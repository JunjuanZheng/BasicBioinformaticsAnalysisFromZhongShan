


#' @keywords GetMEs
#' @title Get new module eigengenes matrix from a given dataset
#' @description Get New module eigengenes matrix from a given dataset
#' @param object the result of \code{\link{FastWGCNA}}
#' @inheritParams ModulePreservation
#' @param log.convert whether to do log scale for expression matrix
#' @importFrom WGCNA labels2colors moduleEigengenes
#' @details Only support wgcna object with 1 Block.
#' @return LuckyWGCNA Object
#' @seealso \code{\link{FastWGCNA}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' expr.matrix = rna.fpkm.tumor
#' design = rna.design.tumor
#' log.convert =T
#' newge <- GetMEs(object,expr.matrix)
#' View(newge$Data$MEs)
#' @export
GetMEs <- function(object,
                   valid.matrix,
                   log.convert = T){
  ## 自定义函数
  GetMEs_1 <- function(object,
                       expr.matrix,
                       log.convert = T){
    ## package
    # nd <- c("WGCNA");Plus.library(nd)

    ## get training module information
    moduleLabels <-  object$net$colors
    moduleColors <-  labels2colors(moduleLabels)
    genes <- colnames(object$dataExpr)
    names(genes) <- moduleColors

    ## get new dataExpr
    dataExpr <- t(expr.matrix[genes,])
    if(log.convert == T){dataExpr <- log2(dataExpr + 1)}

    ## create new module Eigengenes
    ME_cols <- moduleEigengenes(expr = dataExpr,
                                colors = moduleColors)
    ME_cols <- ME_cols$eigengenes
    ME_cols <- orderMEs(ME_cols)
    rownames(ME_cols) <- rownames(dataExpr)

    ## 输出结果
    result <- list(
      Repeat = list(log.convert = log.convert),
      Data = list(MEs = ME_cols,
                  dataExpr = dataExpr,
                  genes = genes),
      Plot = NULL
    )
    return(result)
  }

  ## 基本信息
  moduleLabels <-  object$net$colors
  moduleColors <-  labels2colors(moduleLabels)
  genes <- colnames(object$dataExpr)
  names(genes) <- moduleColors

  ## Valid.matrix
  valid <- NULL
  for(i in 1:length(valid.matrix)){ #i =1
    expr.matrix <- valid.matrix[[i]]
    r.i <- GetMEs_1(object = object,
                    expr.matrix = expr.matrix,
                    log.convert = log.convert)
    ME.i <- r.i$Data$MEs
    dataExpr.i <- r.i$Data$dataExpr
    valid.i <- list(
      MEs = ME.i,
      dataExpr = dataExpr.i
    )
    valid <- c(valid,list(valid.i))
    names(valid)[i] <- names(valid.matrix)[i]
  }

  ## output
  result <- list(
    Repeat = list(log.convert = log.convert),
    Data = valid,
    Plot = NULL
  )
  return(result)

}













