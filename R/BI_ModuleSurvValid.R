

#' @keywords ModuleSurvValid
#' @title Get module-Survival relationship given new validation cohorts
#' @description Get module-Survival relationship given new validation cohorts
#' @inheritParams GetMEs
#' @inheritParams ModuleSurv
#' @inheritParams ModuleTraitValid
#' @param valid.time a named list of time colnames
#' @param valid.status a named list of status colnames
#' @return LuckyWGCNA object
#' @seealso \code{\link{ModuleSurv}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' load("E:/RCloud/RFactory/lucky/love/WGCNA-test/love_wgcna.rda");object = wgcna;rm(wgcna);gc()
#'  log.convert=T
#' valid.matrix = list(Valid1 = rna.fpkm.tumor,
#'                     Valid2 = rna.fpkm.tumor)
#' valid.design = list(Valid1 = rna.design.tumor,
#'                     Valid2 = rna.design.tumor)
#' time = c("OS.time","DFS.time")
#' status = c("OS.status","DFS.status")
#' valid.time = list(Valid1 = time,
#'                   Valid2 = time)
#' valid.status = list(Valid1 = status,
#'                     Valid2 = status)
#' save.path = "WGCNA"
#' names = "love"
#' digits = 3
#' result_MSV <- ModuleSurvValid(object,
#'                               valid.matrix,
#'                               log.convert=T,
#'                               valid.design,
#'                               valid.time,
#'                               valid.status,
#'                               digits = 3,
#'                               save.path = "WGCNA",
#'                               names = "love")
#' @export
ModuleSurvValid <- function(object,
                            valid.matrix,log.convert=T,
                            valid.design,
                            valid.time,
                            valid.status,
                            digits = 3,
                            save.path = "WGCNA",
                            names = "love"){
  ## package
  #nd <-c("survival","WGCNA","broom","plyr")
  #Plus.library(nd)

  ## get multiple MEs object
  result_GMEs <- GetMEs(object = object,
                        valid.matrix = valid.matrix,
                        log.convert = log.convert)

  ## ModuleSurv pipeline for every validation cohort
  valid <- NULL
  for(i in 1:length(valid.design)){ # i=1
    name.i <- names(valid.design)[i]
    ME.i <- result_GMEs$Data[[i]]$MEs
    design.i <- valid.design[[i]]
    time <- valid.time[[i]]
    status <- valid.status[[i]]
    result_MS.i <- ModuleSurv(object = ME.i,
                              design = design.i,
                              time = time,
                              status = status,
                              digits = digits,
                              save.path = paste0(save.path,"/",name.i,"/ModuleSurvValid"),
                              names = names)
    valid <- c(valid,list(result_MS.i))
    names(valid)[i] <- name.i
  }

  ## Data
  Data <- NULL
  for(i in 1:length(valid)){ # i=1
    valid.i <- valid[[i]]
    name.i <- names(valid)[i]
    Data.i <- valid.i$Data
    Data <- c(Data,list(Data.i))
    names(Data)[i] <- name.i
  }

  ## result output
  result <- list(
    Repeat = list(
      log.convert=log.convert,
      valid.design = valid.design,
      valid.time = valid.time,
      valid.status = valid.status,
      save.path = save.path,
      names = names
    ),
    Data = Data,
    Plot = NULL
  )
  return(result)
}



