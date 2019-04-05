

#' @keywords ModuleTraitValid
#' @title Get module-trait relationship given new validation cohorts
#' @description Get module-trait relationship given new validation cohorts. \code{ModuleTraitValid} supports multiple expression matrix.
#' @inheritParams GetMEs
#' @inheritParams ModuleTrait
#' @param valid.design a named list of design object
#' @param valid.variable a named list of trait vector
#' @return LuckyWGCNA object
#' @seealso \code{\link{ModuleTrait}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' load("E:/RCloud/RFactory/lucky/love/WGCNA-test/love_wgcna.rda")
#' object = wgcna;rm(wgcna);gc()
#' log.convert=T
#' valid.matrix = list(Valid1 = rna.fpkm.tumor,
#'                     Valid2 = rna.fpkm.tumor)
#' valid.design = list(Valid1 = rna.design.tumor,
#'                     Valid2 = rna.design.tumor)
#' variable = c("age","his1","gender","N.status","T.status")
#' valid.variable = list(Valid1 = variable,
#'                       Valid2 = variable)
#' corType = c("pearson","bicor")[1]
#' report=F
#' height = 10;width = 8
#' save.path = "WGCNA"
#' names = "love"
#' result_MTV <- ModuleTraitValid(object,
#'                                valid.matrix,log.convert=T,
#'                                valid.design,
#'                                valid.variable,
#'                                corType = c("pearson","bicor")[1],
#'                                height = 10,width = 8,
#'                                save.path = "WGCNA",
#'                                names = "love")
#' @export
ModuleTraitValid <- function(object,
                             valid.matrix,log.convert=T,
                             valid.design,
                             valid.variable,
                             corType = c("pearson","bicor")[1],
                             report=T,
                             height = 10,width = 8,
                             save.path = "WGCNA",
                             names = "love"){
  ## load needed package
  #nd <- c("plyr","WGCNA")
  #Plus.library(nd)

  ## get multiple MEs object
  if(report == T) cat("Get MEs from valiadation cohort...","\n")
  result_GMEs <- GetMEs(object = object,
                        valid.matrix = valid.matrix,
                        log.convert = log.convert)

  ## ModuleTrait pipeline for every validation cohort
  valid <- NULL
  for(i in 1:length(valid.design)){ # i=1
    name.i <- names(valid.design)[i]
    ME.i <- result_GMEs$Data[[i]]$MEs
    design.i <- valid.design[[i]]
    variable <- valid.variable[[i]]
    result_MT.i <- ModuleTrait(object = ME.i,
                               design = design.i,
                               variable = variable,
                               corType = corType,
                               report = F,
                               height = height,
                               width = width,
                               save.path = paste0(save.path,"/",name.i,"/ModuleTraitValid"),
                               names = names)
    valid <- c(valid,list(result_MT.i))
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
      valid.variable = valid.variable,
      corType = corType,
      height = height,
      width = width,
      save.path = save.path,
      names = names
    ),
    Data = Data,
    Plot = NULL
  )
  if(report == T) cat("All done!","\n")
  return(result)

}










