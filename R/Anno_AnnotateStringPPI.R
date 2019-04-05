

#' @title Re-annotation for ppi network from STRING web
#' @description  Re-annotation for ppi network from STRING web
#' @param ppi a data frame from "... as simple tabular text output" in string PPI network analysis
#' @param anno.abnormal a list with wrong annotation name and right annotation elements
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' anno.abnormal <- list(
#'"ENSG00000248993" = "XXbac-BPG181M17.5",
#' "PTPN20A" = "PTPN20",
#' "LRRC33" = "NRROS",
#' "C10orf54"  = "VSIR",
#' "EMR1" = "ADGRE1",
#' "1-Mar"= "MARCH1",
#' "FAIM3"="FCMR",
#' "ADORA3"= "TMIGD3"
#' )
#' ppi2 <- AnnotateStringPPI(ppi,anno.abnormal)
#' @export
AnnotateStringPPI <- function(ppi,anno.abnormal){

  ## 提取蛋白id
  ppi_1 <- ppi
  for(i in 1:length(anno.abnormal)){ # i=1
    wrong <- names(anno.abnormal)[i]
    right <- anno.abnormal[[i]]
    ppi_1$node1[as.character(ppi_1$node1) %in% wrong] <- right
    ppi_1$node2[as.character(ppi_1$node2) %in% wrong] <- right
  }
  return(ppi_1)
}



