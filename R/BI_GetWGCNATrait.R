

#' @title create trait matrix for WGCNA
#' @description \code{GetWGCNATrait} creates trait matrix for WGCNA.If there are some factor variables in your design object,they must be converted into numeric vector before go into \code{\link{ModuleTrait}} function,when \code{GetWGCNATrait} would complete this job.
#' @keywords  GetWGCNATrait
#' @param design design object
#' @param convert.list a list of convert information.See examples.
#' @return a new design object
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#'library(lucky)
#'data(rna.design.tumor)
#' stad.convert.list <- {list(
#'condition = list("tumor" = 1,"normal" = 0),
#'Mol.subtype = list("EMT" = 3,"MSI" = 2,"MSS/TP53-" = 0,"MSS/TP53+" = 1),
#'Mol.subtype2 = list("MSI-H" = 2,"MSI-L"=1,"MSS" = 0),
#'TCGA.subtype = list("CIN"=1,"EBV"=2,"GS"=3,"MSI"=4),
#'gender = list("FEMALE"=0,"MALE"=1),
#'his1 = list("Intestinal Adenocarcinoma"=0,"Adenocarcinoma" =1,"Signet Ring Type"=2),
#'his2 = list("Not Otherwise Specified (NOS)"=0,"Papillary Type"=1,"Tubular Type"=2,"Mucinous Type"=3,"Diffuse Type"=4),
#'his.grade = list("GX"=0,"G1"=1, "G2"=2, "G3"=3),
#'pStage = list("Stage I"=1,"Stage IA"=1,"Stage IB"=1,
#'              "Stage II"=2,  "Stage IIA"=2,  "Stage IIB" =2,
#'              "Stage III" =3,"Stage IIIA"=3, "Stage IIIB"=3,
#'              "Stage IIIC"=3, "Stage IV"=4),
#'pT = list("T1"=1, "T2"=2, "T3"=3, "T4"=4),
#'T.status = list("T_early"=0,"T_later"=1),
#'N.status = list("N0"=0, "Np"=1),
#'M.status  = list("M0"=0, "M1" =1,"MX"=NA)
#')}
#' stad.Trait <- GetWGCNATrait(design,stad.convert.list)
#' @export
GetWGCNATrait <- function(design,convert.list){

  ## 对某一列表中的
  # c.i <- convert.list[[1]]
  # var.i <- names(convert.list)[1]
  # get1(design,c.i,var.i)
  get1 <- function(design,c.i,var.i){
    variable.i <- names(c.i)
    a_1 <- as.character(design[,var.i])
    for(i in 1:length(variable.i)){ # i=2
      v.i <- variable.i[i] #某变量
      a.i <- c.i[[i]]
      a_1[a_1 %in% v.i] <- a.i
    }
    a_1 <- as.numeric(a_1)
    return(a_1)
  }

  ## variables
  design_1 <- design
  for(i in 1:length(convert.list)){ # i=1
    var.i <- names(convert.list)[i]
    c.i <- convert.list[[i]]
    design_1[,var.i] <- get1(design,c.i,var.i)
  }
  return(design_1)
}


# save(stad.Trait,file = "stad.Trait.rda")






