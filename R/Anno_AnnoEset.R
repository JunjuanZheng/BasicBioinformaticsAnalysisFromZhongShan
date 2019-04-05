

#' @title Annotate EsetList via bioconductor
#' @description  Annotate EsetList via bioconductor
#' @param eset a EsetList object
#' @param probeid the colnames of probe id in fData.Default is "ID".
#' @importFrom Biobase fData protocolData
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr left_join
#' @importFrom readxl read_xlsx
#' @return an annotated EsetList object(fData)
#' @seealso \code{\link[clusterProfiler]{bitr}};\code{\link[dplyr]{left_join}};\code{\link[readxl]{read_xlsx}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' eset2 <- AnnoEset(eset)
#' View(Biobase::fData(eset2))
#' @export
AnnoEset <- function(eset,probeid="ID"){
  ## 产生储存文件夹
  old <- options()

  ## get GPL number
  GPL.NO <- Biobase::protocolData(eset)@data
  GPL.NO <- as.character(GPL.NO$platform_id[1])
  LuckyVerbose("Get GPL platform: ",GPL.NO)

  ## get package names of GPL platform
  LuckyVerbose("Get package names of GPL platform...")
  path_1 <- system.file("extdata", "GPLData.xlsx", package = "lucky")
  #path_1 <- "E:/RCloud/RFactory/lucky/inst/extdata/GPLData.xlsx"
  db1 <- readxl::read_xlsx(path_1)
  db1 <- as.data.frame(db1,stringsAsFactors = F)
  package.name <- convert(GPL.NO,"gpl","bioc_package",db1)
  if(is.na(package.name)){
    LuckyVerbose("Not available GPL annotation in Bioconductor →_→ ",levels= 2)
    LuckyVerbose("END!")
  } else {
    LuckyVerbose("Available GPL platform...",levels= 2)
    package.name <- paste0(package.name,".db")
    Plus.library(package.name)

    ## Annotate the fData of eset object
    LuckyVerbose("Annotate fData of eset object...")
    f <- Biobase::fData(eset)
    columns(hgu133plus2.db)
    probe <- as.character(f[,probeid])
    anno1 <- clusterProfiler::bitr(
      probe,
      fromType = "PROBEID",
      toType = c("SYMBOL","ENTREZID","ENSEMBL"),
      OrgDb = package.name
    )
    anno2 <- anno1
    colnames(anno2) <- c(probeid,"GENE_SYMBOL","ENTREZ_ID","ENSEMBL")

    ## merge repeated records
    LuckyVerbose("Merge by probe id...")
    getMerge <- function(x){
      x2 <- paste(unique(as.character(x)),collapse = " // ")
      return(x2)
    }
    geneidfactor <- anno2$ID
    anno3 <- apply(anno2,2,function(x)tapply(x,geneidfactor,getMerge))
    anno3 <- as.data.frame(anno3,stringAsFactor = F)
    options(warn = -1)
    anno4 <- dplyr::left_join(f,anno3,by = probeid)
    rownames(anno4) <- as.character(anno4[,probeid])
    options(warn = old$warn)

    ## replace old fData
    fData(eset) <- anno4
    LuckyVerbose("All done!")
    return(eset)
  }

}

















