


#' @title DEGs screening via common t-test
#' @description DEGs screening via common t-test
#' @param object DesignList or expression matrix
#' @param design design object.If \code{object} is DesignList,it can be NULL
#' @inheritParams DESeq1
#' @param verbose whether report the process
#' @importFrom plyr adply
#' @importFrom Biobase exprs
#' @return data frame of DEGs information
#' @seealso \code{\link{DESeq1}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' res <- ScreenGenes(DesignEset,
#'                    contrast.level = c("tumor","normal"),
#'                    contrast.control = c("normal"),
#'                    cutoff.p=0.05)
#' @export
ScreenGenes <- function(object,
                        design,
                        contrast.col,
                        contrast.level = c("tumor","normal"),
                        contrast.control = c("normal"),
                        cutoff.p=0.05,
                        verbose = F){
  ## 判断是否为DesignEset对象。
  if(all(c("ESet","DesignList") %in% names(object))){
    if(verbose == T) LuckyVerbose("This is a DesignList object...")
    ## 构建contrast对象
    c1 <- paste(setdiff(contrast.level,contrast.control),contrast.control,sep= "-")
    expr <- Biobase::exprs(object$ESet)
    design <- object$DesignList$condition[[c1]]
    contrast.col <- colnames(design)[1]
  } else {
    expr <- object
    design <- design
    contrast.col <- contrast.col
  }
  #expr design contrast.col

  ## 对齐
  expr <- expr[,rownames(design)]

  ## 得到差异性基因
  if(verbose == T) LuckyVerbose("Get DEGs...")
  getSigGenes <- function(x){
    x1 <- as.numeric(as.character(x))
    l <- c(contrast.control,setdiff(contrast.level,contrast.control))
    f1 <- factor(design[,contrast.col],levels = l)
    df1 <- data.frame(
      group = f1,
      expression = x1,
      stringsAsFactors = F
    )
    control <- x1[f1 == l[1]]
    treat <- x1[f1 == l[2]]
    a <- t.test(x = control, y = treat, alternative = "two.sided")
    result <- data.frame(
      mean.control = as.numeric(a$estimate[1]),
      mean.treat = as.numeric(a$estimate[2]),
      logFC = log2(a$estimate[2]/a$estimate[1]),
      statistic = a$statistic,
      pvalue = a$p.value,
      stringsAsFactors = F
    )
    rownames(result) <- 1
    return(result)
  }
  result <- plyr::adply(expr,1,getSigGenes)
  rownames(result) <- as.character(result[,1])
  result <- result[-1]

  ## filter
  result2 <- result[result$pvalue < cutoff.p,]
  if(verbose == T) LuckyVerbose("All done!")
  return(result2)
}






