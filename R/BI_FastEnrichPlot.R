


#' @title Fast way to draw enrich plot after clusterProfilter enrichfunction
#' @description Fast way to draw enrich plot after clusterProfilter enrichfunction
#' @param object result from \code{\link[clusterProfilter]{enrichGO}} or \code{\link[clusterProfilter]{enrichKEGG}} or \code{\link[clusterProfilter]{gseGO}} or \code{\link[clusterProfilter]{enrichKEGG}}
#' @param pvalueCutoff p value cutoff
#' @param qvalueCutoff q value cutoff
#' @param id.col ID you want to show in final plot
#' @param select if \code{select = NULL},the number of \code{topshow} would be shown.If \code{select != NULL},selected terms would be showed and thus \code{topshow} is not effecient.
#' @param x.title the title of x axis
#' @param y.title the title of y axis
#' @param size the size of whole plot
#' @param short.cutoff if the number of character of terms is larger than \code{short.cutoff},it would automatically change row in the showed plot.Defaut is 46.If you don't want change-row,just set \code{short.cutoff = Inf}
#' @param topshow the number of top iterms ordered by FDR you want to show.It works only \code{select = NULL}
#' @importFrom ggplot2 ggplot geom_point scale_fill_gradient coord_flip labs theme_bw theme
#' @return LuckyList Object
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' load("E:/iProjects/LM_WGCNA/data/LM_WGCNA_02A/rda/goList.rda")
#' object <- goList$GOEnrichment$BP
#' res_EP <- FastEnrichPlot(object = object,
#'                          pvalueCutoff = NULL,
#'                          qvalueCutoff = NULL,
#'                          id.col = "Description",
#'                          select = NULL,
#'                          x.title = "Biological Process terms",
#'                          y.title = "-log10(FDR)",
#'                          size = 15,
#'                          short.cutoff = 46,
#'                          topshow = 10)
#' win.graph(12,9);res_EP$Plot
#'
#' pvalueCutoff = 0.05
#' qvalueCutoff = 0.05
#' select = res2$ID[1:10]
#' id.col = "Description"
#' x.title = "Biological Process terms"
#' y.title = "-log10(FDR)"
#' @export
FastEnrichPlot <- function(object,
                           pvalueCutoff = NULL,
                           qvalueCutoff = NULL,
                           id.col = "Description",
                           select = NULL,
                           x.title = "Biological Process terms",
                           y.title = "-log10(FDR)",
                           size = 15,
                           short.cutoff = 46,
                           topshow = 10){
  ## Package
  nd <- c("ggplot2")
  Plus.library(nd)

  ## get result from object
  res <- object@result

  ## get significant data
  if(is.null(pvalueCutoff) & is.null(qvalueCutoff)){
    pvalueCutoff = object@pvalueCutoff
    qvalueCutoff = object@qvalueCutoff
  }
  colnames(res)
  res1 <- subset(res,pvalue < pvalueCutoff & qvalue < qvalueCutoff)
  if(!is.null(select)){
    if(all(is.numeric(select))){
      res2 <- res1[select,]
    } else {
      test.res <- as.matrix(res1)
      a <- test.res %in% select
      a1 <- matrix(a,nrow = nrow(test.res),byrow = F)
      lg1.a1 <- apply(a1,1,is.one.true)
      res2 <- res1[lg1.a1,]
    }
  } else {
    res2 <- res1[order(res1$p.adjust,decreasing = F),]
    res2 <- res2[1:topshow,]
  }

  ## data processing
  res3 <- res2
  res3$qvalue <- -log10(res3$qvalue)
  res3 <- res3[order(res3$qvalue,decreasing = F),]
  res4 <- subset(res3,select = c(id.col,"Count","qvalue"))
  colnames(res4)[1] <- "ID"

  ## deal with too-long names
  getShort <- function(vt){
    vt <- as.character(vt)
    x <- NULL
    for(i in 1:length(vt)){# i=1
      vt.i <- vt[i]
      if(nchar(vt.i) > short.cutoff){
        r <- floor(nchar(vt.i)/short.cutoff)+1 #分成3行
        vt.i2 <- Fastextra(vt.i," ")
        new.name <- cut_vector(vt.i2,nsplit = r)
        nn <- NULL
        for(j in 1:length(new.name)){ #j=1
          nn.j <- new.name[[j]]
          nn.j <- paste(nn.j,collapse = " ")
          nn.j <- paste0(nn.j," \n ")
          nn <- paste0(nn,nn.j)
        }
        nn2 <- substring(nn,1,(nchar(nn)-3)) #[1] "char"
      } else {
        nn2 <- vt.i
      }
      x <- c(x,nn2)
    }
    return(x)
  }
  res4$ID <- getShort(res4$ID)
  res4$ID <- factor(res4$ID,levels = res4$ID)
  res4$Count <- as.integer(res4$Count)

  ## ggplot style
  p <- ggplot(res4, aes(x=ID,y=qvalue,size=Count)) +
    geom_point(aes(colour = qvalue)) +
    scale_fill_gradient(low = mycolor[23],high = mycolor[21],aesthetics = "colour") +
    coord_flip() +
    labs(x = x.title,y = y.title,colour = "-log10(FDR)") +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold",size = size),
      axis.text = element_text(face = "bold",size = size/15*12),
      legend.title = element_text(face = "bold",size = size/15*12),
      legend.text =element_text(face = "bold",size = size/15*12)
    )

  ## 输出结果
  l <- list(
    Repeat = list(pvalueCutoff = pvalueCutoff,
                  qvalueCutoff = qvalueCutoff,
                  id.col = id.col,
                  topshow = topshow,
                  select = select,
                  x.title = x.title,
                  y.title = y.title,
                  short.cutoff = short.cutoff,
                  size = size),
    Data = res3,
    Plot = p
  )
  return(l)
}











