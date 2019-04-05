


#' @title Survival K-M plot base on expression matrix and design object
#' @description Survplot.genedata help draw  KM plot base on expression matrix and design object
#' @param expr.matrix expression matrix
#' @param design design object
#' @param select the genes or markers you want to select
#' @param event.status the colname of the status of event.
#' @param event.timethe colname of the time of event.
#' @param event.lower the lower cutoff of event time.Default is 0.
#' @param mode the mode of cut-off strategy for continuous value like gene expression data.Default is \code{mode = "maxstat"} provided by \code{\link[maxstat]{maxstat.test}}.If you want median value as cut-off,you can just set \code{mode = "median"}
#' @param pval.position the postion of p value in the plot
#' @param size the size of saved plot
#' @param save.file whether to save pdf for plot.If you select more than 1 genes,it would set \code{save.file = T} automatically.
#' @param show.music whether to show music when job was completed.Only available when you select more than 1 genes.
#' @param width the width of saved plot
#' @param height the height of saved plot
#' @param part name of saved pdf file
#' @details There are some different between the output of multi- and uni- variate plot.
#' @seealso \code{\link[maxstat]{maxstat.test}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky);
#' data("rna.tpm",package = "lucky")
#' data("rna.design",package = "lucky")
#'
#' ## Muti-variable Plot.here I focus on patients with more follow-up time than 89 days with event.lower = 89
#' p <- Survplot_genedata(expr.matrix = log2(rna.tpm + 1),
#'                        design = rna.design,
#'                        select = c("ENSG00000004478",
#'                                   "ENSG00000000457"),
#'                        event.status = "OS.status",
#'                        event.time = "OS.time",
#'                        event.lower = 89,
#'                        mode = c("median","maxstat")[2],
#'                        pval.position = c(900,1),# x=900,y=1
#'                        size = 12,
#'                        save.file = F,
#'                        show.music = F,
#'                        width = 12,height = 12,
#'                        names = "love")
#'
#' ## univariable Plot
#' p <- Survplot_genedata(expr.matrix = log2(rna.tpm + 1),
#'                        design = rna.design,
#'                        select = "ENSG00000004478",
#'                        event.status = "OS.status",
#'                        event.time = "OS.time",
#'                        event.lower = 89,
#'                        mode = c("median","maxstat")[2],
#'                        pval.position = c(900,1))
#' @export
Survplot_genedata <- function(expr.matrix,
                              design,
                              select,
                              event.status,
                              event.time,
                              event.lower = 0,
                              mode = c("median","maxstat")[2],
                              pval.position =NULL,
                              size = 12,
                              save.file = F,
                              show.music = F,
                              width = 12,height = 12,
                              names = "love"){
  ## 是否是多个基因
  if(length(select) > 1){
    ## 多个基因
    p <- Survplot_genedata2(expr.matrix=expr.matrix[select,],
                            design=design,
                            event.status = event.status ,
                            event.time = event.time,
                            event.lower = event.lower,
                            mode = mode,
                            pval.position = pval.position,
                            size = size,
                            width = width,height = height,
                            names = names,
                            show.music = show.music)
  } else {
    p <- Survplot_genedata1(expr.matrix = expr.matrix,
                            design = design,
                            select = select,
                            event.status = event.status,
                            event.time = event.time,
                            event.lower = event.lower,
                            mode = mode,
                            pval.position = pval.position,
                            size =  size,
                            save.file = save.file,
                            width = width,height = height,
                            names = names,
                            print = T)
  }

  ## 输出结果
  return(p)
}

## univariate plot
Survplot_genedata1 <- function(expr.matrix,
                              design,
                              select,
                              event.status,
                              event.time,
                              event.lower,
                              mode = c("median","maxstat")[1],
                              pval.position =NULL,
                              size = 12,
                              save.file = T,
                              width = 12,height = 12,
                              names = "love",
                              print = T){
  ## 加载必要的包
  nd <- c("survival","survminer","ggplot2")
  Plus.library(nd)

  ## 对齐
  expr.matrix <- expr.matrix[,rownames(design)]

  ## 提取信息
  list.surv <- extra.surv(expr.matrix =expr.matrix,
                          design = design,
                          select = select,
                          event.status = event.status,
                          event.time = event.time,
                          event.lower = event.lower,
                          mode=mode)
  df1 <- list.surv$Metadata


  ## 生存曲线绘制
  if(is.null(pval.position)){pval.position <- c(max(df1$time),1)}
  p <- FastSurvplot(data = df1,
                    time = "time",status ="status",
                    marker = "gene",
                    title = NULL,
                    color = NULL,
                    legend.title = select,
                    legend = "top",
                    size = size,
                    pval = TRUE,
                    pval.position = pval.position)
  if(print == T) {print(p)}

  ## 保存图片
  if(save.file == T){
    ggsave(paste0(names,"_",select,"_KM.plot.pdf"),p,width = width,height = height)
    print("Plot has been saved in present work space!")

  }

  ## 输出结果
  df2 <- df1;colnames(df2)[3] <- select
  l <- list(
    Survplot = p,
    CutOff = list.surv$CutOff,
    Metadata = df2
  )
  return(l)
}


## Mutiple plot
Survplot_genedata2 <- function(expr.matrix,
                               design,
                               event.status,
                               event.time,
                               event.lower,
                               mode = c("median","maxstat")[1],
                               pval.position =NULL,
                               size = 12,
                               width = 12,height = 12,
                               names = "love",
                               show.music = F){
  ## 对齐
  expr.matrix <- expr.matrix[,rownames(design)]
  select <- rownames(expr.matrix)

  ## 画图
  pdf(paste0(names,"_multiple KMCurves.pdf"),width = width,height = height)
  m <- NULL;cutoff <- NULL
  for(i in 1:length(select)){
    select.i <- select[i]
    List.Surv <- Survplot_genedata1(expr.matrix = expr.matrix,
                                   design = design,
                                   select = select.i,
                                   event.status = event.status,
                                   event.time = event.time,
                                   event.lower = event.lower,
                                   size = size,
                                   mode = mode,
                                   pval.position =pval.position,
                                   save.file = F,
                                   print = F)
    p1 <- List.Surv$Survplot + ggtitle("")
    print(List.Surv$Survplot)
    ## 保存数据
    m <- c(m,list(List.Surv$Metadata))
    cutoff <- c(cutoff,list(List.Surv$CutOff))
    names(m)[i] <- select.i;names(cutoff)[i] <- select.i;
    print(paste0(i," : The Survplot of ",select.i," has been printed!"))
  }
  dev.off()

  ## 结束
  if(show.music == T){tuneR::play(music)}
  print("All done!")
  return(list(Metadata = m,
              CutOff = cutoff))
}

