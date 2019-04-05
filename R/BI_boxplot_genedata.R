



#' @title Fast way to draw boxplot for gene expression
#' @description boxplot_genedata provide fast way to draw boxplot for gene expression based on gene expression matrix and design object
#' @param select.genes the genes you want to plot
#' @param design the design object(data of clinical features)
#' @param expr.matrix expression matrix of genes
#' @param convert Whether to convert Ensembl id as symbols.Default is \code{convert = T}.Note that if you use common names as genes(like miRNAs names),please set \code{convert = F}.
#' @param type the style of boxplot.One of "facet" and "normal".
#' @param contrast the colnames of contrast value.Like "N.status".
#' @param contrast.list the components of contrast.Like c("N0","Np")
#' @param genes.levels the level of genes.Default is NULL,and it's recommanded that you set a self-defined levels for better visualization.
#' @param palette the color of genes ploted.Default is NULL.
#' @param point.alpha the color density of point
#' @param box.alpha the color density of box
#' @param method the method of P value for comparision.Default is "wilcox.test".See also \code{method} of \code{\link[ggpubr]{compare_means}} function.
#' @param x.title the names of x axis of plot
#' @param y.title the names of y axis of plot
#' @param cut whether to cut plot.If the number of genes is very much,it is recommanded to set \code{cut = T} for better visualization.
#' @param ncut the number of cut strategy.Only available when \code{cut = T}.Default is 2.
#' @param width the width of plot
#' @param height the height of plot
#' @param names part name of the saved plot
#' @details the names of genes shold be Ensembl id.
#' @seealso \code{\link[ggpubr]{compare_means}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' data("rna.tpm",package = "lucky")
#' data("rna.design",package = "lucky")
#' p <- boxplot_genedata(select.genes = c("ENSG00000004478",
#'                                        "ENSG00000000457"),
#'                       design = rna.design,
#'                       expr.matrix = log2(rna.tpm + 1),
#'                       contrast = "condition",
#'                       contrast.list = c("normal","tumor"),
#'                       genes.levels = NULL,
#'                       palette=NULL,
#'                       point.alpha = 0.5,
#'                       box.alpha = 1,
#'                       method = "wilcox.test",
#'                       x.title="Genes",
#'                       y.title="The Expression of Genes(log TPM)",
#'                       cut = F,
#'                       width = 12,height = 10,
#'                       names = "love")
#' ## See plot of genes expression
#' print(p$plot)
#'
#' ## See meta data of genes expression
#' View(p$metadata)
#' @export
boxplot_genedata <- function(select.genes,
                             design,
                             expr.matrix,
                             convert = T,
                             type = c("facet","normal")[1],
                             contrast = "N.status",
                             contrast.list = c("N0","Np"),
                             genes.levels = NULL,
                             palette=NULL,
                             point.alpha = 0.5,box.alpha = 1,
                             method = "wilcox.test",
                             x.title=NULL,
                             y.title=NULL,
                             legend.position=c(0.9,0.9),
                             cut = F,ncut = 2,
                             width = 12,height = 10,
                             names = "love"){
  ## 原始函数
  boxplot.x <- function(select.genes,
                        design,
                        expr.matrix,
                        convert = T,
                        type = c("facet","normal")[1],
                        contrast = "N.status",
                        contrast.list = c("N0","Np"),
                        genes.levels = NULL,
                        palette=NULL,
                        point.alpha = 0.5,box.alpha = 1,
                        method = "wilcox.test",
                        x.title=NULL,
                        y.title=NULL,
                        legend.position="right",
                        print = T){
    ## 加载包
    nd <- c("ggplot2","ggpubr")
    Plus.library(nd)

    ## 生成boxplot.data
    expr.matrix <- expr.matrix[,rownames(design)]
    x <- data.frame(expr=NULL,genes =NULL);
    #i=select.genes[1]
    for(i in select.genes){
      x.i <- as.numeric(as.character(expr.matrix[i,]))
      x.i <- as.data.frame(x.i)
      colnames(x.i) <- "expr"
      x.i$genes <- i
      if(convert == T){x.i$symbols <- convert(i)} else {
        x.i$symbols <- i
      }#是否转换id
      p <- match(contrast,colnames(design))
      x.i$col1 <- design[,p]
      colnames(x.i)[match("col1",colnames(x.i))] <- contrast
      x <- rbind(x,x.i)
    }

    ## 提取和contrast.list有关的数据
    logic1 <- x[,contrast] %in% contrast.list
    x <- x[logic1,]

    ## levels
    x[,contrast] <- factor(x[,contrast],levels = contrast.list)
    if(convert == F){
      #不转换id
      x$symbols <- factor(x$symbols,levels = genes.levels)
    } else {
      if(!is.null(genes.levels)){
        if(all(genes.levels %in% select.genes)) {
          #ensembl类
          x$symbols <- factor(x$symbols,levels = convert(genes.levels))
        } else {
          #symbol类
          x$symbols <- factor(x$symbols,levels = genes.levels)
        }
      }
    }

    #method
    if(is.null(method)){
      #不进行比较
      compare <- NULL
    } else {
      compare <- stat_compare_means(comparisons = list(contrast.list),label = "p.signif",method = method)
    }

    ## 颜色
    if(!is.null(palette)){
      palette1 <- palette
    } else {
      palette = rep(mycolor,100)
      palette1 = palette[1:length(select.genes)]
    }

    ## 坐标轴标题
    if(is.null(x.title)){x.t = contrast} else {x.t = x.title}
    if(is.null(y.title)){y.t = "The Expression Level of Genes"} else {y.t = y.title}

    # 画图
    if(type == "facet"){
      p <- ggboxplot(x,
                     x = contrast, y = "expr",
                     fill = "symbols",#箱体填充颜色
                     color = "black",
                     palette = palette1,#按JCO杂志风格来配色
                     add = "jitter",
                     alpha = box.alpha,
                     add.params = list(alpha = point.alpha)) + #加上jitter点
        #scale_fill_manual(values = c(brewer.pal(12, "Set3"),brewer.pal(3, "Dark2"))) +
        labs(x = x.t,y = y.t) +
        theme(axis.text.x = element_text(size = 12,colour = "black",face = "bold")) + #调整X轴变量的属性。titile是图上方，caption在图下方。
        compare +
        facet_grid(. ~ symbols) +
        theme(axis.title.x = element_text(size = 15,colour = "black",face = "bold"),
              axis.title.y = element_text(size = 15,colour = "black",face = "bold"),
              legend.position='none',
              strip.background = element_rect(fill="white")) + #调整x轴和y轴title的属性。
        rotate_x_text(angle = 45) #x轴的text以逆时针旋转45度角
    } else {
      if(type == "normal"){
        ##  普通画法
        x1 <- x
        colnames(x1)[match(contrast,colnames(x1))] <- "condition"
        p <- ggplot(data = x1,aes(x = symbols,y = expr,fill = condition)) +
          geom_boxplot() +
          scale_fill_manual(values = palette1) +
          labs(x = x.t,y = y.t,fill = contrast) +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 12,colour = "black",face = "bold"),
            axis.text.y = element_text(size = 12,colour = "black",face = "bold"),
            axis.title.x = element_text(size = 15,colour = "black",face = "bold"),
            axis.title.y = element_text(size = 15,colour = "black",face = "bold"),
            legend.position = legend.position,
            strip.background = element_rect(fill="white")
          ) +
        rotate_x_text(angle = 45)
      } else {
        cat('Wrong type.It should be one of "facet" or "normal"')
      }
    }
    if(print == T){print(p)}

    # 输出数据框和图片
    l <- list(
      plot = p,
      metadata = x
    )
    return(l)
  }

  ## cut的应用
  if(cut == F){
    a <- boxplot.x(select.genes = select.genes,
                  design = design,
                  expr.matrix = expr.matrix,
                  convert = convert,
                  type = type,
                  contrast = contrast,
                  contrast.list = contrast.list,
                  genes.levels = genes.levels,
                  palette=palette,
                  point.alpha = point.alpha,
                  box.alpha = box.alpha,
                  method = method,
                  x.title=x.title,
                  y.title=y.title,
                  print=T)
    return(a)
  } else {
    print("Use cut.vector strategy...")
    ## 分割为多个图
    l1 <- cut_vector(1:length(select.genes),ncut)
    pdf(paste0(names,"_boxplot_multiple cut.pdf"),width,height)
    for(i in 1:length(l1)){ #i =1
      s.i <- l1[[i]];s.i <- select.genes[s.i]
      a <- boxplot.x(select.genes = s.i,
                     design = design,
                     expr.matrix = expr.matrix,
                     convert = convert,
                     type = type,
                     contrast = contrast,
                     contrast.list = contrast.list,
                     genes.levels = s.i,
                     palette=palette,
                     point.alpha = point.alpha,
                     box.alpha = box.alpha,
                     method = method,
                     x.title=x.title,
                     y.title=y.title,
                     print=F)
     print(a$plot)
    }
    dev.off()
    print("The boxplot had been saved in present work space.")
  }
}






