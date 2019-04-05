

#' Visualize specified module Network after FastWGCNA pipeline
#'
#' @description Visualize specified module Network after FastWGCNA pipeline
#'
#' @param object the result from \code{\link{FastWGCNA}}
#' @param module module name.It must be a color oject like "blue"/"red"
#' @param hubGenes Genes you want to highlight.Default is NULL,which means no gene would be highlighted.
#' @param hubGenes.color the color of hubGenes.Default is mycolor[21].
#' @inheritParams FastWGCNA
#' @inheritParams Fastcornet2
#' @examples
#' ## NOT RUN
#' object = wgcna;rm(wgcna);
#' module = "blue"
#' hubGenes = c("ENSG00000009790","ENSG00000185811","ENSG00000134516","ENSG00000167613","ENSG00000122122","ENSG00000010610","ENSG00000066336","ENSG00000183813","ENSG00000102245","ENSG00000135426")
#'
#' CoNet_moduleBlue <- FastWGCNA_network(object,
#'                                       module,
#'                                       hubGenes,
#'                                       ratio = 0.8,
#'                                       size = 15,
#'                                       linkDistance = 300,
#'                                       linkWidth=6,
#'                                       linkColour.type = NULL,
#'                                       opacityNoHover = 0,
#'                                       save.path = "WGCNA",
#'                                       names = "love")
#' @export
FastWGCNA_network <- function(object,
                              module,
                              hubGenes=NULL,
                              hubGenes.color = mycolor[21],
                              ratio = 0.8,
                              size = 15,
                              linkDistance = 300,
                              linkWidth=6,
                              linkColour.type = NULL,
                              opacityNoHover = 0,
                              save.path = "WGCNA",
                              names = "love"){

  ## 加载必要的包
  nd <- c("networkD3","plyr","grDevices")
  Plus.library(nd)

  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)

  ## 获取数据
  hub <- object$Hub

  ## 进行module观察
  net <- NULL
  for(i in 1:length(object$Hub)){ # i=1
    ## cyt object
    block <- object$Hub[[i]]
    block_cyt <- block$cyt

    ## 找到模块基因
    LuckyVerbose(paste0("Block",i,": Step1-find max mudule and hub genes..."))
    node.i <- block_cyt$nodeData
    node.i <- node.i[node.i$nodeAttr.nodesPresent... %in% module,]
    nrow(node.i)

    ## 选择模块相关数据
    LuckyVerbose(paste0("Block",i,": Step2-edge data preparing..."))
    data.i <- block_cyt$edgeData
    data.i <- data.i[order(data.i$weight,decreasing = T),]#colnames(data.i)
    lg1 <- (data.i$fromNode %in% node.i$nodeName) & (data.i$toNode %in% node.i$nodeName) #table(lg1)
    data.i <- data.i[lg1,]

    ## 根据ratio选择data
    LuckyVerbose(paste0("Block",i,": Step3-select edges via ratio..."))
    portion <- floor(nrow(data.i)*ratio)
    data.i2 <- data.i[1:portion,]

    ## 确定基因分组
    LuckyVerbose(paste0("Block",i,": Step4-decide genes category..."))
    u <- unique(c(data.i2$fromNode,data.i2$toNode))
    if(is.null(hubGenes)){
      color.group <- list(u);names(color.group)[1] <- module
      nodeColour <- list(module);names(nodeColour)[1] <- module
    } else {
      nonhub_genes <- setdiff(u,hubGenes)
      color.group <- list(
        HubGenes = hubGenes,
        OtherGenes = nonhub_genes
      )
      #确定颜色列表
      nodeColour <- list(
        HubGenes = hubGenes.color,
        OtherGenes = module
      )
    }

    #设置边线的颜色
    if(is.null(linkColour.type)){
      linkColour.type = list(
        name =  "YlGnBu",
        select.lower = 1,
        select.upper = 8
      )
    }

    ## 绘制网络图
    LuckyVerbose(paste0("Block",i,": Step5-draw D3 network..."))
    colnames(data.i2)
    nw2 <- Fastcornet2(data= data.i2,
                       source.col= "fromNode",
                       target.col= "toNode",
                       value = "weight",
                       color.group = color.group,
                       control.group = NULL,
                       ratio = 1,
                       size = size,
                       linkColour = NULL,
                       linkColour.type = linkColour.type,
                       control.linkColour = NULL,
                       colourScale.type = NULL,
                       nodeColour = nodeColour,
                       linkDistance = linkDistance,
                       linkWidth = linkWidth,
                       opacityNoHover = opacityNoHover);
    LuckyVerbose(paste0("Block",i,": Step6-save network..."))
    save.Fastcornet(nw2,names = paste0(names,"_",module))

    ## 输出结果
    l.i <- list(
      network = nw2,
      module = list(
        genes = color.group,
        color = nodeColour)
    )
    net <- c(net,list(l.i));
    names(net)[i] <- module
  }

  ## 输出结果
  LuckyVerbose("All done!")
  return(net)

}











