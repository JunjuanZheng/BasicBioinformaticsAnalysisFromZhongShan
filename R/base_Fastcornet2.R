

#' @title Draw Network plot of D3 style
#' @description Draw Network plot of D3 style
#' @param data a data frame with source,target and value cols
#' @param source.col the colnames of source
#' @param target.col the colnames of target
#' @param value the colnames of value
#' @param color.group  a list of elements goup.
#' @param control.group Default is NULL.The group you want highlight.
#' @param ratio The portion of data you select
#' @param legend Whether show lengend
#' @param size a numeric to adjust the size of nodes
#' @param linkColour the link color.If it's NULL,the link color would be provided by \code{linkColour.type} parameter.
#' @param linkColour.type a list of link color for \code{value} mapping.
#' @param control.linkColour the color of highlight group.It only works when \code{control.group} exists.
#' @param colourScale.type d3-scale-chromatic style object.\cr
#' 1="d3.schemeCategory20",\cr
#' 2="d3.schemeCategory20b",\cr
#' 3="d3.schemeCategory20c",\cr
#' 4="d3.schemeCategory10".\cr
#' Default is 1
#' @param  nodeColour the named list of group color.It is unavailable when colourScale.type exists
#' @param linkDistance link distance
#' @param linkWidth the width of link.If it is not a numeric,the default \code{"JS(function(d) { return Math.sqrt(d.value); })" } would be used.
#' @param opacity numeric value of the proportion opaque you would like the graph elements to be
#' @param arrows Whether show direction between two nodes
#' @param opacityNoHover the legend intensity.From 0-1.Default is 0
#' @param zoom logical value to enable (TRUE) or disable (FALSE)
#' @details linkColour.type make full use of RColorBrewer and would be a good color strategy provider.The parameters of linkColour.type is including:\cr
#' 1.names: One of "Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds" ,"YlGn","YlGnBu","YlOrBr" and "YlOrRd";\cr
#' 2.select.lower: the lower color order;\cr
#' 3.select.upper: the upper color order.
#' @return a list contain a network plot information and edge/node data.
#' @seealso \code{\link[networkD3]{forceNetwork}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' data=state.x77;colnames(data)
#' ## get MIC cormatrix
#' result1 <- FastMINE(data,
#'                     transposition = F,
#'                     method = "all.pairs",
#'                     control.markers="Income",
#'                     target.markers=NULL)
#'
#' ## long data frame network
#' data <- result1[["MINE.result"]][,1:3]
#' colnames(data) #[1] "X var"  "Y var" "MIC (strength)"
#' source.col="X var";target.col="Y var";value="MIC (strength)"
#'
#' ## parameters of Fastcornet2
#' color.group = list(
#'   love = c("Murder","Frost","Area","HS Grad"),
#'   hate = c("Life Exp","Illiteracy"),
#'   control = c("Population","Income"))
#' control.group = "control"
#' linkColour.type = list(
#'   name =  "Blues",
#'   select.lower = 1,
#'   select.upper = 3
#' )
#'
#' # control.linkColour = "red"
#' control.linkColour = NULL
#' colourScale.type = NULL
#' nodeColour=list(
#'   love = mycolor[6],
#'   hate = mycolor[8],
#'   control = mycolor[35])
#'
#' ## Quick Start
#' nw2 <- Fastcornet2(data,
#'                    source.col,
#'                    target.col,
#'                    value,
#'                    color.group,
#'                    control.group,
#'                    ratio=0.8,
#'                    legend=T,
#'                    size = 15,
#'                    linkColour=NULL,
#'                    linkColour.type,
#'                    control.linkColour,
#'                    colourScale.type = NULL,
#'                    nodeColour);nw2$netplot
#' @export
Fastcornet2 <- function(data,
                        source.col=NULL,
                        target.col=NULL,
                        value=NULL,
                        color.group = list(),
                        control.group = NULL,
                        ratio=1,
                        legend=T,
                        size = 15,
                        linkColour=NULL,
                        linkColour.type = list(
                          name = "Blues",
                          select.lower = 1,
                          select.upper = 9
                        ),
                        control.linkColour = NULL,
                        colourScale.type = NULL,
                        nodeColour=list(),
                        linkDistance = 300,
                        linkWidth=6,
                        opacity=1,
                        arrows=F,
                        opacityNoHover = 0,
                        zoom = T){

  ## 加载必要的包
  nd <- c("networkD3","stringr")
  Plus.library(nd)

  ## 长型数据
  a <- c(value,source.col,target.col)
  newcolnames <- Fastmatch(a,colnames(data))
  colnames(data)[newcolnames] <- c("weight","source","target")
  cm3 <- data[order(data$weight,decreasing = T),]
  select.nrow <- floor(ratio*nrow(cm3))
  cm3 <- cm3[1:select.nrow,]

  ## 构建link和node对象
  {
    names = unique(c(as.character(cm3$source),as.character(cm3$target)))
    names = as.character(names)
    if(!is.null(control.group)){
      #说明选择了对照组，并将其排序至最底层
      names <- c(setdiff(names,control.group),control.group)
    }
    source <- as.character(cm3$source)
    target <- as.character(cm3$target)
    value <- cm3$weight
    for(i in 1:length(names)){source[source %in% names[i]] <- i-1}
    for(i in 1:length(names)){target[target %in% names[i]] <- i-1}
    source <- as.integer(as.character(source))
    target <- as.integer(as.character( target))

    ## link文件
    link <- data.frame(source = source,
                       target = target,
                       value = value)
    link <- link[order(link$source,link$target,decreasing = F),];
    rownames(link) <- 1:nrow(link)

    ## 计算连接个数
    degree <- c(as.character(link$source),as.character(link$target))
    degree1 <- as.data.frame(table(degree))
    degree1$degree <- names[as.numeric(as.character(degree1$degree))+1]
    rownames(degree1) <- degree1$degree
    degree1 <- degree1[names,]
    size0 = degree1$Freq

    ##自定义分组
    group <- names
    for(z in 1:length(color.group)){ #z=1
      c.i <- color.group[[z]]
      p.i <- Fastmatch(c.i,names)
      group[p.i] <- names(color.group)[z]
    }
    node <- data.frame(name=names,
                       group = group,
                       size = size0*size,
                       row.names = 1:length(names),
                       stringsAsFactors = F)
  }

  ## 连接线颜色是否按相关系数的大小进行渐变
  if(is.null(linkColour)){
    ## 未提供特定的颜色，说明按相关系数进行渐变色
    Plus.library("RColorBrewer")
    color.r <- brewer.pal(linkColour.type$select.upper,linkColour.type$name)[c(linkColour.type$select.lower,linkColour.type$select.upper)] # color range via RColorBrewer
    ol <- colorRamp(c(color.r[1], color.r[2]))
    test.link <- link[order(link$value,decreasing = F),]
    ##变成0-1
    test.link$scales.values <- scales::rescale(test.link$value)
    test.link$linkColour <- rgb(ol(test.link$scales.values),max = 255)
    test.link <- test.link[as.character(1:nrow(link)),]
    ValjeanCols <- test.link$linkColour
  }


  ##连接线的颜色：有对照marker时请高亮
  if(!is.null(control.group)){
    if(is.null(control.linkColour)){
      ##说明未特别对照连接的颜色。用默认颜色
      print("未指定对照，节点连线采用默认颜色。")
      ValjeanCols <- test.link$linkColour
    } else {
      ##说明指定了对照颜色。用此特定颜色高亮对照
      print("指定了对照，对照节点将高亮，其余节点继续采用默认色。")
      id.c <- match(control.group,node$name) - 1
      lg1 <- test.link$source %in% id.c|test.link$target %in% id.c
      test.link$linkColour[lg1] <- control.linkColour
      ValjeanCols <- test.link$linkColour
    }
  }

  ## linkWidth 线宽
  if(is.numeric(linkWidth)){
    linkWidth = JS(paste0("function(d) { return Math.abs(d.value)*",linkWidth,"; }"))
  } else{
    linkWidth = JS("function(d) { return Math.sqrt(d.value); }")#提取value的值求平方根
  }

  ## 节点颜色：colorScale
  if(is.null(colourScale.type)){
    # 说明使用自定义的颜色
    color1 <- unlist(nodeColour)
    # 给color1排序:按node中第一次出现的位置为正确顺序
    p1 <- Fastmatch(names(color1),node$group)
    p1 <- sort(p1)
    p2 <- as.character(node$group[p1])
    color1 <- color1[Fastmatch(p2,names(color1))]
    # 生成JavaScript对象
    color2 <- paste(str_c("\"",color1,"\""),collapse = ",")
    color3 <- paste0('d3.scaleOrdinal([',color2 ,"])")
    colourScale <- JS(color3)
  } else {
    ## 使用系统颜色
    type <- c("d3.schemeCategory20","d3.schemeCategory20b","d3.schemeCategory20c","d3.schemeCategory10")
    s1 <- type[colourScale.type]
    print(paste0("Use ",s1," color scale..."))
    colourScale <-  JS(paste0("d3.scaleOrdinal(",s1,");"))
  }

  ## networkD3绘制动态网络图
  np <- forceNetwork(Links = link,#线性质数据框
                     Nodes = node,#节点性质数据框
                     Source = "source",#连线的源变量
                     Target = "target",#连线的目标变量
                     Value = "value",#连线的粗细值
                     NodeID = "name",#节点名称
                     Group = "group",#节点的分组
                     Nodesize = "size" ,#节点大小，节点数据框中
                     linkDistance = linkDistance,
                     opacity = opacity,
                     ###美化部分
                     fontFamily="Arial",#字体设置如"华文行楷" 等
                     fontSize = 20, #节点文本标签的数字字体大小（以像素为单位）。
                     linkColour=ValjeanCols,#连线颜色,black,red,blue,
                     colourScale = colourScale,
                     linkWidth = linkWidth,#连线的粗细
                     charge = -100,#数值表示节点排斥强度（负值）或吸引力（正值）
                     legend=legend,#显示节点分组的颜色标签
                     arrows=arrows,#是否带方向
                     bounded=F,#是否启用限制图像的边框
                     opacityNoHover = opacityNoHover,#当鼠标悬停在其上时，节点标签文本的不透明度比例的数值
                     zoom = zoom);
  # np
  ## 输出结果
  node$size <- as.integer(as.character(node$size))/size
  node <- node[order(node$size,decreasing = T),]
  l2 <- list(
    edge = link,
    node = node,
    netplot = np
  )
  return(l2)
}


## Save the result of Fastcornet series function
#' @export
save.Fastcornet <- function(object,names = "love"){
  Plus.library("networkD3")
  metadata <- object
  save(metadata,file = paste0(names,"_networkD3_metadata.rda"))
  networkD3.plot <- object$netplot
  saveNetwork(networkD3.plot,file = paste0(names,"_networkD3_netplot.html"))
  LuckyVerbose(paste0(names,": The metadata and plotdata of networkD3 had been saved!"))
}







