

## D3 Java Script 配色库
# https://hijiangtao.github.io/2018/03/23/D3-5.0-is-out-3/
# https://github.com/d3/d3-scale-chromatic

## Usage:
# cormatrix #a matrix of correlation index.
# color.group #a named list of markers categery
# control.group # the names of control group
# ratio # select the most significant portion of correlations.it range 0 to 1.
# legend # whether show legend in the netplot output
# size# the index to enhance the degree difference among nodes.If you consider the size of nodes too small,try make size larger.
# linkColour # link colour
# linkColour.type # a list of parameters to RColorBrewer::brewer.pal
# control.linkColour # the color of links to control group nodes.If NULL,use default colours,which is recommanded
# colourScale.type # d3-scale-chromatic style object.1="d3.schemeCategory20",2="d3.schemeCategory20b",3="d3.schemeCategory20c",4="d3.schemeCategory10".Default is 1
# nodeColour# the named list of group color.It is unavailable when colourScale.type exists
# linkDistance# link distance
# linkWidth # the width of link.If it is not a numeric,the default JS("function(d) { return Math.sqrt(d.value); }") would be used.
# opacity # numeric value of the proportion opaque you would like the graph elements to be
# arrows=F#whether show direction between two nodes
# zoom #logical value to enable (TRUE) or disable (FALSE) zooming
#' @export
Fastcornet <- function(cormatrix,
                       color.group = list(),
                       control.group = "control",
                       ratio=0.8,
                       legend=T,
                       size = 15,
                       linkColour=NULL,
                       linkColour.type = list(
                         name = "YlGn",
                         select.lower = 1,
                         select.upper = 9
                       ),
                       control.linkColour = "red",
                       colourScale.type = NULL,
                       nodeColour=list(),
                       linkDistance = 300,
                       linkWidth=6,
                       opacity=1,
                       arrows=F,
                       zoom = T){
  ## 加载必要的包
  nd <- c("networkD3","stringr")
  Plus.library(nd)

  ## 将cormatrix转变为长型数据,内容作为权重
  l1 <- list()
  for(i in 1:ncol(cormatrix)){
    l.i <- cormatrix[,i]
    l1 <- c(l1,list(l.i))
  }
  names(l1) <- colnames(cormatrix)
  cm1 <- stack(l1)
  cm1$target <- rep(rownames(cormatrix),ncol(cormatrix))
  colnames(cm1)[1:2] <- c("weight","source")
  cm2 <- cm1[order(cm1$weight,decreasing = T),]
  cm2 <- cm2[cm2$weight != 1,]
  cm2 <- cm2[seq(2,nrow(cm2),by=2),]

  ## 按ratio选择显著者
  {
  cm3 <- cm2[1:round(nrow(cm2)*ratio),]
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
                     row.names = 1:length(names))
}

  ## 连接线颜色是否按相关系数的大小进行渐变
  if(is.null(linkColour)){
    ## 未提供特定的颜色，说明按相关系数进行渐变色
    Plus.library("RColorBrewer")
    color.r <- brewer.pal(linkColour.type$select.upper,linkColour.type$name)[c(linkColour.type$select.lower,linkColour.type$select.upper)] # color range via RColorBrewer
    ol <- colorRamp(c(color.r[1], color.r[2]))
    test.link <- link[order(link$value,decreasing = F),]
    test.link$linkColour <- rgb(ol(test.link$value),max = 255)
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
               legend=T,#显示节点分组的颜色标签
               arrows=arrows,#是否带方向
               bounded=F,#是否启用限制图像的边框
               #opacityNoHover=0.8,#当鼠标悬停在其上时，节点标签文本的不透明度比例的数值
               zoom = zoom);# np 允许放缩，双击放大


   ## 输出结果
  l2 <- list(
    edge = link,
    node = node,
    netplot = np
  )
  print("All done!")
  return(l2)
}
















