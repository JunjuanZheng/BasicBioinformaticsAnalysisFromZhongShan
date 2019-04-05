

####======================heatmap.dds2==========================####
## a plus version of heatmap.dds

# dds# dds object.IF NULL,expr.matrix should provide expression matrix.
# transformation #  normTransform or vst
# expr.matrix #If NULL,dds object gives the expression matrix.
# log.convert #whether do a log scale
# select#select markers or genes
# design#a design object
# contrast.col#the colnames of contrast
# contrast.level#If NULL,use default level.
# contrast.list#a specified list.see example.
# rowscale#whether scaling and centering of matrix-like objects
# expr.name #a legend name of expression data
# cluster_rows,cluster_cols #whether cluster in row or col direction
# row.k # the number of annotated color in row direction
# show_column_names,show_row_names #whether to show column or row names
# clustering_distance_row #see ComplexHeatmap::Heatmap
# clustering_distance_column #see ComplexHeatmap::Heatmap
# column_names_gp#see ComplexHeatmap::Heatmap and the example
# row_names_gp #see ComplexHeatmap::Heatmap and the example
#' @export
heatmap.dds2 <- function(dds,
                        transformation = "normTransform",
                        expr.matrix = NULL,log.convert = T,
                        select,
                        design,
                        contrast.col,#对比所在列名
                        contrast.level=NULL,
                        contrast.list,
                        rowscale=T,
                        expr.name = "expression",
                        cluster_rows = T,cluster_cols=T,
                        row.k = 4,
                        show_column_names = F,
                        show_row_names = T,
                        clustering_distance_row = "euclidean",
                        clustering_distance_column = "euclidean",
                        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                        row_names_gp = gpar(fontsize = 16, fontface = "bold")){
  ##获得矩阵
  if(is.null(expr.matrix)==T){
    #依赖dds提供表达矩阵
    Plus.library("DESeq2");
    if(transformation == "normTransform"){mat1 <- assay(normTransform(dds))} else {if(transformation == "vst"){mat1 <- assay(vst(dds))} else print("error:transformation is not exist!")}
  } else {
    #依赖现有的表达矩阵
    mat1 <- expr.matrix
  }

  ##归一化
  mat1 <- mat1[select,]
  if(rowscale==T){mat1 <- rowScale(mat1,log.convert)} else {mat1 <- mat1}
  #把有NaN值的行去除（代表所有的数值相同，没有意义）
  l1 <- mat1[,1] %in% NaN
  mat1 <- mat1[!l1,]

  ##按contrast排序
  if(is.null(contrast.level)){
    unique.contrast <- unique(as.character((design[,contrast.col])))
  } else {
    unique.contrast <- contrast.level
  }
  u.p <- NULL
  for(i in 1:length(unique.contrast)){
    u.i <- unique.contrast[i]
    u.p.i <- design[,contrast.col] %in% u.i
    u.p.i <- Fastgrep(T,u.p.i)
    u.p <- c(u.p,u.p.i)
  }
  design <- design[u.p,]

  ##heatmap参数
  # rownames(mat1) <- convert(rownames(mat1))
  Plus.library(c("circlize","RColorBrewer","dendextend","ComplexHeatmap","grid"))
  mat1 <- mat1[,rownames(design)]#按design的排列对mat1的病例进行排序
  row_dend = hclust(dist(mat1))
  col_dend = hclust(dist(t(mat1)))
  design1 <- subset(design,select = names(contrast.list))
  ha <- HeatmapAnnotation(design1, col = contrast.list)
  row.cluster <- if(cluster_rows == T){color_branches(row_dend, k = row.k)} else {F}
  col.cluster <- if(cluster_cols == T){color_branches(col_dend, k = 2)} else {F}

  #heatmap绘制
  p <- Heatmap(mat1, name = expr.name,
               cluster_rows = row.cluster,
               cluster_columns = col.cluster,
               show_column_names = show_column_names,
               show_row_names = show_row_names,
               clustering_distance_rows =  clustering_distance_row,
               clustering_distance_columns = clustering_distance_column,
               column_names_gp = column_names_gp,
               row_names_gp = row_names_gp,
               top_annotation = ha)
  return(p)
}






