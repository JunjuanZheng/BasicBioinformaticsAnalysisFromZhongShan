

####========================heatmap.dds============================####
##根据选择的基因和dds文件，按一定的参数绘制热图。
# dds来自DESeq2的结果
# select：待画热图的基因
# design：含有分组信息的数据框。画热图的数据框将按design的分组信息来排序！
# contrast：分组的分类
# contrast.list：颜色与对比的对应关系,
# rowscale=T：T代表按行对矩阵进行scale的归一化运算
# transformation：仅支持"normTransform"和"vst"两种，其它会报道。
# cluster_rows = T,cluster_cols=T：是否对行/列进行聚类。
# column_names_gp = gpar(fontsize = 16, fontface = "bold"), row_names_gp = gpar(fontsize = 16, fontface = "bold")：对行/列的labels进行自定义。
#' @export
heatmap.dds <- function(dds,
                        transformation = "normTransform",
                        expr.matrix = NULL,log.convert = T,
                        select,
                        design,
                        contrast,#对比所在列名
                        contrast.list=c(list(c("N0","Np")),list(c("green","red"))),
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
    library(DESeq2);
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

  ##heatmap参数
  # rownames(mat1) <- convert(rownames(mat1))
  library(circlize);library(RColorBrewer);library(dendextend);library(ComplexHeatmap);library(grid)
  mat1 <- mat1[,rownames(design)]#按design的排列对mat1的病例进行排序
  row_dend = hclust(dist(mat1))
  col_dend = hclust(dist(t(mat1)))
  ascol=function(contrast.list){
    col=contrast.list[[2]]
    names(col) <- contrast.list[[1]]
    col <- list(col)
    names(col) <- contrast
    return(col)
  }
  col = ascol(contrast.list)
  ha <- HeatmapAnnotation(subset(design,select = contrast), col = col)
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






