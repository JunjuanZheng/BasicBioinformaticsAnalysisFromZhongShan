
\name{heatmap.dds}
\alias{heatmap.dds}
\title{heatmap.dds}
\description{
根据选择的基因和dds文件，按一定的参数绘制热图。
}
\usage{
  heatmap.dds(dds,
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
              row_names_gp = gpar(fontsize = 16, fontface = "bold"))
}
\arguments{
  \item{dds}{来自DESeq2的结果}
  \item{select}{待画热图的基因}
  \item{design}{含有分组信息的数据框。画热图的数据框将按design的分组信息来排序}
  \item{contrast}{分组的分类。如"N.status"}
  \item{contrast.list}{颜色与对比的对应关系}
  \item{rowscale}{T代表按行对矩阵进行scale的归一化运算}
  \item{transformation}{仅支持"normTransform"和"vst"两种，其它会报道。}
  \item{cluster_rows/cluster_cols}{是否对行/列进行聚类。}
}
\details{

}
\value{
 heatmap plot
}
\references{

}
\author{
  Weibin Huang
}
\note{

}
\seealso{

}
\examples{

}
\keyword{heatmap.dds}

