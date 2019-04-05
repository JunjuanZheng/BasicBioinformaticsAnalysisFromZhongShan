
\name{heatmap.dds2}
\alias{heatmap.dds2}
\title{a quick way to draw heatmap}
\description{
a quick way to draw heatmap.heatmap.dds2 is a plus version of heatmap.dds.
}
\usage{
heatmap.dds2(dds,
             transformation = "normTransform",
             expr.matrix = NULL,log.convert = T,
             select,
             design,
             contrast.col,
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
             row_names_gp = gpar(fontsize = 16, fontface = "bold"))
}
\arguments{
\item{dds}{dds object.IF NULL,expr.matrix should provide expression matrix.}
\item{transformation}{normTransform or vst}
\item{expr.matrix}{If NULL,dds object gives the expression matrix.}
\item{log.convert}{whether do a log scale}
\item{select}{select markers or genes}
\item{design}{a design object}
\item{contrast.col}{the colnames of contrast}
\item{contrast.level}{If NULL,use default level.}
\item{contrast.list}{a specified list.see example.}
\item{rowscale}{whether scaling and centering of matrix-like objects}
\item{expr.name}{a legend name of expression data}
\item{cluster_rows,cluster_cols}{whether cluster in row or col direction}
\item{row.k}{the number of annotated color in row direction}
\item{show_column_names,show_row_names}{whether to show column or row names}
\item{clustering_distance_row}{see ComplexHeatmap::Heatmap}
\item{clustering_distance_column}{see ComplexHeatmap::Heatmap}
\item{column_names_gp}{see ComplexHeatmap::Heatmap and the example}
\item{row_names_gp}{see ComplexHeatmap::Heatmap and the example}
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
ComplexHeatmap::Heatmap()
lucky::heatmap.dds()
}
\examples{
data(expr) #use internal data "expr" in lucky package.
dds=NULL
transformation = "normTransform"
expr.matrix=expr[,1:(ncol(expr)-3)]
select=c("gene1","gene2","gene3","gene4","gene5","gene6")
design=data.frame(condition=Fastextra(colnames(exprs.matirx),"_",2),
                  condition1=Fastextra(colnames(exprs.matirx),"_",2),
                  row.names = colnames(exprs.matirx))
contrast.col="condition"
contrast.list=list(condition = c("cell01" = "green",
                                 "cell02" = "gray",
                                 "cell03" = "darkred"),
                   condition1 = c("cell01" = mycolor[1],
                                  "cell02" = mycolor[2],
                                  "cell03" = mycolor[3]))
expr.name = "exprs"
heatmap.dds2(dds=NULL,
             transformation=NULL,
             expr.matrix,
             log.convert=F,
             select,design,
             contrast.col,contrast.level=NULL,
             contrast.list,
             rowscale=T,
             expr.name,
             cluster_rows=T,cluster_cols=T,
             row.k=2,
             show_column_names=T,show_row_names=T)

## example:the setting of contrast.list of continuous value
library(circlize)
col = list(
 cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
 am = c("0" = "yellow","1" = "orange"),
 mpg = colorRamp2(c(17, 25), c("lightblue", "purple"))
)
}
\keyword{heatmap.dds2}

