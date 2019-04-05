
\name{co.matrix2}
\alias{co.matrix2}
\title{co.matrix2}
\description{
以某些基因为参照，对目的矩阵选择可以注释成功的基因。例如，拟通过valid集对train集限定基因进行差异性分析，而valid和train分别属于两个不同的annotation系统，此时co.matrix2()函数将十分快速、有效地找到valid与train的交集基因，并输出统一注释的基因矩阵。co.matrix2()将ENSEMBL作为唯一性标识进行基因的识别。推荐co.matrix2()而不是co.matrix()进行芯片数据的预处理。
}
\usage{
  co.matrix2(list.exprs.matrix,
             list.annotation,
             list.annotation.type,
             list.co_colname,
             control.co_colname,
             control.genes,
             control.genes.type,
             control.annotation,
             output.annotation.type)
}
\arguments{
  \item{list.exprs.matrix}{list(x1=gene.mt,x2=gene.mt)}
  \item{list.annotation}{list(x1=annotation,x2=annotation)}
  \item{list.annotation.type}{list("probeid","probeid").芯片矩阵的原标注类型。}
  \item{list.co_colname}{list(x1 = c("Symbol","ENTREZID"),x2=c("Symbol","ENTREZID"))}
  \item{control.co_colname}{c("SYMBOL","ENTREZID").参考annotation对应的列}
  \item{control.genes}{rownames(dds.valid).参考基因名称向量，保证输出的矩阵包含的基因在此限定内。}
  \item{control.genes.type}{"ENSEMBL".参考基因对应的类型。}
  \item{control.annotation}{db.anotation.总注释库。有最全的注释信息}
  \item{output.annotation.type}{ "ENSEMBL".通常为最唯一的注释类型。一般为ENSEMBL。}
  \item{geneidfactor <- 0}{保证并行不出错的假向量}
}
\details{

}
\value{

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
\keyword{co.matrix2}

