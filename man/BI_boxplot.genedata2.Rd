
\name{boxplot.genedata2}
\alias{boxplot.genedata2}
\title{quick way to draw boxplot for multiple markers.}
\description{
boxplot.genedata2 is a plus version of boxplot.genedata.It use the exogenous p values from p.val.matrix from limma package or other differentially analysis tools.
}
\usage{
  boxplot.genedata2(select.genes,
                  design,
                  expr.matrix,
                  p.val.matrix,
                  p.val.col="P.Value",
                  p.val.position=NULL,
                  contrast = "N.status",
                  contrast.list = c("N0","Np"))
}
\arguments{
  \item{select.genes}{marker names}
  \item{design}{a design object of a expression matrix.}
  \item{expr.matrix}{a expression matrix.}
  \item{p.val.matrix}{a result object containing p.val}
  \item{p.val.col}{colnames of p value in the p.val.matrix}
  \item{p.val.position}{the position of significant symbol in the plot}
  \item{contrast}{the colname of the group in every facet.}
  \item{contrast.list}{the group in every facet}
  \item{returndata}{whether output the data of boxplot based.}
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
## example
select.genes=co.genes
design=design.x
expr.matrix=expr.x
contrast = "condition"
contrast.list = c("N0","Np")
method = c("t.test", "wilcox.test", "anova", "kruskal.test")[2]
p.val.matrix= p.72305
p.val.col = "P.Value"
boxplot.genedata2(select.genes,
                  design,
                  expr.matrix,
                  p.val.matrix,
                  p.val.col,
                  p.val.position=NULL,
                  contrast,
                  contrast.list)
}
\keyword{boxplot.genedata2}

