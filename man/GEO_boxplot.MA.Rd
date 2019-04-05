
\name{boxplot.MA}
\alias{boxplot.MA}
\title{quick way to draw a boxplot for specified contrast via MAList object}
\description{
 boxplot.MA is based on boxplot.genedata2 and a quick way to draw a boxplot for specified contrast via MAList object.boxplot.MA2() is based on boxplot.MA and appropriate for multiple significant markers plotting.
}
\usage{
boxplot.MA(MA1,
           design,
           matrix.type="M",
           contrast,
           select,
           p.val.matrix,
           p.val.col="P.Value",
           p.val.position=NULL)

boxplot.MA2(MA1,
            matrix.type="M",
            cut.n,
            design,
            contrast,
            p.val.matrix,
            p.val.col="P.Value",
            p.val.position=NULL,
            names="test1")
}
\arguments{
 \item{MA1}{MA object after lucky::MApipeline2()}
 \item{design}{a design matrix from model.matrix()}
 \item{matrix.type}{M or A matrix}
 \item{contrast}{the group in every facet}
 \item{select}{marker names}
 \item{p.val.matrix}{a result object containing p.val}
 \item{p.val.col}{colnames of p value in the p.val.matrix}
 \item{p.val.position}{the position of significant symbol in the plot}
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
boxplot.MA(MA1,
           design,
           matrix.type="M",
           contrast=c("Np","N0"),
           select=sig.genes,
           p.val.matrix=sig.expr,
           p.val.col="P.Value",
           p.val.position=NULL)

boxplot.MA(MA1,
           design,
           matrix.type="A",
           contrast=c("Np","N0"),
           select=sig.genes,
           p.val.matrix=sig.expr,
           p.val.col="P.Value",
           p.val.position=NULL)
}
\keyword{boxplot.MA}

