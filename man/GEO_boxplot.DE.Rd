
\name{boxplot.DE}
\alias{boxplot.DE}
\title{quick way to draw a boxplot for specified contrast via DesignEset object}
\description{
 boxplot.DE() is based on boxplot.genedata2 and a quick way to draw a boxplot for specified contrast via DesignEset object.
}
\usage{
boxplot.DE(DesignEset1,
           contrast,
           select,
           p.val.matrix=sig.exprs,
           p.val.col="P.Value",
           p.val.position=NULL)

boxplot.DE2(DesignEset1,
            cut.n=11,
            contrast,
            p.val.matrix=sig.exprs,
            p.val.col="P.Value",
            p.val.position=NULL,
            names="test1",
            width = 12,height = 8)
}
\arguments{
 \item{DesignEset1}{a DesignEset object.It is always be one after lucky::AnnotateDesignEset()}
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
## example
library(lucky)
data(test-DesignEset1)
contrast=c("Np","N0")
data(test-sig.exprs)
p.val.matrix=sig.exprs
select <- rownames(sig.exprs)[1:10]
p.val.col="P.Value"
boxplot.DE(test-DesignEset1,
           contrast,
           select,
           p.val.matrix=sig.exprs,
           p.val.col,
           p.val.position=NULL)

boxplot.DE2(test-DesignEset1,
            cut.n=11,
            contrast,
            p.val.matrix=sig.exprs,
            p.val.col="P.Value",
            p.val.position=NULL)
}
\keyword{boxplot.DE}

