
\name{sametrend}
\alias{sametrend}
\title{find out genes that can be plot appropriately in boxplot with student test.}
\description{
 find out genes that can be plot appropriately in boxplot with student test.
}
\usage{
  sametrend(expr.matrix.list,
          res.list,
          logFC.col = "logFC",
          design,
          contrast="N.status",
          sigdif.co.id,
          count.test = T,
          log.convert = T,
          method="wilcox.test",
          p.cutoff = 0.05)
}
\arguments{
  \item{expr.matrix.list}{a list of different gene expression matrix.}
  \item{res.list}{a list of different results which contain logFC columns.}
  \item{design}{the design of gene expression matrix}
  \item{contrast}{the colname of design representing contrast.}
  \item{sigdif.co.id}{the range of test genes}
  \item{logFC.col}{the colnames of results}
  \item{count.test}{whether to make test like t test and so on.Always T.}
  \item{log.convert}{whether to make log scale in gene expression}
  \item{method}{method to do test}\\item{p.cutoff}{cut-off of P value}
}
\details{

}
\value{

}
\references{

}
\author{
  Weibin Huang<654751191@qq.com>
}
\note{

}
\seealso{

}
\examples{
## example
expr.matrix.list=list(test1 = expr,test2=expr)
res.list = list(test1 = sig.exprs,test2 = sig.exprs)
design = condition
contrast="condition"
sigdif.co.id = rownames(sig.exprs)
colnames(sig.exprs)
logFC.col = "logFC"
count.test =T
log.convert = F
method="wilcox.test"
p.cutoff=0.05
test1 <- sametrend(expr.matrix.list,
                   res.list,
                   logFC.col,
                   design,
                   contrast,
                   sigdif.co.id,
                   count.test,
                   log.convert,
                   method,
                   p.cutoff)
}
\keyword{sametrend}




