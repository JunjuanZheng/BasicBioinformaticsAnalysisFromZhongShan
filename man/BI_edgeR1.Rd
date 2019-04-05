
\name{edgeR1}
\alias{edgeR1}
\title{Differencially expression analysis via edgeR package}
\description{
    edgeR1 help do rawcounts-based differencially expression analysis via edgeR package.
}
\usage{
edgeRList <- edgeR1(counts,
                    design,
                    contrast.col,
                    contrast.level = c("Treat","Control"),
                    contrast.control = "Control",
                    cutoff.lFC = 1,
                    cutoff.padj = 0.1,
                    show.report = F,
                    save.file = T,
                    names = "love")
}
\arguments{
\item{counts}{raw counts}
\item{design}{design object}
\item{contrast.col}{contrast colname in the design object}
\item{contrast.level}{the order should be treat then control}
\item{contrast.control}{the control of the contrast level}
\item{cutoff.lFC}{the cut-off of log fold change}
\item{cutoff.padj}{the cut-off of FDR}
\item{save.report}{whether to show quality control report}
\item{save.file}{whether to save edgeRList}
\item{names}{part name of saved file}
}
\details{

}
\value{
edgeRList object that containing significant genes,significant information and edgeR fit object
}
\references{

}
\author{
  Weibin Huang<654751191@qq.com>
}
\note{

}
\seealso{
 edgeR::calcNormFactors
 edgeR::estimateDisp
 edgeR::glmQLFit
 edgeR::glmQLFTest
}
\examples{
## This is a simulative process and NOT RUN
counts = data.counts;rm(data.counts)
design = design.train2
contrast.col = "N.status"
contrast.level = c("Np","N0")
contrast.contol = "N0"
count.filter=10
cutoff.lFC = 1
cutoff.padj = 0.1
save.file = T
names = "love"

## Quick Start
edgeRList <- edgeR1(counts,
                    design,
                    contrast.col,
                    contrast.level,
                    contrast.control,
                    names = "love")

## Complex Setting
edgeRList <- edgeR1(counts,
                    design,
                    contrast.col,
                    contrast.level = c("Treat","Control"),
                     contrast.level = "Control",
                    cutoff.lFC = 1,
                    cutoff.padj = 0.1,
                    save.file = T,
                    names = "love")
}
\keyword{edgeR1}




