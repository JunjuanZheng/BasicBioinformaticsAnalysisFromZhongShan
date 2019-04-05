
\name{MApipeline}
\alias{MApipeline}
\title{transform a raw MA object to a comparable MA which is needed in limma:lmFit}
\description{
    Raw MA objects always contains probes with the same gene,or lots of control probes which are considered not differentially expressed in every arrays,and all of this need to deal with before MA go to the lmFit pipeline.I believe MApipeline() give a fast way to do this job.
    MApipeline2 is the plus version of MApipeline().It annotates MA object via sequence enhanced annotation strategy derectly..
}
\usage{
  MApipeline(MA,
             MA.probeid,
             MA.colnames,
             control.probes,
             array.annotation,
             eset,
             enhanced.annotation=T,
             gpl.path=NULL,
             control.weight = F)

  MApipeline2(MA,
              MA.probeid,
              MA.colnames,
              control.probes,
              array.annotation,
              eset,
              gpl.path=NULL,
              control.weight = F)
}
\arguments{
  \item{MA}{a MA List}
  \item{MA.probeid}{the probeid of M matrix or A matrix in MA List.Note that it was the same compose with eset feature data.}
  \item{MA.colnames}{the new colnames of M matrix or A matrix in MA List.Note that it was the subset of colnames of eset exprs data.}
  \item{control.probes}{a list containing control probes.See AddDesignList()}
  \item{array.annotation}{a list containing array annotation information.See AddDesignList()}
  \item{eset}{a project from getElist()}
  \item{control.weight}{whether to gather the expression of control probes and weight them with 2.In most time we recommand control.weight.F.You can also compared two results and decide which one is the best choice.}
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
library(lucky)
load("E:/RCloud/database/DataDownload/GEOqueryDownload/GSE17154/GSE17154_rawEset.rda")
load("E:/RCloud/database/DataDownload/GEOqueryDownload/GSE17154/GSE17154_MAList.rda")
load("E:/RCloud/database/DataDownload/GEOqueryDownload/GSE17154/targets.rda")

## parameters
MA.probeid=1:44544
MA.colnames = Fastextra(targets$FileName,"_",1)
control.probes <- list(
  l1=list(control.col="Reporter Group[role]",control.symbol = "Control"))
array.annotation = list(
  probeid.col="ID",
  symbol.cols=c("Symbol"),
  sequence.col="SEQUENCE",
  anno.cols=c("SYMBOL"),
  db.anno=common.annot)

## Runing
MA1 <- MApipeline(MA,
                  MA.probeid,
                  MA.colnames,
                  control.probes,
                  array.annotation,
                  eset,
                  control.weight = F)

MA2 <- MApipeline(MA,
                  MA.probeid,
                  MA.colnames,
                  control.probes,
                  array.annotation,
                  eset,
                  control.weight = T)

ma1 <- MApipeline2(MA,
                   MA.probeid,
                   MA.colnames,
                   control.probes,
                   array.annotation,
                   eset,
                   gpl.path=NULL,
                   control.weight = T)

ma1 <- MApipeline2(MA,
                   MA.probeid,
                   MA.colnames,
                   control.probes,
                   array.annotation,
                   eset,
                   gpl.path=NULL,
                   control.weight = F)

}
\keyword{MApipeline}




