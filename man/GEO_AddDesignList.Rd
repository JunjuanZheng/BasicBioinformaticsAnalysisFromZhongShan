
\name{AddDesignList}
\alias{AddDesignList}
\title{AddDesignList}
\description{
 After getElist and make a DesignList.
 }
\usage{
  AddDesignList(eset,
              factor,
              levels,
              control.probes,
              array.annotation,
              NA.method = c("one","all")[1],
              normalizeBetweenArrays=F)
}
\arguments{
\item{factor}{a group defined vector.}
\item{levels}{the levels of factor.Note that Control should be placed last.}
\item{control.probes}{a list of control probe information including control.col and control.symbol}
\item{array.annotation}{array annotation information including probeid.col,symbol.cols,anno.cols,sequence.col and db.anno.}
\item{NA.method}{one of "one" and "all".one means the row should be deleted just as one NA value exists,and all means the row should be deleted only when all elements were NA value.Default is "one".}
}
\details{

}
\value{
 DesignList Object
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
## parameters
factor = ifelse(pdata$pnstage == "N0","N0","Np")
levels = c("Np","N0") # factor的level
control.probes = list(
   l1=list(control.col="Probe_Type",control.symbol = c("A","I")),
   l1=list(control.col="Source",control.symbol = "ILMN_Controls"))
array.annotation = list(
     probeid.col="ID",
     symbol.cols=c("Symbol","Unigene_ID","Entrez_Gene_ID"),
     anno.cols=c("SYMBOL","UNIGENE","ENTREZID"),
     sequence.col = "SEQUENCE",
     db.anno=common.annot)
## Running
AddDesignList(eset,
              factor,
              levels,
              control.probes,
              array.annotation,
              NA.method = c("one","all")[1],
              normalizeBetweenArrays=F)
}
\keyword{AddDesignList}

