
\name{getexpr}
\alias{getexpr}
\title{get long data frame from a gene expression matrix.}
\description{
getexpr help get long data frame from a gene expression matrix.
}
\usage{
getexpr(exprs.matirx,
        design,
        select,
        contrast.col,
        contrast.list=NULL)
}
\arguments{
  \item{exprs.matirx}{genes expression matrix or a matrix with individual colnames and marker rownames.}
  \item{design}{design object}
  \item{select}{selected markers or genes}
  \item{contrast.col}{the colnames of the contrast}
  \item{contrast.list}{the elements of a contrast like Np and N0. If NULL,use all the contrast element.}
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

data(expr)
exprs.matirx <- expr[,1:(ncol(expr)-3)]
design=data.frame(condition=Fastextra(colnames(exprs.matirx),"_",2))
select = c("gene1","gene2","gene3","gene4","gene5","gene6")
contrast.col = "condition"
contrast.list = as.character(unique(design$condition))

df1 <- getexpr(exprs.matirx,
               design,
               select,
               contrast.col,
               contrast.list=NULL)

df2 <- getexpr(exprs.matirx,
               design,
               select,
               contrast.col,
               contrast.list)

df3 <- getexpr(exprs.matirx,
               design,
               select,
               contrast.col,
               contrast.list=contrast.list[1:2])
}
\keyword{getexpr}

