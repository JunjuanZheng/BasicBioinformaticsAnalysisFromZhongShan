
\name{filter.low}
\alias{filter.low}
\title{filtering lots of different dataset via lots of condition based on lower cut off}
\description{
 filter.low() help quickly filtering lots of different dataset via lots of condition based on lower cut off number.It returns,by every dataset,a targets vetor and the related sub-data via list style.
}
\usage{
  filter.low(list.data,
             list.cutoff.col,
             list.cut.off,
             list.target.col,
             range)
}
\arguments{
  \item{list.data}{a list of data frame}
  \item{list.target.col}{a list of character vectors containing target marker colnames of related data.}
  \item{list.cut.off}{a list of numeric vectors containing target marker cut-off value.}
  \item{range}{the range you want to do the filter.}
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
library(lucky)
state1 <- state.x77;state1 <- as.data.frame(state1)
state1$city <- rownames(state1)
list.data = list(data1 = state1,data2 = state1);colnames(state1)
list.cutoff.col = list(c("Income","Life Exp"),
                       c("Income","Life Exp")) # select people with long life and high income
list.cutoff = list(c(3380,70),
                   c(3380,70))
list.target.col = list("city","city")
range = state1$city[state1$city != "Alabama"] # I don't care Alabama State
l1 <- filter.low(list.data,
                 list.cutoff.col,
                 list.cut.off,
                 list.target.col,
                 range)
show(l1[["data1"]][["select.marker"]])
View(l1[["data1"]][["select.data"]])
}
\keyword{filter.low}




