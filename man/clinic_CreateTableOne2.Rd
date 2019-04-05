
\name{CreateTableOne2}
\alias{CreateTableOne2}
\title{Optimized baselined information production from tableone package::CreateTableOne}
\description{
    CreateTableOne2 is a simplified and enhanced version of table::CreateTableOne.It can refactor the data frame and make visualization and file saving easier.
}
\usage{
CreateTableOne2(cluster,
                control = NULL,
                strata,
                data,
                save.file = T,
                names = "love")
}
\arguments{
  \item{cluster}{the row markers like Gender}
  \item{control}{the control of row marker like FEMALE}
  \item{strata}{the col marker like N.status}
  \item{data}{a data frame}
  \item{save.file}{whether to save file}
  \item{names}{part name of saved file}
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
## This is a simulative process and NOT RUN
cluster = c("age","gender","his1","his2","T.status","N.status","Group1","Np.count","M.status","Stage","risk.score")
t1 <- CreateTableOne2(cluster,
                      control = NULL,
                      strata = "score.status",
                      data = design.train1,
                      save.file = T,
                      names = "love")

}
\keyword{CreateTableOne2}




