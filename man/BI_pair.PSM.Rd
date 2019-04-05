
\name{pair.PSM}
\alias{pair.PSM}
\title{pair.PSM}
\description{
通过倾向性评分对样本进行1：1的matched。
支持"nonrandom","MatchIt"两种策略。
}
\usage{
  pair.random(data,
              contrast,
              cluster,
              seed = 2018,
              ratio=0.2,
              pair = TRUE,
              summary = T)
}
\arguments{
  \item{contrast}{string.colnames.}
  \item{seed}{seed used in sample()}
  \item{pair}{logic.If T,give a random-matched pair data by giving a contrast}
  \item{summary}{logic.If T,give a simple data.}

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
data(data.selfprocess2.3)
contrast="N.status"
treat = "Np"
cluster = c("age","gender","pT")

# strategy = "nonrandom"
strategy = c("nonrandom","MatchIt")[1]
pr1 <- pair.PSM(data,
                contrast,treat,cluster,
                strategy,
                seed = 2018,
                ratio = 1,
                summary = F)

## strategy = "MatchIt"
strategy = c("nonrandom","MatchIt")[2]
method = c("nearest","exact", "full", "genetic" , "optimal","subclass")[1]
pr2 <- pair.PSM(data,
                contrast,treat,cluster,
                strategy,method,#strategy = "nonrandom"时可忽略
                seed = 2018,
                ratio=1,
                summary = F)
}
\keyword{pair.PSM}

