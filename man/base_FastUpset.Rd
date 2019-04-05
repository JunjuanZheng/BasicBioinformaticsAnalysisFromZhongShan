
\name{FastUpset}
\alias{FastUpset}
\title{Give a fast way to draw an Upset plot base on upsetR}
\description{
    FastUpset is based on upsetR and give a fast way to draw an Upset plot
}
\usage{
FastUpset(data,
          nsets = NULL,nintersects =NA,
          matrix.color = NULL,
          main.bar.color = NULL,
          mainbar.y.label = "Intersection Size",
          sets.bar.color = NULL,
          sets.x.label = "Set Size",
          order.by = c("freq","degree")[1],
          point.size = 3,line.size =1,
          decreasing = c(T,F),
          shade.color = mycolor[30], shade.alpha = 0.25,
          text.scale = 2)
}
\arguments{
\item{data}{a list of character vector like markers or genes.The names of list was recommanded to be set}
\item{nsets}{If NULL,all combs was printed}
\item{nintersects}{Number of intersections to plot. If set to NA, all intersections will be plotted}
\item{matrix.color}{Color of the intersection points}
\item{main.bar.color}{Color of the main bar plot}
\item{mainbar.y.label}{The y-axis label of the intersection size bar plot}
\item{sets.bar.color}{Color of set size bar plot}
\item{sets.x.label}{The x-axis label of the set size bar plot}
\item{order.by}{How the intersections in the matrix should be ordered by. Options include frequency (entered as "freq"), degree, or both in any order}
\item{point.size}{Size of points in matrix plot}
\item{line.size}{Width of lines in matrix plot}
\item{decreasing}{How the variables in order.by should be ordered. "freq" is decreasing (greatest to least) and "degree" is increasing (least to greatest)}
\item{shade.color}{Color of row shading in matrix}
\item{shade.alpha}{Transparency of shading in matrix}
\item{text.scale}{Numeric, value to scale the text sizes, applies to all axis labels, tick labels, and numbers above bar plot. Can be a universal scale, or a vector containing individual scales in the following format: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)}"
}
\details{

}
\value{
 an Upset Plot
}
\references{

}
\author{
  Weibin Huang<654751191@qq.com>
}
\note{

}
\seealso{
upsetR::upset()
}
\examples{
## make a data list
data <- list(
  Cohort1 = sample(1:100,50),
  Cohort2 = sample(1:120,60),
  Cohort3 = sample(1:150,75)
  )

## A simple setting
win.graph(width = 12,height = 8);
FastUpset(data,
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          order.by = c("freq","degree")[1])

## More complex setting
FastUpset(data,
          nsets = NULL,nintersects =NA,
          matrix.color = NULL,
          main.bar.color = NULL,
          mainbar.y.label = "Intersection Size",
          sets.bar.color = NULL,
          sets.x.label = "Set Size",
          order.by = c("freq","degree")[1],
          point.size = 3,line.size =1,
          decreasing = c(T,F),
          shade.color = mycolor[30], shade.alpha = 0.25,
          text.scale = 2)

}
\keyword{FastUpset}




