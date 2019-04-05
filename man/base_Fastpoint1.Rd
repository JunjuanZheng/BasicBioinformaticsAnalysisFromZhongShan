
\name{Fastpoint1}
\alias{Fastpoint1}
\title{fast way to make a point plot based on ggplot2}
\description{
 fast way to make a point plot
}
\usage{
Fastpoint1(data,
           select,
           color="#FB8072",
           select.color=1:2,
           legend.position = "right",
           half.border = T,
           alpha = 0.5,
           width = 10,
           height = 10,
           save.file = T,
           names = "love")
}
\arguments{
  \item{data}{a data frame}
  \item{select}{selected marker colnames}
  \item{color}{a color vector or a colname of data frame}
  \item{select.color}{if color is a colname of data frame,points are colored by groups.Thus,select.color set which colors should be used.}
  \item{legend.position }{ "right".legend position}
  \item{half.border }{ T.if show half of black border of the plot}
  \item{width }{ 10.new showing window width}
  \item{height }{ 10.new showing window height}
  \item{save.file }{ T.whether save PDF file}
  \item{names }{ "love".part of PDF file name}
}
\details{

}
\value{
 a point plot
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
data=mtcars
select = c("disp","hp")
data$gear <- factor(data$gear,levels = c(3,4,5))
#all border
p <- Fastpoint1(data,
                select,
                color = "gear",
                select.color=1:3,
                legend.position = "right",
                half.border = F)
#half border
p <- Fastpoint1(data,
                select,
                color = "gear",
                select.color=1:3,
                legend.position = "right",
                half.border = T)

#top legend
p <-Fastpoint1(data,
               select,
               color = "gear",
               select.color=1:3,
               legend.position = "top",
               half.border = T)

#new window
p <- Fastpoint1(data,
                select,
                color = "gear",
                select.color=1:3,
                legend.position = "right",
                half.border = F,
                width = 12,
                height = 12)

#new color.Here legend.position would not work.
p <- Fastpoint1(data,
                select,
                color = "#8DD3C7",
                legend.position = "right",
                half.border = F,
                width = 12,
                height = 12)

Fastpoint1(data = data.i1,
                  select = colnames(data.i1),
                  color = mycolor[13],
                  legend.position = "right",
                  half.border = T,
                  save.file = F) +
       stat_smooth(method = "auto",se=T,
                   color=mycolor[4],
                   fill = mycolor[5],
                   formula = y ~ x)

}
\keyword{Fastpoint1}




