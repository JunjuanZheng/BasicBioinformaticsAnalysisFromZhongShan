% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BI_CreateBioProject.R
\name{CreateBioProject}
\alias{CreateBioProject}
\title{Create a bioproject with constant format}
\usage{
CreateBioProject(project = "..EXAMPLE", path = NULL, subNO = "01A")
}
\arguments{
\item{project}{character;The name of a project you want to create.Default is "..EXAMPLE".}

\item{path}{character;project location.If NULL,it uses default "E:/iProjects".}

\item{subNO}{character;The id of sub project like "01A".}
}
\description{
Create a bioporject with constant format
}
\details{
A bio-project contains "constant","data","experiment","learn" and "report".
}
\examples{
CreateBioProject(subNO="01B") # create 01B sub project
CreateBioProject(subNO="01B") # 01B would not be covered
CreateBioProject(subNO="01C") # create new sub project
## More normal usage
CreateBioProject(project="myproject",
                 path="E:/..Example",
                 subNO="01A")
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
\keyword{CreateBioProject}
