
\name{AnnotateDesignEset}
\alias{AnnotateDesignEset}
\title{AnnotateDesignEset}
\description{
  Annotate DesignEset object.
}
\usage{
 AnnotateDesignEset(DesignEset,
                    enhanced.annotation=T,
                    gpl.path=NULL)

 AnnotateDesignEset2(DesignEset,
                     gpl.path=NULL)
}
\arguments{

}
\details{
    AnnotateDesignEset() help annotate chip via basic annotation and also probe enhanced annotation.It base on co.matrix3() so that merge repeated annotation probe via mean strategy.
    AnnotateDesignEset2() is the plus version of AnnotateDesignEset().AnnotateDesignEset2() annotation chip via sequece enhanced stragegy.Before using AnnotateDesignEset2(),you must prepare a file like "GPLXXXX_information.rda",which was definited by "gpl.path" parameter.
    We recommand that AnnotateDesignEset2() is the prior choice when probe sequences are available.If probe sequences were not available,then use AnnotateDesignEset and set enhanced.annotation=F.
}
\value{
 A New DesignEset object.
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

}
\keyword{AnnotateDesignEset}

