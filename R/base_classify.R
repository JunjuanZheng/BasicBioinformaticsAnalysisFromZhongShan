

#' @title Classify a vector by specified classifier
#' @description classify help classify a vector by specified classifier
#' @param vector a vector
#' @param classifier  a list of classifier
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## example
#' vector <- 1:100
#'
#' ## classifier exist
#' classifier1 <- list(
#'   lower = c(1:49),
#'   upper = c(50:90)
#' )
#' v1 <- classify(vector,classifier1)
#' table(v1)
#'
#' ## upper not exist
#' classifier2 <- list(
#'   lower = c(1:49),
#'   upper = c(101:900)
#' )
#' v1 <- classify(vector,classifier2)
#' table(v1)
#' @export
classify <- function(vector,classifier){
  vector <- as.character(vector)
  vector1 <- rep("Not Available",length(vector))
  for(i in 1:length(classifier)){ # i = 1
    n.i <- names(classifier)[i]
    c.i <- as.character(classifier[[i]])
    vector1[vector %in% c.i] <- n.i
  }
  return(vector1)
}















