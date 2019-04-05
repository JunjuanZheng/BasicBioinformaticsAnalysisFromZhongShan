

## Usage:cut a vector into lots of parts.
#' @title cut a vector into lots of parts
#' @description cut a vector into lots of parts
#' @param vt a vector
#' @param nsplit the number of parts you want to split
#' @return a list contain lots of vectors
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' l1 <- cut.vector(1:54800,nsplit=98)
#' View(l1)
#' @export
cut_vector <- function(vt,nsplit=100){
  if(nsplit==1){
    v2 <- list(vt);
    names(v2) <- paste0("1-",length(v2))
  } else {
    ##转化成位置向量
    len.vt=1:length(vt)

    ##间隔
    len1 <- floor(length(len.vt)/nsplit)

    ## low.ci and upper.ci
    low.ci <- 1;
    for(i in 1:(nsplit-1)){
      low.ci[i+1] <- low.ci[i] + len1
    }
    upper.ci <- len1
    for(i in 1:(nsplit-1)){
      upper.ci[i+1] <- upper.ci[i] + len1
    }

    ## upper.ci的最后一个值
    upper.ci[length(upper.ci)] <- length(vt)

    ## 形成列表
    v2 <- NULL
    for(i in 1:length(low.ci)){
      v2.i <- low.ci[i]:upper.ci[i]
      v2.i <- vt[v2.i]
      v2 <- c(v2,list(v2.i))
      names(v2)[i] <- paste0(low.ci[i],"-",upper.ci[i])
    }
  }
  ## 输出结果
  return(v2)
}



cut.vector <- function(vt,nsplit=100){
  if(nsplit==1){
    v2 <- list(vt);
    names(v2) <- paste0("1-",length(v2))
  } else {
    ##转化成位置向量
    len.vt=1:length(vt)

    ##间隔
    len1 <- floor(length(len.vt)/nsplit)

    ## low.ci and upper.ci
    low.ci <- 1;
    for(i in 1:(nsplit-1)){
      low.ci[i+1] <- low.ci[i] + len1
    }
    upper.ci <- len1
    for(i in 1:(nsplit-1)){
      upper.ci[i+1] <- upper.ci[i] + len1
    }

    ## upper.ci的最后一个值
    upper.ci[length(upper.ci)] <- length(vt)

    ## 形成列表
    v2 <- NULL
    for(i in 1:length(low.ci)){
      v2.i <- low.ci[i]:upper.ci[i]
      v2.i <- vt[v2.i]
      v2 <- c(v2,list(v2.i))
      names(v2)[i] <- paste0(low.ci[i],"-",upper.ci[i])
    }
  }
  ## 输出结果
  return(v2)
}
