



# round2 is based on lucky::round.plus.I recommand use round2 to do an optimized show.
#' @keywords round2
#' @title publish-level round
#' @description round.plus is a specific round method for print decimals in publication.round2 is based on lucky::round.plus.round2 is recommanded to do an optimized show.
#' @param x a numeric vector
#' @param digits decimal places
#' @details round2 would not deal with value like NA and Inf.round2 supports multiple vector.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' round2(200.4526)
#' round2(1)
#' round2(0.9)
#' round2(0.00258)
#' round2(0.002)
#' round2(0.001)
#' round2(0.0004)
#' a <- c(200.4526,1,0.9,0.00258,0.002,0.001,0.0004,NA,Inf,-Inf)
#' round2(a)
#' @export
round2 <- function(x,digits=3){
  test <- c(NA,Inf,-Inf)
  out <- NULL
  for(i in x){
    l1 <- i %in% test
    if(l1){
      #特殊值
      out.i <- as.character(i)
    } else {
      #普通值
      out.i <- round.plus(i,digits)
    }
    out <- c(out,out.i)
  }
  return(out)
}


####=========================round.plus()=======================####
### 在实际文章中，把小数按特定位数保存是很重要的。round.plus在round的基础上进行改良，使之符合文章的P值输出习惯。输出类型可能是字符型。
##判断一个数字的小数点前面和后面有多少位数。
#' @export
round.plus <- function(x,digits=3){
  AfterZero <- function(a){
    if(as.integer(a)==a){return(0)} else {
      #a是小数
      if(a<0){a = -a + 1} else {
        if(a<=1){a=a+1}else{
          a = a
        }
      }
      n1 <- floor(log(a, 10) + 1)
      total.len <- nchar(as.character(a))
      n2 <- total.len - 1 - n1
      return(n2)
    }
  }
  if(x>1){round(x,digits)} else {
    if(x==1){
      x1 <- paste0("1.",paste0(rep(0,digits),collapse = ""),collapse = "")
      return(x1)
    } else {
      #数字小于1
      #判断数的小数点后几位
      after0 <- AfterZero(x)#小数点后有几位数。
      if(after0<digits){
        #直接在数的后方加够0
        x <- format(x, scientific=F)
        x1 <- paste0(x,paste0(rep(0,digits-after0),collapse = ""),collapse = "");return(x1)
      } else {
        if(after0==digits){return(x)} else{
          ##小数点后位数较多。
          #digit指定后的最小小数。
          l.vt <- paste0("0.",paste0(rep(0,digits-1),collapse = ""),"1",collapse = "");l.vt <- as.numeric(l.vt)
          if(x<l.vt){return(paste0("<",l.vt,collapse = ""))} else {
            if(x==l.vt){return(l.vt)} else {
              return(round(x,digits = digits))
            }
          }
        }
      }

    }
  }

}
# round.plus(200.4526)
# round.plus(1)
# round.plus(0.9)
# round.plus(0.00258)
# round.plus(0.002)
# round.plus(0.001)
# round.plus(0.0004)




