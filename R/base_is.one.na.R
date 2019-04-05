
## 向量是否至少包含一个空值
#' @export
is.one.na <- function(vt){
  l2 <- table(is.na(vt))
  l2.i <- T %in% names(l2) #FALSE则说明没有F,即全为空NA值。
  if(l2.i){
    #说明至少有一个NA值在向量里
    return(T)
  } else {
    return(F)
  }
}
# a <- c(NA,NA);is.one.na(a)
# b <- c(NA,2);is.one.na(b)
# c <- c(1,3);is.one.na(c)
