

## 判断一个向量是否全为T值
#vt <- c(T,"2","1","3");is.one.true(vt)
#vt <- c(F,"2","1","3");is.one.true(vt)
#vt <- c(T,2,1,3);is.one.true(vt)
#vt <- c(T,F);is.one.true(vt)
#vt <- c(F,F);is.one.true(vt)
#' @export
is.one.true <- function(vt){
 if(!is.logical(vt)){
   #vt不全是逻辑值
   if(!is.numeric(vt)){
     #vt不是数值型变量
     l1 <- T %in% vt
     return(l1)
   } else {
     #vt是数值型向量
     print("向量是数值型，无法区分。")
   }
 } else {
   #vt是逻辑型向量
   l1 <- T %in% vt
   return(l1)
 }

}

