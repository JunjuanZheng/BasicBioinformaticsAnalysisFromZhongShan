

####=========================set.default()====================####
# Usage: help self-defined function to set default values

# input# a list contain input value
# default# a list contain defalut value
# output#  a list contain value that you want to output when the input was as the same as default.
#' @export
set.default <- function(input,default,output){
  # 函数：对于某个变量，输出默认变量
  set.default1 <- function(value,default,default.output){
    if(is.null(default)){l1=is.null(value)} else {
      if(is.na(default)){l1=is.na(value)} else {
        l1 <- value == default
      }
    }
    if(l1){
      return(default.output)
    } else {
      return(value)
    }
    #End
  }
  # 对于多个变量，输出对应的默认变量
  for(i in 1:length(input)){
    input[[i]] <- set.default1(input[[i]],default[[i]],output[[i]])
  }
  return(input)
}




