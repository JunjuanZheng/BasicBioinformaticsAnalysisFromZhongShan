

#' @title convertStatus
#' @description convert time and status to a binary variable
#' @param data the data containing time and status columns
#' @param time colname of time
#' @param status colname of status
#' @param event the symbol of event happening
#' @param cutoff the cut off of time to create a binary variable
#' @param name the colnames of new binary variables
#' @param filter if \code{filter = T},then the sample that with less time and event not happening would be consider as \code{NA},which means that we could not do a correct judge.
#' @return a data frame of binary variables
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' data("rna.design")
#' df <- convertStatus(data = rna.design,
#'                     time = "OS.time",status = "OS.status",
#'                     event = "1",
#'                     cutoff = c(365,365*2),
#'                     names = c("Year_1","Year_2"))
#' View(df)
#' @export
convertStatus <- function(data,
                          time,status,
                          event = "1",
                          cutoff = c(365,365*2),
                          names = c("Year_1","Year_2"),
                          filter = T){
  ## 数据转换
  df_1 <- data[c(time,status)]
  colnames(df_1) <- c("time","status")
  df_1$status <- as.character(df_1$status)
  df_1$time <- as.numeric(as.character(df_1$time))

  ## 当时间小于给定值，且事件发生（event=1），状态记为1，其余情况为0.
  ## filter=T,对于未发生事件（status=0）,且时间 < 给定的时间，说明无法进行判断
  df_2 <- NULL
  if(filter != T){
    for(i in 1:length(cutoff)){ # i=1
      c.i <- cutoff[i]
      a.i <- ifelse(df_1$time <= c.i & df_1$status == event,1,0)
      df_2 <- cbind(df_2,a.i)
      colnames(df_2)[i] <- names[i]
    }
  } else {
    for(i in 1:length(cutoff)){
      c.i <- cutoff[i]
      a.i <- ifelse(df_1$time <= c.i & df_1$status == event,1,ifelse(df_1$time <= c.i & df_1$status != event,NA,0))
      df_2 <- cbind(df_2,a.i)
      colnames(df_2)[i] <- names[i]
    }
  }
  df_2 <- as.data.frame(df_2)
  rownames(df_2) <- rownames(df_1)
  #df_3 <- cbind(df_1,df_2)

  ## 输出结果
  return(df_2)
}

















