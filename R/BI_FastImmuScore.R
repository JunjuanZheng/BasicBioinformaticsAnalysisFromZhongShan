
#' @title Fast way to calculate scores for multiple markers
#' @description Fast way to calculate scores for multiple markers
#' @param data a data frame
#' @param marker the colnames of markers
#' @param time the colname of time value for survival data
#' @param status the colname of status value for survival data
#' @param cutoff a list of cutoff for classification of each marker.See example please
#' @param score a list of score for classification of each marker.See example please
#' @importFrom maxstat maxstat.test
#' @importFrom survival Surv
#' @return a LuckyList
#' @seealso \code{\link{FastSurvplot}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' data <- readxl::read_xlsx("E:/RCloud/Experimental meterials/others/lijin/immu_lijin.xlsx")
#' data <- as.data.frame(data,stringAsFactor = F)
#' data$OS <- as.numeric(as.character(data$OS))
#' data$Osoutcome <- as.numeric(as.character(data$Osoutcome))
#' data$DFS <- as.numeric(as.character(data$DFS))
#' data$DFSoutcome <- as.numeric(as.character(data$DFSoutcome))
#' colnames(data)[1:4] <- c("I3","M3","I8","M8")
#' marker = c("I3","M3","I8","M8")
#' time = "OS"
#' status = "Osoutcome"
#' ## 3 cutoff(n)
#' cutoff <- list(I3 = c(25,50,75)/100,
#'                M3 = c(25,50,75)/100,
#'                I8 = c(25,50,75)/100,
#'                M8 = c(25,50,75)/100)
#' ## 4 scores(n+1)
#' score <- list(I3 = c(1,2,3,4),
#'               M3 = c(1,2,3,4),
#'               I8 = c(1,2,3,4),
#'               M8 = c(1,2,3,4))
#' @export
FastScore <- function(data,
                      marker,
                      time,status,
                      cutoff,score){
  ## 将连续变量变成分类评分
  data.marker <- data[marker]
  data.marker2 <- data.marker
  for(i in 1:ncol(data.marker)){ # i=1
    m.i <- data.marker2[,i]
    name.i <- colnames(data.marker2)[i]
    data.marker2[,i] <- quantile_2(vector = m.i,
                                   probs = cutoff[[name.i]],
                                   score = score[[name.i]])
  }

  ## 计算总分
  sum.score <- rowSums(data.marker2)
  new <- cbind(data[c(time,status)],score = sum.score)
  colnames(new)[1:2] <- c("time","status")
  new1 <- cbind(data.marker2,score = sum.score)

  ## maxstat计算cutoff值
  #nd <- c("survival","maxstat");Plus.library(nd)
  mtHL <- maxstat.test(Surv(time, status) ~ score,
                                data=new,
                                smethod="LogRank",
                                pmethod="none")
  all.cutoff <- as.numeric(mtHL[["estimate"]])
  new$class <- ifelse(new$score >= all.cutoff,paste0("score>=",all.cutoff),paste0("score<",all.cutoff))

  ## 生存曲线
  test <- apply(new,1,is.one.na)
  p <- FastSurvplot(data = new,
                    time = "time",
                    status = "status",
                    marker = "class",
                    color = NULL,
                    size = 10,
                    legend.title = "",
                    pval.position = c(max(new$time),1),
                    saveplot = F ,
                    legend = "top")

  ## 输出结果
  l <- list(
    Repeat = list(
      rawData = data,
      marker = marker,
      time = time,
      status = status,
      cutoff = cutoff
    ),
    Data = list(survData = new,
                metaData = new1,
                cutoff = all.cutoff),
    Plot = p
  )
  return(l)
}

#' @export
quantile_2 <- function(vector,
                       probs=c(25,50,75)/100,
                       score = c(1,2,3,4)){
  a <- quantile(vector,probs = probs)
  probs1 <- c(0,a,Inf);probs1 <- unique(probs1)
  vector2 <- cut(vector,breaks = as.numeric(probs1),labels = score,right = F) # vector2
  vector3 <- as.numeric(as.character(vector2))
  return(vector3)
}
# quantile_2(vector,probs,score)





