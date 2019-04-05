


#' @keywords ModuleSurv
#' @title Explore relationship between Modules and survival in FastWGCNA pipeline
#' @description ModuleSurv explores relationship between Modules and survival in FastWGCNA pipeline.
#' @param time the time colname in design object
#' @param status the status colname in design object
#' @inheritParams ModuleTrait
#' @inheritParams base::round
#' @importFrom WGCNA labels2colors
#' @importFrom survival coxph Surv
#' @importFrom plyr ddply
#' @importFrom broom tidy
#' @details I.Multiple time and status is supported.    II.If object is a ME matrix,then you must make the row of ME matrix as the same as design.Pay more attention to \code{\link{ModuleSurvValid}}
#' @return LuckyWGCNA object
#' @seealso \code{\link{FastWGCNA}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' library(lucky)
#' object = wgcna;rm(wgcna);gc()
#' design = rna.design.tumor
#' time = c("OS.time","DFS.time")
#' status = c("OS.status","DFS.status")
#' result_MS <- ModuleSurv(object,
#'                         design,
#'                         time,status,
#'                         digits = 3,
#'                         height = 10,width = 8,
#'                         save.path = "WGCNA-test",
#'                         names = "love")
#' @export
ModuleSurv <- function(object,
                       design,
                       time,status,
                       digits = 3,
                       save.path = "WGCNA",
                       names = "love"){
  ## 包
  #nd <-c("survival","WGCNA","broom","plyr")
  #Plus.library(nd)

  ## 产生储存文件夹
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = 0)

  ## get design matrix from every variable
  traitData <- design[,c(time,status)]

  ## get ME matrix from wcgna object
  object.names <- names(object)
  lg1 <- all(c("dataExpr","sft","net","Hub") %in% object.names)
  if(lg1){
    MEs_col <- object$net$MEs
    nSamples <- nrow(MEs_col)
    # gene counts of MEs
    moduleLabels <-  object$net$colors
    moduleColors <-  labels2colors(moduleLabels)
    df_gc <- as.data.frame(table(moduleColors))
    # select traitData
    select <- rownames(object$dataExpr)
    traitData <- traitData[select,]
  } else {
    #这里可能会出错：trait未进行选择
    MEs_col <- object
    nSamples <- nrow(MEs_col)
    df_gc <- NULL
  }

  ## univariate cox regression
  df_result_all <- NULL
  for(i in 1:length(time)){ # i = 1
    time.i <- time[i];status.i <- status[i]
    type.i <- Fastextra(time.i,"[.]",1)
    tD.i <- subset(traitData,select = c(time.i,status.i))
    df_surv <- cbind(tD.i,MEs_col)
    colnames(df_surv)[1:2] <- c("time","status")

    # 计算单变量的univariate cox
    # df_surv_i <- df_surv[1:3];ME = "MEblack"
    get1 <- function(ME){
      df_surv_i <- subset(df_surv,select = c("time","status",as.character(ME)))
      a <- coxph(Surv(time,status) ~ ., data = df_surv_i)
      a_1 <- tidy(a)
      a_2 <- as.data.frame(summary(a)[["coefficients"]])
      a_2 <- a_2["exp(coef)"]
      a_3 <- cbind(a_1,a_2)
      # exp convert
      #x <- a_3[,colnames(a_3) %in% c("conf.low","conf.high")]
      #a_3[,colnames(a_3) %in% c("conf.low","conf.high")] <- exp(x)
      s <- colnames(a_3)[c(2,8,6,7,4,5)]
      a_4 <- subset(a_3,select = s)
      return(a_4)
    }
    test1 <- data.frame(Module = colnames(MEs_col),
                        stringsAsFactors = F)
    df_result <- ddply(.data = test1,
                       .variables = "Module",
                       .fun = get1)
    df_result$type <- type.i
    df_result$FDR <- p.adjust(df_result$p.value, method = "fdr")

    ## 输出表格
    df_result_2 <- data.frame(
      Module = df_result$Module,
      HR = round(df_result$`exp(coef)`,digits),
      CI = paste(round(df_result$conf.low,digits),round(df_result$conf.high,digits),sep = "~"),
      `p-value` = round2(df_result$p.value,digits),
      FDR = round2(df_result$FDR,digits),
      stringsAsFactors = F
    )
    write.csv(df_result_2,paste0(dir,names,"_Module-",type.i," relationship.csv"),row.names = F)

    ## 汇总
    df_result_all <- rbind(df_result_all,df_result)
  }

  ## output
  result <- list(
    Repeat = list(
      time = time,
      status = status,
      digits = digits,
      save.path = save.path,
      names = names
    ),
    Data = list(HR = df_result_all,
                Freq = df_gc),
    Plot = NULL
  )
  return(result)
}










