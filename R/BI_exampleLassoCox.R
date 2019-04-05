



#' @title Visualize and analysis of cox.optimized2 result by giving a specified seed
#' @description Visualize and analysis of cox.optimized2 result by giving a specified seed
#' @inheritParams cox.optimized2
#' @param seed a good example from the result of \code{uniqueModel}
#' @param plot.label whether to show labels in the plot of glmnet object
#' @param plot.lwd the line width of the plot of glmnet object
#' @param line.lty the line type of the cutoff line in the plot of glmnet object
#' @param line.lwd the line width of the cutoff line in the plot of glmnet object
#' @param line.col the line color of the cutoff line in the plot of glmnet object
#' @param verbose whether to do a plot report
#' @param save.file whether to save plot
#' @param names part name of the saved file
#' @importFrom glmnet cv.glmnet glmnet predict.cv.glmnet
#' @seealso \code{\link[glmnet]{glmnet}};\code{\link[glmnet]{cv.glmnet}};
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
exampleLassoCox <- function(expr.matrix,
                            design,
                            select,
                            event.status,
                            event.time,
                            event.lower,
                            k=10,
                            seed,
                            optimize.method = "min",
                            plot.label = T,
                            plot.lwd = 2,
                            line.lty=3,
                            line.lwd=2,
                            line.col="blue",
                            verbose = T,
                            save.file = T,
                            names = "love"){
  ## loading packages
  #nd <- c("glmnet");Plus.library(nd)

  ## 矩阵对齐
  expr.matrix1 <- expr.matrix[select,rownames(design)]
  expr.matrix1 <- t(expr.matrix1)

  ##形成cox模型用的数据
  cox.data <- cbind(design,expr.matrix1)

  ### 形成建模数据
  l <- list();
  for(i in 1:length(event.status)){ # i=1
    prognosis <- c(event.time[i],event.status[i])

    ## 去除空值
    cox.data1 <- cox.data[!is.na(cox.data[,prognosis[1]]),]
    cox.data1 <- cox.data1[!is.na(cox.data1[,prognosis[2]]),]

    ## 空值报警系统
    report1 <- sum(is.na(cox.data[,prognosis[1]]))
    report2 <- sum(is.na(cox.data1[,prognosis[2]]))
    if(report1 != 0){
      LuckyVerbose(paste0("NOTE! There are ",report1," NA values in ",prognosis[1])," and would be removed.")
    }
    if(report2 != 0){
      LuckyVerbose(paste0("NOTE! There are ",report2," NA values in ",prognosis[2])," and would be removed.")
    }

    ## 数值型
    time <- as.numeric(as.character(cox.data1[,prognosis[1]]))
    status <- as.numeric(as.character(cox.data1[,prognosis[2]]))
    cox.data1 <- subset(cox.data1,select=select)
    cox.data1 <- cbind(time,status,cox.data1)

    ## 过滤数据
    cox.data1 <- cox.data1[time > event.lower[i],]
    ## 建立formula
    y <- subset(cox.data1,select = c("time","status"))
    x <- subset(cox.data1,select = setdiff(colnames(cox.data1),c("time","status")))
    y <- as.matrix(y)
    x <- as.matrix(x)

    ## 选择最优lamda值。
    set.seed(seed);foldid <- sample(rep(seq(k),length=nrow(x)))#用以进行交叉验证的k向量。
    cv.fit <- cv.glmnet(x,y,family="cox",foldid=foldid);
    if(optimize.method == "min"){
      s = cv.fit$lambda.min
    } else {
      if(optimize.method == "1se"){
        s = cv.fit$lambda.1se
      } else {
        LuckyVerbose("Not a right set of 'optimize.method'.Use 'min' as default.")
        s = cv.fit$lambda.min
      }
    }

    if(verbose == T){
      win.graph(8,8);plot(cv.fit) # 图1
    }

    ## 根据某个seed建立lasso Cox模型
    fit <- glmnet(x,y,family="cox")
    p <- predict(fit,s = s,type = "coefficients")
    p <- as.matrix(p)
    p <- p[which(p[,1] != 0),]
    if(verbose == T){
      win.graph(8,8);# 图2
      plot(fit,xvar = "lambda",label = plot.label,lwd = plot.lwd)
      abline(v=log(s),lty=line.lty,lwd=line.lwd,col=line.col)
    }

    ## 判断排除次序
    if(T){
      o1 <- colnames(x);names(o1) <- 1:length(o1)
      c <- coef(fit)
      c <- as.matrix(c);c.logic <- c == 0;
      o <- rowSums(c.logic);o <- sort(o,decreasing = F)
      o2 <- names(o);names(o2) <- 1:length(o)
    }

    ## 输出图像
    if(save.file == T){
      pdf(paste0(names,"_",event.time,"_Plot of lasso Cox model.pdf"),8,8)
      plot(cv.fit)
      plot(fit,xvar = "lambda",label = plot.label,lwd = plot.lwd)
      abline(v=log(s),lty=line.lty,lwd=line.lwd,col=line.col)
      dev.off()
    }

    ## 输出结果
    l1 <- list(
      Repeat = list(
        seed = seed,
        event.status=event.status,
        event.time=event.time,
        event.lower = event.lower,
        k= k,
        seed = seed,
        optimize.method = optimize.method,
        plot.label = plot.label,
        plot.lwd =  plot.lwd,
        line.lty= line.lty,
        line.lwd=line.lwd,
        line.col=line.col
      ),
      Data = list(
        cv.glmnet = cv.fit,
        glmnet = fit,
        coef = p,
        priority = o2,
        colnames = o1
      )
    )

    ## 循环
    l <- c(l,list(l1))
    names(l)[i] <- event.time
  }
  return(l)
}












