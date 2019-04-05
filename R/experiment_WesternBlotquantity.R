
#' @title Tool for protein quantity in Western Blot
#' @description Tool for protein quantity in Western Blot
#' @param name the name of saved .rda
#' @param path Default is NULL.You have to set up an available.Please see the example file in lucky package
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme_bw ggtitle labs annotate theme
#' @importFrom tidyr %>%
#' @importFrom DT datatable
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## See the example data
#' example <- system.file("extdata", "wb quantity.xlsx", package = "lucky")
#' View(readxl::read_xlsx(example))
#'
#' ## A common and useful running
#' res <- WesternBlotQuantity()
#' View(res$Data$metadata)
#' View(res$Data$result)
#'
#' ## If you want a backup,just set the 'name' parameter.The .rda file would be saved in the specified path.
#' res <- WesternBlotQuantity("I love you")
#' @export
WesternBlotQuantity <- function(name = NULL,
                                path = NULL){

  ## Package
  #nd <- c("ggplot2");Plus.library(nd)

  ## data input
  if(is.null(path)){
    path <- "E:/Common/实验技术/Western Blot/R语言算蛋白定量.xlsx"
  }
  data <- readxl::read_xlsx(path,na = c(""))
  data1 <- na.omit(data)

  ## BSA Concentration
  bsa.con <- as.numeric(data1$Value[grepl("BSA concentration",data1$Type)]) #BSA concerntration
  s.v <- as.numeric(data1$Value[grepl("standard volume",data1$Type)]) #standard volume
  sc.volum <- data1$Type[grepl("BSA reference",data1$Type)];
  sc.vol <- as.numeric(Fastextra(sc.volum,"_",2)) #标准品体积
  sc.con <- sc.vol*bsa.con/s.v  #标准品浓度

  ## standard curve
  y <- sc.con
  pattern <- paste("BSA reference_",sc.vol,sep = "")
  x <- as.numeric(data1$Value[Fastmatch(pattern,data1$Type)])
  fit <- lm(y~x)
  a <- as.numeric(fit$coefficients[2])
  b <- as.numeric(fit$coefficients[1])
  sy <- summary(fit);print(sy)


  ## Plot of standard curve
  d <- data.frame(x = x,y = y)
  label.formula <- paste0("Cont == ",round(a,5),"*Abs + ",round(b,5),collapse = "")
  label.r <- paste0("R ^ 2 == ",round(sy$r.squared,5))
  label.adj.r <- paste0("adj.R ^ 2 == ",round(sy$adj.r.squared,5))
  #label <- paste(label.formula,label.r,label.adj.r,collapse = "\n")

  # ggplot strategy
  p <- ggplot(aes(x = x,y = y),data = d) +
   geom_point(color = mycolor[5]) +
   geom_smooth(method='lm',color = mycolor[5],fill = mycolor[5]) +
    theme_bw() +
    ggtitle("Standard Curve") +
    labs(x  = "Absorbance",y = "Concentration") +
    #annotate("text",x = 0.5*max(x),y=max(y),label = label,size = 5,parse = T) +
    annotate("text",x = 0.5*max(x),y=max(y),label = label.formula,size = 5,parse = T) + annotate("text",x = 0.5*max(x),y=max(y)*0.945,label =label.r,size = 5,parse = T) + annotate("text",x = 0.5*max(x),y=max(y)*0.88,label =label.adj.r,size = 5,parse = T)  + # 文字
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold",size = 18,hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(face = "bold",size = 15),
      axis.text = element_text(face = "bold",size = 12),
      legend.title = element_text(face = "bold",size = 12),
      legend.text =element_text(face = "bold",size = 12)
    )
  win.graph(10,10);print(p)

  ## Calculate protein concentration of samples
  logic1 <- unique(Fastgrep(c("BSA concentration","standard volume","protein volume","loading protein","reference","western volume"),data1$Type))
  logic1 <- setdiff(1:nrow(data1),logic1)
  data2 <- data1[logic1,]
  # result
  Sample = data2$Type
  Abs = as.numeric(data2$Value)
  xsb = (as.numeric(data1$Value[data1$Type %in% "standard volume"])/as.numeric(data1$Value[data1$Type %in% "protein volume"])) #稀释比
  Cont = (a*Abs + b) * xsb
  Volume = as.numeric(data1$Value[data1$Type %in% "loading protein"])/Cont
  loadingBuffer = 0.2 * as.numeric(data1$Value[data1$Type %in% "western volume"])
  ComplementBuffer = as.numeric(data1$Value[data1$Type %in% "western volume"]) - loadingBuffer - Volume

  d2 <- data.frame(
    Sample = Sample,
    Absorbance = Abs,
    Concentration = Cont,
    Volume = round(Volume,1),
    ComplementBuffer =round(ComplementBuffer,1),
    loadingBuffer =round(loadingBuffer,1)
  )
  colnames(d2) <- c("样品","吸光度","浓度(ug/ul)","蛋白液体积(ul)","裂解液体积(ul)","loadingBuffer(ul)")
  print(datatable(d2)) #DT Package

  ## OutPut
  l <- list(
    Data = list(metadata = data1,
                result = d2,
                lmfit = fit),
    Plot = p
  )

  if(!is.null(name)){
    # library(tidyr)
    save.path <- Fastextra(path,"[/]") %>% .[-length(.)] %>% paste0(.,collapse = "/") %>% paste0(.,"/","WesternBlotQuantity_",Sys.Date(),"_",name,".rda",collapse = "")
    save(l,file = save.path)
  }
  return(l)
}


















