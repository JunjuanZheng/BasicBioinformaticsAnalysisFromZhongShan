

## FastUpset is based on upsetR and give a fast way to draw an Upset plot

# data # a list of character vector like markers or genes.The names of list was recommanded to be set
# nsets # If NULL,all combs was printed
# nintersects #Number of intersections to plot. If set to NA, all intersections will be plotted
# matrix.color#Color of the intersection points
# main.bar.color #Color of the main bar plot
# mainbar.y.label#The y-axis label of the intersection size bar plot
# sets.bar.color #Color of set size bar plot
# sets.x.label #The x-axis label of the set size bar plot
# order.by #How the intersections in the matrix should be ordered by. Options include frequency (entered as "freq"), degree, or both in any order
# point.size#Size of points in matrix plot
# line.size #Width of lines in matrix plot
# decreasing #How the variables in order.by should be ordered. "freq" is decreasing (greatest to least) and "degree" is increasing (least to greatest)
# shade.color#Color of row shading in matrix
# shade.alpha#Transparency of shading in matrix
# text.scale#Numeric, value to scale the text sizes, applies to all axis labels, tick labels, and numbers above bar plot. Can be a universal scale, or a vector containing individual scales in the following format: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)

#' @export
FastUpset <- function(data,
                      nsets = NULL,nintersects =NA,
                      matrix.color = NULL,
                      main.bar.color = NULL,
                      mainbar.y.label = "Intersection Size",
                      sets.bar.color = NULL,
                      sets.x.label = "Set Size",
                      order.by = c("freq","degree")[1],
                      point.size = 3,line.size =1,
                      decreasing = c(T,F),
                      show.numbers = "yes",
                      shade.color = mycolor[30], shade.alpha = 0.25,
                      text.scale = 2
){

  ## 加载必要的包
  need <- c("UpSetR","dplyr")
  Plus.library(need)

  ## 取唯一值
  data2 <- data
  for(i in 1:length(data)){
    t.i <- length(data[[i]])
    data2[[i]] <- unique(data[[i]])
    if(length(data2[[i]]) < t.i){print(paste0(names(data)[i]," is not a unique vector."))}
  }

  ## 进行计算
  fromList1 <- function(data){

    data1 <- data
    ## 整合
    name.marker <- NULL
    for(i in 1:length(data1)){
      name.marker <- c(name.marker,data1[[i]])
    }
    name.marker <- unique(name.marker)
    name.marker <- data.frame(name.marker)

    ## 算数
    options(warn = -1) #隐藏warning
    for(i in 1:length(data1)){ #i=1
      n.i <- data1[[i]]
      n.i2 <- data.frame(name.marker = n.i,id = 1)
      colnames(n.i2)[2] <- names(data1)[i]
      name.marker <- left_join(name.marker,n.i2,by = "name.marker")
    }
    name.marker[is.na(name.marker)] <- 0
    rownames(name.marker) <- name.marker$name.marker
    name.marker <- subset(name.marker,select = names(data1))
    options(warn = 0) #改回默认值
    return(name.marker)
  }
  data1 <- fromList1(data2)

  ## 计算当前数据的组合类型数
  calculate.combn <- function(data1){
    ## 计算每一行的组合类型
    #data1.i = data1[1,]
    getcombn1 <- function(data1.i){
      a <- as.character(data1.i)
      a1 <- names(data1.i)[a %in% "1"]
      a2 <- paste(a1,collapse = "-")
      return(a2)
    }
    a <- apply(data1,1,getcombn1)
    a1 <- as.data.frame(table(a))
    a1 <- a1[order(a1$Freq,decreasing = T),]
    colnames(a1)[1] <- "Comb"

    ## 对某种分类，获得对应的symbol
    uni <- unique(a);l1 <- NULL
    for(x in 1:length(uni)){ # x =1
      l1.i <- names(a[a %in% uni[x]])
      l1 <- c(l1,list(l1.i))
      names(l1)[x] <- uni[x]
    }

    ## 输出结果
    l2 = list(
      freq = a1,
      category = l1
    )
    return(l2)
  }
  comb <- calculate.combn(data1)
  inter <- nrow(comb$freq)

  ## test
  test.comb <- function(comb){
    category <- comb$category
    for(i in 1:length(category)){ #i=2
      c.i <- category[[i]]
      n.i <- Fastextra(names(category)[i],"-")
      lg <- NULL
      for(n.j in n.i){
        lg.j <- all(c.i %in% data2[[n.j]])
        lg <- c(lg,lg.j)
      }
     lg1 <- all(lg)
     if(lg1){
       print(paste0(names(category)[i]," correct!"))
     } else {
       print(paste0(names(category)[i]," error!"))
     }
    }
  }
  print("correction testing...")
  test.comb(comb)

  ## nsets
  if(is.null(nsets)){
    nsets <- inter
  } else {
    if(nsets > inter){
      print("The nsets you set is too large and had been corrected.")
      nsets <- inter
    } else {
      nsets <- nsets
    }
  }

  ## color
  #矩阵颜色
  if(is.null(matrix.color)){
    matrix.color <- mycolor[4]
  }
  #bar颜色
  if(is.null(main.bar.color)){
    c0 <- c(1,3:20)
    c <- mycolor[c(c0,setdiff(1:length(mycolor),c0))]
    c <- rep(c,10)
    main.bar.color <- c[1:nsets]
  }
  #set颜色
  if(is.null(sets.bar.color)){
    c1 <- mycolor[46:70]
    c1 <- rep(c1,100)
    sets.bar.color <- c1[1:ncol(data1)]
  }

  ## Upset图
  upset(data1,
        nsets = nsets,
        nintersects = nintersects,
        matrix.color = matrix.color,
        main.bar.color = main.bar.color,
        #main.bar.color = "#B3CDE3",
        mainbar.y.label = mainbar.y.label,
        sets.bar.color = sets.bar.color,
        sets.x.label = sets.x.label,
        order.by = order.by,
        #order.by = "freq",
        point.size = point.size,line.size = line.size ,
        decreasing = decreasing,
        show.numbers = show.numbers,
        shade.color = shade.color, shade.alpha = shade.alpha,
        text.scale = text.scale)

  ## 输出
  return(comb)
}






