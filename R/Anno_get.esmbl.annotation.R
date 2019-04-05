


#'@export
get.esmbl.annotation <- function(ENSEMBL){

  #加载包
  need <- c("rvest","stringr","stringi")
  Plus.library(need)

  #report1
  print(str_c("共有",length(journal.names),"本杂志的ISSN号待搜索..."))

  #按URL获得ISSN号
  df <- NULL;
  for(n in 1:length(journal.names)){
    #转换格式
    journal.i <- journal.names[n]
    special1 = c("[ ]","-","&","_");
    special2 = c("+","-","%26","_")
    for(i in 1:length(special1)){
      journal.i <- gsub(special1[i], special2[i], journal.i)
    }

    gurl <- str_c("https://portal.issn.org/api/search?search[]=MUST=keyproper,keyqualinf,keytitle,notcanc,notinc,notissn,notissnl,unirsrc=",journal.i)
    a <- tryCatch(md <- gurl %>%  read_html(encoding="GBK"), error = function(e)e, finally = NULL)
    if(is.null(a[["message"]])){
      ##程序没有报错
      #提取ISSN号
      ISSN <- md %>%
        html_nodes("div.squaredFour") %>%
        html_nodes("input") %>%
        html_attr("value")

      if(length(ISSN) == 0){
        print(str_c("无法正确提取第",n,"个杂志的ISSN编码。"))
        #说明url虽然没错，但是无法获取ISSN
        df.i <- data.frame(journal =journal.names[n],search.names = "ISSN.NULL",ISSN.code ="ISSN.NULL")
      } else {
        #提取名称
        jn <- md %>%
          html_nodes("div.item-result-block-title") %>%
          html_nodes("h5") %>%
          html_text()

        #提取国家
        #c <- md %>%
        # html_nodes("div.flex-zero") %>%
        # html_nodes("p") %>%
        #html_text()
        #c1 <- c[c(2,7,11)]
        #c2 <- unlist(strsplit(c1,"[:]"))[c(2,4,6)]

        ##组合
        p <- grep("Online",jn)
        if(length(p) == 0){
          print(str_c("第",n,"个杂志并没有网络版"))
          #说明并无Online的版本。不适合用于网络检索。
          df.i <- data.frame(journal =journal.names[n],search.names = "NON.Online",ISSN.code =paste0(ISSN,collapse = "_"))
        } else {
          #说明有Online的版本，适合用于网络检索。
          df.i <- data.frame(journal =journal.names[n],search.names = jn[p],ISSN.code =ISSN[p])
          print(str_c("完成第",n,"个杂志的ISSN编码搜索。"))
        }
      }

    } else {
      ##程序有报错。
      print(str_c("地址解析错误，无法完成第",n,"个杂志的ISSN编码搜索。"))
      ##组合
      df.i <- data.frame(journal =journal.names[n],search.names = "error",ISSN.code ="error")
    }
    df <- rbind(df,df.i)
  }

  #输出结果
  if(show.music == T){
    tuneR::play(music)
    return(df)
  } else {
    return(df)
  }

}

