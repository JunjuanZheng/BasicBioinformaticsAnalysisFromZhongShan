

#' @export
ensemblSummary <- function(genes){

  ##加载必要的包
  need <- c("rvest","stringr","plyr")
  Plus.library(need)


  gurl <- "http://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000186092;r=1:65419-71585"
  md <- gurl %>%  read_html(encoding="UTF-8")
  a <- md %>% html_nodes("div.lhs") %>% html_text();a
  b <- md %>% html_nodes("div.rhs") %>% html_text();b
  c <- md %>% html_nodes("div.twocol") %>% html_nodes("div.row")

}


