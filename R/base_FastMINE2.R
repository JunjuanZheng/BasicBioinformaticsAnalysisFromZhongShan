

#' @export
FastMINE2 <- function(data,
                      transposition = F,
                      control.markers,
                      target.markers=NULL,
                      nsplit = 55){
  ##




  H19 <- cut.vector(1:length(p1),nsplit=55)
  #process2
  for(i in 19:37){ #i=1 length(H19)
    position.i <- H19[[i]]
    data.x <- data1[,c(p,position.i)]
    result.x <- FastMINE(data=data.x,
                         transposition = transposition,
                         method = "one.pair",
                         control.markers=control.markers,
                         target.markers=target.markers)
    H19[[i]] <- result.x
    save(result.x,file = paste0("./correlation/result_",names(H19)[i],".rda"))
  }




}



