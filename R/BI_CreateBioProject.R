

#' @keywords CreateBioProject
#' @title Create a bioproject with constant format
#' @description Create a bioporject with constant format
#' @param project character;The name of a project you want to create.Default is "..EXAMPLE".
#' @param path character;project location.If NULL,it uses default "E:/iProjects".
#' @param subNO character;The id of sub project like "01A".
#' @details A bio-project contains "constant","data","experiment","learn" and "report".
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' CreateBioProject(subNO="01B") # create 01B sub project
#' CreateBioProject(subNO="01B") # 01B would not be covered
#' CreateBioProject(subNO="01C") # create new sub project
#' ## More normal usage
#' CreateBioProject(project="myproject",
#'                  path="E:/..Example",
#'                  subNO="01A")
#' @export
CreateBioProject <- function(project="..EXAMPLE",
                             path=NULL,
                             subNO="01A"){

  ## default values
  if(is.null(path)){
    path <- "E:/iProjects"
    #path <- "E:/RCloud/RFactory/lucky assistant/bioproject"
  }

  ## project path
  path_1 <- paste0(path,"/",project)
  project_1 <- paste0(project,"_",subNO)

  ## close warning
  old <- options() # save old options
  options(warn = -1) #隐藏warning

  ## create project directory
  dir.create(path_1,recursive=T)

  ## internal function
  #get1 create new directory
  get1 <- function(subdir.i){
    path_2 <- paste0(path_1,"/",subdir.i)
    dir.create(path_2,recursive=T)
  }
  #get2 create biomedical object files
  get2 <- function(newobjects.i){
    path.newob <- paste0(path_1,newobjects.i)
    gc_2 <- file.create(path.newob)
  }

  ## create sub-directory
  dirs_exist <- list.dirs(path = path_1)

  # cost exist
  lg1 <- grep("constant",dirs_exist)
  lg1 <-  length(lg1) == 0
  if(lg1){
    get1("constant")
    constant.files <- c(
      paste0("/constant/",project,"_cost record.docx"),
      paste0("/constant/",project,"_search items.docx"))
    get2(constant.files)
  } else {
    cat("Directory 'constant' exists and would not be created.","\n")
  }

  ## other exist
  lg1 <- grep(subNO,dirs_exist)
  lg1 <- length(lg1) != 0
  if(lg1){
    cat("Sub project",subNO,"exists! Don't be fool to creat it.","\n")
  } else {
    subdir_1 <- c("data","discussion","experiment","learn","report")
    subdir_1 <- paste(subdir_1,"/",project_1,sep = "")
    gc_1 <- sapply(subdir_1,get1)

    ## name biomedical objects
    discussion_QA <- paste0("/discussion/",project_1,"/",project_1,"_discussion_question & answer.docx")
    discussion_read <- paste0("/discussion/",project_1,"/",project_1,"_discussion_paper reading.docx")
    experiment_discription <-  paste0("/experiment/",project_1,"/",project_1,"_experiment_discription.docx")
    experiment_discription2 <-  paste0("/experiment/",project_1,"/",project_1,"_experiment_discription.pptx")
    learn_read <- paste0("/learn/",project_1,"/",project_1,"_learn_paper reading.docx")
    learn_display <- paste0("/learn/",project_1,"/",project_1,"_learn_paper display.pptx")
    report_display <- paste0("/report/",project_1,"/",project_1,"_report_result & display.pptx")
    newobjects <- c(discussion_QA,discussion_read,
                    experiment_discription,
                    experiment_discription2,
                    learn_read,learn_display,
                    report_display)

    ## create biomedical object files
    gc_2 <- sapply(newobjects,get2)
    cat("All done!","\n")
  }
  on.exit(options(old), add = TRUE) # back to old options

}














