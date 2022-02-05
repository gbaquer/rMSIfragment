#rMSIAnnot Class
#' @import rMSI
#' @import rMSIproc
#' @import ggnet
#' @import network
#' @import sna
#' @import ggplot2
#' @import gridExtra
#' @import PRROC
#' @import magic
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @import data.table

#'
print.rMSILipidAnnot<-function(x){
  title=as.character(substitute(x))
  #class(x)<-"data.frame"
  if(is.null(x$correlation))
    View((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","adduct_score")]),title)
  else
    View((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","lipid_occurences","correlation","adduct_score")]),title)

}
