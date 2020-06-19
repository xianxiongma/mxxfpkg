#' merge multiple files with the more than one intersected rownames
#'
#' merge multiple files with the more than one intersected rownames
#'
#' @param pattern character string containing a regular expression to be matched in the given character vector.
#' @param header a logical value indicating whether the file contains the names of the variables as its first line.
#' @param sep the field separator character of the separate file.
#' @param stringsAsFactors logical: should character vectors be converted to factors?
#' @param namesSep colnames separator.
#' @param colSelect column index to merge.
#' @param nameSelect name index of separated colnames to be new colnames.
#' @return generate a expression matrix (eSet).
#' @export
#' @keywords downGSE
#' @examples
#'  mergeStarOutput(pattern =  ".*count.tab",header = F,sep = "\t",stringsAsFactors = F,namesSep="\\.",colSelect=c(1,2),nameSelect=1)


mergeStarOutput <- function(pattern =  ".*count.tab",header = F,sep = "\t",stringsAsFactors = F,namesSep="\\.",colSelect=c(1,2),nameSelect=1){
  nameList <- list.files(pattern =  pattern)
  matrix <- read.table(nameList[1],header = header,sep = sep,stringsAsFactors = stringsAsFactors)[,colSelect]
  for (i in 2:length(nameList)){
    matrix <- dplyr::inner_join(matrix,read.table(nameList[i],header = header,sep = sep,stringsAsFactors = stringsAsFactors)[,colSelect],by=colnames(matrix)[1])
  }
  colnames(matrix)[2:ncol(matrix)]<-unlist(lapply(nameList,function(nameList) strsplit(as.character(nameList),namesSep)[[1]][nameSelect]))
  mat=as.matrix(matrix[2:ncol(matrix)])
  rownames(mat) <- matrix[,1]
  return(mat)
}