#' using limma package to do differential expression analysis
#'
#' using limma package to do differential expression analysis
#'
#' @param GSE a GSE id or a data.frame. if a data.frame, first column must be gene symbol.
#' @param columnGSE which column should be used to do differential expression analysis.
#' @param GPL the matched GPL platform of GSE id.
#' @param method max or mean.
#' @param columnGPL which column to convert probe id.
#' @param group if GSE is a data.frame, a group vector should be provided. Otherwise omit.
#' @param sepTitle separator of the title.
#' @param groupIndex which vector should be taken after the title was split.
#' @param logFC log fold change.
#' @param pvalue cutoff to screen result.
#' @return generate a differential expression analysis results.
#' @export
#' @keywords limmaResult
#' @examples
#' createLimmaOutput(GSE = fpkm, group = group)
#' createLimmaOutput(GSE="GSE19136",columnGSE = c(1,2,4,5,7,8,10,11),GPL='GPL570', method="max", columnGPL=c(1,11),logFC=1.5,pvalue=0.05)


createLimmaOutput <- function(GSE="GSE65858", columnGSE=1:ncol(eset), GPL='GPL10558', method="max", columnGPL=NULL, group=NULL, sepTitle="_", groupIndex=5, logFC=1.5, pvalue=0.05){
  if (is.data.frame(GSE) | is.matrix(GSE)){
    if (method=="mean"){
      aggr1 <- aggregate(GSE[,2:ncol(GSE)],by=list(GSE[,1]),mean)
      rownames(aggr1) <- aggr1[,1]
      aggr1[,1] <- NULL
	  group_list <- factor(group)
    } else {
      aggr1 <- aggregate(GSE[,2:ncol(GSE)],by=list(GSE[,1]),max)
      rownames(aggr1) <- aggr1[,1]
      aggr1[,1] <- NULL
	  group_list <- factor(group)
    }
  } else {
    eSet<-GEOquery::getGEO(GSE, destdir='./', getGPL=F)
    eset <- as.data.frame(exprs(eSet[[1]]))[,columnGSE]
    metadata <- pData(eSet[[1]]) 
    x<-metadata$title[columnGSE]
    gpl <- GEOquery::getGEO(GPL, destdir=".") 
    if (is.null(columnGPL)){
      merg1 <- merge(Table(gpl)[,c(1,which(grepl("ymbol|YMBOL",colnames(Table(gpl)))))],eset,by.x=colnames(Table(gpl))[1],by.y="row.names")
    } else {
      merg1 <- merge(Table(gpl)[,columnGPL],eset,by.x=colnames(Table(gpl))[1],by.y="row.names")
    }
    
    if (method=="mean"){
      aggr1 <- aggregate(merg1[,3:ncol(merg1)],by=list(merg1[,2]),mean)
    } else {
      aggr1 <- aggregate(merg1[,3:ncol(merg1)],by=list(merg1[,2]),max)
    }
    rownames(aggr1) <- aggr1[,1]
    aggr1[,1] <- NULL
	group_list <- factor(unlist(lapply(x,function(x) strsplit(as.character(x),sepTitle)[[1]][groupIndex])))
  }
  design<-model.matrix(~-1+group_list) 
  contrast.matrix <- makeContrasts(contrasts = paste0("group_list",sort(unique(as.character(group_list)))[1],"-group_list",sort(unique(as.character(group_list)))[2]), levels = design)
  if (min(aggr1) < 0 | max(aggr1) < 50) {
    fit <- lmFit(aggr1,design)  
    fit1 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit1)  
    tempOutput <- topTable(fit2, coef=paste0("group_list",sort(unique(as.character(group_list)))[1],"-group_list",sort(unique(as.character(group_list)))[2]), n=nrow(fit2),lfc=log2(logFC))
    dif <- subset(tempOutput,tempOutput[,4] < pvalue)
    return(dif)
  } else {
    fit <- lmFit(log2(aggr1+0.000001),design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)
    tempOutput <- topTable(fit2, coef=paste0("group_list",sort(unique(as.character(group_list)))[1],"-group_list",sort(unique(as.character(group_list)))[2]), n=nrow(fit2),lfc=log2(logFC))
    dif <- subset(tempOutput,tempOutput[,4] < pvalue)
    return(dif)
  }
}