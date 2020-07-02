# THIS FILE CONTAINS ONE FUNCTION (GetGenesHTML) WHICH IS CALLED BY THE CLASS COMPARISON,                   ##
# CLASS PREDICTION, SURVIVAL, AND QUANTITATIVE TRAIT ANALYSES.  THE OUTPUT OF THIS FUNCTION                 ##
# IS A LIST CONTAINING STRINGS WHICH CAN BE INSERTED INTO HTML OUTPUT FILES TO CREATE                       ##
# THE TABLE OF SIGNIFICANT GENESETS AND TABLE OF GENES WITHIN THOSE GENESETS (FOR GENE SET                  ##
# COMPARISON) OR THE TABLE OF 'SIGNIFICANT GENES' (FOR ALL OTHER ANALYSIS TYPES), AND THE GENE              ##
# ONTOLOGY OBSERVED v. EXPECTED TABLE (IF APPLICABLE).  THE CHROMOSOMAL DISTRIBUTION BARPLOT                ##
# IS CREATED AS A SIDE EFFECT (BUT NOT A RETURNED OBJECT) OF THIS FUNCTION.                                 ##
#                                                                                                           ##
#  Notes:                                                                                                   ##
#  - 'classifiers' and 'GenesFiltered' are defined in <GetXXXXResults.R> files.                             ##
#  - 'AnalysisType', 'IsSingleChannel', 'GeneListFileName' defined in <RunXXXX.R> scripts.                  ##
#  - RunPermP' seems to be used only for 'ClassComparison', and is defined elsewhere ??                     ##
#  - codeSurv' is defined in 'GetComparisonResults.R' file.                                                 ##
#  - CVsupport', 'methods' and 'SameN' are defined in 'GetPredictionResults.R' file.                        ##
#  - SampleCode', 'SampleGroup' and 'GeneExpData' are defined in 'GetBinaryTreePredictionResults.R' file.   ##
#  - CVsupport', 'coxcoef' are defined in 'GetRiskResults.R' file.                                          ##
#   'ncomp' and 'pHotelling' required for Hotelling test are defined in 'RunGOComparison.R',                ##
#                 'GetGOcomparison.R' and other gene set comparison file.                                   ##
#  - 'StatusVariable' is defined in <RunXXXX.R> script.                                                     ##
#  - 'genepairIndex' is defined in 'GetPredictionResults.R> file.                                           ##
#  - 'xikprime' is defined in 'GetPAMResults.R> file.                                                       ##
#  - genelistfilename parameter is used in analysis which outputs more than one gene list.                  ##
#    If it is not specified, it is defined as GeneListFileName.                                             ##

GetGenesHTML <- function(
    index
  , classlabels=NULL
  , paired=F
  , TruncLevel=10^(-7)
  , RefClass = NULL
  , univtest = FALSE
  , puniv = NULL
  , DoMasterAnnotations=T
  , genelistfilename = NULL
  , WorkPath
  , Hotelling = FALSE
  , GeneSetGlobalTest = FALSE
  , GeneSetGSA = FALSE
  , pGeneSet = NULL
  , npermGSA = NULL
  , GSAscore = NULL
  , GeneSetUnit = NULL
  , GeneTableCaption = NULL
#  , IsSingleChannel = NULL
  , AnalysisType = NULL
  , SameN = NULL
  , classifiers = NULL
  , LogBase = NULL
  , methods = NULL
  , CVsupport = NULL
  , GeneListFileName = NULL
  , ProjectPath = NULL
  , AnalysisName = NULL
  , ArrayToolsPath = NULL
  , genepairIndex = NULL
  , uselasso = NULL
) {

  Trunc <- function(x, cutoff=NULL, pos=T, TruncLevel=10^(-7), npermGSA=NULL,
                    GSAscore=NULL, classlabels=NULL) {
     # Unify the way a numerical number is printed out in a HTML table.
     # Try not to use round, signif individually in the code unless there is a reason.
	 # Input:
	 #   pos: is the value always positive? Useful for printing correlation, t-statistic and coefficients from a regression.
     TruncDigits <- as.integer(log(TruncLevel, 10)*(-1))
     if (is.null(cutoff)) {
       if (pos) {
         ifelse(x >= TruncLevel, round(x, TruncDigits), paste("<", TruncLevel))
       } else {
         round(x, 3) # should compare with the header 't-value' width. Too more can result in word wrap.
       }
     } else if (!is.null(npermGSA)) {
      	options(scipen=4) #use fix notation for decimal value.
        TruncDigits <- 5
        TruncLevel <- 1/npermGSA  # override TruncLevel for GSA method!
        sgn <- if (length(unique(classlabels)) == 2 | exists("StatusVariable")) {
                  ifelse(GSAscore >=0, "(+)", "(-)")
               } else rep("", length(x))
        ret.x <- ifelse(x >= TruncLevel,
                		ifelse(x <= cutoff, paste("<font color=\"red\">", round(x, TruncDigits), sgn, "</font>"),
                           paste(round(x, TruncDigits), sgn)),
     		            paste("<font color=\"red\">", "<", round(TruncLevel, TruncDigits), sgn, "</font>"))
      	options(scipen=0) # restore default
      	return(ret.x)
     } else {
      	options(scipen=4) #use fix notation for decimal value.
        ret.x <- ifelse(x >= TruncLevel,
                		ifelse(x < cutoff, paste("<font color=\"red\">", round(x, TruncDigits), "</font>"),
                           round(x, TruncDigits)),
     		            paste("<font color=\"red\">", "<", TruncLevel, "</font>"))
      	options(scipen=0) # restore default
      	return(ret.x)
     }
  }

  #Qian 5/7/09
  signifNew <- function(x, digits=3) {
    return(ifelse(x >= 0,
      ifelse(x >= 1, round(x, digits), ifelse(x >= 1e-7, signif(x, digits), "< 1e-07")),
      ifelse(x <= -1, round(x, digits), ifelse(x <= -1e-7, signif(x, digits), "> -1e-07")))
    )
  }
  #EndQian

  # DEFINE PARAMETERS FOR DIFFERENT ANALYSIS TYPES.

  GeneSetComparison <- F
  HasClasses <- F
#  ValueDef <- ifelse(IsSingleChannel,"intensities","ratios")
  ValueDef <- "values"

  if (AnalysisType == "ClassComparison" | AnalysisType == "BlockClassComparison") {
    HasClasses <- T
  }

  if (AnalysisType == "SAM" ) {
    #Note that SAM does not have the 'classifiers' object.
    RunPermP <- F
    HasClasses <- T
  }

 if (AnalysisType == "ClassPrediction" | AnalysisType == "BinaryTreePrediction" | AnalysisType == "PAM") {
    RunPermP <- F
    HasClasses <- T
    DoGOOvE <- F; MinObserved <- MinRatio <- NA
  }

  if (AnalysisType == "GOComparison" | AnalysisType == "PathwayComparison" | AnalysisType == "GeneListComparison") {
    RunPermP <- F
    GeneSetComparison <- T
    HasClasses <- if (!exists("StatusVariable")) T else F
    DoGOOvE <- F; MinObserved <- MinRatio <- NA
  }
  GONumTotCat <- GONumSigCat <- NULL

  if (AnalysisType == "SurvivalRisk") {
    paired <- F
    univtest <- F
    RunPermP <- F
    DoGOOvE <- F; MinObserved <- MinRatio <- NA
  }

  if (AnalysisType == "Survival") {
    paired <- F
  }

  if (AnalysisType == "Correlation") {
    paired <- F
    # RunPermP <- F MLI 2/17/2010
  }

  # CREATE CODE FOR THE HTML TABLE OF SIGNIFICANT GENE SETS, FOR GENE SET COMPARISON.

  if (GeneSetComparison) {

    # First, create a table of gene set categories ('GeneSetIndex' is the list of significant
    # categories in the increasing p-value order.

    # NOTE: For Leonid's original code for case when AddGeomeanColumns=T, see
    #       'PrepareGOComparisonTableHTML.R' file from v3.4.  AddGeomeanColumns has been
    #       hard-coded as FALSE for some time now, thus GeneSetTable2 is also no longer used.

    GeneSetIndex <- index

    # Get the R object with gene set categories as a list.
    load(paste(WorkPath, "/GeneSetList.rda", sep=""))
    GeneSetList.raw <- GeneSetList[GeneSetIndex] #note GeneSetList is the object name stored in GeneSetList.rda

    GeneSetList <- list(2) #note now GeneSetList is a list !

    for(i in 1:length(GeneSetList.raw) ) {
       GeneSetListElement <- strsplit(GeneSetList.raw[i],split=" ")[[1]]
       GeneSetList[[i]] <- as.double(GeneSetListElement[!is.element(GeneSetListElement,c(""," "))])
    }

    load(paste(WorkPath, "/GeneSetName.rda", sep=""))

    if (AnalysisType == "GOComparison") {
      GeneSetIdName <- "GO category"
      GeneSetId <- GeneSetName[GeneSetIndex]
      GeneSetnames <- GeneSetId
      GeneSetCode <- GeneSetId
      load(paste(WorkPath, "/GODesc.rda", sep=""))
      GOOntologyName <- "GO ontology"
      GOOntologySig <- GOOntology[match(GeneSetnames, names(GOOntology))]
      GeneSetDescriptionName <- "GO term"
      GeneSetdescriptions <- GOTerm[match(GeneSetnames, names(GOTerm))]
      GONumTotCat <- c(length(which(GOOntology=="CC")),length(which(GOOntology=="MF")),length(which(GOOntology=="BP")))
      GONumSigCat <- c(length(which(GOOntologySig=="CC")),length(which(GOOntologySig=="MF")),length(which(GOOntologySig=="BP")))

    } else {
      GOOntologyName <- NULL
      GOOntologySig <- NULL
    }

    if (AnalysisType == "PathwayComparison") {
      GeneSetId <- GeneSetName[GeneSetIndex]
      GeneSetnames <- GeneSetId
      GeneSetCode <- GeneSetId
      GeneSetIdName <- paste(PathwayType, "Pathway")
      GeneSetDescriptionName <- "Pathway description"
      if (PathwayType == "Biocarta") {
        BiocartaTXT <- read.table(paste(ArrayToolsPath, "/Pathways/Output/Names_Biocarta.txt", sep=""), sep="\t", header=F)
        GeneSetdescriptions <- as.character(BiocartaTXT[[2]]) [match(GeneSetnames,as.character(BiocartaTXT[[1]]))]
        GeneSetdescriptions <- paste("<A HREF=\"http://cgap.nci.nih.gov/Pathways/BioCarta/", GeneSetId, "\">", GeneSetdescriptions, "</A>", sep="")
      }
      if (PathwayType == "Kegg") {
        KeggTXT <- read.table(paste(ArrayToolsPath, "/Pathways/Output/Names_KEGG.txt", sep=""), sep="\t", header=F)
        GeneSetdescriptions <- as.character(KeggTXT[[2]]) [match(GeneSetnames,as.character(KeggTXT[[1]]))]
		    GeneSetdescriptions <- paste("<A HREF=\"http://www.kegg.jp/kegg-bin/show_pathway?", GeneSetId, "\">", GeneSetdescriptions, "</A>", sep="")
      }
    }

    if (AnalysisType == "GeneListComparison") {

      GeneSetName <- unlist(lapply(strsplit(GeneSetName,"/"), function(x) rev(x)[1]))
      GeneSetName <- unlist(strsplit(GeneSetName,".txt"))
      GeneSetnames <- GeneSetName[GeneSetIndex]

      if (GeneSetGroup == "TranscriptionFactor") {
        GeneSetCode <- unlist(lapply(strsplit(GeneSetnames, "_"), function(x) rev(x)[1]))
        GeneSetCode <- gsub(" ", "", GeneSetCode)
      	GeneSetnames <- paste("<a href=\"http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=", GeneSetCode, "\">", GeneSetnames, "</a>", sep="")
      }

  	  # Ting added code to support predicted TF target genelists 8/11/2014
  	  if (GeneSetGroup == "TranscriptionFactorPredict") {
          GeneSetCode <- unlist(lapply(strsplit(GeneSetnames, "_"), function(x) rev(x)[1]))
          GeneSetCode <- gsub(" ", "", GeneSetCode)
  		GeneSetnames <- paste("<a href=\"http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=", GeneSetCode, "&rm=present&collection=CORE", "\">", GeneSetnames, "</a>",sep = "")
  	  }

      if (GeneSetGroup == "MicroRNA") {
        GeneSetCode <- unlist(lapply(strsplit(GeneSetnames, "_"), function(x) rev(x)[1]))
        GeneSetCode <- gsub(" ", "", GeneSetCode)
      	GeneSetnames <- paste("<a href=\"http://microrna.sanger.ac.uk/cgi-bin/sequences/query.pl?terms=", GeneSetCode, "\">", GeneSetnames, "</a>", sep="")
      }
      if (GeneSetGroup == "ProteinDomain") {
        GeneSetCode <- unlist(lapply(strsplit(GeneSetnames, "_"), function(x) x[1]))
        GeneSetCode <- gsub(" ", "", GeneSetCode)
      	if (tolower(substr(GeneSetCode[1],1,2))=="pf")
          GeneSetnames <- paste("<a href=\"http://pfam.sanger.ac.uk/family?acc=", GeneSetCode, "\">", GeneSetnames, "</a>", sep="")
        if (tolower(substr(GeneSetCode[1],1,2))=="sm")
          GeneSetnames <- paste("<a href=\"http://smart.embl.de/smart/do_annotation.pl?DOMAIN=", GeneSetCode, "\">", GeneSetnames, "</a>", sep="")
      }
      if (GeneSetGroup == "Broad") {
      	GeneSetCode <- GeneSetnames
        GeneSetnames <- paste("<a href=\"http://www.broad.mit.edu/gsea/msigdb/cards/", GeneSetnames, ".html\"  TARGET=_blank>", GeneSetnames, "</a>", sep="")
      }
      #Qian 3/13/09
      if (GeneSetGroup == "Lymphoid") {
      	GeneSetCode <- GeneSetnames
      	GeneSetnames <- paste("<a href=\"http://lymphochip.nih.gov/cgi-bin/signaturedb/sig_id_query.cgi?Signature=", GeneSetCode, "\">", GeneSetnames, "</a>", sep="")
      }
      #EndQian
      if (GeneSetGroup == "GeneList") {
      	GeneSetCode <- GeneSetnames
      }

      GeneSetId <- GeneSetnames
      GeneSetIdName <- paste(GeneSetGroup, " GeneSets")
      GeneSetdescriptions <- NULL
      GeneSetDescriptionName <- NULL
    }

    n <- length(GeneSetIndex)
    order.pvalue <- order(pGeneSet$pLS)
    pGeneSet <- pGeneSet[order.pvalue,,drop=F]
    GSAscore <- GSAscore[order.pvalue] # MLI 11/10/2008

    GeneSetIndex <- GeneSetIndex[order.pvalue]
    GeneSetList <- GeneSetList[order.pvalue]
    GeneSetId <- GeneSetId[order.pvalue]
    GeneSetnames <- GeneSetnames[order.pvalue]
    GeneSetdescriptions <- GeneSetdescriptions[order.pvalue]
    GeneSetCode <- GeneSetCode[order.pvalue]
    if (AnalysisType == "GOComparison") { GOOntologySig  <- GOOntologySig[order.pvalue] }

    # Compute number of genes in the gene set catagory:
    GeneTableFileName <- "GeneSetGenesTable"
    if (HasClasses) {
      NgenesGeneSet <- rep("", n)
      HeatmapGeneSet <- rep("", n)
      for(i in 1:n) {
          #6/8/07 Yong modified to separate gene tables out from gene-set table html file.
          #NgenesGeneSet[i] <- paste("<a href=#", GeneSetCode[i], ">", length(GeneSetList[[i]]), "</a>", sep="")
          FileNum <- ((i-1) %/% GeneSetUnit + 1)
          NgenesGeneSet[i] <- paste("<a href=\"", GeneTableFileName, FileNum, ".html#", GeneSetCode[i], "\">", length(GeneSetList[[i]]), "</a>", sep="")
          #add gene set heatmap link 2/14/08
          HeatmapGeneSet[i] <- paste("<a href=\"", "heatmap/geneset", i, ".html", "\" TARGET=_blank>", "heatmap", "</a>", sep="")
      }
    } else {
      NgenesGeneSet <- rep("", n)
      for(i in 1:n) {
          FileNum <- ((i-1) %/% GeneSetUnit + 1)
          NgenesGeneSet[i] <- paste("<a href=\"", GeneTableFileName, FileNum, ".html#", GeneSetCode[i], "\">", length(GeneSetList[[i]]), "</a>", sep="")
      }
    }
    # Create HTML table of gene set categories.
    fillNULL <- NULL
    DoKStest <- "pKS" %in% colnames(pGeneSet)
    GeneSetTable <- cbind(GeneSetId, GOOntologySig, GeneSetdescriptions, NgenesGeneSet, if (HasClasses) HeatmapGeneSet,
                      Trunc(pGeneSet$pLS, cutoff=alpha, TruncLevel=10^(-5)),
                      if(DoKStest) Trunc(pGeneSet$pKS, cutoff=alpha, TruncLevel=10^(-5)),
                      if(!GeneSetGSA) NULL else Trunc(pGeneSet$pGSA, cutoff=alpha, npermGSA=npermGSA, GSAscore=GSAscore, classlabels=classlabels),
                      if(!Hotelling) NULL else Trunc(pGeneSet$pHotelling, cutoff=alpha),
                      if(!GeneSetGlobalTest) NULL else Trunc(pGeneSet$pGlobal, cutoff=alpha))

    GeneSetTableNames <- c(GeneSetIdName, GOOntologyName, GeneSetDescriptionName, "Number of <NOBR>genes</NOBR>", if (HasClasses) "Heatmap link",
                      "LS permutation <NOBR>p-value</NOBR>",
                      if(DoKStest) "KS permutation <NOBR>p-value</NOBR>",
                      if(!GeneSetGSA) NULL else "Efron-Tibshirani's GSA test <NOBR>p-value</NORB>",
                      if(!Hotelling) NULL else "Hotelling test <NOBR>p-value</NOBR>",
                      if(!GeneSetGlobalTest) NULL else "Goeman's global test <NOBR>p-value</NOBR>")

    dimnames(GeneSetTable)[[2]] <- GeneSetTableNames
    GeneSetTable1 <- CreateHtmlTable(GeneSetTable)

    # Leonid's original code created GeneSetTable2 only when data is paired AND AddGeomeanColumns=T.
    GeneSetTable2 <- NULL

    # Create index of genes in gene set categories
    # sort genes in the gene table by the pvalues in each gene set
    # and also define 'Pvalues' as the ordered pvalues for signficiant genes.
    AllPvalues <- read.table(paste(WorkPath,"/Pvalues.txt", sep=""), header=T)[,2]
    GeneTablePartition <- rep(0, (n %/% GeneSetUnit + 1))
    GeneRowNum <- c(0)
    prevTableNum <- TableNum <- 0
    for(i in 1:n) {
      prevTableNum <- TableNum
      TableNum <- i %/% GeneSetUnit
      if(i == 1) {
        tmp <- GeneSetList[[i]]
        index <- tmp[order(AllPvalues[tmp])]
        genesGeneSetnames <- c(paste("<a name=", GeneSetCode[i], "></a>", GeneSetnames[i]), rep("", length(GeneSetList[[i]])-1) )
        genesGOTerms <- rep(GOOntologySig[i],length(GeneSetList[[i]]) )
        genesGeneSetdescriptions <- rep(GeneSetdescriptions[i],length(GeneSetList[[i]]) )
        GeneRowNum <- rep(1:length(GeneSetList[[i]]))
      } else {
        tmp <- GeneSetList[[i]]
        index <- c(index, tmp[order(AllPvalues[tmp])])
        genesGeneSetnames <- c(genesGeneSetnames,c(paste("<a name=", GeneSetCode[i], "></a>", GeneSetnames[i]),rep("", length(GeneSetList[[i]])-1) ))
        genesGOTerms <- c(genesGOTerms,rep(GOOntologySig[i],length(GeneSetList[[i]]) ))
        genesGeneSetdescriptions <- c(genesGeneSetdescriptions,rep(GeneSetdescriptions[i],length(GeneSetList[[i]]) ))
        GeneRowNum <- c(GeneRowNum, rep(1:length(GeneSetList[[i]])))
      }
      if (TableNum > prevTableNum)  GeneTablePartition[TableNum] <- length(index)
    }
    GeneTablePartition[TableNum + 1] <- length(index)
    Pvalues <- AllPvalues[index]
  } else {  # if not GeneSetComparison.
    GeneSetTable1 <- NULL
    GeneSetTable2 <- NULL
  }

  # THE REST OF THIS CODE CREATES CODE FOR THE HTML TABLE OF SIGNIFICANT GENES
  # (OR GENES CONTAINED WITHIN SIGNIFICANT GENE SETS, FOR GENE SET COMPARISON).

  # FOR BINARY TREE PREDICTION, CREATE GENE EXPRESSION DATA TABLE FOR EACH SUBTREE ???

  if (AnalysisType == "BinaryTreePrediction") {
    GeneExpSampleCode <- as.double(unlist(GeneExpData[1,])[-1])
    GeneExpGroupCode <- SampleCode[GeneExpSampleCode]
    GeneExpData <- GeneExpData[-1,]
    if( any(index != GeneExpData[[1]]) ) stop("Incorrect input: index != GeneExpData[[1]] ")
    GeneExpData <- GeneExpData[,-1]
  }

  # READ GENE ID TABLE FROM FILES.

  load(paste(WorkPath, "/GeneIdNames.rda", sep=""))
  NumGeneIds <- length(GeneIdNames)
  load(paste(WorkPath, "/GeneId1.rda", sep=""))
  NumGenesInAnalysis <- length(GeneId1)
  GeneId <- matrix("", nrow=NumGenesInAnalysis, ncol=NumGeneIds)
  GeneId[,1] <- GeneId1
  GeneId[,1] <- gsub("@#$%&", "'", GeneId[,1], fixed=T)
  GeneId[,1] <- gsub("@#$%!", "\"", GeneId[,1], fixed=T)
  if (NumGeneIds > 1) {
    for (i in 2:NumGeneIds) {
      load(paste(WorkPath, "/GeneId", i, ".rda", sep=""))
      GeneId[,i] <- get(paste("GeneId",i,sep=""))
      GeneId[,i] <- gsub("@#$%&", "'", GeneId[,i], fixed=T)
      GeneId[,i] <- gsub("@#$%!", "\"", GeneId[,i], fixed=T)
    }
  }

  ### MLI added for pairwise class comparison 12/24/08
  GeneId <- GeneId[index,,drop=F ]
  if (file.exists(paste(WorkPath, "/pairwiseCC.rda", sep=""))) {
     load(paste(WorkPath, "/pairwiseCC.rda", sep=""))
     GeneId <- cbind(GeneId, sapply(pairwiseCC, function(x) paste(x, collapse=", ")))
     GeneIdNames <- c(GeneIdNames, "Pairwise significant")
  }

  if (AnalysisType == "Survival") {
    # add KM curve hyperlinks to the end of GeneId
    GeneId <- cbind(GeneId, paste("<a href=kmcurve/gene", 1:length(index), ".html>Curve</a>", sep=""))
    GeneIdNames <- c(GeneIdNames, "Kaplan-Meier-<BR>curve")
    NumGeneIds <- NumGeneIds + 1
  }

  # IDENTIFY PROBESETS COLUMN, IF ANY.

  ProbeSets <- NA
  GeneIdColumn <- match("Probe set", GeneIdNames, nomatch=NA)
  if (!is.na(GeneIdColumn)) ProbeSets <- GeneId[,GeneIdColumn]

  # DETERMINE NUMBER OF CLASSES, IF APPLICABLE.

  if (HasClasses) {
    nclasses <- if (paired) 2 else length(unique(classlabels))
  }

  # SORT GENE ID TABLE (BY P- OR T-VALUES).

  ordercode <- NULL

  if ((AnalysisType != "ClassPrediction") & (AnalysisType != "SAM") & (AnalysisType != "PAM") & (AnalysisType != "SurvivalRisk") & (!GeneSetComparison)) {
    ordercode <- order(classifiers$pvalue)
    if (length(index) > 1) {
      GeneId <- GeneId[ordercode,,drop=F]
      classifiers <- classifiers[ordercode,]
    }
  } else if (GeneSetComparison) {
    ordercode <- 1:length(index)
  }
  if (AnalysisType == "SurvivalRisk") {
    if (!uselasso) {
		ordercode <- order(puniv)
		if (length(index) > 1) {
		   GeneId <- GeneId[ordercode,,drop=F]
		   puniv <- puniv[ordercode]
		   CVsupport <- CVsupport[ordercode]
		}
	}
  }

  if (AnalysisType == "ClassPrediction") {
    if(nclasses == 2 & SameN) {
      ordercode  <-  order(classifiers$t)
    } else {
      ordercode  <-  order(classifiers$pvalue)
    }
    if (length(index) > 1) {
      GeneId <- GeneId[ordercode,,drop=F]
      classifiers <- classifiers[ordercode,]
    }
  }

  # READ OTHER DATA COLUMNS FROM FILES.
  if (HasClasses) {
    if (paired) {
      load(paste(WorkPath, "/meandiff.rda", sep=""))
      meandiff <- meandiff[index]
      meandiff <- LogBase^meandiff
      if (AnalysisType != "SAM") meandiff <- meandiff[ordercode]
    } else {
      if (AnalysisType == "BinaryTreePrediction") {
        #if(IsSingleChannel) {
        #  load(paste(WorkPath, "/medianLR.rda", sep=""))
        #  medianLR <- medianLR[index]
        #} else {
        #  medianLR <- rep(0,length(GeneExpData[[1]]))
        #}
        #medianLR <- medianLR[ordercode]
        for( i in 1:2) {
          temp <- GeneExpData[,GeneExpGroupCode==i]
          if(any(temp < -9999998)) temp[temp < -9999998] <- NA
          #meanLR <- LogBase^(rowMeans(temp, na.rm=T)+medianLR)
          meanLR <- LogBase^(rowMeans(temp, na.rm=T))
          if(i == 1) {meanLRall <- meanLR} else { meanLRall <- cbind(meanLRall,meanLR)}
        }
      } else {
        for( i in 1:nclasses) {
          load(paste(WorkPath, "/meanLR",i,".rda", sep=""))
          meanLR <- LogBase^meanLR[index]
          if(i == 1) {meanLRall <- meanLR} else { meanLRall <- cbind(meanLRall,meanLR)}
        }
      }
      if (length(ordercode)!=0) {
         meanLRall <- meanLRall[ordercode,,drop=F]
      }
    }
  }
  # FOR CLASS PREDICTION PROCEDURES, CREATE CV SUPPORT COLUMNS.

  if (AnalysisType == "ClassPrediction") {
    nMethods <- length(methods[methods==1])
    CVsupport.selected <- round(100*CVsupport[,c(methods==1),drop=F])
    SameSupport <- T
    if(nMethods > 1) {
      for(i in 1:nMethods) {
        SameSupport <- SameSupport & all(CVsupport[,1] == CVsupport[,i])
      }
    }
    if (length(index) > 1) {
      CVsupport   <- CVsupport[ordercode,]
      CVsupport.selected   <- CVsupport.selected[ordercode,]
    }
    if(SameSupport) {
       CVsupport.selected <- CVsupport.selected[[1]]
       CVsupport.names <- "% CV <NOBR>support</NOBR>"
    } else {
      CVsupport.names <- c("% CV <NOBR>support</NOBR><NOBR>CCP</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>LDA</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>K=1</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>K=3</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>Centroid</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>SVM</NOBR>",
                      "% CV <NOBR>support</NOBR><NOBR>BCPP</NOBR>")
      CVsupport.names <- CVsupport.names[methods==1]
    }
  }

  if (AnalysisType == "BinaryTreePrediction") {
    CVsupport.selected <- round(100*classifiers$CVsupport,0)
    CVsupport.names <- "% CV <NOBR>support</NOBR>"
    CVsupport.selected <- CVsupport.selected[ordercode]
  }

  # OUTPUT SIGNIFICANT GENES TO GENELIST.

  # GeneListFileName is defined in RunXXXX.R (comment out calling OutputGeneList in batr package)
  # if (missing(genelistfilename)) genelistfilename <- GeneListFileName
  # if (!GeneSetComparison) {
  #    OutputGeneList(GenesVec=c(GeneId[,1]), GeneIdType=GeneIdNames[1], FileName=genelistfilename,
  #                   LocalToProject=T, ProjectPath = ProjectPath)
  #}

  # ADD GENE SET INFORMATION TO 'GeneId' TABLE, FOR GENE SET COMPARISON.

  # NOTE: If the number of gene set column(s) appended before the unique id in 'GeneId' is changed,
  #       then the GeneIdColumnToMatch parameter should also be modified correspondingly.

  if (GeneSetComparison) {
    if (AnalysisType == "GeneListComparison") {
      GeneId <- cbind(genesGeneSetnames,GeneId)
      NumGeneIds <- NumGeneIds + 1
      GeneIdNames <- c(GeneSetIdName, GeneIdNames)
      GeneIdColumnToMatch <- 2
    } else {
      GeneId <- cbind(genesGeneSetnames, genesGOTerms, genesGeneSetdescriptions, GeneId)
      if (is.null(genesGOTerms)) {
         NumGeneIds <- NumGeneIds + 2
      } else {
         NumGeneIds <- NumGeneIds + 3
      }
      GeneIdNames <- c(GeneSetIdName, GOOntologyName, GeneSetDescriptionName, GeneIdNames)
      if (is.null(genesGOTerms)) {
         GeneIdColumnToMatch <- 3
      } else {
         GeneIdColumnToMatch <- 4
      }
    }
  }

  # DETERMINE n (=NUMBER OF GENES).

  if (length(index) > 1) {
    n <- if (is.matrix(GeneId)) nrow(GeneId) else length(GeneId)
  } else {
    n <- 1
  }

  # CONCATENATE COLUMNS THAT BELONG IN HTML TABLE OF SIGNIFICANT GENES.
  # NOTE: ALL COLUMNS SHOULD ALREADY BE SORTED BY 'ordercode' EXCEPT SAM, PAM AND GENESET COMPARISON.

  GeneIdTemp <- GeneIdTempPair <- NULL
  GeneIdNamesTemp <- GeneIdNamesTempPair <- NULL

  if (AnalysisType == "Correlation") {
    GeneIdTemp <- cbind(GeneIdTemp, Trunc(classifiers$correlation, pos=F, TruncLevel=TruncLevel))
    GeneIdNamesTemp <- c(GeneIdNamesTemp, "Correlation <NOBR>coefficient</NOBR>")
  }

  if ((AnalysisType != "SAM") & (AnalysisType != "PAM") & (AnalysisType != "SurvivalRisk")) {
      if (!GeneSetComparison) Pvalues <- classifiers$pvalue
      GeneIdTemp <- cbind(GeneIdTemp, Trunc(Pvalues, TruncLevel=TruncLevel))
      GeneIdNamesTemp <- c(GeneIdNamesTemp, "Parametric <NOBR>p-value</NOBR>")
  }

  if (AnalysisType == "SAM") {
	# MLI added for v4.2 10/18/2010
      GeneIdTemp <- cbind(GeneIdTemp, Trunc(classifiers$di, pos=F, TruncLevel=TruncLevel))
      GeneIdNamesTemp <- c(GeneIdNamesTemp, "d(i)")
  }

  if (AnalysisType == "ClassPrediction" | AnalysisType == "BinaryTreePrediction") {
    GeneIdTemp <- cbind(GeneIdTemp, Trunc(classifiers$t, pos=F, TruncLevel=TruncLevel), CVsupport.selected)
    GeneIdNamesTemp <- c(GeneIdNamesTemp, "t-value", CVsupport.names)
  }

  if (AnalysisType == "SurvivalRisk") {
	if (!uselasso) {
		GeneIdTemp <- cbind(GeneIdTemp, Trunc(puniv, TruncLevel=TruncLevel), round(CVsupport,2))
		GeneIdNamesTemp <- c(GeneIdNamesTemp, "p-value", "% CV Support")
	} else {
		GeneIdTemp <- cbind(GeneIdTemp, Trunc(coxcoef, pos=F, TruncLevel=10^(-5)), round(CVsupport,2))
		GeneIdNamesTemp <- c(GeneIdNamesTemp, "Coefficient", "% CV Support")
	}
  }

  if (univtest) {
    longp <- PAllgenes
    fdr <- p.adjust(longp, "fdr")[classifiers$Index]
    GeneIdTemp <- cbind(GeneIdTemp, signifNew(fdr, digits=3)) #Qian 5/7/09
    GeneIdNamesTemp <- c(GeneIdNamesTemp,"FDR")
  }
  if (AnalysisType == "ClassComparison" && RunLFDR) {
    GeneIdTemp <- cbind(GeneIdTemp, signifNew(classifiers$lfdr, digits=3)) #MLI 2/19/2015
    GeneIdNamesTemp <- c(GeneIdNamesTemp,"<nobr>Local FDR</nobr>")
  }

  if(RunPermP) {
    GeneIdTemp <- cbind(GeneIdTemp, Trunc(classifiers$permpvalue, TruncLevel=TruncLevel))
    GeneIdNamesTemp  <- c(GeneIdNamesTemp,"Permutation <NOBR>p-value</NOBR>")
  }

  if (HasClasses) {
    if (paired) {
      DoSwap <- T
      if(is.null(RefClass)) {
        DoSwap <- F
      } else {
        if(as.character(RefClass) == as.character(classlabels[2]) ) {
          DoSwap <- F
        }
      }
      if(DoSwap) {
        GeneIdTemp  <- cbind(GeneIdTemp, signifNew(1/meandiff, digits=2)) #Qian 5/7/09
        GeneIdNamesTemp <- c(GeneIdNamesTemp,paste("Geometric mean of ",ValueDef,"<BR>(class ",classlabels[2],
                         "/class ",classlabels[1],")",collapse=""))
      } else {
        GeneIdTemp  <- cbind(GeneIdTemp, signifNew(meandiff, digits=2)) #Qian 5/7/09
        GeneIdNamesTemp <- c(GeneIdNamesTemp,paste("Geometric mean of ",ValueDef,"<BR>(class ",classlabels[1],
                         "/class ",classlabels[2],")",collapse=""))
      }
    } else {
      if (AnalysisType == "BinaryTreePrediction") {
        for(i in 1:2) {  #Print geometric means for the two branches of the node.
          GeneIdTemp <- cbind(GeneIdTemp, signifNew(meanLRall[,i], digits=2)) #Qian 5/7/09
          GeneIdNamesTemp <- c(GeneIdNamesTemp,paste("Geom mean of ",ValueDef," in group ", i, collapse=""))
        }
      } else {
        for(i in 1:nclasses) {
          GeneIdTemp <- cbind(GeneIdTemp, signifNew(meanLRall[,i], digits=2)) #Qian 5/7/09
          GeneIdNamesTemp <- c(GeneIdNamesTemp,paste("Geom mean of ",ValueDef," in class ", i, collapse="")) #modified by Yong 12/15/06
        }
      }
      if(nclasses==2) {
        GeneIdTemp <- cbind(GeneIdTemp, signifNew(meanLRall[,1]/meanLRall[,2], digits=2)) #Qian 5/7/09
        GeneIdNamesTemp <- c(GeneIdNamesTemp,"Fold-change") # changed from "Ratio of geom means" 1/15/08
      }
	  if (AnalysisType == "PAM") { # Ming-chung 8/5/2011
	     for(i in 1:nclasses) {
	       GeneIdTemp <- cbind(GeneIdTemp,
		                       ifelse(round(xikprime[,i], 2)==0, format(xikprime[, i], digits=3), round(xikprime[,i], 2)))
		   GeneIdNamesTemp <- c(GeneIdNamesTemp, paste("Shrunken centroid in class ", i, collapse=""))
		 }
		 GeneIdTemp <- cbind(GeneIdTemp, signifNew(si, digits=2))
		 GeneIdNamesTemp <- c(GeneIdNamesTemp, "Standard deviation s<sub>i</sub>+s<sub>0</sub>")
	  }
	}
    #Not sure why this is necesary, but Amy copied it over from original PreparePredictionTableHTML.R file when merging.  6/29/06.
    if (AnalysisType == "ClassPrediction") {for(ii in 1:length(GeneId)) GeneId[[ii]] <- as.character(GeneId[[ii]])}
  }

  if (AnalysisType == "Survival") {
    LogTruncLevel <- log(TruncLevel)
    HazardRatio <- classifiers$b
    HazardRatio <- ifelse(HazardRatio < LogTruncLevel,paste("< ",TruncLevel),
                          ifelse(HazardRatio > -LogTruncLevel,paste("> ",1/TruncLevel),
                                 round(exp(HazardRatio),3) ))
    GeneIdTemp <- cbind(GeneIdTemp, HazardRatio, round(classifiers$SD,3))
    GeneIdNamesTemp  <- c(GeneIdNamesTemp,"Hazard <NOBR>Ratio</NOBR>",paste("SD of log<NOBR>",ValueDef,"</NOBR>"))
  }

  if (GeneSetComparison) {
    GeneId  <- cbind(GeneId, GeneIdTemp)
    GeneIdNamesTemp <- c(GeneIdNames, GeneIdNamesTemp)
  } else {
    GeneId  <- cbind(GeneIdTemp, GeneId)
    GeneIdNamesTemp <- c(GeneIdNamesTemp, GeneIdNames)
    GeneIdColumnToMatch <- ncol(GeneIdTemp) + 1
  }

  GeneId <- as.data.frame(GeneId)
  names(GeneId) <- GeneIdNamesTemp

  # CALL FUNCTION TO CREATE HTML TABLE OF SIGNIFICANT GENES (INCLUDING "INFO" LINK),
  # CHROMOSOMAL DISTRIBUTION BARPLOT, AND GO ANALYSIS (IF APPLICABLE).

  #6/8/07 Yong modified to probihit showing Chromsome bar plot when doing gene set comparison.
  #RunChromDistBarPlot <- if (Hotelling) F else T
  # RunChromDistBarPlot <- if (GeneSetComparison) F else T
  RunChromDistBarPlot <- FALSE # batr package
  GeneId1B <- GeneId1 # debug use
  TableID <- if(file.exists(paste(WorkPath, "/pairwiseCC.rda", sep=""))) "id_of_table" else NULL
  if (!GeneSetComparison) {
    GeneId1=CreateGeneTable(GeneTable=GeneId
          , AllGenesInAnalysis=GeneId1B
          , InsertAfterColumn=GeneIdColumnToMatch
          , RunChromDistBarPlot=RunChromDistBarPlot
          , RunGOOvE=DoGOOvE
          , MinObserved=MinObserved
          , MinRatio=MinRatio
          , RedoExpectedTable=RedoExpectedTable
          , FolderPath=paste(ProjectPath, "/Output/", AnalysisName, sep="")
          , ReturnHtmlCode=T
          , GeneTableFile=""
          , AnnotationsTableFile="Annotations.html"
          , DoMasterAnnotations=DoMasterAnnotations
          , ChromDistBarPlotFile=paste(AnalysisName, "-ChromDist.jpg")
          , HeaderStringsVec=""
          , TableCaption=""
          , GeneIdColumnToMatch=GeneIdColumnToMatch
          , TableID=TableID
          , ProjectPath = ProjectPath)

    # Specifically designed for pairwise comparison.
    if (file.exists(paste(WorkPath, "/pairwiseCC.rda", sep=""))) {
       ncolgenetable <- ncol(GeneId)+1*GeneId1$DoAnnotations
       save(ncolgenetable, file=paste(WorkPath, "/ncolgenetable.rda", sep=""))
       # the final number of columns may be one more (annotation column) if annotation has been done successfully
    }

    OEhtml <- GeneId1$GOTableHtmlCode
    GeneId1 <- GeneId1$GeneTableHtmlCode
  } else {
    #6/8/07 generate separate gene table html files for gene set comparison analysis.
    GeneId1 <- NULL
    for (i in 1:length(GeneTablePartition)) {
      if (i==1) {
        startRow <- 1
        endRow <- GeneTablePartition[i]
      } else {
        startRow <- GeneTablePartition[i-1] + 1
        endRow <- GeneTablePartition[i]
      }
      GeneIdUnit <- GeneId[startRow:endRow,]
      GeneRowNumUnit <- GeneRowNum[startRow:endRow]
      GeneTableFile <- paste(GeneTableFileName, i, ".html", sep="")
      CreateGeneTable(GeneTable=GeneIdUnit
            , AllGenesInAnalysis=NULL
            , InsertAfterColumn=GeneIdColumnToMatch
            , RunChromDistBarPlot=RunChromDistBarPlot
            , RunGOOvE=DoGOOvE
            , MinObserved=MinObserved
            , MinRatio=MinRatio
            , RedoExpectedTable=RedoExpectedTable
            , FolderPath=paste(ProjectPath, "/Output/", AnalysisName, sep="")
            , ReturnHtmlCode=F
            , GeneTableFile=GeneTableFile
            , AnnotationsTableFile=""
            , DoMasterAnnotations=F
            , ChromDistBarPlotFile=""
            , HeaderStringsVec=""
            , TableCaption=GeneTableCaption
            , GeneIdColumnToMatch=GeneIdColumnToMatch
            , GeneRowNum=GeneRowNumUnit)

      GeneId1 <- c(GeneId1, ifelse(is.null(GeneId1), "", " ... "),
                  paste("<a href=\"", GeneTableFile, "\">",
                  "Gene Sets ", (i-1)*GeneSetUnit+1, " - ", ifelse(i<length(GeneTablePartition), i*GeneSetUnit, length(GeneSetIndex)), "</a>", sep=""))
    }
    OEhtml <- NULL

	  # Color code KEGG pathway
	  #Generate URL to color code up/down regulated genes for KEGG pathway by using KEGG API  Qihao 8/15/2013
	  ColorCodeUseKEGGAPI <- function(KEGGpathwayID,GeneIdentifersList,GeneExpressionToComputeFoldChange=NULL,Pvalues) {
  		#save(list=c('KEGGpathwayID','GeneIdentifersList','GeneExpressionToComputeFoldChange','Pvalues'),file='c:/temp/test')
  		#Qihao modified codes 7/22/2015.
  		if (!(is.null(GeneExpressionToComputeFoldChange))){
  		GeneExpressionToComputeFoldChange<-sapply(GeneExpressionToComputeFoldChange,function(x){as.double(levels(x)[x])})
  		}
  		#Qihao added codes to handle genes without gene symbol.
  		if (all((GeneIdentifersList==''))|all(is.na(GeneIdentifersList))){
  			return('No graph available')
  		}
  		if(any((GeneIdentifersList==''|is.na(GeneIdentifersList)))){
  			MisssymbolGenes<-which(((GeneIdentifersList!='')&(!(is.na(GeneIdentifersList)))))
  			GeneIdentifersList<-GeneIdentifersList[MisssymbolGenes]
  			if (!(is.null(GeneExpressionToComputeFoldChange))){GeneExpressionToComputeFoldChange<-GeneExpressionToComputeFoldChange[MisssymbolGenes,,drop=FALSE]}
  			Pvalues<-Pvalues[MisssymbolGenes]
  		}
  		Pvalues<-sapply(as.character(levels(Pvalues)[Pvalues]),function(x){if(x=='< 1e-07'){as.numeric(x=0.0000000)}else{as.numeric(x)}})
  		index<-(Pvalues<=0.01)
  		GeneIdentifersList<-(GeneIdentifersList[index])

  		if (sum(index)==0){
  			#URLstring<-paste('<A href=\"http://www.kegg.jp/kegg-bin/show_pathway?map=',KEGGpathwayID,'\">','View','</A>',sep='')
  			return('No graph available')
  		} else {
    		if (is.null(GeneExpressionToComputeFoldChange)) {
    			URLstring<-paste(lapply(GeneIdentifersList,function(x){paste(x,'+%23D74B4B,',sep='')}),collapse='/')
    			#GeneSetdescriptions <- paste("<A HREF=\"http://www.kegg.jp/kegg-bin/show_pathway?", GeneSetId, "\">", GeneSetdescriptions, "</A>", sep="")
    		} else {
    			#Assume GeneExpressionToComputeFoldChange must be a matrix having only two columns.
    			GeneExpressionToComputeFoldChange<-GeneExpressionToComputeFoldChange[index,,drop=FALSE]
    			#Qihao modified codes. 7/22/2015
    			UpDownRegulatedIndex<-apply(GeneExpressionToComputeFoldChange,1,function(x){if(x[1]==x[2]){c(1)}else{if(x[1]>x[2]){c(2)}else{c(3)}}})

    			if (any(UpDownRegulatedIndex==2)){
    			URLstring<-paste(lapply(GeneIdentifersList[(UpDownRegulatedIndex==2)],function(x){paste(x,'+%23D74B4B,',sep='')}),collapse='/')
    			}
    			if (any(UpDownRegulatedIndex==3)){
    			DownURLstring<-paste(lapply(GeneIdentifersList[(UpDownRegulatedIndex==3)],function(x){paste(x,'+%233399CC,',sep='')}),collapse='/')
    			if(exists('URLstring')){URLstring<-paste(URLstring,DownURLstring,sep='/')}else{URLstring<-DownURLstring}
    			}
    			if (any(UpDownRegulatedIndex==1)){
    			EqualURLstring<-paste(URLstring,lapply(GeneIdentifersList[(UpDownRegulatedIndex==1)],function(x){paste(x,'+%23CDB99C,',sep='')}),collapse='/')
    			if(exists('URLstring')){URLstring<-paste(URLstring,EqualURLstring,sep='/')}else{URLstring<-EqualURLstring}
    			}
    			#save(list=c('URLstringUp','URLstringDown','URLstringNormal','UpDownRegulatedIndex'),file='c:/temp/test')
    			#URLstring<-paste(URLstringUp,URLstringDown,URLstringNormal,collapse='/')
    		}

    		URLstring<-paste('<A href=\"http://www.kegg.jp/kegg-bin/show_pathway?map=',KEGGpathwayID,'&multi_query=',URLstring,'\">','View','</A>',sep='')
    		return(URLstring)
  		} # end of if (sum(index) == 0)
	  } # end of ColorCodeUseKEGGAPI()

    #save(list=c('GeneId1','GeneSetTable','GeneId','GeneSetList','GeneSetId','GeneSetdescriptions','GeneIdNames'),file='c:/temp/qqk')
    # save.image(file='c:/temp/qqy')
	  if (AnalysisType == "PathwayComparison") {
	    #Qihao fixed #515
	    if((PathwayType == "Kegg")&&(!(paired))) {
		    KEGG.colorcode.pathwayName<-GeneSetDescriptionName
	    	j<-0
		    KEGG.color.Position<-unlist(lapply(GeneSetList,length))

		    if (length(unique(classlabels)) == 2){

			    for (i in 1:length(GeneSetId)){
    				#List<-GeneId[1+j:length(GeneSetList[[i]]),grep('*symbol',GeneIdNames,ignore.case=T)]
    				KEGG.color.Begin<-1+j
    				KEGG.color.End<-j+KEGG.color.Position[i]
    				#KEGG.colorcode.pathwayName[i]<-ColorCodeUseKEGGAPI(GeneSetId[i],GeneId[KEGG.color.Begin:KEGG.color.End,grepl('*symbol',colnames(GeneId),ignore.case=T)],subset(GeneId[KEGG.color.Begin:KEGG.color.End,(grepl('Geom mean',colnames(GeneId)))]),GeneId[KEGG.color.Begin:KEGG.color.End,grepl('p-value',colnames(GeneId),ignore.case=T)])
    				KEGG.colorcode.pathwayName[i]<-ColorCodeUseKEGGAPI(GeneSetId[i],GeneId[KEGG.color.Begin:KEGG.color.End,grepl('*EntrezID',colnames(GeneId),ignore.case=T)],subset(GeneId[KEGG.color.Begin:KEGG.color.End,(grepl('Geom mean',colnames(GeneId)))]),GeneId[KEGG.color.Begin:KEGG.color.End,grepl('p-value',colnames(GeneId),ignore.case=T)]) ## Ting used EntrezID to generate KEGG pathway web links with map coloring 4/18/2017
    				j<-KEGG.color.End
    			}
		    } else {
			    for (i in 1:length(GeneSetId)){
    				#List<-GeneId[1+j:length(GeneSetList[[i]]),grep('*symbol',GeneIdNames,ignore.case=T)]
    				KEGG.color.Begin<-1+j
    				KEGG.color.End<-j+KEGG.color.Position[i]
    				#KEGG.colorcode.pathwayName[i]<-ColorCodeUseKEGGAPI(GeneSetId[i],GeneId[KEGG.color.Begin:KEGG.color.End,grepl('*symbol',colnames(GeneId),ignore.case=T)],as.null(),GeneId[KEGG.color.Begin:KEGG.color.End,grepl('p-value',colnames(GeneId),ignore.case=T)])
    				KEGG.colorcode.pathwayName[i]<-ColorCodeUseKEGGAPI(GeneSetId[i],GeneId[KEGG.color.Begin:KEGG.color.End,grepl('*EntrezID',colnames(GeneId),ignore.case=T)],as.null(),GeneId[KEGG.color.Begin:KEGG.color.End,grepl('p-value',colnames(GeneId),ignore.case=T)]) ## Ting used EntrezID to generate KEGG pathway web links with map coloring 4/18/2017
    				j<-KEGG.color.End
    			}
		    }
		    if (!(is.null(GeneSetTable1))) {
			    GeneSetTable <- cbind(GeneSetId, GOOntologySig, GeneSetdescriptions,KEGG.colorcode.pathwayName,NgenesGeneSet, if (HasClasses) HeatmapGeneSet,
							  Trunc(pGeneSet$pLS, cutoff=alpha, TruncLevel=10^(-5)),
							  if(DoKStest) Trunc(pGeneSet$pKS, cutoff=alpha, TruncLevel=10^(-5)),
							  if(!GeneSetGSA) NULL else Trunc(pGeneSet$pGSA, cutoff=alpha, npermGSA=npermGSA, GSAscore=GSAscore, classlabels=classlabels),
							  if(!Hotelling) NULL else Trunc(pGeneSet$pHotelling, cutoff=alpha),
							  if(!GeneSetGlobalTest) NULL else Trunc(pGeneSet$pGlobal, cutoff=alpha))
			    GeneSetTableNames <- c(GeneSetIdName, GOOntologyName, GeneSetDescriptionName
							  ,'Color Coded KEGG Pathway Graphs'
							  ,"Number of <NOBR>genes</NOBR>", if (HasClasses) "Heatmap link",
							  "LS permutation <NOBR>p-value</NOBR>",
							  if(DoKStest) "KS permutation <NOBR>p-value</NOBR>",
							  if(!GeneSetGSA) NULL else "Efron-Tibshirani's GSA test <NOBR>p-value</NORB>",
							  if(!Hotelling) NULL else "Hotelling test <NOBR>p-value</NOBR>",
							  if(!GeneSetGlobalTest) NULL else "Goeman's global test <NOBR>p-value</NOBR>")
			    dimnames(GeneSetTable)[[2]] <- GeneSetTableNames
			    GeneSetTable1 <- CreateHtmlTable(GeneSetTable)
		    }
	    } # end of if((PathwayType == "Kegg")&&(!(paired)))
	  } # end of if (AnalysisType == "PathwayComparison")
  } # end of if (!GeneSetComparison)

  # For greedy pair method in class prediction we output gene table sorted by gene pairs             ##
  if (AnalysisType == "ClassPrediction" && !is.null(genepairIndex)) {
        tmp <- GeneId[match(genepairIndex, index[ordercode]),,drop = F]
        tmp <- cbind(data.frame(Pair = rep(1:(nrow(tmp)/2), each =2)), tmp)
        GeneIdPair=CreateGeneTable(GeneTable=tmp
                   , InsertAfterColumn=GeneIdColumnToMatch+1
                   , RunChromDistBarPlot=F
                   , RunGOOvE=F
                   , MinObserved=MinObserved
                   , MinRatio=MinRatio
                   , RedoExpectedTable=RedoExpectedTable
                   , FolderPath=paste(ProjectPath, "/Output/", AnalysisName, sep="")
                   , ReturnHtmlCode=T
                   , GeneTableFile=""
                   , AnnotationsTableFile="Annotations.html"
                   , DoMasterAnnotations=DoMasterAnnotations
                   , ChromDistBarPlotFile=""
                   , HeaderStringsVec=""
                   , TableCaption=""
                   , GeneIdColumnToMatch=GeneIdColumnToMatch+1)$GeneTableHtmlCode
  } else GeneIdPair <- NULL

  # FOR QUANTITATIVE TRAITS OR ANY PAIRED ANALYSIS, RE-SORT GENES AND CREATE A SECOND TABLE OF SIGNIFICANT GENES.
  if (paired) {
    if(length(meandiff) > 1) {
      GeneId <- GeneId[rev(order(meandiff)), ,drop=F]
    }
  }

  if (AnalysisType == "Correlation") {
    if(length(classifiers$correlation) > 1) {
      GeneId <- GeneId[rev(order(classifiers$correlation)), ,drop=F]
    }
  }

#  if ((paired)|(AnalysisType == "Correlation")) {  change for batr package
  if (AnalysisType == "Correlation") {
    GeneId2=CreateGeneTable(GeneTable=GeneId
      , AllGenesInAnalysis=GeneId1B
      , InsertAfterColumn=GeneIdColumnToMatch
      , RunChromDistBarPlot=T
      , RunGOOvE=DoGOOvE
      , MinObserved=MinObserved
      , MinRatio=MinRatio
      , RedoExpectedTable=RedoExpectedTable
      , FolderPath=paste(ProjectPath, "/Output/", AnalysisName, sep="")
      , ReturnHtmlCode=T
      , GeneTableFile=""
      , AnnotationsTableFile="Annotations.html"
      , DoMasterAnnotations=DoMasterAnnotations
      , ChromDistBarPlotFile=paste(AnalysisName, "-ChromDist.jpg")
      , HeaderStringsVec=""
      , TableCaption=""
      , GeneIdColumnToMatch=GeneIdColumnToMatch
    )
    GeneId2 <- GeneId2$GeneTableHtmlCode
  } else {
    GeneId2 <- NULL
  }

  # RETURN STRINGS CONTAINING HTML CODE TO BE INSERTED INTO HTML OUTPUT FILES, AND OTHER PARAMETERS.
  out <- list(GeneSetTable1=GeneSetTable1, GeneSetTable2=GeneSetTable2,
       GeneId1=GeneId1, GeneId2=GeneId2, OEhtml=OEhtml, GeneIdPair=GeneIdPair,
       NumGenesInAnalysis=NumGenesInAnalysis, ProbeSets=ProbeSets
       , GeneIndex=index # return significant gene index, Yong 8/23/06
       , GONumTotCat=GONumTotCat, GONumSigCat=GONumSigCat # return number of total investigated and significant GO categories, Yong 4/2/07
       ,GeneIdForIngenuity=GeneId #Qihao modified this.
	   )
  invisible(out)
} # end of GetGenesHTML()


CreateGeneTable <- function(
  GeneTable
  , AllGenesInAnalysis=NULL
  , InsertAfterColumn=1
  , RunChromDistBarPlot=T
  , RunGOOvE=T
  , MinObserved=5
  , MinRatio=2
  , RedoExpectedTable=T
  , FolderPath=""
  , ReturnHtmlCode=T
  , GeneTableFile=""
  , AnnotationsTableFile=""
  , DoMasterAnnotations=T
  , ChromDistBarPlotFile=""
  , HeaderStringsVec="Output of plugin analysis"
  , TableCaption=""
  , ColumnNum
  , GeneIdColumnToMatch=1
  , GeneRowNum=NULL
  , TableID=NULL
  , ProjectPath = NULL
) {
## INPUT:
##
## (REQUIRED PARAMETERS:)
##   GeneTable          = data frame containing genes, where first column is the gene id used in first column of 'Gene identifiers' worksheet.
##                Columns with names matching specific keywords ('Clone', 'Probe set', 'GB acc', 'UG cluster', or 'Gene symbol')
##                will also be hyperlinked in the resulting HTML gene table.  Periods in columns names will be replaced by spaces before matching.
##   AllGenesInAnalysis     = REQUIRED ONLY IF RunGOOvE IS TRUE.  Vector of gene ids for ALL genes used in analysis.
##
## (OPTIONAL PARAMETERS:)
##   InsertAfterColumn      = number giving column after which annotations link will be inserted.
##   RunChromDistBarPlot    = Boolean denoting whether or not to run the chromosomal distribution barplot.
##   RunGOOvE           = Boolean denoting whether or not to run the Gene Ontology observed versus expected table.
##   MinObserved        = Used only if RunGOOvE is TRUE.  Minimum observed value required in order for GO term to be printed in table.
##   MinRatio           = Used only if RunGOOvE is TRUE.  Minimum observed to expected ratio required in order for GO term to be printed in table.
##   RedoExpectedTable  = Used only if RunGOOvE is TRUE.  If RedoExpectedTable is TRUE, then compute expected values from scratch, instead of using saved values from file.
##   ReturnHtmlCode     = Boolean denoting whether function should return the HTML code, or write the HTML output to GeneTableFile.
##                NOTE: The HTML code returned by this function must be written to an HTML file contained in the folder specified
##                by FolderPath.  Otherwise the links to AnnotationsTableFile and ChromDistBarPlotFile will be broken!
##   FolderPath         = string denoting path of folder in which all output files will be written.
##   GeneTableFile      = string denoting filename for user-provided gene table into which annotations link will be inserted.
##   AnnotationsTableFile   = string denoting filename for annotations table in a separate file.
##   DoMasterAnnotations    = Boolean denoting whether or not to create one master annotations table for multiple gene tables in the HTML output. Yong 8/21/06
##   ChromDistBarPlotFile   = string denoting filename for chromosomal distribution barplot, if RunChromDistBarPlot is TRUE.
##   HeaderStringsVec       = vector of strings which contain lines of text to be inserted into GeneTableFile before the actual gene table.
##   TableCaption       = string to use as a caption for gene table
##   ColumnNum          = number giving column of numerical values. So numerical values will not be wrapped in html output. It is assumed columns from ColumnNum to ncol(GeneTable) are all numeric.
##   GeneIdColumnToMatch    = Column number of GeneTable which contains the column to match (e.g., primary id).
##
## OUTPUT:
##
##   If ReturnHtmlCode is TRUE, then the returned value of the function is a list containing components named
##   GeneTableHtmlCode, GOTableHtmlCode, and ChromDistHtmlCode.  Each component is a vector of strings which can
##   be printed to an HTML file using the 'cat' function in R.  If RunGOOvE is FALSE, then the GOTableHtmlCode
##   component is an empty string ("").  Similarly, if RunChromDistBarPlot is FALSE, then the ChromDistHtmlCode
##   is an empty string ("").  The HTML code in the GeneTableHtmlCode component contains links to AnnotationsTableFile,
##   and the HTML code in the ChromDistHtmlCode contains a link to ChromDistBarPlotFile if RunChromDistBarPlot is TRUE.
##   These links are relative links, so it is important that the HTML code returned by this function is written to an
##   HTML file which is contained in the folder specified by the FolderPath parameter.
##
##   If ReturnHtmlCode is FALSE, then GeneTableFile will be created in the folder specified by FolderPath, and will
##   contain relative links to AnnotationsTableFile and ChromDistBarPlotFile (if RunChromDistBarPlot
##   is TRUE).
##
## Called by GetGenesHTML().

  ArrayToolsPath <- path.package("classpredict")

  #Set default path and filenames if user did not specify any.
  if (FolderPath=="") FolderPath <- paste(ProjectPath, "/Output/Plugins", sep="")
  if (GeneTableFile=="") GeneTableFile <- "GeneTableFile.html"
  if (AnnotationsTableFile=="") AnnotationsTableFile <- "Annotations.html"
  if (ChromDistBarPlotFile=="") ChromDistBarPlotFile <- "ChromDist.jpg"

  AnnotationsFileName <- AnnotationsTableFile
  ChromDistFileName <- ChromDistBarPlotFile

  GeneTableFile <- paste(FolderPath, GeneTableFile, sep="/")
  AnnotationsTableFile <- paste(FolderPath, AnnotationsTableFile, sep="/")
  ChromDistBarPlotFile <- paste(FolderPath, ChromDistBarPlotFile, sep="/")

  #Make sure the user did not specify duplicate filenames.
  if (GeneTableFile==AnnotationsTableFile) {
    ErrMsg <- "ERROR: GeneTableFile and AnnotationsTableFile have the same filename."
    require(tcltk)
    tkmessageBox(message=ErrMsg,icon="info",type="ok")
    cat(ErrMsg)
  }
  if (RunChromDistBarPlot & (GeneTableFile==ChromDistBarPlotFile)) {
    ErrMsg <- "ERROR: GeneTableFile and ChromDistBarPlotFile have the same filename."
    require(tcltk)
    tkmessageBox(message=ErrMsg,icon="info",type="ok")
    cat(ErrMsg)
  }
  if (RunChromDistBarPlot & (AnnotationsTableFile==ChromDistBarPlotFile)) {
    ErrMsg <- "ERROR: AnnotationsTableFile and ChromDistBarPlotFile have the same filename."
    require(tcltk)
    tkmessageBox(message=ErrMsg,icon="info",type="ok")
    cat(ErrMsg)
  }

  #Create FolderPath if it doesn't yet exist.
  CreatePath(FolderPath)

  GeneId <- as.matrix(GeneTable)
  GeneIdsMatch <- c(as.character(GeneId[,GeneIdColumnToMatch]))

  ## (Need to check if annotations functions will automatically convert to character during matching??)
  ## Also need to test on case where GeneTable has only one row or one column.
  ##  ==> One col is OK, and one row with one col is OK, but one row with more than one col does not hyperlink.

  ColNames <- dimnames(GeneId)[[2]]
  ColNames <- gsub(pattern="\\.", replacement=" ", x=ColNames)
  NumGeneIds <- ncol(GeneId)

  ## TRY TO CREATE ANNOTATIONS:

  options(warn=2)  #Converts warnings to errors in try function.
#  cat("Attempting to create annotations ...\n")
  try.out <- try ( annotations <- CreateAnnotations(GeneIds=GeneIdsMatch, PathName=FolderPath, AnnotationsFileName=AnnotationsFileName,
                                                    CallFromPlugins=T, DoMasterAnnotations=DoMasterAnnotations,
                                                    ProjectPath = ProjectPath) ) # add DoMasterAnnotations, Yong 8/21/06
  DoAnnotations <- F
  if ("try.out" %in% ls()) {   #If "try.out" object exists.
    if (!is.null(class(try.out))) {  #Necessary only for R1.6.1, which evaluates class(NULL) and class(NA) as NULL whereas R1.7.0 evaluates it as "NULL" and "logical" respectively.
      if (class(try.out) != "try-error") {
        if (!is.null(annotations)) {
          DoAnnotations <- T
        }
      } else { #Yong added testing try-error 1/4/07
        require(tcltk)
        tkmessageBox(message=paste('Error occurred in Annotations.R:', geterrmessage(), sep='\n'),icon="info",type="ok")
        cat(geterrmessage());
        stop();
      }
    } else {  #In R1.6.1, class(try.out) evaluates as NULL even if function evaluates correctly and try.out is a valid non-null object.
      if ("annotations" %in% ls()) {
        if (!is.null(annotations)) {
          DoAnnotations <- T
        }
      }
    }
    rm(try.out)
  }

  if (DoAnnotations) {
    NumRows <- nrow(GeneId)
    NumCols <- ncol(GeneId)
    ColNames <- dimnames(GeneId)[[2]]
    ColNames <- gsub(pattern="\\.", replacement=" ", x=ColNames)

    if (InsertAfterColumn < 1) {
      cat("Error: InsertAfterColumn parameter must be >= 1.")
    } else {
      if (InsertAfterColumn > NumCols) {
        cat(paste("Error: InsertAfterColumn parameter must be <= ", NumCols, ".", sep=""))
      }
    }

    if (!missing(ColumnNum)) ColumnNum <- ColumnNum + 1 # because of annotations column is added to GeneId
    if (InsertAfterColumn == NumCols) {  #INSERT ANNOTATIONS INTO LAST COLUMN.
      ColNames <- c(ColNames, "Anno-<BR>tations")
      if (NumRows > 1) {
        GeneId <- cbind(GeneId, annotations)
      } else {    #For some reason, cbind() doesn't work properly when there is only one row.
        GeneId <- c(GeneId, annotations)
        GeneId <- matrix(GeneId, nrow=1)
      }
      dimnames(GeneId)[[2]] <- ColNames
    } else {
      ColNames <- c(ColNames[1:InsertAfterColumn], "Anno-<BR>tations", ColNames[(InsertAfterColumn+1):NumCols])
      if (NumRows > 1) {
        GeneId <- cbind(GeneId[,1:InsertAfterColumn], annotations, GeneId[,(InsertAfterColumn+1):NumCols])
      } else {    #For some reason, cbind() doesn't work properly when there is only one row.
        GeneId <- c(GeneId[1:InsertAfterColumn], annotations, GeneId[(InsertAfterColumn+1):NumCols])
        GeneId <- matrix(GeneId, nrow=1)
      }
      dimnames(GeneId)[[2]] <- ColNames
    }
    NumCols <- NumCols+1
    NumGeneIds <- ncol(GeneId)
  } else {
    # cat("==> 'CreateAnnotations' function could not be run.\n")  modified for batr package
  } # end of if (DoAnnotations)

  ## TRY TO CREATE CHROMOSOMAL DISTRIBUTION BARPLOT:

  if (RunChromDistBarPlot) {
    AnalysisName <- ""
    cat("Attempting to create chromosomal distribution barplot ...\n")
    try.out <- try ( ChromDistBarPlot(GeneIdsMatch, ProjectPath, AnalysisName, ChromDistBarPlotFile, AllGenesInAnalysis) )
    DoChromDistBarPlot <- F
    if ("try.out" %in% ls()) {   #If "try.out" object exists.
      if (!is.null(class(try.out))) {  #Necessary only for R1.6.1, which evaluates class(NULL) and class(NA) as NULL whereas R1.7.0 evaluates it as "NULL" and "logical" respectively.
        if (class(try.out) != "try-error") {
          if (file.access(ChromDistBarPlotFile) == 0) {
            DoChromDistBarPlot <- T
          }
        }
      } else {  #In R1.6.1, class(try.out) evaluates as NULL even if function evaluates correctly and try.out is a valid non-null object.
        if (file.access(ChromDistBarPlotFile) == 0) {
          DoChromDistBarPlot <- T
        }
      }
      rm(try.out)
    }
    if (!DoChromDistBarPlot) {
      cat("==> 'ChromDistBarPlot' function could not be run.\n")
    }
  }
  ## TRY TO RUN GENE ONTOLOGY ANALYSIS:
  if(RunGOOvE){
    cat("Attempting to run Gene Ontology analysis ...\n")
    try.out <- try ( OEhtml <- GOanalysis(GeneIdsMatch, ProjectPath, MinObserved, MinRatio) )
    DoGOanalysis <- F
    if ("try.out" %in% ls()) {   #If "try.out" object exists.
      if (!is.null(class(try.out))) {  #Necessary only for R1.6.1, which evaluates class(NULL) and class(NA) as NULL whereas R1.7.0 evaluates it as "NULL" and "logical" respectively.
        if (class(try.out) != "try-error") {
          DoGOanalysis <- T
          if (length(OEhtml)==1) {
            if (!is.na(OEhtml)) {
              if (OEhtml=="") {
                DoGOanalysis <- F
              }
            }
          }
        }
      } else {  #In R1.6.1, class(try.out) evaluates as NULL even if function evaluates correctly and try.out is a valid non-null object.
        if ("OEhtml" %in% ls()) {
          DoGOanalysis <- T
          if (length(OEhtml)==1) {
            if (!is.na(OEhtml)) {
              if (OEhtml=="") {
                DoGOanalysis <- F
              }
            }
          }
        }
      }
      rm(try.out)
    }
    if (!DoGOanalysis) {
      cat("==> 'GOanalysis' function could not be run.\n")
      OEhtml <- ""
    }
  }else{
    OEhtml <- ""
  }

  options(warn=0)  #Set back to default.
#  cat("Preparing output tables ...\n")

  ## HYPERTEXT CLONE IDS.
  indic <- (ColNames=="CloneID")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://nciarray.nci.nih.gov/cgi-bin/clone_report.cgi?CRITERIA=clone&PARAMETER=",
    ifelse(regexpr(":",temp) > 0, temp, paste("IMAGE:", temp, sep="")), "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT PROBE SET IDS.
  indic <- (ColNames=="ProbeSet")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    #temp <- ifelse(temp=="","",paste("<A HREF=\"https://www.affymetrix.com/LinkServlet?probeset=", temp, "\" TARGET=_blank>", temp, "</A>", sep=""))
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT GENBANK ACCESSION NUMBERS.
  indic <- (ColNames=="Accession")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Nucleotide&term=", temp, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT MicroRNA ACCESSION NUMBERS.
  indic <- (ColNames=="MicroRNAID")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?id=", temp, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT Entrez ID.
  indic <- (ColNames=="EntrezID")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=search&term=", temp, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT UNIGENE CLUSTER IDS.
  indic <- (ColNames=="UGCluster")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    PeriodIndex <- regexpr("\\.",temp)
    Org <- ifelse(PeriodIndex==-1, "", substring(temp, 1, PeriodIndex-1))
    Number <- ifelse(PeriodIndex==-1, "", substring(temp, PeriodIndex+1))
    temp <- paste("<A HREF=\"http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=", Org, "&CID=", Number, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT Gene Symbol.
  indic <- (ColNames %in% c("Gene symbol", "Symbol"))
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=gene&term=", temp, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## HYPERTEXT Ensembl ID.
  indic <- (ColNames == "Ensembl")
  if (sum(indic) > 0) {
    temp <- GeneId[ , (1:NumGeneIds)[indic],drop=F ]
    temp <- paste("<A HREF=\"http://www.ncbi.nlm.nih.gov/gene/?term=", temp, "\" TARGET=_blank>", temp, "</A>", sep="")
    GeneId[ , (1:NumGeneIds)[indic] ] <- temp
  }

  ## CREATE GENE IDS AS A VECTOR OF STRINGS.
  n <- nrow(GeneId)
  temp <- rep( "", n+1 )  #Each row in vector will contain HTML code for one row in resulting table, plus header row.

  # Create HTML code for row numbers.
  if (!missing(ColumnNum)) ColumnNum <- ColumnNum + 1
  GeneId <- if (is.null(GeneRowNum)) {
              cbind( paste("<TR><TD>", 1:n), GeneId )
            } else {
              cbind( paste("<TR><TD>", GeneRowNum), GeneId )
            }

  # Create HTML code for header row.
  # force unwrap in "t-value"
  temp[1] <- ""
  for(i in 1:length(ColNames))
      temp[1] <- paste(temp[1], ifelse(ColNames[i]=="t-value", "<TH NOWRAP>", "<TH>"), ColNames[i])
  temp[1] <- paste("<TR><TH>&nbsp;", temp[1], collapse="") #Qian 6/1/09
  NGeneIds <- ncol(GeneId)
  if (missing(ColumnNum)) {
     for (i in 1:n)
       temp[i+1] <- paste(paste(GeneId[i,], collapse="<BR><TD>"),"<BR>",sep="")
  } else {
     for (i in 1:n) {
       temp[i+1] <- paste(paste(GeneId[i,1:(ColumnNum-1)], collapse="<TD>"),
                    paste(GeneId[i,ColumnNum:NGeneIds], collapse="<TD NOWRAP ALIGN=\"RIGHT\">"),
                    sep="<TD NOWRAP ALIGN=\"RIGHT\">")
    }
  }

  ## GENERATE THE HTML CODE.

  GeneTableHtmlCode <- ""
  ChromDistHtmlCode <- ""

  GeneTableHtmlCode <- c(
      ifelse(length(grep("sorted and grouped", TableCaption)) > 0 #Qian 6/2/09 sortable
    , "<table border='1'>" #Qian 6/2/09 sortable
    , "<table border='1' class='sortable'>") #Qian 6/2/09 sortable
    , "<CAPTION style='text-align:left'>" #use style='text-align:left' to replace ALIGN=Left because the former is both supported in IE and Mozilla Firefox.
    , TableCaption
    , "</CAPTION>"
    , temp
    , "</TABLE>"
  )

  if (RunChromDistBarPlot){
    if (file.access(ChromDistBarPlotFile, mode=4)==0) {  #If file is read-accessible.
      ChromDistHtmlCode <- c("<IMG SRC=\"", ChromDistFileName, "\">")
    }
  }


 ## EITHER RETURN THE HTML CODE, OR WRITE THE HTML CODE TO A FILE.

  if (ReturnHtmlCode) {

    return(
      list(GeneTableHtmlCode = GeneTableHtmlCode
        , GOTableHtmlCode = OEhtml
        , ChromDistHtmlCode = ChromDistHtmlCode
        , DoAnnotations = DoAnnotations
      )
    )

  } else {

    cat("<!-- This is an HTML document.  Please make sure the filename has either"
      , "     an .htm or .html extension, and view this file in a web browser. -->"
      , ""
      , "<HTML>"
      , "<HEAD>"
      , "<STYLE>"
      , "P {font-size: 9pt}"
      , "TABLE {font-size: 9pt}"
      , "</STYLE>"
      , paste("<script src='file:///", ArrayToolsPath, "/Misc/sorttable.js'></script>", sep="") #Qian 5/29/09 sortable
      , "</HEAD>"
      , "<BODY>"
      ," ", file=GeneTableFile, sep="\n", append=F
    )
    HeaderStringsVec <- paste(HeaderStringsVec, collapse="<BR>")
    cat(HeaderStringsVec, "<BR><BR>", file=GeneTableFile, sep="\n", append=T)
    cat( "<HR SIZE=5 WIDTH=\"100%\" NOSHADE>"
      , "<BR><B><U>Gene table:</U></B><BR><BR>"
      , GeneTableHtmlCode
      , file=GeneTableFile, sep="\n", append=T
    )
    if (RunGOOvE) {
      cat( "<P></P>"
        , "<HR SIZE=5 WIDTH=\"100%\" NOSHADE>"
        , "<P></P><B><U>"
        , paste("'Observed v. Expected' table of GO classes and parent classes,",
            "in list of", nrow(GeneTable), "genes shown above:")
        , "</U></B><P></P>"
        , OEhtml
        , file=GeneTableFile, sep="\n", append=T
      )
    }
    if (RunChromDistBarPlot){
      if (file.access(ChromDistBarPlotFile, mode=4)==0) {  #If file is read-accessible.
        cat( "<P></P>"
          , "<HR SIZE=5 WIDTH=\"100%\" NOSHADE>"
          , "<P></P>"
          , ChromDistHtmlCode
          , file=GeneTableFile, sep="\n", append=T
        )
      }
    }
    cat( "<P></P>"
      , "</BODY>"
      , "</HTML>"
      , file=GeneTableFile, sep="\n", append=T
    )
  }

}

CreateAnnotations <- function(GeneIds, PathName="", AnnotationsFileName="", CallFromPlugins=F, DoMasterAnnotations=T,
                              ProjectPath = "") {
  # Called by CreateGeneTable() and CreateMasterAnnotationsTable()

  ArrayToolsPath <- path.package("classpredict")

  PrimaryIdFile <- paste(ProjectPath, "/Annotations/PrimaryId.txt", sep="")
  if (file.access(PrimaryIdFile) < 0) return(NULL)
  PrimaryIdName <- scan(PrimaryIdFile, what="", sep="\n", quiet = T)
  PrimaryIdFile <- paste(ProjectPath, "/Annotations/", PrimaryIdName, ".txt", sep="")
  if (file.access(PrimaryIdFile) < 0) return(NULL)

  if (DoMasterAnnotations==T) { # Yong 8/21/06
    NameFile <- paste(ProjectPath, "/Annotations/Name.txt", sep="")
    GBaccFile <- paste(ProjectPath, "/Annotations/Accession.txt", sep="")
    UniGeneFile <- paste(ProjectPath, "/Annotations/UGCluster.txt", sep="")
    SymbolFile <- paste(ProjectPath, "/Annotations/Symbol.txt", sep="")
    LLIDFile <- paste(ProjectPath, "/Annotations/LLID.txt", sep="")
    if (!file.exists(LLIDFile)) {
          LLIDFile <- paste(ProjectPath, "/Annotations/GeneID.txt", sep="")
    }
    if (!file.exists(LLIDFile)) {
          LLIDFile <- paste(ProjectPath, "/Annotations/EntrezID.txt", sep="")     #Lori 12/22/09
    }
    GOFile <- paste(ProjectPath, "/Annotations/GO.txt", sep="")
    ChromosomeFile <- paste(ProjectPath, "/Annotations/Chromosome.txt", sep="")
    CytobandFile <- paste(ProjectPath, "/Annotations/Cytoband.txt", sep="")
    SumFuncFile <- paste(ProjectPath, "/Annotations/SumFunc.txt", sep="")
    SPFunctionFile <- paste(ProjectPath, "/Annotations/SPFunction.txt", sep="")
    SPLocalFile <- paste(ProjectPath, "/Annotations/SPLocal.txt", sep="")
    PathwayFile <- paste(ArrayToolsPath, "/Pathways/Output/Genes.txt", sep="")
    BiocartaFile <- paste(ArrayToolsPath, "/Pathways/Output/Names_Biocarta.txt", sep="")
    KeggFile <- paste(ArrayToolsPath, "/Pathways/Output/Names_KEGG.txt", sep="")
    OrganismFile <- paste(ProjectPath, "/Annotations/Organism.txt", sep="")
    BroadFile <- paste(ArrayToolsPath, "/Pathways/Output/Names_Broad.txt", sep="")

    NameExists <- ifelse(file.access(NameFile) < 0, F, T)
    GBaccExists <- ifelse(file.access(GBaccFile) < 0, F, T)
    UniGeneExists <- ifelse(file.access(UniGeneFile) < 0, F, T)
    SymbolExists <- ifelse(file.access(SymbolFile) < 0, F, T)
    LLIDExists <- ifelse(file.access(LLIDFile) < 0, F, T)
    GOExists <- ifelse(file.access(GOFile) < 0, F, T)
    ChromosomeExists <- ifelse(file.access(ChromosomeFile) < 0, F, T)
    CytobandExists <- ifelse(file.access(CytobandFile) < 0, F, T)
    SumFuncExists <- ifelse(file.access(SumFuncFile) < 0, F, T)
    SPFunctionExists <- ifelse(file.access(SPFunctionFile) < 0, F, T)
    SPLocalExists <- ifelse(file.access(SPLocalFile) < 0, F, T)
    OrganismExists <- ifelse(file.access(OrganismFile) < 0, F, T)

    #PathwayExists <- ifelse((file.access(PathwayFile) < 0)|(file.access(BiocartaFile) < 0)|(file.access(KeggFile) < 0), F, T)
    #Deyun modified 12/1/05 for adding Broad pathway; Amy modified 1/19/06 to make sure downloaded "Genes.txt" file has 5 columns.
    FileFind <- file.exists(BroadFile) &&
                file.exists(PathwayFile) &&
                (length(scan(PathwayFile, sep="\t", what="", nlines=1, quiet=T))==5)
    if (FileFind){
    	  PathwayExists <- ifelse((file.access(PathwayFile) < 0)|(file.access(BiocartaFile) < 0)|(file.access(KeggFile) < 0)|(file.access(BroadFile) < 0), F, T)
      } else {
  	    PathwayExists <- ifelse((file.access(PathwayFile) < 0)|(file.access(BiocartaFile) < 0)|(file.access(KeggFile) < 0), F, T)
    }


    PrimaryId <- read.table(file=PrimaryIdFile, sep="\t", header=T, fill=T, as.is=T, comment.char = "", quote="")  #Deyun modified 4/20/06, MLI added quote="" 5/28/08
    index1 <- match(GeneIds,PrimaryId[,2])

    PrimaryId <- PrimaryId[index1,2]
    if (NameExists) {
      Name <- scan(file=NameFile, sep="\n", skip=1, what="", quiet=T)
      Name <- matrix(unlist(lapply(strsplit(Name, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      Name <- gsub("@#$%&", "\'", Name, fixed=T)                    #Replace by backslash-apostrophe combination.
      Name <- gsub("@#$%!", "\"", Name, fixed=T)                    #Replace by backslash-double-quote combination.
    }
    if (GBaccExists) {
      GBacc <- scan(file=GBaccFile, sep="\n", skip=1, what="", quiet=T)
      GBacc <- matrix(unlist(lapply(strsplit(GBacc, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      GBacc <- gsub("@#$%&", "'", GBacc, fixed=T)              #Replace by backslash-apostrophe combination.
      GBacc <- gsub("@#$%!", "\"", GBacc, fixed=T)             #Replace by backslash-double-quote combination.
    }
    if (UniGeneExists) {
      UniGene <- scan(file=UniGeneFile, sep="\n", skip=1, what="", quiet=T)
      UniGene <- matrix(unlist(lapply(strsplit(UniGene, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      UniGene <- gsub("@#$%&", "'", UniGene, fixed=T)              #Replace by backslash-apostrophe combination.
      UniGene <- gsub("@#$%!", "\"", UniGene, fixed=T)              #Replace by backslash-double-quote combination.
    }
    if (SymbolExists) {
      Symbol <- scan(file=SymbolFile, sep="\n", skip=1, what="", quiet=T)
      Symbol <- matrix(unlist(lapply(strsplit(Symbol, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      Symbol <- gsub("@#$%&", "'", Symbol, fixed=T)                #Replace by backslash-apostrophe combination.
      Symbol <- gsub("@#$%!", "\"", Symbol, fixed=T)                #Replace by backslash-double-quote combination.
    }
    if (LLIDExists) {
      LLID <- scan(file=LLIDFile, sep="\n", skip=1, what="", quiet=T)
      LLID <- matrix(unlist(lapply(strsplit(LLID, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      LLID <- gsub("@#$%&", "'", LLID, fixed=T)                    #Replace by backslash-apostrophe combination.
      LLID <- gsub("@#$%!", "\"", LLID, fixed=T)                    #Replace by backslash-double-quote combination.
    }
    if (GOExists) {
      GO <- scan(file=GOFile, sep="\n", skip=1, what="", quiet=T)
      GO <- matrix(unlist(lapply(strsplit(GO, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      GO <- gsub("@#$%&", "'", GO, fixed=T)                        #Replace by backslash-apostrophe combination.
      GO <- gsub("@#$%!", "\"", GO, fixed=T)                        #Replace by backslash-double-quote combination.
    }
    if (ChromosomeExists) {
      Chromosome <- scan(file=ChromosomeFile, sep="\n", skip=1, what="", quiet=T)
      Chromosome <- matrix(unlist(lapply(strsplit(Chromosome, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      Chromosome <- gsub("@#$%&", "'", Chromosome, fixed=T)        #Replace by backslash-apostrophe combination.
      Chromosome <- gsub("@#$%!", "\"", Chromosome, fixed=T)        #Replace by backslash-double-quote combination.
    }
    if (CytobandExists) {
      Cytoband <- scan(file=CytobandFile, sep="\n", skip=1, what="", quiet=T)
      Cytoband <- matrix(unlist(lapply(strsplit(Cytoband, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      Cytoband <- gsub("@#$%&", "'", Cytoband, fixed=T)            #Replace by backslash-apostrophe combination.
      Cytoband <- gsub("@#$%!", "\"", Cytoband, fixed=T)            #Replace by backslash-double-quote combination.
    }
    if (SumFuncExists) {
      SumFunc <- scan(file=SumFuncFile, sep="\n", skip=1, what="", quiet=T)
      SumFunc <- matrix(unlist(lapply(strsplit(SumFunc, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      SumFunc <- gsub("@#$%&", "'", SumFunc, fixed=T)              #Replace by backslash-apostrophe combination.
      SumFunc <- gsub("@#$%!", "\"", SumFunc, fixed=T)              #Replace by backslash-double-quote combination.
    }
    if (SPFunctionExists) {
      SPFunction <- scan(file=SPFunctionFile, sep="\n", skip=1, what="", quiet=T)
      SPFunction <- matrix(unlist(lapply(strsplit(SPFunction, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      SPFunction <- gsub("@#$%&", "'", SPFunction, fixed=T)        #Replace by backslash-apostrophe combination.
      SPFunction <- gsub("@#$%!", "\"", SPFunction, fixed=T)        #Replace by backslash-double-quote combination.
    }
    if (SPLocalExists) {
      SPLocal <- scan(file=SPLocalFile, sep="\n", skip=1, what="", quiet=T)
      SPLocal <- matrix(unlist(lapply(strsplit(SPLocal, "\t"), function(x) if (length(x)==1) c(x, "") else x)), ncol=2, byrow=T)[index1,2]
      SPLocal <- gsub("@#$%&", "'", SPLocal, fixed=T)              #Replace by backslash-apostrophe combination.
      SPLocal <- gsub("@#$%!", "\"", SPLocal, fixed=T)              #Replace by backslash-double-quote combination.
    }

    BlankStr <- rep("", length(PrimaryId))
    TDStr <- rep("<TD>", length(PrimaryId))

    GeneIds1 <- sub(" ", "", GeneIds)  #Get rid of all blanks in GeneId names to be used as anchors.
    GeneIds1 <- sub("\\.", "", GeneIds1)  #Get rid of all periods in GeneId names to be used as anchors.

    GeneIdStr <- paste("<TR><TD VALIGN=TOP><a name=", GeneIds1, "><b>", PrimaryIdName, "</b>: ", GeneIds, "<BR>", sep="")

    if (NameExists) {
      NameStr <- ifelse(Name=="" | is.na(Name) , "", paste("<BR><b>Name: </b>", Name, ""))
    } else {
      NameStr <- BlankStr
    }

    if (GBaccExists) {
      GBaccStr <- ifelse(GBacc=="" | is.na(GBacc) , "", paste("<BR><b>Accession: </b><a href=http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Nucleotide&term=", GBacc, ">", GBacc, "</a>", sep=""))
    } else {
      GBaccStr <- BlankStr
    }

    if (UniGeneExists & OrganismExists) {
      UniGenePrefix <- scan(OrganismFile, what="", sep="\n", quiet=T)
      if(toupper(UniGenePrefix) %in% c("HUMAN", "HS", "HOMO SAPIENS")) UniGenePrefix <- "Hs"
      if(toupper(UniGenePrefix) %in% c("MOUSE", "MM", "MUS MUSCULUS")) UniGenePrefix  <-"Mm"
      if(toupper(UniGenePrefix) %in% c("RAT", "RN", "RATTUS NORVEGICUS")) UniGenePrefix <- "Rn"
      UniGeneStr <- ifelse(is.na(UniGene) | (UniGene %in% c(""," ","No Valid ID","Data not found")), "", ifelse(UniGene=="In multiple clusters", "<BR><b>UniGene: </b>In multiple clusters", paste("<BR><b>UniGene: </b><a href=http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?CMD=search&DB=unigene&term=", UniGene, ">", UniGene, "</a>", sep="")))
      Urlfile <- read.delim(paste(ArrayToolsPath, "/Misc/url.txt", sep=""), sep="\t", header=F)
    	SOURCEInd <- which(Urlfile[,1]=="sourcemain")
    	SOURCEUrl <- as.character(Urlfile[SOURCEInd, 2])
    	SOURCEStr <- ifelse(is.na(UniGene) | (UniGene %in% c(""," ","No Valid ID","Data not found","In multiple clusters")), "", paste("<BR><a href=", SOURCEUrl, "/cgi-bin/source/sourceResult?choice=Gene&option=CLUSTER&organism=", UniGenePrefix, "&criteria=", UniGene, ">SOURCE</a>", sep="")) #Qian 3/17/09
    } else {
      UniGeneStr <- BlankStr
      SOURCEStr <- BlankStr   #Lori 12/14/09
    }

    if (SymbolExists) {
      SymbolStr <- ifelse(Symbol=="" | is.na(Symbol), "", paste("<BR><b>Symbol: </b><a href=http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=gene&term=", Symbol, ">", Symbol, "</a>", sep=""))
      GENECARDStr <- ifelse(Symbol=="" | is.na(Symbol), "", paste("<BR><a href=http://www.genecards.org/cgi-bin/carddisp?", Symbol, ">GENECARD</a>", sep="")) #Lori 10/12/2015
      DRUGBANKStr <- ifelse(Symbol=="" | is.na(Symbol), "", paste("<BR><b>DrugBank: </b><a href=http://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&query=", Symbol, "&searcher=targets&button=>", Symbol, "</a>", sep="")) #Qian 8/27/2015 DrugBank 4.2 and 4.3.
  	} else {
      SymbolStr <- BlankStr
      GENECARDStr <- BlankStr #Lori 12/14/09
      DRUGBANKStr <- BlankStr #Lori 12/14/09
    }

    if (LLIDExists) {
      LLIDStr <- ifelse(LLID=="" | is.na(LLID), "", paste("<BR><b>EntrezID: </b><a href=http://www.ncbi.nlm.nih.gov/entrez?db=gene&cmd=Retrieve&dopt=summary&list_uids=", LLID, ">", LLID, "</a>", sep=""))
    } else {
      LLIDStr <- BlankStr
    }

    if (ChromosomeExists) {
      ChromosomeStr <- ifelse(Chromosome=="" | is.na(Chromosome), "", paste("<BR><b>Chromosome: </b>", Chromosome, ""))
    } else {
      ChromosomeStr <- BlankStr
    }

    if (CytobandExists) {
      CytobandStr <- ifelse(Cytoband=="" | is.na(Cytoband), "", paste("<BR><b>Cytoband: </b>", Cytoband, ""))
    } else {
      CytobandStr <- BlankStr
    }

    if (SumFuncExists) {
      SumFuncStr <- ifelse(SumFunc=="" | is.na(SumFunc), "<TD VALIGN=TOP><b>SumFunc: </b>", paste("<TD VALIGN=TOP><b>SumFunc: </b>", SumFunc, ""))
    } else {
      SumFuncStr <- BlankStr
    }

    if (SPFunctionExists) {
      SPFunctionStr <- ifelse(SPFunction=="" | is.na(SPFunction), "<TD VALIGN=TOP><b>SPFunction: </b>", paste("<TD VALIGN=TOP><b>SPFunction: </b>", SPFunction, ""))
    } else {
      SPFunctionStr <- BlankStr
    }

    if (SPLocalExists) {
      SPLocalStr <- ifelse(SPLocal=="" | is.na(SPLocal), "<TD VALIGN=TOP><b>SPLocal: </b>", paste("<TD VALIGN=TOP><b>SPLocal: </b>", SPLocal, ""))
    } else {
      SPLocalStr <- BlankStr
    }

    BiocartaStr <- BlankStr
    KeggStr <- BlankStr
    #Deyun added 12/1/05 for Broad pathway
    if (FileFind){
    	BroadStr <- BlankStr
    }
    if (UniGeneExists) {
      if (PathwayExists) {
        ## READ 'GENES.TXT' FILE.
        PathwayTable <- scan(PathwayFile, sep="\t", what="", quiet=T)
        PathwayTable <- gsub("@#$%&", "'", PathwayTable, fixed=T)  #Replace by backslash-apostrophe combination.
        PathwayTable <- gsub("@#$%!", "\"", PathwayTable, fixed=T)  #Replace by backslash-double-quote combination.
        #Deyun modified 12/02/05 for Broad pathway
        if (FileFind){
          PathwayTable <- matrix(PathwayTable, ncol=5, byrow=T)   #Note extra EOL column at end.
        } else {
          PathwayTable <- matrix(PathwayTable, ncol=4, byrow=T)   #Note extra EOL column at end.
        }
        dimnames(PathwayTable) <- list(PathwayTable[,1], PathwayTable[1,])
        PathwayTable <- as.data.frame(PathwayTable[2:nrow(PathwayTable),2:(ncol(PathwayTable)-1)])
        ## READ 'NAMES_BIOCARTA.TXT' FILE.
        BiocartaTable <- scan(BiocartaFile, sep="\t", what="", quiet=T)
        BiocartaTable <- gsub("@#$%&", "'", BiocartaTable, fixed=T)  #Replace by backslash-apostrophe combination.
        BiocartaTable <- gsub("@#$%!", "\"", BiocartaTable, fixed=T)  #Replace by backslash-double-quote combination.
        BiocartaTable <- matrix(BiocartaTable, ncol=2, byrow=T)   #Note extra EOL column at end.
        ## READ 'NAMES_KEGG.TXT' FILE.
        KeggTable <- scan(KeggFile, sep="\t", what="", quiet=T)
        KeggTable <- gsub("@#$%&", "'", KeggTable, fixed=T)  #Replace by backslash-apostrophe combination.
        KeggTable <- gsub("@#$%!", "\"", KeggTable, fixed=T)  #Replace by backslash-double-quote combination.
        KeggTable <- matrix(KeggTable, ncol=2, byrow=T)   #Note extra EOL column at end.
        #Deyun added 12/1/05 for Broad pathway
        if (FileFind){
          BroadTable <- scan(BroadFile, sep="\t", what="", quiet=T)
          BroadTable <- gsub("@#$%&", "'", BroadTable, fixed=T)  #Replace by backslash-apostrophe combination.
          BroadTable <- gsub("@#$%!", "\"", BroadTable, fixed=T)  #Replace by backslash-double-quote combination.
          BroadTable <- matrix(BroadTable, ncol=2, byrow=T)   #Note extra EOL column at end.
        }

        ## CREATE PATHWAYS.
        index1 <- match(UniGene, rownames(PathwayTable))
        Biocarta <- as.character(PathwayTable$Biocarta.Pathways)[index1]
        Kegg <- as.character(PathwayTable$KEGG.Pathways)[index1]
        #Deyun added 12/1/05 for Broad pathway
        if (FileFind){
                Broad <- as.character(PathwayTable$Broad.Pathways)[index1]
        }

        for (i in 1:length(index1)) {
          BiocartaStr[i] <- "<TD VALIGN=TOP><b>BioCarta Pathways:</b><BR>"
          KeggStr[i] <- "<TD VALIGN=TOP><b>KEGG Pathways:</b><BR>"
          #Deyun added 12/1/05 for Broad pathway
          if (FileFind){
        		BroadStr[i] <- "<TD VALIGN=TOP><b>Broad/MIT Pathways:</b><BR>"
          }
          if (!(is.na(index1[i]))) {  #index1 may be NA if UniGene id was not found in 'GENES.TXT' file.
            if (Biocarta[i] != "") {  #Biocarta may be blank if that gene did not have a BioCarta pathway.
              entries <- unlist(strsplit(Biocarta[i], split="\\|"))
              CurrentPathways <- match(entries, BiocartaTable[,1])
              CurrentPathways <- BiocartaTable[CurrentPathways,2]
              for (j in 1:length(entries)) {
                BiocartaStr[i] <- paste(BiocartaStr[i], "<BR>", j, ": <a href=http://cgap.nci.nih.gov/Pathways/BioCarta/", entries[j], ">", CurrentPathways[j], "</a><BR>", sep="")
              }
            }
            if (Kegg[i] != "") {  #Kegg may be blank if that gene did not have a KEGG pathway.
              entries <- unlist(strsplit(Kegg[i], split="\\|"))
              CurrentPathways <- match(entries, KeggTable[,1])
              CurrentPathways <- KeggTable[CurrentPathways,2]
              for (j in 1:length(entries)) {
              #KeggStr[i] <- paste(KeggStr[i], "<BR>", j, ": <a href=http://cgap.nci.nih.gov/Pathways/Kegg/", entries[j], ">", CurrentPathways[j], "</a><BR>", sep="")
              KeggStr[i] <- paste(KeggStr[i], "<BR>", j, ": <a href=http://www.kegg.jp/kegg-bin/show_pathway?", entries[j], ">", CurrentPathways[j], "</a><BR>", sep="")
              }
            }

            #Deyun added 12/1/05 for Broad pathway
        	  if (FileFind){
        		  if (Broad[i] != "") {  #Broad may be blank if that gene did not have a Broad pathway.
        		    entries <- unlist(strsplit(Broad[i], split="\\|"))
        		    CurrentPathways <- match(entries, BroadTable[,1])
        		    CurrentPathways <- BroadTable[CurrentPathways,2]
        		    for (j in 1:length(entries)) {
        		      BroadStr[i] <- paste(BroadStr[i], "<BR>", j, ": ", entries[j], ifelse(is.na(CurrentPathways[j]),"",CurrentPathways[j]), "<BR>", sep="")
        		    }
        		  }
        	  } # end of if (FileFind)

          } # end of if (!(is.na(index1[i])))
        } # end of for (i in 1:length(index1))
      } # end of if (PathwayExists)
    } # end of if (UniGeneExists)

    #	 AnnotationStr <- ifelse(NameStr=="" & GBaccStr=="" & UniGeneStr=="" & SymbolStr=="" & LLIDStr=="" & ChromosomeStr=="" & CytobandStr=="" & SOURCEStr=="" & GENECARDStr=="" & SumFuncStr=="" & SPFunctionStr=="" & SPLocalStr=="" & GOStr=="" & BiocartaStr=="" & KeggStr=="", "",
    #	 paste(GeneIdStr,NameStr,GBaccStr,UniGeneStr,SymbolStr,LLIDStr,ChromosomeStr,CytobandStr,SOURCEStr,GENECARDStr,SumFuncStr,SPFunctionStr,SPLocalStr,GOStr,BiocartaStr,KeggStr))
    #Deyun modified 12/1/05 for Broad pathway
     if (FileFind){
       AnnotationStr <- ifelse(NameStr=="" & GBaccStr=="" & UniGeneStr=="" & SymbolStr=="" & LLIDStr=="" & ChromosomeStr=="" & CytobandStr=="" & SOURCEStr=="" & GENECARDStr=="" & DRUGBANKStr=="" & SumFuncStr=="" & SPFunctionStr=="" & SPLocalStr=="" & #GOStr=="" & BiocartaStr=="" & KeggStr=="" & BroadStr=="", "",
  	 BiocartaStr=="" & KeggStr=="" & BroadStr=="", "", ## Ting removed GO column, 4/11/2017
       paste(GeneIdStr,NameStr,GBaccStr,UniGeneStr,SymbolStr,LLIDStr,ChromosomeStr,CytobandStr,SOURCEStr,GENECARDStr,DRUGBANKStr,#SumFuncStr,SPFunctionStr,SPLocalStr,GOStr,BiocartaStr,KeggStr,BroadStr))
  	 SumFuncStr,SPFunctionStr,SPLocalStr,BiocartaStr,KeggStr,BroadStr)) ## Ting removed GO column, 4/11/2017
      } else {
    	 AnnotationStr <- ifelse(NameStr=="" & GBaccStr=="" & UniGeneStr=="" & SymbolStr=="" & LLIDStr=="" & ChromosomeStr=="" & CytobandStr=="" & SOURCEStr=="" & GENECARDStr=="" & DRUGBANKStr=="" & SumFuncStr=="" & SPFunctionStr=="" & SPLocalStr=="" & #GOStr=="" & BiocartaStr=="" & KeggStr=="", "",
  	 BiocartaStr=="" & KeggStr=="", "", ## Ting removed GO column, 4/11/2017
    	 paste(GeneIdStr,NameStr,GBaccStr,UniGeneStr,SymbolStr,LLIDStr,ChromosomeStr,CytobandStr,SOURCEStr,GENECARDStr,DRUGBANKStr,#SumFuncStr,SPFunctionStr,SPLocalStr,GOStr,BiocartaStr,KeggStr))
  	 SumFuncStr,SPFunctionStr,SPLocalStr,BiocartaStr,KeggStr)) ## Ting removed GO column, 4/11/2017
      }

    # Amy modified 4/19/04.
    if (CallFromPlugins==F) {
      AnnotationsFileName <- "Annotations.html"
      OutFile <- paste(ProjectPath, "/Output/", AnalysisName, "/", AnnotationsFileName, sep="")
    } else {
      if (PathName=="") PathName <- paste(ProjectPath, "/Output/Plugins", sep="")
      if (AnnotationsFileName=="") AnnotationsFileName <- "Annotations.html"
      OutFile <- paste(PathName, "/", AnnotationsFileName, sep="")
    }

    unlink(OutFile)

    cat("<HTML><HEAD><STYLE>TABLE {font-size: 9pt}</STYLE>"
      , paste("<script src='file:///", ArrayToolsPath, "/Misc/sorttable.js'></script>", sep="") #Qian 5/29/09 sortable
      , "</HEAD><BODY><table border='1' class='sortable'>" #Qian 5/29/09 sortable
      , ""
      , sep="\n", file=OutFile, append=F)

    #Deyun modified 12/2/05  for Broad pathway###
    if (FileFind){
    cat("<TR><TH > Gene Info </TH>", ifelse(SumFuncExists,"<TH > Sum Func </TH>",""), ifelse(SPFunctionExists,"<TH > Function </TH>"#,""), ifelse(SPLocalExists,"<TH > SP Local </TH>",""), ifelse(GOExists,"<TH > Gene Ontology </TH>","")
    ,""), ifelse(SPLocalExists,"<TH > SP Local </TH>","") ## Ting removed GO column, 4/11/2017
      , ifelse((UniGeneExists & PathwayExists),"<TH >  Biocarta Pathways </TH><TH >    Kegg Pathways </TH><TH >   Broad Pathways   </TH></TR>","")
      , sep="\n", file=OutFile, append=T)
    } else {
    cat("<TR><TH > Gene Info </TH>", ifelse(SumFuncExists,"<TH > Sum Func </TH>",""), ifelse(SPFunctionExists,"<TH > Function </TH>",#""), ifelse(SPLocalExists,"<TH > SP Local </TH>",""), ifelse(GOExists,"<TH > Gene Ontology </TH>","")
    ""), ifelse(SPLocalExists,"<TH > SP Local </TH>","") ## Ting removed GO column, 4/11/2017
      , ifelse((UniGeneExists & PathwayExists),"<TH >  Biocarta Pathways </TH><TH >    Kegg Pathways </TH>","")
      , sep="\n", file=OutFile, append=T)
    }

    cat(AnnotationStr, "</TABLE></BODY></HTML>"
      , sep="\n", file=OutFile, append=T)

    AnnotationStr <- paste("<a href=", AnnotationsFileName, "#", GeneIds1, ">Info</a>", sep="")
    return(AnnotationStr)
  } else {
 # DoMasterAnnotations is False, Yong  8/21/06
  GeneIds1 <- sub(" ", "", GeneIds)  #Get rid of all blanks in GeneId names to be used as anchors.
  GeneIds1 <- sub("\\.", "", GeneIds1)  #Get rid of all periods in GeneId names to be used as anchors.
  AnnotationStr <- paste("<a href=", AnnotationsFileName, "#", GeneIds1, ">Info</a>", sep="")
  return(AnnotationStr)
 }
}

#######################################################################
## Create.Class.Description.Table                                    ##
#######################################################################

## FUNCTION TO WRITE SENSITIVITY TABLES FOR CLASS PREDICTION ANALYSIS.

Create.Class.Description.Table <- function(SensitivityData,ClassLabels,groupid,digits=3) {
  Classes <- sort(unique(SensitivityData[[1]]))
  nClass  <- length(Classes)
  table <- data.frame(Class = Classes, Sensitivity = double(nClass),
                      Specificity = double(nClass), PPV = double(nClass),
                      NPV = double(nClass))
  for(i in 1:nClass) {
    iClass <- Classes[i]
    SensitivityData.Class <- SensitivityData[[2]][SensitivityData[[1]] == iClass]

    table$Sensitivity[i] <- round(SensitivityData.Class[1],digits)
    table$Specificity[i] <- round(SensitivityData.Class[2],digits)
    table$PPV[i] <- ifelse(SensitivityData.Class[3] > -0.5,
                           round(SensitivityData.Class[3],digits),NA)
    table$NPV[i] <- ifelse(SensitivityData.Class[4] > -0.5,
                           round(SensitivityData.Class[4],digits),NA)
  }
  HTMLtable <- character(nClass+1)
  HTMLtable[1] <- paste("<TR><TH><BR>",paste(names(table), collapse="<TH>"), "</BR>",collapse="")
  for (i in 1:nClass) {
    HTMLtable[i+1] <- paste( paste("<TR><TH>",ClassLabels[table$Class[i] ] ), "<TD>",
                             paste(table[i,-1], collapse="<TD>") )
  }
  return(list(table=table,HTMLtable=HTMLtable))
}

outputClassname <- function(classname, classsize) { # a helper function, Yong 12/15/06
  # used in <GetPredictionResults.R>
  str1 <- ""
  bclasssize <- !missing(classsize)
  for (i in 1:length(classname)) {
   str2 <- paste("<font size=2>Class ", i, ": ", "<I>", classname[i], "</I>",
                 if (bclasssize) paste0(" (", classsize[i], " samples)"), "</font>", sep="") # add classsize, MLI 5/16/17
   str1 <- paste(str1, str2, sep="")
   if (i < length(classname)) str1 <- paste(str1, "; ", sep="")
  }
  str1 <- paste(str1, ".", sep="")
}
