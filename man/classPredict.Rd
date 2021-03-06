% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classPredict.R
\name{classPredict}
\alias{classPredict}
\title{Class prediction}
\usage{
classPredict(
  exprTrain,
  exprTest = NULL,
  isPaired = FALSE,
  pairVar.train = NULL,
  pairVar.test = NULL,
  geneId,
  cls,
  pmethod = c("ccp", "bcc", "dlda", "knn", "nc", "svm"),
  geneSelect = "igenes.univAlpha",
  univAlpha = 0.001,
  univMcr = 0.2,
  foldDiff = 2,
  rvm = FALSE,
  filter = NULL,
  ngenePairs = 25,
  nfrvm = 10,
  cvMethod = 1,
  kfoldValue = 10,
  bccPrior = 1,
  bccThresh = 0.8,
  nperm = 0,
  svmCost = 1,
  svmWeight = 1,
  fixseed = 1,
  prevalence = NULL,
  projectPath,
  outputName = "ClassPrediction",
  generateHTML = FALSE
)
}
\arguments{
\item{exprTrain}{matrix of gene expression data for training samples. Rows are genes and columns are arrays. Its column names must be provided.}

\item{exprTest}{matrix of gene expression data for new samples. Its column names must be provided.}

\item{isPaired}{logical. If \code{TRUE}, samples are paired.}

\item{pairVar.train}{vector of pairing variables for training samples.}

\item{pairVar.test}{vector of pairing variables for new samples.}

\item{geneId}{matrix/data frame of gene IDs.}

\item{cls}{vector of training sample classes.}

\item{pmethod}{character string vector of prediction methods to be employed.
\itemize{
  \item "ccp":  Compound Covariate Predictor
  \item "bcc":  Bayesian Compound Covariate Predictor
  \item "dlda": Diagonal Linear Discriminant Analysis
  \item "knn":  1-Nearest Neighbor/ 3-Nearest Neighbor
  \item "nc":   Nearest Centroid
  \item "svm":  Support Vector Machine
}}

\item{geneSelect}{character string for gene selection method.
\itemize{
  \item "igenes.univAlpha": select individual genes univariately significantly differentially expressed between the classes at the specified threshold significance level.
  \item "igenes.grid":      select individual genes that optimize over the grid of alpha levels.
  \item "igenes.univMcr":   select individual genes with univariate misclassification rate below a specified value.
  \item "gpairs":           select gene pairs bye the "greedy pairs" method.
  \item "rfe":              select genes by recursive feature elimination.
 }}

\item{univAlpha}{numeric for a significance level. Default is 0.001.}

\item{univMcr}{numeric for univariate misclassification rate. Default is 0.2.}

\item{foldDiff}{numeric for fold ratio of geometric means between two classes exceeding. 0 means not to enable this option. Default is 2.}

\item{rvm}{logical. If \code{TRUE}, random variance model will be employed. Default is \code{FALSE}.}

\item{filter}{vector of 1/0's of the same length as genes. 1 means to keep the gene while 0 means to exclude genes
from class comparison analysis. If \code{rvm = TRUE}, all genes will be used in random variance model estimation. Default is \code{FALSE}.}

\item{nfrvm}{numeric specifying the number of features selected by the support vector machine recursive feature elimination method. Default is 10.}

\item{cvMethod}{numeric for the cross validation method. Default is 1.
\itemize{
  \item 1:  leave-one-out CV,
  \item 2:  k-fold CV,
  \item 3:  0.632+ bootstrap.
  }}

\item{kfoldValue}{numeric specifying the number of folds if K-fold method is selected. Default is 10.}

\item{bccPrior}{numeric specifying the prior probability option for the Baysian compound covariate prediction.
If \code{bccPrior == 1}, equal prior probabilities will be applied.
If \code{bccPrior == 2}, prior probabilities based on the proportions in training data are applied.
Default is 1.}

\item{bccThresh}{numeric specifying the uncertainty threshold for the Bayesian compound covariate prediction. Default is 0.8.}

\item{nperm}{numeric specifying the number of permutations for the significance test of cross-validated mis-classification rate.
It should be equal to zero or greater than 50. Default is 0.}

\item{svmCost}{numeric specifying the cost values for SVM. Default is 1.}

\item{svmWeight}{numeric specifying the weight values for SVM. Default is 1.}

\item{fixseed}{numeric. \code{fixseed == 1} if a fixed seed is used; otherwise, \code{fixseed == 0}. Default is 1.}

\item{prevalence}{vector for class prevalences. When prevalence is \code{NULL}, the proportional of samples in each class will be the estimate of class prevalence. Default is \code{NULL}.
Names of vector should be provided and consistent with classes in \code{cls}.}

\item{projectPath}{character string specifying the full project path.}

\item{outputName}{character string specifying the output folder name. Default is "ClassPrediction".}

\item{generateHTML}{logical. If \code{TRUE}, an HTML page will be generated with detailed class prediction results saved in <projectPath>/Output/<outputName>/<outputName>.html.}

\item{ngenePairs:}{numeric specifying the number of gene pairs selected by the greedy pairs method. Default is 25.}
}
\value{
A list that may include the following objects:
\itemize{
   \item \code{performClass}: a data frame with the performance of classifiers during cross-validation:
   \item \code{percentCorrectClass}: a data frame with the mean percent of correct classification for each sample using
             different prediction methods.
   \item \code{predNewSamples}:s a data frame with predicted class for each
             new sample. `NC` means that a sample is not classified. In this example, there are four new samples.
   \item \code{probNew}: a data frame with the predicted probability of each new sample belonginG to the class (BRCA1) from the the Bayesian Compound Covariate method.
   \item \code{classifierTable}: a data frame with composition of classifiers such as geometric means of values in each class, p-values and Gene IDs.
   \item \code{probInClass}: a data frame with predicted probability of each training sample belonging to
             aclass during cross-validation from the Bayesian Compound Covariate
   \item \code{CCPSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                  negative prediction value) of the Compound Covariate Predictor Classifier.
   \item \code{LDASenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                                negative prediction value) of the Diagonal Linear Discriminant Analysis Classifier.
   \item \code{K1NNSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                   negative prediction value) of the 1-Nearest Neighbor Classifier.
   \item \code{K3NNSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                                 negative prediction value) of the 3-Nearest Neighbor Classifier.
   \item \code{CentroidSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                                     negative prediction value) of the Nearest Centroid Classifierr.
   \item \code{SVMSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                                negative prediction value) of the Support Vector Machine Classifier.
   \item \code{BCPPSenSpec}: a data frame with performance (i.e., sensitivity, specificity, positive prediction value,
                                                                 negative prediction value) of the Bayesian Compound Covariate Classifier.
   \item \code{weightLinearPred}: a data frame with gene weights for linear predictors such as Compound Covariate Predictor,
             Diagonal Linear Discriminat Analysis and Support Vector Machine.
   \item \code{thresholdLinearPred}: a numeric vector of the thresholds for the linear prediction rules related with \code{weightLinearPred}.
             Each prediction rule is defined by the inner sum of the weights (\eqn{w_i})
             and log expression values (\eqn{x_i}) of significant genes.
             In this case, a sample is classified to the class BRCA1 if
             the sum is greater than the threshold; that is, \eqn{\sum_i w_i x_i > threshold}.
   \item \code{GRPCentroid}: a data frame with centroid of each class for each predictor gene.
   \item \code{ppval}: a vector of permutation p-values of statistical significance tests of cross-validated estimate of misclassification rate from specified #'    prediction methods.
   \item \code{pmethod}: a vector of prediction methods that are specified.
   \item \code{workPath}: the path for fortran and other intermediate outputs.
   }
}
\description{
This function calculates multiple classifiers that are used to predict the class of a new sample.
It implements the class prediction tool with multiple methods in BRB-ArrayTools.
}
\details{
Please see the BRB-ArrayTools manual (\url{https://brb.nci.nih.gov/BRB-ArrayTools/Documentation.html}) for details.
}
\examples{
dataset<-"Brca"
# gene IDs
geneId <- read.delim(system.file("extdata", paste0(dataset, "_GENEID.txt"),
                     package = "classpredict"), as.is = TRUE, colClasses = "character")
# expression data
x <- read.delim(system.file("extdata", paste0(dataset, "_LOGRAT.TXT"),
                package = "classpredict"), header = FALSE)
# filter information, 1 - pass the filter, 0 - filtered
filter <- scan(system.file("extdata", paste0(dataset, "_FILTER.TXT"),
               package = "classpredict"), quiet = TRUE)
# class information
expdesign <- read.delim(system.file("extdata", paste0(dataset, "_EXPDESIGN.txt"),
                        package = "classpredict"), as.is = TRUE)
# training/test information
testSet <- expdesign[, 10]
trainingInd <- which(testSet == "training")
predictInd <- which(testSet == "predict")
ind1 <- which(expdesign[trainingInd, 4] == "BRCA1")
ind2 <- which(expdesign[trainingInd, 4] == "BRCA2")
ind <- c(ind1, ind2)
exprTrain <- x[, ind]
colnames(exprTrain) <- expdesign[ind, 1]
exprTest <- x[, predictInd]
colnames(exprTest) <- expdesign[predictInd, 1]
projectPath <- file.path(Sys.getenv("HOME"),"Brca")
outputName <- "ClassPrediction"
generateHTML <- TRUE
resList <- classPredict(exprTrain = exprTrain, exprTest = exprTest, isPaired = FALSE,
                        pairVar.train = NULL, pairVar.test = NULL, geneId,
                        cls = c(rep("BRCA1", length(ind1)), rep("BRCA2", length(ind2))),
                        pmethod = c("ccp", "bcc", "dlda", "knn", "nc", "svm"),
                        geneSelect = "igenes.univAlpha",
                        univAlpha = 0.001, univMcr = 0, foldDiff = 0, rvm = TRUE,
                        filter = filter, ngenePairs = 25, nfrvm = 10, cvMethod = 1,
                        kfoldValue = 10, bccPrior = 1, bccThresh = 0.8, nperm = 0,
                        svmCost = 1, svmWeight =1, fixseed = 1, prevalence = NULL,
                        projectPath = projectPath, outputName = outputName, generateHTML)
if (generateHTML)
  browseURL(file.path(projectPath, "Output", outputName,
            paste0(outputName, ".html")))


}
