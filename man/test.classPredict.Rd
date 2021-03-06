% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classPredict.R
\name{test.classPredict}
\alias{test.classPredict}
\title{Test classpredict() function}
\usage{
test.classPredict(
  dataset = c("Brca", "Perou", "Pomeroy"),
  projectPath,
  outputName = "ClassPrediction",
  generateHTML = FALSE
)
}
\arguments{
\item{dataset}{character string specifying one of "Brca", "Perou" or "Pomeroy" datasets.}

\item{projectPath}{characteer string specifying the project path. Default is C:/Users/UserName/Documents/$dataset.}

\item{outputName}{character string for the output folder name.}

\item{generateHTML}{logical. If \code{TRUE}, an html page will be generated with detailed prediction results.}
}
\value{
A list as returned by \code{classPredict}.
}
\description{
This function will load a test dataset to run \code{classPredict} function.
}
\details{
If the random variance model is enabled, all genes will be used for the model estimation.
}
\examples{
test.classPredict("Pomeroy")
}
\seealso{
\code{\link{classPredict}}
}
