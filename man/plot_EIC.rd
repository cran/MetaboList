\name{plot_EIC}
\alias{plot_EIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot an Extracted Ion Chromatogram (EIC)
}
\description{
Plot an Extracted Ion Chromatogram (EIC) either from processed MS1 or MS/MS file with \link[MetaboList]{FullMS} or \link[MetaboList]{AIF}.
}
\usage{
plot_EIC(fullms,peakID=333,ms=1,CE=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullms}{
   A fullms or DIA file processed by \link[MetaboList]{FullMS} or \link[MetaboList]{AIF}.
}
  \item{peakID}{
Identity of the EIC desired to plot. The peakID is indicated in the ouput obtained with  \link[MetaboList]{FullMS} or \link[MetaboList]{AIF}.
  }
  \item{ms}{
  numeric. MS level of EIC desired to plot.
  }
  \item{CE}{
  numeric. Collision energy for the file processed.
  }



}


\author{Manuel D Peris Diaz}
\references{

1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.


}
\examples{
library(MetaboList)
#Reading the file.mzXML
# file<-fullMS.mzXML

#Reading the database.csv file:
# database<- read.csv("C:/FullMS1.csv")

#Processing peak-picking and annotation with default parameters
# FullMS_results<-FullMS(file,database, ion_mode = "positive",)

#Output:
#FullMS_results$ms

#plot_EIC(fullmsposH,peakID=413)


}
