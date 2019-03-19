\name{ScoresDIA}
\alias{ScoresDIA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical analysis for a pair of peak grouped metabolites from LC-MS/MS DIA analysis.
}
\description{
Peak-to-peak Pearson correlation coefficient, peak-to-peak shape ratio and product/precursor ion intensity ratios are calculated for a product and precursor metabolites from LC-MS/MS DIA experiment.
}
\usage{
ScoresDIA(input,file,ID1,ID2,CE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{input}{
   Peak grouped for a particular metabolite obtained with the \link[MetaboList]{PeakGroup}.
}
  \item{file}{
  LC-MS/MS DIA file processed by the \link[MetaboList]{AIF}.
}
  \item{ID1}{
     PeakID of the precursor ion metabolite.}
  \item{ID2}{
     PeakID of the product ion metabolite.}
  \item{CE}{
  numeric. Collision energy for the file processed.}
}


\value{
  \item{Score}{Peak-to-peak Pearson correlation coefficient for a pair of EIC peaks.}
  \item{IntensityRatio}{ Peak intensity ratio between product and precursor ion metabolite.}
  \item{AssymetriRatio}{Score for the chromatogram peak shape based on assymmetry factor.}

}


\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.

}
\examples{
library(MetaboList)

#CE.isolation("AIFpos1000-AIF.mzXML","fileposB")

#Reading the database.csv file:
# database<- read.csv("C:/database.csv")

#Processing peak-picking and annotation with default parameters

#aif5<-AIF(fileMS,fileMS2CE5,database,CE=5, ion_mode = "positive")
#aif10<-AIF(fileMS,fileMS2CE10,database,CE=10, ion_mode = "positive")

#Peakgroup<-PeakGroup(aif5,aif10)

#Scores5<-ScoresDIA(Peakgroup$Glutamine,aif5,ID1=90, ID2 = 95,CE=5)

}
