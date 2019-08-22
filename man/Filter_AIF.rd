\name{Filter_AIF}
\alias{Filter_AIF}
\title{Automatic Metabolite Annotation from LC-MS DIA experiments.}
\description{Analysis and annotation of LC-MS/MS DIA data with the use of in-house mass spectral libraries.}
\usage{Filter_AIF(aif5,a=1,database,rttol=0.2)}
\arguments{
  \item{aif5}{
 DIA file processed
}

  \item{a}{
   Rule to restrict number of product ions.
}
  \item{database}{
  database employed for targeted annotation.
  }
   \item{rttol}{
  numeric. Retention time tolerance for aligning fragments and precursors.Default value=0.2min
  }

}


\author{Manuel D Peris Diaz}
\references{
1. R-MetaboList: a flexible tool for metabolite extraction from high-resolution data-independent acquisition mass spectrometry analysis. Metabolites. Soon

2. A Survey of Orbitrap All Ion Fragmentation Analysis Assessed by an R MetaboList Package to Study Small-Molecule Metabolites. Chromatographia. 2018, 81, 981-994.
}
\examples{
library(MetaboList)

CE.isolation("AIFpos1000-AIF.mzXML","fileposB")

#Reading the database.csv file:
# database<- read.csv("C:/database.csv")

#Processing peak-picking and annotation with default parameters

#aif5<-AIF(fileMS,fileMS2CE5,database,CE=5, ion_mode = "positive")

#Filter_AIF(aif5,full=TRUE,a=0,database)

}
